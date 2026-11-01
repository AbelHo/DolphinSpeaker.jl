# Included by config.jl. Definitions only: calibrated parameters load on request.
if !isdefined(@__MODULE__, :ACTIVE_CLICKER_CALIBRATION)
    ACTIVE_CLICKER_CALIBRATION = nothing
end

active_clicker_calibration() = ACTIVE_CLICKER_CALIBRATION

# A preset may be activated and used inside the same compiled caller. Resolve
# both the newly loaded module binding and its methods in the latest world.
function _clicker_call(name::Symbol,args...;kwargs...)
    Base.invokelatest() do
        module_value=getfield(@__MODULE__,:ClickerCalibration)
        getfield(module_value,name)(args...;kwargs...)
    end
end

function _clicker_minmax_tdoas(data, peaks; window,ref_channel=3)
    size(data,2)==3 && ref_channel==3 || throw(ArgumentError("Calibrated D3 expects selected channels 3–5 and selected reference 3"))
    result=zeros(Float64,length(peaks),3)
    for (i,peak) in enumerate(peaks)
        lo=max(1,peak+first(window));hi=min(size(data,1),peak+last(window))
        result[i,:]=_clicker_call(:prep_minmax,@view(data[lo:hi,:]);corrected=true)
    end
    return result
end

# D3/Zoom F6 runs at 192 kHz, the same default as _mr_estimator_spec and prep_gcc.
_clicker_sample_rate(cal)=Float64(get(get(cal,"estimator_spec",Dict()),"sample_rate_hz",192000.0))

function _clicker_method_tdoas(data,peaks;window,ref_channel=3)
    calibration=active_clicker_calibration()
    calibration===nothing && error("No active method-matched calibration")
    return _clicker_call(:estimate_method_tdoas,calibration,data,peaks;fs=_clicker_sample_rate(calibration),window=window,ref_channel=ref_channel)
end

# The only lazy loader for ClickerCalibration.jl; also used by the DolphinSpeaker wrappers.
function _ensure_clicker_module()
    isdefined(@__MODULE__,:ClickerCalibration) ||
        Base.include(@__MODULE__,joinpath(@__DIR__,"ClickerCalibration.jl"))
    return nothing
end

# A function, not a const: config.jl and this file are included many times.
_clicker_default_path(device)=joinpath(@__DIR__,"calibrations","clicker_20260209_$(device).json")

# References (not copies) to every non-const global of the including module.
function _clicker_global_snapshot()
    mod=@__MODULE__
    return Dict{Symbol,Any}(n=>getglobal(mod,n) for n in names(mod;all=true)
        if !startswith(String(n),"#") && isdefined(mod,n) && !isconst(mod,n))
end

_clicker_changed_names(before)=Symbol[n for (n,v) in before if !(isdefined(@__MODULE__,n) && getglobal(@__MODULE__,n)===v)]

function _activate_clicker_calibration!(device, calibration)
    set_device__ophk_acoustic_D3_clicker()
    global ACTIVE_CLICKER_CALIBRATION=calibration
    global c=calibration["sound_speed_m_s"]
    global default_c=c
    global ref_channel=3
    global default_window=window_impulsive
    global imsize=(Int(calibration["height"]),Int(calibration["width"]))
    camera=calibration["camera"]
    global fov_angle=[rad2deg(2atan(imsize[2]/(2camera["fx"]))),
                      rad2deg(2atan(imsize[1]/(2camera["fy"]))),missing]
    global default_getTDOA_func=calibration["method"]=="minmax2" ? _clicker_minmax_tdoas : _clicker_method_tdoas
    global FILE_device_ID="ophk_acoustic_D3_clicker_$(device)"
    return calibration
end

"""
    with_clicker_calibration(f, device; calibration_path=_clicker_default_path(device))

Run `f()` with clicker `device` ("A" or "B") active, then restore every global the activation
changed. Returns `f()`'s value. Blocks nest and are exception safe; a bad path changes nothing.
Run calibrated work inside the block, leave no tasks running after it, and do not call `set_device__*` inside it.
"""
function with_clicker_calibration(f,device::AbstractString;calibration_path=_clicker_default_path(device))
    _ensure_clicker_module()
    calibration=_clicker_call(:load_calibration,calibration_path;device=device)  # validate before touching any global
    before=_clicker_global_snapshot()
    changed=Symbol[]
    try
        try
            _activate_clicker_calibration!(device,calibration)
        finally
            append!(changed,_clicker_changed_names(before))
        end
        return f()
    finally
        # Names first defined during activation are not in `before`; they are left as they are.
        for n in changed
            setglobal!(@__MODULE__,n,before[n])
        end
    end
end

function calibrated_tdoa2dir(tdoas,receivers,fs,calibration;return_residual=false)
    Float64(fs)==_clicker_sample_rate(calibration) || error("sample rate differs from active clicker calibration")
    geometry=reduce(vcat,[permutedims(Float64.(row)) for row in calibration["geometry_m"]])
    size(receivers)==(3,3) && maximum(abs.(receivers.-permutedims(geometry)))<1e-12 ||
        throw(ArgumentError("Receiver geometry differs from the active calibration"))
    size(tdoas,2)==3 || throw(ArgumentError("Expected three selected-channel delays"))
    relative=Float64.(tdoas).-Float64.(tdoas[:,3:3])
    projected=_clicker_call(:project_tdoas,calibration,relative;fs=fs)
    angles=fill(NaN,size(tdoas,1),2);residuals=fill(NaN,size(tdoas,1))
    offsets=Float64.(calibration["relative_delay_offsets_s"])
    for (i,result) in enumerate(projected)
        u=result.direction_array
        all(isfinite,u) || continue
        angles[i,:]=[atan(u[1],u[3]),asin(clamp(u[2],-1,1))]
        expected=-(geometry.-permutedims(geometry[3,:]))*u/calibration["sound_speed_m_s"]+offsets
        residuals[i]=sum(abs2,expected.*fs.-relative[i,:])
    end
    return return_residual ? (angles,residuals) : angles
end

function calibrated_angle2px(angles,calibration,shape;degrees=false)
    size(angles,2)==2 || throw(ArgumentError("Expected N×2 Julia azimuth/elevation angles"))
    records=_clicker_call(:project_angles,calibration,angles;units=degrees ? "degrees" : "radians")
    result=fill(NaN,size(angles,1),2)
    for (i,record) in enumerate(records)
        record.valid || continue
        result[i,1]=record.x*shape[2]/calibration["width"]
        result[i,2]=record.y*shape[1]/calibration["height"]
    end
    return result
end
