# Included by config.jl. Definitions only: calibrated parameters load on request.
# Defined once, so re-including config.jl keeps the current state.
isdefined(@__MODULE__, :ACTIVE_CLICKER_CALIBRATION) || (ACTIVE_CLICKER_CALIBRATION = nothing)
isdefined(@__MODULE__, :ACTIVE_CLICKER_TOKEN) || (ACTIVE_CLICKER_TOKEN = nothing)            # (rx_vect, fov_angle) set by the activation
isdefined(@__MODULE__, :CLICKER_UNCALIBRATED_TDOA) || (CLICKER_UNCALIBRATED_TDOA = nothing)  # default_getTDOA_func before activation

# The calibration counts only while the arrays its activation set are still in place. Any later
# set_device__* (or re-included config.jl) rebinds rx_vect or fov_angle, which switches it off.
function active_clicker_calibration()
    calibration=ACTIVE_CLICKER_CALIBRATION
    calibration===nothing && return nothing
    token=ACTIVE_CLICKER_TOKEN
    return token isa Tuple && token[1]===rx_vect && token[2]===fov_angle ? calibration : nothing
end

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

# default_getTDOA_func while a calibration is set: the calibrated estimator, or the previous function once switched off.
function _clicker_tdoas(data,peaks;kwargs...)
    calibration=active_clicker_calibration()
    calibration===nothing && return CLICKER_UNCALIBRATED_TDOA(data,peaks;kwargs...)
    return calibration["method"]=="minmax2" ? _clicker_minmax_tdoas(data,peaks;kwargs...) : _clicker_method_tdoas(data,peaks;kwargs...)
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

# c, default_c and default_window are left alone: the calibrated tdoa2dir/angle2px read sound speed
# and geometry from the JSON, and most setters never reset those globals after a switch.
function _activate_clicker_calibration!(device, calibration)
    set_device__ophk_acoustic_D3_clicker()
    global ref_channel=3
    global imsize=(Int(calibration["height"]),Int(calibration["width"]))
    camera=calibration["camera"]
    global fov_angle=[rad2deg(2atan(imsize[2]/(2camera["fx"]))),
                      rad2deg(2atan(imsize[1]/(2camera["fy"]))),missing]
    if default_getTDOA_func!==_clicker_tdoas
        global CLICKER_UNCALIBRATED_TDOA=default_getTDOA_func
    end
    global default_getTDOA_func=_clicker_tdoas
    global FILE_device_ID="ophk_acoustic_D3_clicker_$(device)"
    global ACTIVE_CLICKER_CALIBRATION=calibration
    global ACTIVE_CLICKER_TOKEN=(rx_vect,fov_angle)
    return calibration
end

function _set_clicker_calibration!(device,calibration_path)
    _ensure_clicker_module()
    calibration=_clicker_call(:load_calibration,calibration_path;device=device)  # validate before touching any global
    _activate_clicker_calibration!(device,calibration)
    return nothing
end

"""
    set_device__ophk_acoustic_D3_clicker_A(; calibration_path)

D3 clicker settings plus the fitted device-A calibration, until the next `set_device__*` call switches it off.
Only `tdoa2dir`, `angle2px` and the default TDOA estimator are calibrated; see docs/clicker_calibration.md.
"""
set_device__ophk_acoustic_D3_clicker_A(;calibration_path=_clicker_default_path("A"))=_set_clicker_calibration!("A",calibration_path)

"""Same as `set_device__ophk_acoustic_D3_clicker_A`, for device B."""
set_device__ophk_acoustic_D3_clicker_B(;calibration_path=_clicker_default_path("B"))=_set_clicker_calibration!("B",calibration_path)

"""
    with_clicker_calibration(f, device; calibration_path=_clicker_default_path(device))

Run `f()` with clicker `device` ("A"/"B" or :A/:B) active, then restore every global the activation
changed. Returns `f()`'s value. Blocks nest and are exception safe; a bad path changes nothing.
Run calibrated work inside the block and leave no tasks running after it. A `set_device__*` call
inside the block switches the calibration off for the rest of it.
"""
function with_clicker_calibration(f,device::Union{AbstractString,Symbol};calibration_path=_clicker_default_path(device))
    device=string(device)  # accept :A/:B as used by run_full_clicker_calibration
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
