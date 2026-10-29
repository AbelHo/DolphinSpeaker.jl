# Included by config.jl. Definitions only: calibrated parameters load on request.
if !isdefined(@__MODULE__, :ACTIVE_CLICKER_CALIBRATION)
    ACTIVE_CLICKER_CALIBRATION = nothing
    CLICKER_CALIBRATION_RESTORE = nothing
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

function clear_clicker_calibration!()
    if CLICKER_CALIBRATION_RESTORE !== nothing
        old=CLICKER_CALIBRATION_RESTORE
        global c=old.c
        old.default_c === nothing || (global default_c=old.default_c)
        old.tdoa === nothing || (global default_getTDOA_func=old.tdoa)
    end
    global ACTIVE_CLICKER_CALIBRATION=nothing
    global CLICKER_CALIBRATION_RESTORE=nothing
    return nothing
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

function _set_clicker_calibration(device, path)
    if !isdefined(@__MODULE__,:ClickerCalibration)
        Base.include(@__MODULE__,joinpath(@__DIR__,"ClickerCalibration.jl"))
    end
    # Validate before changing the existing active device.
    calibration=_clicker_call(:load_calibration,path;device=device)
    set_device__ophk_acoustic_D3_clicker()
    global CLICKER_CALIBRATION_RESTORE=(c=c,
        default_c=isdefined(@__MODULE__,:default_c) ? default_c : nothing,
        tdoa=isdefined(@__MODULE__,:default_getTDOA_func) ? default_getTDOA_func : nothing)
    global ACTIVE_CLICKER_CALIBRATION=calibration
    global c=calibration["sound_speed_m_s"]
    global default_c=c
    global ref_channel=3
    global imsize=(Int(calibration["height"]),Int(calibration["width"]))
    camera=calibration["camera"]
    global fov_angle=[rad2deg(2atan(imsize[2]/(2camera["fx"]))),
                      rad2deg(2atan(imsize[1]/(2camera["fy"]))),missing]
    global default_getTDOA_func=_clicker_minmax_tdoas
    global FILE_device_ID="ophk_acoustic_D3_clicker_$(device)"
    return calibration
end

function calibrated_tdoa2dir(tdoas,receivers,fs,calibration;return_residual=false)
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
