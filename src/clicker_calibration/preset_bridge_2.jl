# Included by the last line of preset_bridge.jl (so by config.jl, many times): definitions only, guarded
# where state is involved. v2 profiles (clicker-physical-calibration/2.0): estimator-agnostic, explicit
# device, runtime sound speed, scoped activation that keeps the detection preset's TDOA estimator.
# Nothing here edits the v1 functions in preset_bridge.jl. See docs/clicker_calibration_v2_contract.md.
isdefined(@__MODULE__, :CLICKER_SITE_SOUND_SPEED_M_S) || (CLICKER_SITE_SOUND_SPEED_M_S = 1530.0)

# A function, not a const: this file is included many times.
default_v2_path(device)=joinpath(@__DIR__,"calibrations","clicker_20260209_$(device)_v2.json")

_v2_estimator_name(f)=f isa Function ? string(nameof(f)) : string(typeof(f))

# Everything the v2 activation installs, in one place. Called only inside with_clicker_calibration_2.
function _activate_clicker_calibration_2!(device,calibration,c_m_s,detection,tdoa_estimator)
    detection===:dolphin ? set_device__ophk_acoustic_D3() : set_device__ophk_acoustic_D3_clicker()
    global ref_channel=3
    global imsize=(Int(calibration["height"]),Int(calibration["width"]))
    camera=calibration["camera"]
    global fov_angle=[rad2deg(2atan(imsize[2]/(2camera["fx"]))),
                      rad2deg(2atan(imsize[1]/(2camera["fy"]))),missing]
    if tdoa_estimator!==nothing
        global default_getTDOA_func=tdoa_estimator
    elseif default_getTDOA_func===_clicker_tdoas   # left by a v1 activation: v2 never uses the v1 estimator
        CLICKER_UNCALIBRATED_TDOA===nothing && error("v1 clicker estimator is active but its uncalibrated estimator is unknown")
        global default_getTDOA_func=CLICKER_UNCALIBRATED_TDOA
    end
    global FILE_device_ID="ophk_acoustic_D3_$(device)2"
    global c=c_m_s
    global default_c=c_m_s
    global ACTIVE_CLICKER_CALIBRATION=calibration
    global ACTIVE_CLICKER_TOKEN=(rx_vect,fov_angle)
    return calibration
end

"""
    with_clicker_calibration_2(f, device; c_m_s=CLICKER_SITE_SOUND_SPEED_M_S, detection=:dolphin,
                               calibration_path=default_v2_path(device), tdoa_estimator=nothing)

Run `f()` with the v2 profile of `device` ("A"/"B"/:A/:B) active, then restore every global the activation
changed. `detection=:dolphin` applies `set_device__ophk_acoustic_D3()`, `:clicker` the clicker preset.
`c_m_s` (1400..1600) replaces the profile's sound speed and also sets the globals `c` and `default_c`.
The TDOA estimator is `tdoa_estimator` when given; otherwise the preset's own, except that the v1
`_clicker_tdoas` (left by a v1 setter) is swapped for the uncalibrated one it wrapped. Exception safe;
invalid device, detection, `c_m_s` or profile fail before any global changes.
"""
function with_clicker_calibration_2(f,device::Union{AbstractString,Symbol};c_m_s=CLICKER_SITE_SOUND_SPEED_M_S,
        detection=:dolphin,calibration_path=default_v2_path(device),tdoa_estimator=nothing)
    device=string(device)
    device in ("A","B") || throw(ArgumentError("device must be \"A\" or \"B\""))
    detection in (:dolphin,:clicker) || throw(ArgumentError("detection must be :dolphin or :clicker"))
    c_m_s isa Real && isfinite(c_m_s) && 1400<=c_m_s<=1600 || throw(ArgumentError("c_m_s must be within 1400..1600 m/s"))
    _ensure_clicker_module()
    # ClickerCalibration may have just been loaded: everything that touches a CalibrationV2 runs in the latest world.
    return Base.invokelatest(_with_clicker_calibration_2,f,device,Float64(c_m_s),detection,calibration_path,tdoa_estimator)
end

function _with_clicker_calibration_2(f,device,c_m_s,detection,calibration_path,tdoa_estimator)
    profile=_clicker_call(:load_calibration_2,calibration_path;device=device)   # validate before touching any global
    calibration=_clicker_call(:_v2_with_sound_speed,profile,c_m_s)
    before=_clicker_global_snapshot()
    changed=Symbol[]
    try
        try
            _activate_clicker_calibration_2!(device,calibration,c_m_s,detection,tdoa_estimator)
        finally
            append!(changed,_clicker_changed_names(before))
        end
        return f()
    finally
        for n in changed   # names first defined during activation are not in `before`; they are left as they are
            setglobal!(@__MODULE__,n,before[n])
        end
    end
end

# True only while this device's v2 profile, the runtime c, and the preset's estimator are all in place.
function _v2_active_status(device,c_m_s)
    calibration=active_clicker_calibration()
    calibration===nothing && return (ok=false,why="no active clicker calibration")
    get(calibration,"schema_version","")=="clicker-physical-calibration/2.0" || return (ok=false,why="active calibration is not schema 2.0")
    string(calibration["device"])==string(device) || return (ok=false,why="active calibration is for device $(calibration["device"])")
    calibration["sound_speed_m_s"]==Float64(c_m_s) && c==Float64(c_m_s) && default_c==Float64(c_m_s) || return (ok=false,why="sound speed differs from the requested $(c_m_s) m/s")
    default_getTDOA_func===_clicker_tdoas && return (ok=false,why="the v1 clicker estimator is installed")
    return (ok=true,why="")
end

function _v2_assert_active(device,c_m_s)
    status=_v2_active_status(device,c_m_s)
    status.ok || error("v2 clicker calibration is not active: $(status.why)")
    return active_clicker_calibration()
end

function _v2_video_geometry(path)
    info=Base.invokelatest(get_media_info,path)
    stream=first(filter(s->get(s,"codec_type","")=="video",info["streams"]))
    tags=get(get(info,"format",Dict()),"tags",Dict())
    return (width=Int(stream["width"]),height=Int(stream["height"]),creation_time=get(tags,"creation_time",nothing))
end

"""
    run_contiguous_folders_calibrated(folder, device; res_dir, c_m_s=CLICKER_SITE_SOUND_SPEED_M_S,
        calibration_path=default_v2_path(device), detection=:dolphin, tdoa_estimator=nothing, sync=:affine, kwargs...)

Run `run_contiguous_folders_2` (exact-time frames, drift-aware sync, NaN-safe overlay CSV) on `folder` inside
`with_clicker_calibration_2`. Before any heavy work it refuses a video whose resolution or GPMF lens mode
disagrees with the profile (`python3 -I tools/gpmf_lens_info.py`: VFOV, ZFOV within 0.01°, CASN of every chapter; unavailable → `@warn`) and warns outside the profile's
validity window. It asserts the calibration is active before and after the run and writes
`calibration_used.json` into the results folder (`res_dir/<folder name>`). Returns what
`run_contiguous_folders_2` returns with `flag_return=true`; remaining `kwargs` go to it.
"""
function run_contiguous_folders_calibrated(folder::AbstractString,device::Union{AbstractString,Symbol};res_dir,
        c_m_s=CLICKER_SITE_SOUND_SPEED_M_S,calibration_path=default_v2_path(device),detection=:dolphin,
        tdoa_estimator=nothing,sync=:affine,kwargs...)
    device=string(device)
    isempty(res_dir) && throw(ArgumentError("res_dir is required so calibration_used.json has a home"))
    _ensure_clicker_module()
    return Base.invokelatest(_run_contiguous_folders_calibrated,folder,device,res_dir,c_m_s,calibration_path,detection,tdoa_estimator,sync,kwargs)
end

function _run_contiguous_folders_calibrated(folder,device,res_dir,c_m_s,calibration_path,detection,tdoa_estimator,sync,kwargs)
    profile=_clicker_call(:load_calibration_2,calibration_path;device=device)   # fail before any heavy work
    isdefined(@__MODULE__,:run_contiguous_folders_2) || Base.include(@__MODULE__,joinpath(@__DIR__,"..","run_example_2.jl"))
    vidlist,audlist=Base.invokelatest(split_vid_au,folder;readdir_func=readdir_all)
    isempty(vidlist) && error("No video files found in $folder")

    geometries=[_v2_video_geometry(v) for v in vidlist]
    lens_infos=_clicker_call(:_v2_gpmf_lens_infos,vidlist)   # python3 -I tools/gpmf_lens_info.py, all chapters in one call
    date=_clicker_call(:_v2_recording_date,vidlist[1];creation_time=geometries[1].creation_time)
    for (v,g) in zip(vidlist,geometries)   # every chapter must match the profile; lens mode and date come from the first
        (g.width,g.height)==(Int(profile["width"]),Int(profile["height"])) ||
            error("Video $(basename(v)) is $(g.width)x$(g.height); the calibration profile is $(profile["width"])x$(profile["height"]) (device $device); refusing to run")
    end
    checks=_clicker_call(:_v2_validate_video,profile;width=geometries[1].width,height=geometries[1].height,lens_infos=lens_infos,recording_date=date)

    results_dir=joinpath(res_dir,basename(folder))
    used_path=joinpath(results_dir,"calibration_used.json")
    run_kwargs=[kv for kv in kwargs if first(kv)!==:flag_return]
    return with_clicker_calibration_2(device;c_m_s=c_m_s,detection=detection,calibration_path=calibration_path,tdoa_estimator=tdoa_estimator) do
        _v2_assert_active(device,c_m_s)
        record=Dict{String,Any}("status"=>"started","device"=>device,"folder"=>abspath(folder),
            "calibration_path"=>abspath(calibration_path),"profile_sha256"=>_clicker_call(:sha256_file,calibration_path),
            "c_m_s"=>Float64(c_m_s),"detection_preset"=>string(detection),
            "tdoa_estimator"=>_v2_estimator_name(default_getTDOA_func),"sync"=>string(sync),
            "lens_mode"=>checks.lens.lens_mode,"zfov_deg"=>checks.lens.zfov_deg,"casn"=>checks.lens.casn,
            "profile_lens_mode"=>checks.lens.profile_lens_mode,"profile_zfov_deg"=>checks.lens.profile_zfov_deg,"profile_casn"=>checks.lens.profile_casn,
            "gpmf_lens_info"=>lens_infos,
            "video_resolution"=>[geometries[1].width,geometries[1].height],"recording_date"=>date===nothing ? nothing : string(date),
            "valid_from"=>get(profile,"valid_from",nothing),"valid_until"=>get(profile,"valid_until",nothing),
            "in_validity_window"=>checks.in_validity_window,"warnings"=>checks.warnings,
            "n_videos"=>length(vidlist),"n_audio"=>length(audlist),
            "git_commit"=>_clicker_call(:_v2_git_commit),"started_utc"=>string(_clicker_call(:_v2_now_utc)))
        _clicker_call(:write_json,used_path,record)
        out=nothing
        try
            out=Base.invokelatest(run_contiguous_folders_2,folder;res_dir=res_dir,sync=sync,flag_return=true,run_kwargs...)
            after=_v2_active_status(device,c_m_s)
            record["active_after_run"]=after.ok
            after.ok || error("v2 clicker calibration was switched off during the run: $(after.why)")
            record["status"]="completed"
        catch err
            record["status"]=get(record,"active_after_run",true)===false ? "calibration_switched_off" : "failed"
            record["error"]=first(sprint(showerror,err),2000)
            rethrow()
        finally
            record["finished_utc"]=string(_clicker_call(:_v2_now_utc))
            _clicker_call(:write_json,used_path,record)
            out!==nothing && hasproperty(out,:res_dir2) && out.res_dir2!=results_dir && _clicker_call(:write_json,joinpath(out.res_dir2,"calibration_used.json"),record)
        end
        return out
    end
end
