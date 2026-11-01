module DolphinSpeaker
    version = "2025-10-08T10:00"
    @info "version v$version =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-"
    include("utils.jl")
    export showall, skiphiddenfiles, process_files
    
    include("audio.jl")
    export mat2flac, mat2wav, mat2flac_check, bin2flac, bin2flac_check

    include("synchronization.jl")
    export findBlip_bothVidAudio, flac2signal, multisync

    include("run_example.jl")
    export process_one_set, process_folder, process_dir, run_contiguous_folders

    # from dsp.jl
    export extrema_in_file

    include("beampattern.jl")
    export stack_audio_videos, run_func_fileauto

    # from config.jl
    export set_device__ophk_acoustic_D3, set_device__calf_hk, set_device__hk_clicker, set_device__ophk_acoustic_D3_clicker
    export with_clicker_calibration

    # Thin wrappers: ClickerCalibration.jl loads lazily via preset_bridge.jl's _ensure_clicker_module.
    """Run the independent clicker calibration workflow; see docs/clicker_calibration.md."""
    run_clicker_calibration(args...;kwargs...) = (_ensure_clicker_module(); _clicker_call(:run_clicker_calibration,args...;kwargs...))
    """Run the one-call clicker calibration order (prepare through delivery); see docs/clicker_calibration.md."""
    run_full_clicker_calibration(args...;kwargs...) = (_ensure_clicker_module(); _clicker_call(:run_full_clicker_calibration,args...;kwargs...))
    export run_clicker_calibration, run_full_clicker_calibration
end
