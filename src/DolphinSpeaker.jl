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
    export process_one_set, process_folder, process_dir

    # from dsp.jl
    export extrema_in_file

    include("beampattern.jl")
    export stack_audio_videos, run_func_fileauto

end
