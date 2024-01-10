module DolphinSpeaker
@info "version v2024-01-11T01:00 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-"
    include("utils.jl")
    export showall
    
    include("audio.jl")
    export mat2flac, mat2wav, mat2flac_check, bin2flac, bin2flac_check
    export FILE_device_ID, FILE_location_ID, FILE_gain_setting, FILE_device_location

    include("synchronization.jl")
    export findBlip_bothVidAudio, flac2signal

    include("run_example.jl")
    export process_one_set, process_folder, process_dir

end
