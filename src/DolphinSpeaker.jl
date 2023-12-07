module DolphinSpeaker
@info "version v2023-12-03T19:01:00 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-"
    include("utils.jl")
    export showall
    
    include("audio.jl")
    export mat2flac, mat2wav, mat2flac_check

    include("synchronization.jl")
    export findBlip_bothVidAudio, flac2signal

    include("run_example.jl")
    export process_one_set

end
