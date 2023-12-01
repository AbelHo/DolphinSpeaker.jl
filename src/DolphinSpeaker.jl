module DolphinSpeaker
    include("audio.jl")
    export mat2flac, mat2wav

    include("synchronization.jl")
    export findBlip_bothVidAudio, flac2signal

    include("run_example.jl")
    export process_one_set

end
