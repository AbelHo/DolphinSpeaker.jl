include("video.jl")

v1 = "/Users/abel/Documents/data/concretecho/2024-01-24/2024-01-22_10.38.31_T007_topview.mkv"
v2 = "/Users/abel/Documents/data/concretecho/2024-01-24/2024-01-22_10.38.31_T007_uw1.mkv"
au = "/Users/abel/Documents/data/concretecho/2024-01-24/2024-01-22_10.38.31_T007.ogg"

using Plots
plotlyjs()

interval = [-.2 .5]

for fname in readdir("/Users/abel/Documents/data/concretecho/2024-01-24"; join=true)
    if fname[end-3:end] == ".ogg"
        fname_split = splitext(fname)
        thisfiletype = fname_split[2]
        fname_split = fname_split[1]
    end
end
            
findVidAudioBlip(v1; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5)
findVidAudioBlip(v2; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5)
findAudioBlip(au; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5)


include("audio.jl")
include("dsp.jl")
include("plotting.jl")
plotlyjs()
aufname = "/Users/abel/Documents/data/calf/coop/20231120/20231120_14.20.14_log.flac"

aufname = "/Users/abel/Documents/data/calf/sample/20211214_12.12.51_282s.flac"
aufname = "/Users/abel/Documents/data/calf/sample/click_20240124_61s.flac"
aufname = "/Users/abel/Documents/data/calf/sample/20231120_14.20.14_159s_1.flac"
data, fs, _,_, timestamp = readAudio(aufname)

data_filt = filter_simple(data, [1000 Inf]; fs=fs)
plot_norm(data_filt; norm_func=x->maximum(abs.(x),dims=1))
plot_norm!(data_filt |> tkeo; norm_func=x->maximum(abs.(x),dims=1))

plot_fft(truncate_fft_end(data), fs; type=:log)

plotlyjs()
plot(data)
data_filt = filter_simple(data[:,1:end-1], [1000 Inf]; fs=fs)
Plots.plot( filter_simple(abs.(hilbert(data_filt)), [0 10_000]; fs=fs) )
Plots.plot!( abs.(hilbert(data_filt)) )


data_filt = filter_simple(data[:,1:end-1], [1000 Inf]; fs=fs)
data_filt = filter_simple(abs.(hilbert(data_filt)), [0 10_000]; fs=fs)
data_filt = [data_filt @view(data[:,end])] ./ maximum(@view(data_filt[:,1:3]))
data = nothing; GC.gc()
wavwrite(data_filt, "/Users/abel/Documents/data/calf/coop/20231120/20231120_14.20.14_log.wav"; Fs=500_000)

wavwrite([data_filt @view(data[:,end])] ./ maximum(@view(data_filt[:,1:3])), 
    "/Users/abel/Documents/data/calf/coop/20231120/20231120_14.20.14_log.wav"; Fs=500_000)