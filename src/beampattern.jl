using Plots
using Interpolations
using ProgressMeter
include("dsp.jl")
include("detector.jl")
include("audio.jl")
include("video.jl")

# rx_vect = randn(3,16)

# function beampattern(θ, φ, f, c, a, d)
#     k = 2π*f/c
#     r = sqrt(d^2 .+ a^2 .- 2*a*d.*cos.(θ))
#     return 20*log10.(abs.(1 .+ (1im/k)*r.*exp.(1im*k*r)./(4π*r)))
# end



# rx_dist = 0.125
# x = range(-rx_dist*1.5, rx_dist*1.5, length=4)
# y = range(-rx_dist*1.5, rx_dist*1.5, length=4)
# z = @. cos(x.*3) + sin(y')

# using Interpolations, Plots
# Data
# x = range(-2, 3, length=20)
# y = range(3, 4, length=10)
# z = @. cos(x) + sin(y')
# Interpolation object (caches coefficients and such)
function interp_2D(x,y,z; interp_len=1000)
    # itp = LinearInterpolation((x, y), z) # CubicSplineInterpolation
    itp = CubicSplineInterpolation((x, y), z) 
    # Fine grid
    x2 = range(extrema(x)..., length=interp_len)
    y2 = range(extrema(y)..., length=interp_len)
    # Interpolate
    z2 = [itp(x,y) for y in y2, x in x2]
    # @info extrema(z2)
    maxi = maximum(z2)
    # Plot
    # p = heatmap(x2, y2, z2, clim=(-2,2), title="Interpolated heatmap", colorbar=nothing)

    p = heatmap(x2, y2, z2, title="Interpolated heatmap", colorbar=nothing)
    scatter!(p, [x for _ in y for x in x], [y for y in y for _ in x], zcolor=z[:]; lab="original data",
        xlabel="X", ylabel="Y", aspect_ratio=:equal, labels=nothing,
        xlims=extrema(x2), ylims=extrema(y2))

    contour!(p, x2, y2, z2; color=:white, legend=false, clabels=true, levels=[0, -3, -6] .+ maxi)
    # scatter!(p, [x for _ in y for x in x], [y for y in y for _ in x], zcolor=z[:]; lab="original data", clim=(-2,2),
    #     xlabel="X", ylabel="Y", aspect_ratio=:equal, labels=nothing)
end


function plot_beam(pos_x,pos_y,z, x, y)
    z = z .|> Float64
    sorted1 = sortperm(pos_x)
    z2 = z[sorted1]
    pos_y2 = pos_y[sorted1]
    pos_x2 = pos_x[sorted1]

    sorted2 = sortperm(pos_y2)
    pos_x3 = pos_x2[sorted2]
    pos_y3 = pos_y2[sorted2]
    z3 = z2[sorted2]

    interp_2D(x,y,reshape(z3,4,4))
    # title!("$i") #|> display
end

function plot_beams(aufname_data_fs, res_dir, window=-300:300; rx_vect2=rx_vect)
    aufname, data, fs = aufname_data_fs

    rx_vect = rx_vect2
    pos_x = rx_vect[1,:]
    pos_y = rx_vect[2,:]

    res = detect_impulseNtonal((aufname,data,fs,nothing), res_dir; return_datafilt=true, rx_vect=rx_vect)

    # @info (signal(res.res_impulse.data_filt,fs)), map(x-> x.+ window, res.res_impulse.pind_good)#; func=x->( mapslices(extrema,x;dims=1) ) )
    a = funcOnWindows(res.res_impulse.data_filt, map(x-> x.+ window, res.res_impulse.pind_good); func=x->( mapslices(extrema,x;dims=1) ) )
    z = map( x-> 20*log10(x[2]-x[1]), a)

    println("Creating video frame for each click..." )
    progress = Progress(size(a, 1), 1);
    anim = @animate for i in 1:size(a,1)
        plot_beam(pos_x,pos_y,z[i,:], x, y)
        title!("$i") #|> display
        next!(progress)
    end

    outgifname = joinpath(res_dir, splitext(basename(aufname))[1] *"_beam.gif")
    gif(anim,  outgifname; fps = 10)
    outvidname = splitext(outgifname)[1]*".mp4"
    @ffmpeg_env run(`$ffmpeg -i $outgifname -movflags faststart -pix_fmt yuv420p $outvidname -hide_banner -y`)

end

function plot_beams(aufname::String, args...; kwargs...)
    data, fs = readAudio(aufname)
    plot_beams((aufname,data,fs,nothing), args...; kwargs...)
end




# rx_dist = 0.125
# x = range(-rx_dist*1.5, rx_dist*1.5, length=4)
# y = range(-rx_dist*1.5, rx_dist*1.5, length=4)


# aufname = "/Users/abel/Documents/data/concretecho/2023-12-04/acoustic/2023-12-04_12.35.06.ogg"
# aufname = "/Users/abel/Documents/data/concretecho/2023-12-04/acoustic/2023-12-04_11.19.24.ogg"
# aufname = "/Users/abel/Documents/data/concretecho/2023-12-04/acoustic/2023-12-04_12.33.54.ogg"
# aufname = "/Users/abel/Documents/data/concretecho/2023-12-04/acoustic/2023-12-04_11.56.27.ogg"
# res_dir = "/Users/abel/Documents/data_res/concretecho/2023-12-04_2"

# aufname = "/Users/abel/Documents/data/concretecho/Ella/acoustic/2024/01/2024-01-17_12.33.17/2024-01-17_12.36.19_T004_P2_C1_acoustic.ogg"
# res_dir = "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33"
# dist_impulsive = 1500


# rx_dist = 0.125
# pos_x = repeat( range(-rx_dist*1.5, rx_dist*1.5, length=4), 4)
# pos_y = repeat( range(rx_dist*1.5, -rx_dist*1.5, length=4), inner=4)
# rx_vect = [pos_x'; pos_y'; zeros(1,length(pos_x))]
# # z = [ 6 1 1 0; 4 2 2 0; 0 3 3 0; 0 1 1 8]

# data, fs = readAudio(aufname)
# res = detect_impulseNtonal((aufname,data,fs,nothing), res_dir; return_datafilt=true)

# a = funcOnWindows(signal(res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulse.pind_good); func=x->( mapslices(extrema,x;dims=1) ) )
# # a = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulse.pind_good); func=x->( mapslices(extrema,x;dims=1) ) )
# # a = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulsetrain.pind_good); func=x->( mapslices(extrema,x;dims=1) ) )
# a2 = map( x-> 20*log10(x[2]-x[1]), a)

# println("Creating video frame for each click..." )
# progress = Progress(size(a, 1), 1);
# anim = @animate for i in 1:size(a,1)
#     z = a2[i,:] .|> Float64
#     sorted1 = sortperm(pos_x)
#     z2 = z[sorted1]
#     pos_y2 = pos_y[sorted1]
#     pos_x2 = pos_x[sorted1]

#     sorted2 = sortperm(pos_y2)
#     pos_x3 = pos_x2[sorted2]
#     pos_y3 = pos_y2[sorted2]
#     z3 = z2[sorted2]

#     interp_2D(x,y,reshape(z3,4,4))
#     title!("$i") #|> display
#     next!(progress)
# end

# outgifname = joinpath(res_dir, splitext(basename(aufname))[1] *"_beam-time.gif")
# gif(anim,  outgifname; fps = 10)
# outvidname = splitext(outgifname)[1]*".mp4"
# @ffmpeg_env run(`$ffmpeg -i $outgifname -movflags faststart -pix_fmt yuv420p $outvidname -hide_banner -y`)





############################################################################################################
# include("video.jl")
# rx_dist = 0.125
# x = range(-rx_dist*1.5, rx_dist*1.5, length=4)
# y = range(-rx_dist*1.5, rx_dist*1.5, length=4)

# v1 = "/Users/abel/Documents/data/concretecho/Ella/topview/2024/01/2024-01-17_12.33.17/2024-01-17_12.36.19_T004_P2_C1_topview.mkv"
# v2 = "/Users/abel/Documents/data/concretecho/Ella/uw1/2024/01/2024-01-17_12.33.17/2024-01-17_12.36.19_T004_P2_C1_uw1.mkv"

# aufname = "/Users/abel/Documents/data/concretecho/2024-03-06/2024-03-06_11.30.10_T010_acoustic.ogg"
# v1 = "/Users/abel/Documents/data/concretecho/2024-03-06/2024-03-06_11.30.10_T010_topview.mkv"
# v2 = "/Users/abel/Documents/data/concretecho/2024-03-06/2024-03-06_11.30.10_T010_uw1.mkv"

# res_dir = joinpath("/Users/abel/Documents/data_res/concretecho", basename(aufname) |> x -> x[1:findlast("_",x)[1]-1])

function stack_audio_videos(aufname, v1, v2, res_dir; 
    rx_dist = 0.125,
    x = range(-rx_dist*1.5, rx_dist*1.5, length=4),
    y = range(-rx_dist*1.5, rx_dist*1.5, length=4),

    pos_x = repeat( range(-rx_dist*1.5, rx_dist*1.5, length=4), 4),
    pos_y = repeat( range(rx_dist*1.5, -rx_dist*1.5, length=4), inner=4),
    kwargs...)

    @info (aufname, v1, v2, res_dir)

    isdir(res_dir) || mkpath(res_dir)

    data, fs = readAudio(aufname)
    res = detect_impulseNtonal((aufname,data,fs,nothing), res_dir; return_datafilt=true)
    # res_impulse = res.res_impulse
    # data_filt = filter_simple(data, [1000 Inf]; fs=fs)

    fps = get_fps(v1) #25
    a = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulsetrain.pind_good); func=x->( mapslices(extrema,x;dims=1) ) )
    a2 = map( x-> 20*log10(x[2]-x[1]), a)


    println("Creating video frame for each click..." )
    progress = Progress(size(a, 1), 1);
    anim = @animate for t in 0:1/fps:(size(data,1) / fs)
        clicks_in_timeframe = findall( x-> t <= x <= t+1/fps, res.res_impulsetrain.pind_good_inS)
        if isempty(clicks_in_timeframe) 
            interp_2D(x,y,zeros(4,4));
            title!(" ") 
            next!(progress)
        else
            i = clicks_in_timeframe[ maximum( a2[clicks_in_timeframe, :]; dims=2) |> argmax ]

            z = a2[i,:] .|> Float64
            sorted1 = sortperm(pos_x)
            z2 = z[sorted1]
            pos_y2 = pos_y[sorted1]
            pos_x2 = pos_x[sorted1]

            sorted2 = sortperm(pos_y2)
            pos_x3 = pos_x2[sorted2]
            pos_y3 = pos_y2[sorted2]
            z3 = z2[sorted2]

            interp_2D(x,y,reshape(z3,4,4));
            title!("$i") #|> display
            next!(progress)
        end
    end

    outgifname = joinpath(res_dir, splitext(basename(aufname))[1] *"_beam-time.gif")
    gif(anim,  outgifname; fps = fps)
    outvidname = splitext(outgifname)[1]*".mp4"
    @ffmpeg_env run(`$ffmpeg -i $outgifname -movflags faststart -pix_fmt yuv420p $outvidname -hide_banner -y`)

    signal_plot_fname = joinpath(res_dir, splitext(basename(aufname))[1] *"_time-ch1.mp4")
    plot_signal2vid(data[:,1],fs,signal_plot_fname; fps=fps, yaxis=false, size=(900,150))

    vidoutname = joinpath(res_dir, v1[1:findlast('_', v1)] |> basename)
    vid_combine_fname =  vidoutname * "vid-combined.mp4"
    combine_2v1a(v1,v2,aufname, vid_combine_fname)

    vid_combine_beam_fname  = splitext(vid_combine_fname)[1]  * "_beam.mp4"
    @ffmpeg_env run(`ffmpeg -i $vid_combine_fname -i $outvidname -filter_complex "[0:v]scale=-1:400[v0];[1:v][v0]hstack" $vid_combine_beam_fname`)

    vid_fullcombined_fname = vidoutname * "beam-vid-combined_timeplot.mp4"
    @ffmpeg_env run(`$ffmpeg -i $vid_combine_beam_fname -i $signal_plot_fname -filter_complex "[1:v]scale=900:120[v1];[0:v][v1]vstack" -metadata comment="$aufname,$v1,$v2" $vid_fullcombined_fname -hide_banner`)

end

# ffmpeg -i "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_vid-combined.mp4" -i "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_acoustic_beam-time.mp4" -filter_complex "[0:v]scale=-1:400[v0];[1:v][v0]hstack" "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_beam-vid-combined.mp4"
# ffmpeg -i "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_beam-vid-combined.mp4" -i "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_acoustic_time-ch1.mp4" -filter_complex "[1:v]scale=900:120[v1];[0:v][v1]vstack" "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_beam-vid-combined_timeplot7.1.mp4"
# ffmpeg -i "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_beam-vid-combined.mp4" -i "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_acoustic_time-ch1_2.mp4" -filter_complex "[1:v]scale=900:120[v1];[0:v][v1]vstack" "/Users/abel/Documents/data_res/concretecho/2024-01-17_12.33/2024-01-17_12.36.19_T004_P2_C1_beam-vid-combined_timeplot7.2.mp4"






# interval = [-.2 .5]
# include("synchronization.jl")
# findVidAudioBlip(v1; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5)
# findVidAudioBlip(v2; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5)
# findAudioBlip(aufname; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5)

# aufname = "/Users/abel/Documents/data/concretecho/2024-03-06/2024-03-06_11.30.10_T010_acoustic.ogg"
# res = detect_impulseNtonal((aufname,data,fs,nothing), res_dir; return_datafilt=true)
# res_impulse = res.res_impulse
# data_filt = filter_simple(data, [1000 Inf]; fs=fs)