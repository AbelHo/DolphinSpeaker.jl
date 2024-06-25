using Plots
using Interpolations
using ProgressMeter
include("dsp.jl")
include("detector.jl")
include("audio.jl")
include("video.jl")

# set_device__rwsnus()

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

function interp2D(x,y,z; interp_len=1000)
    # itp = LinearInterpolation((x, y), z) # CubicSplineInterpolation
    itp = CubicSplineInterpolation((x, y), z) 
    # Fine grid
    x2 = range(extrema(x)..., length=interp_len)
    y2 = range(extrema(y)..., length=interp_len)
    # Interpolate
    z2 = [itp(x,y) for y in y2, x in x2]

    return z2
end
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
    @ffmpeg_env run(`$ffmpeg -i $outgifname -pix_fmt yuv420p $outvidname -hide_banner -y`)
    # @ffmpeg_env run(`$ffmpeg -i $outgifname -movflags faststart -pix_fmt yuv420p $outvidname -hide_banner -y`)

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
# run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/03";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/03"); func=stack_audio_videos, skipdone=true); 

"""
    stack_audio_videos(aufname, v1, v2, res_dir; rx_dist = 0.125, x = range(-rx_dist*1.5, rx_dist*1.5, length=4), y = range(-rx_dist*1.5, rx_dist*1.5, length=4), pos_x = repeat( range(-rx_dist*1.5, rx_dist*1.5, length=4), 4), pos_y = repeat( range(rx_dist*1.5, -rx_dist*1.5, length=4), inner=4), skipdone = false, kwargs...)

Stacks audio and video files, applies a beamforming algorithm, and generates a combined video output.

# Arguments
- `aufname`: Path to the audio file.
- `v1`: Path to the first video file.
- `v2`: Path to the second video file.
- `res_dir`: Path to the directory where the result will be saved.

# Keyword Arguments
- `rx_dist`: Distance between receivers. Default is 0.125.
- `x`: Range of x-coordinates for the receivers. Default is a range from -1.5*rx_dist to 1.5*rx_dist with 4 points.
- `y`: Range of y-coordinates for the receivers. Default is a range from -1.5*rx_dist to 1.5*rx_dist with 4 points.
- `pos_x`: Positions of the receivers in the x-axis. Default is a repeated range from -1.5*rx_dist to 1.5*rx_dist with 4 points, repeated 4 times.
- `pos_y`: Positions of the receivers in the y-axis. Default is a repeated range from 1.5*rx_dist to -1.5*rx_dist with 4 points, repeated 4 times internally.
- `skipdone`: If true, the function will skip processing if the output file already exists. Default is false.

# Returns
- Nothing. The function saves the output video in the `res_dir` directory.

# Example
```julia
stack_audio_videos("audio.wav", "video1.mp4", "video2.mp4", "results")
run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/03";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/03"); func=stack_audio_videos, skipdone=true);
run_func_fileauto("/Users/abel/Documents/data/concretecho/data/Shakeela/2024/05/2024-05-07_15.29.13", "/Users/abel/Documents/data_res/concretecho/temp/2024-05-07_15.29.13"; func=stack_audio_videos)
run_func_fileauto.(readdir(infol;join=true)|>skiphiddenfiles, Ref(outfol); func=stack_audio_videos, skipdone=true); 
```
"""
function stack_audio_videos(aufname, v1, v2, res_dir; 
    rx_dist = 0.125,
    x = range(-rx_dist*1.5, rx_dist*1.5, length=4),
    y = range(-rx_dist*1.5, rx_dist*1.5, length=4),

    pos_x = repeat( range(-rx_dist*1.5, rx_dist*1.5, length=4), 4),
    pos_y = repeat( range(rx_dist*1.5, -rx_dist*1.5, length=4), inner=4),
    skipdone = false,
    kwargs...)

    @info (aufname, v1, v2, res_dir)

    vidoutname = joinpath(res_dir, v1[1:findlast('_', v1)] |> basename)
    vid_combine_fname =  vidoutname * "vid-combined.mp4"
    vid_combine_beam_fname  = splitext(vid_combine_fname)[1]  * "_beam.mp4"
    vid_fullcombined_fname = vidoutname * "beam-vid-combined_timeplot.mp4"

    skipdone && isfile(vid_fullcombined_fname) && (@info "skip......"; return nothing)

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
    @ffmpeg_env run(`$ffmpeg -i $outgifname -pix_fmt yuv420p $outvidname -hide_banner -y`)

    signal_plot_fname = joinpath(res_dir, splitext(basename(aufname))[1] *"_time-ch1.mp4")
    plot_signal2vid(data[:,1],fs,signal_plot_fname; fps=fps, yaxis=false, size=(900,150))

    # vidoutname = joinpath(res_dir, v1[1:findlast('_', v1)] |> basename)
    # vid_combine_fname =  vidoutname * "vid-combined.mp4"
    combine_2v1a(v1,v2,aufname, vid_combine_fname)

    # vid_combine_beam_fname  = splitext(vid_combine_fname)[1]  * "_beam.mp4"
    @ffmpeg_env run(`ffmpeg -i $vid_combine_fname -i $outvidname -filter_complex "[0:v]scale=-1:400[v0];[1:v][v0]hstack" $vid_combine_beam_fname`)

    # vid_fullcombined_fname = vidoutname * "beam-vid-combined_timeplot.mp4"
    # @ffmpeg_env run(`$ffmpeg -i $vid_combine_beam_fname -i $signal_plot_fname -filter_complex "[1:v]scale=900:120[v1];[0:v][v1]vstack[v2];[v2]drawtext=text='%{n}': x=10: y=35: fontsize=24: fontcolor=black" -metadata comment="$aufname,$v1,$v2" $vid_fullcombined_fname -hide_banner`)
    @ffmpeg_env run(`$ffmpeg -i $vid_combine_beam_fname -i $signal_plot_fname -filter_complex "[1:v]scale=900:120[v1];[0:v][v1]vstack" -metadata comment="$aufname,$v1,$v2" $vid_fullcombined_fname -hide_banner`)
    # rm(vid_combine_beam_fname)@ffmpeg_env run(`$ffmpeg -i $vid_combine_beam_fname -i $signal_plot_fname -filter_complex "[1:v]scale=900:120[v1];[0:v][v1]vstack;[0:v]drawtext=text='%{n}': x=10: y=10: fontsize=24: fontcolor=white[v2];[v2][v1]vstack" -metadata comment="$aufname,$v1,$v2" $vid_fullcombined_fname -hide_banner`)
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
ff(x; band_pass=[0, Inf]) = mapslices(extrema, filter_simple(x, band_pass; fs=fs, mapslices2=mapslices, dims=1);dims=1)


# rx_dist = 0.125
# x = range(-rx_dist*1.5, rx_dist*1.5, length=4)
# y = range(-rx_dist*1.5, rx_dist*1.5, length=4)

# aufname = "/Users/abel/Documents/data/concretecho/Ella/acoustic/2024/01/2024-01-17_12.33.17/2024-01-17_12.36.19_T004_P2_C1_acoustic.ogg"
# # aufname = "/Users/abel/Documents/data/concretecho/2024-03-06/2024-03-06_11.30.10_T010_acoustic.ogg"
# res_dir = "/Users/abel/Documents/data_res/concretecho/temp"
# data, fs = readAudio(aufname)
# res = detect_impulseNtonal((aufname,data,fs,nothing), res_dir; return_datafilt=true)



# a_1 = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulsetrain.pind_good); func=x-> ff(x; band_pass=[1000, 70_000]) )
# a_2 = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulsetrain.pind_good); func=x-> ff(x; band_pass=[70_000, 120_000]) )
# a_3 = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulsetrain.pind_good); func=x-> ff(x; band_pass=[120_000,170_000]) )

# a2_1 = map( x-> 20*log10(x[2]-x[1]), a_1)
# a2_2 = map( x-> 20*log10(x[2]-x[1]), a_2)
# a2_3 = map( x-> 20*log10(x[2]-x[1]), a_3)

# anim = @animate for ind = 1:size(a_1,1)
#     z1 = interp2D(x,y,reshape(a2_1[ind,:],4,4))
#     z2 = interp2D(x,y,reshape(a2_2[ind,:],4,4))
#     z3 = interp2D(x,y,reshape(a2_3[ind,:],4,4))

#     img = RGB.(z1./maximum(z1),z2./maximum(z2),z3./maximum(z3))
#     plot(img)
# end
# gif(anim, joinpath(res_dir, "$(basename(aufname))_color_beam.gif"))
# @ffmpeg_env run(`$ffmpeg -i $(joinpath(res_dir, basename(aufname)*"_color_beam.gif")) -pix_fmt yuv420p $(joinpath(res_dir, "$(basename(aufname))_color_beam.mp4"))`)


# a = funcOnWindows(signal(res.res_impulse.data_filt,fs), map(x-> x.+ (-300:300), res.res_impulsetrain.pind_good); func=x->( mapslices(extrema,x;dims=1) ) )
# a2 = map( x-> 20*log10(x[2]-x[1]), a)
# a2_old = a2
# for 
# anim = @animate for ind = 1:size(a_1,1)
#     interp_2D(x,y,reshape(a2[ind,:],4,4))
# end
# gif(anim, joinpath(res_dir, "$(basename(aufname))_beam.gif"))
# @ffmpeg_env run(`$ffmpeg -i $(joinpath(res_dir, "$(basename(aufname))_beam.gif")) -pix_fmt yuv420p $(joinpath(res_dir, "$(basename(aufname))_beam.mp4"))`)

# GC.gc()



#############################################################################
# run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/03";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/03"); func=stack_audio_video, skipdone=true); 
# run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/04";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/04"); func=stack_audio_video, skipdone=true);
# run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/02";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/02"); func=stack_audio_video, skipdone=true); run_func_fileauto.(readdir("/Volumes/One Touch/data/concretecho/try8/Ella/2024/01";join=true)|>skiphiddenfiles, Ref("/Volumes/One Touch/res/concretecho/outvid/Ella/01"); func=stack_audio_video, skipdone=true)




# infol = "/Users/abel/Documents/data/concretecho/data/Shakeela/2024/05"
# outfol = "/Users/abel/Documents/data_res/concretecho/beam/Shakeela/05_frame"
# run_func_fileauto.(readdir(infol;join=true)|>skiphiddenfiles, Ref(outfol); func=stack_audio_videos, skipdone=true); 

# infol = "/Volumes/data/Concretecho/data/temp/try8/Shiye/2024/06"#/2024-06-20_13.21.12"
# outfol = "/Users/abel/Documents/data_res/concretecho/temp/stack2"
# run_func_fileauto.(readdir(infol;join=true)|>skiphiddenfiles, Ref(outfol); func=stack_audio_videos, skipdone=true); 

