using VideoIO
using FFMPEG
include("media_info.jl")
include("synchronization.jl")
#  @time open_video_out(newvidname, img2, framerate=get_fps(vidfname), encoder_options=encoder_options) do writer
    
# resolution=(1080,720)
# fig = Figure(;resolution=result_resolution)

# v1="/Users/abel/Documents/data/concretecho/2023-12-04/cam_uw1/2023-12-04_12.33.03_uw1.mkv"
# v2="/Users/abel/Documents/data/concretecho/2023-12-04/cam_topview/2023-12-04_12.33.03_topview.mkv"
# au="/Users/abel/Documents/data/concretecho/2023-12-04/acoustic/2023-12-04_12.33.03.ogg"
function combine_2v1a(v1,v2,au,output_file; sync_type=:new)
    @info(v1,v2,au,output_file)
    isdir(output_file) && (output_file = joinpath(output_file, splitext(basename(au))[1]*"_norm.mp4" ))
    # output_file_norm = splitext(output_file)[1]*"_norm.mp4"

    if sync_type == :old
        v1_trigger = findVidAudioBlip(v1; argmax_len=0, plot_window_inS=nothing)
        v2_trigger = findVidAudioBlip(v2; argmax_len=0, plot_window_inS=nothing)
        au_trigger = findAudioBlip(au; argmax_len=0, plot_window_inS=nothing)
    else
        interval = [-.2 .5]
        v1_trigger = findVidAudioBlip(v1; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5, flag_savefig=dirname(output_file))
        v2_trigger = findVidAudioBlip(v2; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5, flag_savefig=dirname(output_file))
        au_trigger = findAudioBlip(au; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5, flag_savefig=dirname(output_file))
    end
    # mini = min(v1_trigger, v2_trigger, au_trigger)
    # @info mini
    v1_trigger_n = v1_trigger - au_trigger
    v2_trigger_n = v2_trigger - au_trigger
    # au_trigger -= mini
    @info "delays: " * string(v1_trigger) * " " * string(v2_trigger) * " " * string(au_trigger)

    # width = get_whatever(v1, 0, "v"; entries_custom="width") |> Int
    try
        println(`ffmpeg -itsoffset $v1_trigger_n -i "$v1" -itsoffset $v2_trigger_n -i "$v2" -i "$au" -filter_complex "[2:a]loudnorm[a];[a]asplit[a0][a1]; [a1]showwaves=s=1920x400:mode=cline:colors=red:rate=25,format=yuv420p[vau]; [0:v][1:v][vau]vstack=inputs=3[v]" -map "[v]" -map "[a0]" $output_file -hide_banner`)
        @ffmpeg_env run(`ffmpeg -itsoffset $v1_trigger_n -i "$v1" -itsoffset $v2_trigger_n -i "$v2" -i "$au" -filter_complex "[2:a]loudnorm[a];[a]asplit[a0][a1]; [a1]showwaves=s=1920x400:mode=cline:colors=red:rate=25,format=yuv420p[vau]; [0:v][1:v][vau]vstack=inputs=3[v]" -map "[v]" -map "[a0]" $output_file -hide_banner`)
    catch err
        @error "FAILED!:    $au"
    end
    # @ffmpeg_env run(`ffmpeg -i "$v1" -i "$v2" -i "$au" -filter_complex "[0:v][1:v]vstack=inputs=2[v];[2:a]loudnorm[a]" -map "[v]" -map "[a]" $output_file`)
    # @ffmpeg_env run(`ffmpeg -i $output_file `) #-af loudnorm=I=-16:LRA=11:TP=-1.5
    # ffmpeg -i "$v1" -i "$v2" -i "$au" -filter_complex "[0:v][1:v]vstack=inputs=2[v];[2:a]anull[a]" -map "[v]" -map "[a]" output.mp4
    return v1_trigger, v2_trigger, au_trigger
end

function combine_2v1a_auto(aufolder, outfolder; filetype=".ogg")
    mkpath(outfolder)
    dname = dirname(aufolder)
    for fname in readdir(aufolder)|>skiphiddenfiles
        fname_split = splitext(fname)
        thisfiletype = fname_split[2]
        fname_split = fname_split[1]
        if thisfiletype != filetype
            continue
        end
        @info fname

        try
            combine_2v1a(joinpath(dname,"cam_topview",fname_split*"_topview.mkv"), joinpath(dname,"cam_uw1",fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
        catch err
            combine_2v1a(joinpath(aufolder,fname_split*"_topview.mkv"), joinpath(aufolder,fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
        end
    end
end

function plot_signal2vid(data,fs, outvidname; fps=25, kwargs...)
    sig = signal(data,fs)
    t_width = 1/fps
    min_max = extrema(data)
    anim = @animate for t in 0:1/fps:(size(data,1) / fs)
        plot(sig; kwargs...)
        plot!([t, t+t_width, t+t_width, t], [min_max[1], min_max[1], min_max[2], min_max[2]], fill=false, c=:red, alpha=0.8)
    end

    outgifname = splitext(outvidname)[1]*".gif"
    gif(anim,  outgifname; fps = fps)
    # outvidname = splitext(outgifname)[1]*".mp4"
    @ffmpeg_env run(`$ffmpeg -i $outgifname -pix_fmt yuv420p $outvidname -hide_banner -y`)
end


# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=ow=1920:oh=ih*ow/iw:mode=line:rate=25,format=yuv420p[bottom]" -map "[top]" -map "[bottom]" -y output.mp4
# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=1920x480:mode=line:rate=25,format=yuv420p[bottom]" -map "[top]" -map "[bottom]" -y output.mp4

# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=ow=1920:oh=ih*ow/iw:mode=line:rate=25[bottom_audio]" -map "[top]" -map "[bottom_audio]" -c:a aac -strict experimental -y output.mp4

# input_resolution=$(ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of csv=s=x:p=0 input1.mp4)
# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=ow=${input_resolution%:*}:oh=${input_resolution#*:}/iw:mode=line:rate=25[bottom_audio]" -map "[top]" -map "[bottom_audio]" -c:a aac -strict experimental -y output.mp4


# "[2:a]showwaves=s=1920x400:mode=cline:colors=blue:rate=25[waves];[0:v][1:v][waves]vstack=inputs=3[v];[2:a]loudnorm[a];"