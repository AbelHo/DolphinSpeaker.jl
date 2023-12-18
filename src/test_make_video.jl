using VideoIO, FFMPEG
using ProgressMeter
using ImageDraw
include("config.jl")
using GLMakie

# res_fol = "/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
# newvidname = process_video(vidfname, res_fol; func=overlay_point!, extra_arg=(pind_vidframes, p_pixels), postfix="_t"*string(thresh)*"_d"*string(dist))

# function get_number_frames(file::AbstractString, streamno::Integer = 0)
#     streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
#     frame_strs = FFMPEG.exe(
#         `-v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 $file`,
#         command = FFMPEG.ffprobe,
#         collect = true,
#     )
#     frame_str = frame_strs[1]
#     if occursin("No such file or directory", frame_str)
#         error("Could not find file $file")
#     elseif occursin("N/A", frame_str)
#         return nothing
#     end
#     return parse(Int, frame_str)
# end
# function get_fps(file::AbstractString, streamno::Integer = 0)
#     streamno >= 0 || throw(ArgumentError("streamno must be non-negative"))
#     fps_strs = FFMPEG.exe(
#         `-v 0 -of compact=p=0 -select_streams 0 -show_entries stream=r_frame_rate $file`,
#         command = FFMPEG.ffprobe,
#         collect = true,
#     )
#     fps = split(fps_strs[1], '=')[2]
#     if occursin("No such file or directory", fps)
#         error("Could not find file $file")
#     elseif occursin("N/A", fps)
#         return nothing
#     end
# 	return reduce(//, parse.(Int, split(fps,'/')) )
#     # return round(reduce(/, parse.(Float64, split(fps,'/')) ), digits=3)
# end

function display_corners!(img, corners_list; dot_size=25, colors=[1,0,0])
    len_cl = length(corners_list)
    for j in 1:len_cl
        corner = corners_list[j]
        for i in 1:size(corner,2)
            draw!(img, Ellipse(CirclePointRadius(Int(round(corner[1,i])), Int(round(corner[2,i])), dot_size)), typeof(img[1])(colors...))
            #@debug [corner[1,i], corner[2,i]]
        end
    end
    return img
end

# encoder_options = (crf=23, preset="ultrafast") #(crf=0, preset="ultrafast")
# vidfname = "/Volumes/dd/Bahamas_2022/2022.06.27/0000/Vid_2022-06-27_10.47.50.mkv"
# res_fol = "/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
# newvidname = joinpath(res_fol, reduce((a,b) -> a*"_overlaid"*b , splitext(basename(vidfname))) )

# res_fol = "/Users/abel/Documents/data_res/aspod/real/bahamas_2022"
# newvidname = process_video(vidfname, res_fol; func=overlay_point!, extra_arg=(pind_vidframes, p_pixels), postfix="_t"*string(thresh)*"_d"*string(dist))

# process_video(vidfname, res_fol; func=overlay_point_auto!, extra_arg=1080) # test

# run(`combineVidAu.sh -n $newvidname $aufname $newvidname`)
func2_I(a;b) = a[4]
func2_I(a) = a[4]

function process_video(vidfname, res_fol; func=(a,b,c)->(a,b,c), extra_arg=nothing, encoder_options=(crf=23, preset="ultrafast"), postfix="", func2=func2_I, extra_arg2=nothing, ref_channel=ref_channel)
    newvidname = joinpath(res_fol, reduce((a,b) -> a*"_overlaid"*postfix*b , splitext(basename(vidfname))))
    
    vid = VideoIO.openvideo(vidfname)
    img = nothing
    try
        img = read(vid)
    catch err
        @warn "Failed to read video first frame, trying again"
        @warn err
        img = read(vid)
    end
    imsize = raw_frame_size(vid)

    img2 = func(img, 1, extra_arg)
    if func2 != func2_I
        @info "extra plot function"
        data, fs, vidau_syncdiff = extra_arg2
        fig = Figure(;resolution=result_resolution)
        aspect_ratio = reduce(//, result_resolution);
        fig[denominator(aspect_ratio), numerator(aspect_ratio)] = GridLayout();
        img2 = func2(fig, @view(data[1:fs,1]), fs, img)
    end
    # event_frames = round.(Int, pind_good_inS * fps) .+ 1
    # @info size(img2) # (newvidname, img2, get_fps(vidfname), encoder_options=encoder_options)
    @time open_video_out(newvidname, img2, framerate=get_fps(vidfname), encoder_options=encoder_options) do writer
        seekstart(vid)
        num_of_frames = get_number_frames(vidfname)
        fps = get_fps(vidfname)

        counter=1
        # while !eof(vid)
        @showprogress "Encoding video frames.." for i in 1:num_of_frames
            # @debug string([counter, num_of_frames])
            try
                read!(vid, img)
            catch err
                @error "Failed to read video frame $counter/$fps, trying again"
                @warn "producing a blank frame and continue...................................."
                img = zeros(eltype(img), size(img))
                try 
                    seekstart(vid); seek(vid, (counter+1)/fps)
                catch err
                    @error "failed seeking to next frame, skip frame again......"
                end
            end
            func(img, counter, extra_arg;)
            # display_corners!(img, [[mod(counter, imsize[1]);200]])
            # Do something with frames
            
            if func2 == func2_I
                write(writer, img)
            else
                vTime = (counter-1)/fps |> Float64
                aTime = vTime - vidau_syncdiff
                aFrame = round(Int, aTime*fs+1)
                maxWindow = fs*5
                plot_window = aFrame .+ (1:maxWindow)

                if aFrame < 1
                    plot_window = 1:plot_window[end]
                end
                if aFrame + maxWindow > size(data,1)
                    plot_window = plot_window[1]:size(data,1)
                end

                if length(plot_window) == maxWindow
                    plot_data =  @view(data[plot_window,ref_channel])
                else
                    @debug (typeof(plot_window), size(data), typeof(data), aFrame)
                    @debug size(plot_window)
                    @debug "TEST"
                    if aFrame < 1
                        plot_data = [zeros(-aFrame); @view(data[plot_window,ref_channel])]
                    else
                        plot_data =  @view(data[plot_window,ref_channel])
                    end
                    if length(plot_data) != maxWindow
                        plot_data = [plot_data; zeros(maxWindow-length(plot_data))]
                    end
                    @debug (size(plot_data), typeof(plot_data))
                end
                # func2(fig, plot_data, fs, img; 
                    # title="vTime: $(vTime)s Frame: $counter aTime: $aTime")
                write(writer, func2(fig, plot_data, fs, img; 
                    title="vTime: $(vTime)s Frame: $counter aTime: $aTime") )
            end
            counter += 1
        end
        close(vid)
    end
    return newvidname
end

overlay_point_auto!(img, counter, extra_arg) = display_corners!(img, [[mod(counter, size(img,2)); extra_arg]])

function overlay_point!(img, counter, extra_arg) 
    pind_vidframes, p_pixels = extra_arg
    
    detections_in_frame_ind = findall(counter .== pind_vidframes)
    isempty(detections_in_frame_ind) && return nothing

    @debug detections_in_frame_ind
    @debug [p_pixels[detections_in_frame_ind,:]]

    try
        display_corners!(img, [p_pixels[detections_in_frame_ind,:]]')
    catch err
        @error "DISPLAY DETECTION FAILED!"
        @error detections_in_frame_ind
        @error [p_pixels[detections_in_frame_ind,:]]'
    end
end

function overlay_points!(img, counter, extra_arg) 
    # pind_vidframes_list, p_pixels_list, colour, ptsize = extra_arg
    
    for element in extra_arg #eachindex(pind_vidframes_list)
        pind_vidframes, p_pixels, colour, ptsize = element

        detections_in_frame_ind = findall(counter .== pind_vidframes)
        isempty(detections_in_frame_ind) && continue #return nothing

        @debug (counter, detections_in_frame_ind, p_pixels[detections_in_frame_ind,:])
        # open("~/Documents/data_res/aspod/temp/vid_out.txt", "a") do f
        #     write(f, string(counter, detections_in_frame_ind, p_pixels[detections_in_frame_ind,:])*"\n")
        # end
        

        try
            display_corners!(img, [p_pixels[detections_in_frame_ind,:]]'; dot_size=ptsize, colors=colour)
        catch err
            @error "DISPLAY DETECTION FAILED!"
            @error detections_in_frame_ind
            @error [p_pixels[detections_in_frame_ind,:]]'
        end
    end
    return img
end

@info "END"

