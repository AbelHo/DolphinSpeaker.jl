include("test_run.jl")
# using Pipe
# using DelimitedFiles
# using ProgressMeter
include("config.jl")
include("test_make_video.jl")
include("tonal_detector.jl")
include("utils.jl")
include("plotting.jl")

#~ ######## temporary fix
include("detector_impulsive.jl")
include("config.jl")
include("synchronization.jl")
include("aspod.jl")
include("video.jl")
band_pass = tonal_band_pass
threshold_tonal = -15

#############################

using Glob
using FFMPEG
using JLD2, FileIO
# folname = "/Users/abel/Documents/data/aspod/field/maui_2022/2022.03.01"
# vids = glob("*.MP4", folname)

# folname = "/Volumes/dd/Bahamas_2022/2022.06.25/0001"
# vidtype=r".mkv|.MP4|.avi|.mp4"; autype=r".wav|.mat|.flac.mp3"
function extract_vid_au(folname; vidtype= "*" .* vidtypes, autype= "*" .* autypes)

    vids = []
    for vt in vidtype
        @info vt
        vids = vcat(vids, glob(vt, folname));
        # if !isempty(vids)
        #     break;
        # end
    end

    audios = []
    for vt in autype
        audios = vcat(audios,glob(vt, folname));
        # if !isempty(audios)
        #     break;
        # end
    end

    return vids, audios
end

function match_recording(vids, aus; timediff_tolerance=3, tolerance_dur_to_skip=3, get_duration=get_duration)
    if length(vids)==0 || length(aus)==0
        return vids, aus
    end

    if occursin("_1.mat", aus[1]) #calf recording
        filter!( x->occursin("_1.mat",x), aus)
        if length(vids) == length(aus)
            return vids, aus
        else
            @warn "Different total number of videos and corresponding audio file!!! shorten audio list"
            return vids, aus[1:length(vids)]
        end
    end

    
    vids_new = Vector{typeof(first(vids))}()
    aus_new = Vector{typeof(first(aus))}()

    # len_aus = length(aus)
    # au_durs = get_duration.(aus)
    aus_copy = deepcopy(aus)

    # au_ind = 1
    for vid ∈ vids
        vid_dur = get_duration(vid)
        vid_dur < tolerance_dur_to_skip && continue
        for ind ∈ eachindex(aus_copy)
            if abs(vid_dur - get_duration(aus_copy[ind])) < timediff_tolerance
                @debug (basename(vid), basename(aus_copy[ind]), (vid_dur - get_duration(aus_copy[ind])))
                push!(vids_new, vid)
                push!(aus_new, aus_copy[ind])
                deleteat!(aus_copy, ind)
                break
            end
        end
    end

    return vids_new, aus_new
end

# vids = glob("*.mkv", folname);
# audios = glob("*.wav", folname);
# [get_duration.(vids) get_duration2.(audios)];
# dif = get_duration.(vids) - get_duration2.(audios);
# map( x -> isnan(x) ? 0 : abs(x), dif) |> maximum

function process_dir(folname; func=(a,b)->x, arg=nothing, no_overwrite_func=nothing)
    for (root, dirs, files) in walkdir(folname)
        # println("Directories in $root")
        for dir in dirs
            println(joinpath(root, dir)) # path to directories
            try
                func(joinpath(root, dir), joinpath(arg, dir); no_overwrite_func=no_overwrite_func)
            catch err
                try 
                    func(joinpath(root, dir), joinpath(arg, dir))
                catch err
                    @error exception=(err, catch_backtrace())
                    @error (joinpath(root, dir), joinpath(arg, dir))
                end
            end
        end
        # println("Files in $root")
        # for file in files
        #     println(joinpath(root, file)) # path to files
        # end
    end
end
function process_vidau_dir(folname, res_dir=nothing; no_overwrite_func=nothing)
    @debug folname
    vids, audios = extract_vid_au(folname)
    vids, audios = match_recording(vids, audios)
    if length(vids)==0 || length(audios)==0
        return
    end
    @debug [length(vids), length(audios)]
    process_one_set.(vids, audios, res_dir; no_overwrite_func=no_overwrite_func)

end

function check_output_exist(vidfname, aufname, res_dir; postfix="_t163.835_d800")
    @debug joinpath(res_dir, reduce((a,b) -> a*"_overlaid"*postfix*b , splitext(basename(vidfname))) *"_normalized-audio.mkv")
    isfile( joinpath(res_dir, reduce((a,b) -> a*"_overlaid"*postfix*b , splitext(basename(vidfname))) *"_normalized-audio.mkv") )
end
function process_one_set(vidfname, aufname, res_dir; skiplist=[], no_overwrite_func=nothing, savejld=true, funcs=[x->x], pt_config=pt_config)
    # skiplist = ["Vid_2022-06-09_10.05.57.mkv"]
    @info (vidfname, aufname, res_dir)
    if !isnothing(no_overwrite_func) # skipoverwrite
        if no_overwrite_func(vidfname, aufname, res_dir) ##FIXME no postfix, so it wouldn't reject
            @warn ("Skipping...... "*vidfname)
            return
        end
    end
    # if basename(vidfname) in skiplist
    #     @warn ("Skipping...... "*vidfname)
    #     return
    # end
    @debug "IN"
    data, fs = readAudio(aufname)
    if aufname isa Array
        aufname = aufname[1]
    end
    # pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect, fs = process_audioVideo( (aufname,data,fs), vidfname, res_dir)
    # detection_b = process_audioVideo_tonal1( (aufname,data,fs), vidfname, res_dir)
    # detector_set = [detection_b,
    #                 (pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect)
    #                 ]
    # pt_config = [((1,1,0),25), ((1,0,0),30), ((0,1,0),20), ((1,0,1),15), ((1,1,1),10)]
    
    fps = get_fps(vidfname)
    vidau_syncdiff = findVidAudioBlip(vidfname; plot_window_inS=nothing, band_pass=[2900 3100]) - findAudioBlip(aufname; plot_window_inS=nothing, band_pass=[2900 3100])
    @info "Audio started later by $(vidau_syncdiff)s"
    # # pixel_related = map( (x,y)->(x[1:2]..., y...), detector_set, pt_config)
    # pixel_related = map( (x,y)->((x[1] .+(vidau_syncdiff*fps) .|>round.|>Int,x[2])..., y...), detector_set, pt_config)

    pixel_estimated_set = process_detections(aufname, vidfname; res_dir=res_dir)
    @debug "--------------------------------"
    @debug size(pixel_estimated_set)
    @debug pixel_estimated_set
    
    pixel_estimated_set = pixel_estimated_set[DETECTION_TYPES]
    # @info "---------------- ------------------"
    # @info pixel_estimated_set
    pt_config = pt_config[1:length(pixel_estimated_set)]
    pixel_related = map( (x,y)->((x[1] .+(vidau_syncdiff*fps) .|>round.|>Int,x[2])..., y...), pixel_estimated_set, pt_config)


    try
        if savejld
            savefname = splitext(basename(aufname))[1] *"_t"*string(impulsive_autothreshold_median_ratio)*"_d"*string(dist_impulsive)*".jld2"
            if !Sys.islinux()
                jldsave(joinpath(res_dir, savefname); pixel_estimated_set, pixel_related)
                # jldsave(joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".jld2"); pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect,  vidfname, aufname, res_dir, detector_set, pixel_related)
            #@error(err)
            #@error("Cant save jld2 file, saving locally and copying instead")
            #rm(joinpath(res_dir, splitext(basename(aufname))[1] *"_t"*string(thresh)*"_d"*string(dist)*".jld2"))
            else
                jldsave(joinpath("", savefname); pixel_estimated_set, pixel_related)
                mv(joinpath("", savefname),
                    joinpath(res_dir, savefname))
            end
            
        end
    catch err
        @error(err)
        savefname = splitext(basename(aufname))[1] *"_t"*string(impulsive_autothreshold_median_ratio)*"_d"*string(dist_impulsive)*".jld2"
        @error joinpath(res_dir,savefname)
        @error("Failed to save JLD file")
    end

    # newvidname = process_video(vidfname, res_dir; func=overlay_points!, extra_arg=pixel_related, postfix="_t"*string(thresh)*"_d"*string(dist),
    #              func2=plot_summary!, extra_arg2=(data, fs, vidau_syncdiff))
    newvidname = process_video(vidfname, res_dir; func=overlay_points!, extra_arg=pixel_related, postfix="_t"*string(impulsive_autothreshold_median_ratio)*"_d"*string(dist_impulsive),
                 extra_arg2=(data, fs, vidau_syncdiff))
    # newvidname = process_video(vidfname, res_dir; func=overlay_points!, extra_arg=pixel_related, postfix="_t"*string(thresh)*"_d"*string(dist),
                #  func2=plot_summary_plots_img!, extra_arg2=(data, fs, vidau_syncdiff))
    try
        aufname_old = nothing;
        if aufname[end-2:end] == "mat"
            aufname_old = aufname;
            aufname = joinpath(res_dir, "temp__" * (splitext(aufname)[1]*".wav" |> basename))
            wavwrite(aufname, data ./ maximum(data), fs)
        end
        println("RAM: ", round(Sys.free_memory() / 1024 / 1024 / 1024, digits=2), "/", round(Sys.total_memory() / 1024 / 1024 / 1024, digits=2), " GB")
        GC.gc()
        println("RAM: ", round(Sys.free_memory() / 1024 / 1024 / 1024, digits=2), "/", round(Sys.total_memory() / 1024 / 1024 / 1024, digits=2), " GB")

        #~ combine video and audio
        if MERGE_VID_AU_DYNAMIC_NORM
            println(`$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 -f matroska "$newvidname""_normalized-audio.mkv"`)
            output = @ffmpeg_env run(`$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 "$newvidname""_normalized-audio.mp4"`)
            # run(`ffmpeg -i "$newvidname" -i "$aufname" -map 0:v -map 1:a -vcodec copy -af loudnorm=I=-16:LRA=11:TP=-1.5 -f matroska "$newvidname""_normalized-audio.mkv"`)
        else
            # Get max volume and normalize
            tempfile = tempname()
            # read(pipeline(`ffmpeg -i $aufname -filter:a volumedetect -f null /dev/null`; stderr = tempfile))
            @ffmpeg_env read(pipeline(`ffmpeg -i $aufname -af astats=metadata=1 -f null /dev/null`; stderr = tempfile))

            ss = read(tempfile, String)
            rm(tempfile)
            m = eachmatch(r"Channel: (\d+)", ss)
            channel_list = [parse(Int, m.captures[1]) for m in m]
            m = eachmatch(r"Peak level dB: (.*)", ss)
            peak_db_list = [parse(Float64, m.captures[1]) for m in m]
            norm_gain = -maximum(peak_db_list[get_relevant_channels(rx_vect)])

            # m = match(r"max_volume: (.*) dB", ss)
            # max_volume = m !== nothing ? parse(Float64, m.captures[1]) : nothing
            # norm_gain = -max_volume
            cmd = `$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af "volume=$(norm_gain)dB" "$newvidname""_normalized-audio.mp4"`
            println(cmd)
            output = @ffmpeg_env run(cmd)
        end
        isfile("$newvidname"*"_normalized-audio.mp4") && rm(newvidname) # delete video without audio
        if aufname_old isa String; rm(aufname); end
    catch err
        @error(err)
        @error("Failed to add audio to overlaid video: " * newvidname)
    end

   
    return aufname, vidfname, vidau_syncdiff, pixel_estimated_set

    # return pind_vidframes, p_pixels, thresh, dist, ang, tdoa_raw, tdoa, window, threshold_indices, pind_good, pind_good_inS, pind, ppeak, ref_channel, c, rx_vect, newvidname
end

# process_dir(folname; func=process_vidau_dir, arg=res_dir, no_overwrite_func=check_output_exist)

process_folder(foldername, outfolder; kwargs...) = process_folder(foldername; outfolder=outfolder, kwargs...)


function process_folder(foldername; outfolder=foldername, skipdone=false, kwargs...)
    mkpath(outfolder)
    a=split_vid_au(foldername)
    # process_one_set.(a[1], a[2], outfolder)

    for fname in a[2] |> skiphiddenfiles #readdir(foldername)|>skiphiddenfiles
        fname_split = splitext(fname)[1]
        
        if skipdone 
            if !isempty(filter( x->occursin(Regex("(?=.*" *basename(fname_split)* ")(?=.*_normalized-audio.mp4)"),x), readdir(outfolder)))
                @info "Done! Skipping: $fname"
                continue
            end
        end

        try
            vidfname = fname_split*vidtypes[1]
            if !isfile(vidfname)
                for vt in vidtypes[2:end]
                    vidfname = fname_split*vt
                    if isfile(vidfname)
                        break
                    end
                end
            end
            process_one_set(joinpath(foldername, vidfname), fname, outfolder; kwargs...)
        catch err
            @error(err)
            @error "Failed to process: "*fname exception=(err, catch_backtrace())
        end
        # combine_2v1a(joinpath(dname,"cam_topview",fname_split*"_topview.mkv"), joinpath(dname,"cam_uw1",fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
    end
end


"""
    run_analysis_split_vidau(folname; res_dir="")

Analyze and synchronize video and audio files within a specified folder.

# Arguments
- `folname::AbstractString`: Path to the folder containing video and audio files.
- `res_dir::AbstractString=""`: Optional directory to store results. Defaults to a subdirectory named after `folname`.

# Description
This function:
1. Identifies video and audio files in `folname` using predefined file type lists (`vidtypes`, `autypes`).
2. Combines all video files into a single video using FFmpeg.
3. Computes the synchronization delay between the first video and audio file using `find_vid_vs_audio_syncdiff_timesegment`.
4. Saves the synchronization delay and confidence score to a CSV file in the results directory.

# Returns
- `delays`: The computed synchronization delay (in seconds).
- `conf`: Confidence score of the synchronization.
"""
function run_analysis_split_vidau(folname; res_dir="",
     flag_verbose=true, kwargs...)

    # occursin.( Ref(Regex(join(vidtypes, '|'))), readdir(folname))
    # vidlist = filter( x -> occursin(Regex(join(vidtypes, "|\\")), x|>lowercase), readdir(folname; join=true))
    # audlist = filter( x -> occursin(Regex(join(autypes, "|\\")), x|>lowercase), readdir(folname; join=true))
    vidlist, audlist = split_vid_au(folname; readdir_func=readdir_all)

    delays, conf = find_vid_vs_audio_syncdiff_timesegment(vidlist[1], audlist[1]; flag_verbose=flag_verbose, flag_return_conf=true, kwargs...)
    # combine all video files into one file
    if !isempty(res_dir) && !isnothing(res_dir)
        res_dir = joinpath(res_dir, basename(folname))
        mkpath(res_dir)
        write(joinpath(res_dir, "sync_delay.csv"), "foldername,delay_s,confidence\n$(basename(folname)),$(delays),$(conf)\n")

        temp_filelist = joinpath(res_dir, "temp_filelist.txt")
        write(temp_filelist, join(["file '$v'" for v in vidlist], "\n"))
        output_vidname = "$res_dir/combined__$(join(basename.(vidlist), '_')).mp4"
        cmd = `ffmpeg -hide_banner -loglevel error -f concat -safe 0 -i $temp_filelist -c copy $output_vidname`
        print(cmd)
        @async try
            @ffmpeg_env run(cmd)
            rm(temp_filelist)
            @info "Succesfully combined videos into $output_vidname !!!"
        catch e
            @error "FFmpeg command failed: $e"
            rm(temp_filelist)
        end
    end

    return delays, conf, output_vidname, vidlist, audlist, res_dir
end

function run_contiguous_folders(folname; res_dir="", overlay_radius=32, detection_types = 1:2,
    flag_overlayvideo=true, flag_overlayimages=false, flag_return = false, kwargs...)
    @info "Processing folder: $folname ............."
    
    delays, conf, output_vidname, vidlist, audlist, res_dir2 = run_analysis_split_vidau(folname; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true, kwargs...)
    # results = process_detections.(audlist, Ref(vidlist[1]); res_dir=res_dir2)
    # run process_detections on each audio file in parallel, preserving order
    results = Vector{Any}(undef, length(audlist))
    Threads.@threads for i in eachindex(audlist)
        try
            results[i] = process_detections(audlist[i], vidlist[1]; res_dir=res_dir2, flag_return=flag_return)
        catch err
            @error "process_detections failed for $(audlist[i])" exception=(err, catch_backtrace())
            results[i] = nothing
        end
    end
    # res[1][1] = res[1][1] .+ (delays*get_fps(vidfname))

    detection_pixels = joinpath(res_dir2, "detection_pixels.csv")
    # detection_types = 1:2  # 1: impulsive, 2: tonal, 3: boat
    # cum_duration = 0.0
    open( detection_pixels, "w") do io
        writedlm(io, ["frame" "px" "py" "radius" "r" "g" "b" "a" "shape"], ',')
        for detection_type in detection_types
            # write_mode = detection_type==1 ? "w" : "a"
            cum_duration = 0.0
            for (ind, res) in enumerate(results)
                writedlm(io, [res[detection_type][1] .+ ( (cum_duration + delays)*get_fps(vidlist[1])) res[detection_type][2]], ',') # add sync delay and cumulative duration
                # write_mode = ind==1 ? "w" : "a"
                # open( detection_pixels, write_mode) do io
                #     (write_mode == "w") && writedlm(io, ["frame" "px" "py" "radius" "r" "g" "b" "a" "shape"], ',')
                #     writedlm(io, [res[detection_type][1] .+ ( (cum_duration + delays)*get_fps(vidlist[1])) res[detection_type][2]], ',') # add sync delay and cumulative duration
                # end
                cum_duration += get_duration(audlist[ind])
                # results[1].res.res_impulsetrain.pind_good
            end
        end
    end
    
    detection_pixels_df = CSV.read(detection_pixels, DataFrame)
    # detection_pixels_df.tag_name .= map(x -> x==1 ? :impulsive : x==2 ? :tonal : :boat, detection_pixels_df.type)
    detection_pixels_df.type = detection_pixels_df.shape
    detection_pixels_df[:, ["r", "g", "b"]] = detection_pixels_df[:, ["r", "g", "b"]] .* 255.0 # convert to 0-255 range
    CSV.write(splitext(detection_pixels) |> x-> x[1]*"__color255_webui"*x[2], detection_pixels_df)

    dfs_impulse = CSV.read.(joinpath.(res_dir2 |> Ref, [results[i].res.res_impulsetrain.outfname * ".txt" for i in 1:length(results)]), DataFrame; header=false)
    dfs_tonal = CSV.read.(joinpath.(res_dir2 |> Ref, [results[i].res.res_tonalsegment.outfname * ".txt" for i in 1:length(results)]), DataFrame; header=false)
    dfs_impulsetrain = CSV.read.(joinpath.(res_dir2 |> Ref, [results[i].res.res_impulsetrain.outfname * "_train-only.txt" for i in 1:length(results)]), DataFrame; header=false)

    cum_duration = 0.0; cum_index_impulse = 0; cum_index_tonal = 0; cum_index_impulsetrain = 0;
    for ind = 1:length(results)
        dfs_impulse[ind][:, 1:2] .+= cum_duration
        dfs_tonal[ind][:, 1:2] .+= cum_duration
        dfs_impulsetrain[ind][:, 1:2] .+= cum_duration

        dfs_impulse[ind][:, 3] .+= cum_index_impulse
        dfs_tonal[ind][:, 3] .+= cum_index_tonal
        dfs_impulsetrain[ind][:, 3] .+= cum_index_impulsetrain

        cum_duration += get_duration(audlist[ind])
        cum_index_impulse += nrow(dfs_impulse[ind])
        cum_index_tonal += nrow(dfs_tonal[ind])
        cum_index_impulsetrain += nrow(dfs_impulsetrain[ind])
    end
    CSV.write(joinpath(res_dir2, "combined_impulse__audacity.txt"), vcat(dfs_impulse...); writeheader=false, delim='\t')
    CSV.write(joinpath(res_dir2, "combined_tonal__audacity.txt"), vcat(dfs_tonal...); writeheader=false, delim='\t')
    CSV.write(joinpath(res_dir2, "combined_impulsetrain__audacity.txt"), vcat(dfs_impulsetrain...); writeheader=false, delim='\t')

    vid_ready = false
    while !vid_ready
        try
            @ffmpeg_env run(`ffprobe $output_vidname`)
            vid_ready = true
        catch err
            @warn "Video file not ready yet: $output_vidname. Retrying in 30 seconds..."
            sleep(30)
        end
    end
    @info "Combined video ready: $output_vidname, proceeds..."
    # flag_overlayvideo && overlay_boxes_on_video(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaid.mp4"; radius=overlay_radius)
    if flag_overlayvideo
        try
            out_vid_path = overlay_annotations_on_video(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaid.mkv"; 
                mode=:VideoIO, radius=overlay_radius, #) #mode=:stream) #
                kwargs...)
                # radius=:in_annotations, default_color=:in_annotations, default_alpha=:in_annotations, default_shape=:in_annotations)
            combine_vidau(out_vid_path, audlist; vidau_syncdiff=delays, MERGE_VID_AU_DYNAMIC_NORM=true, rx_vect=rx_vect, kwargs...)
        catch err
            @error "Failed to overlay boxes on video($output_vidname)" exception=(err, catch_backtrace())
        end
    end
    flag_overlayimages && 
    (overlay_boxes_on_video_imageonly(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaidIMG"; 
        radius=(overlay_radius isa Number ? overlay_radius : OVERLAY_RADIUS));
    pic2vid(splitext(output_vidname)[1]*"_overlaidIMG", splitext(output_vidname)[1]*"_overlaidIMG.mp4"; auto_mode=true)
    )

    # vidpath = "/media/spin/anas2/data_res/dolphin/calf/temp/delete/1/combined__1.GoPro_Clicker.MP4.mp4"

    # task_overlayvideo = @async overlay_boxes_on_video(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaid.mp4"; radius=100)
    # task_overlayvideo_img = @async overlay_boxes_on_video_imageonly(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaidIMG"; radius=100)
    # wait(task_overlayvideo)
    # wait(task_overlayvideo_img)
    # pic2vid(splitext(output_vidname)[1]*"_overlaidIMG", splitext(output_vidname)[1]*"_overlaidIMG.mp4"; auto_mode=true)
    if flag_return
        return (;results, detection_pixels_df, delays, conf, output_vidname, vidlist, audlist, res_dir2)
    end
end


function test()
    println("test 2 ...")
end