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
function extract_vid_au(folname; vidtype=["*.mkv","*.MP4","*.avi","*.mp4"], autype=["*.wav","*.mat","*.flac","*.mp3"])

    vids = []
    for vt in vidtype
        vids = glob(vt, folname);
        if !isempty(vids)
            break;
        end
    end

    audios = []
    for vt in autype
        audios = glob(vt, folname);
        if !isempty(audios)
            break;
        end
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
function process_one_set(vidfname, aufname, res_dir; skiplist=[], no_overwrite_func=nothing, savejld=true, funcs=[x->x])
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
    pt_config = [((1,1,0),25), ((1,0,0),30), ((0,1,0),20), ((1,0,1),15), ((1,1,1),10)]
    
    fps = get_fps(vidfname)
    vidau_syncdiff = findVidAudioBlip(vidfname; plot_window_inS=nothing, band_pass=[2900 3100]) - findAudioBlip(aufname; plot_window_inS=nothing, band_pass=[2900 3100])
    @info "Audio started later by $(vidau_syncdiff)s"
    # # pixel_related = map( (x,y)->(x[1:2]..., y...), detector_set, pt_config)
    # pixel_related = map( (x,y)->((x[1] .+(vidau_syncdiff*fps) .|>round.|>Int,x[2])..., y...), detector_set, pt_config)

    pixel_estimated_set = process_detections(aufname, vidfname; res_dir=res_dir)
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

        # combine video and audio
        println(`$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 -f matroska "$newvidname""_normalized-audio.mkv"`)
        output = @ffmpeg_env run(`$ffmpeg -i "$newvidname" -itsoffset $vidau_syncdiff -i "$aufname" -map 0:v -map 1:a -pix_fmt yuv420p -af loudnorm=I=-16:LRA=11:TP=-1.5 "$newvidname""_normalized-audio.mp4"`)
        isfile("$newvidname"*"_normalized-audio.mp4") && rm(newvidname) # delete video without audio
        # run(`ffmpeg -i "$newvidname" -i "$aufname" -map 0:v -map 1:a -vcodec copy -af loudnorm=I=-16:LRA=11:TP=-1.5 -f matroska "$newvidname""_normalized-audio.mkv"`)
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
            process_one_set(joinpath(foldername, fname_split*".mkv"), fname, outfolder; kwargs...)
        catch err
            @error(err)
            @error "Failed to process: "*fname exception=(err, catch_backtrace())
        end
        # combine_2v1a(joinpath(dname,"cam_topview",fname_split*"_topview.mkv"), joinpath(dname,"cam_uw1",fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
    end
end