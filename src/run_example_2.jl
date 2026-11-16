# run_contiguous_folders_2: a NEW COPY of run_contiguous_folders (run_example.jl, untouched) with
#   (a) frames computed from exact detection times through a sync model (drift-aware affine fit from
#       fit_folder_sync_affine, or the existing constant delay),
#   (b) an optional NaN-safe overlay CSV (detection_pixels_valid.csv).
# Everything else (audacity files, overlay video, kwargs routing) behaves like the original.
# Included from DolphinSpeaker.jl after run_example.jl; or lazily from run_contiguous_folders_calibrated.
# Frame convention (lead decision, T6 finding): the original draws detections 1-2 frames late (aspod adds +1, the
# overlay uses 0-based frame_idx and round). Here frame = floor(t_video*fps), a 0-based Integer, no +1.
include("detection_pixels_valid.jl")

# Exact detection times (s, within one audio file) for pixel_related_* index `detection_type`, or `nothing`
# when they cannot be recovered with the same length as the pixel rows (caller falls back to the rounded frames).
function _rc2_detection_times(rec,detection_type,nrows,fps)
    try
        r=rec.res
        t=if detection_type in (1,3,5)
            r.res_impulsetrain.pind_good_inS
        elseif detection_type==2
            seg=r.res_tonalsegment
            PARAM_TONALSPREAD ? vcat([collect(seg.train_start[i]/seg.fs : 1/fps : seg.train_end[i]/seg.fs) for i in eachindex(seg.train_start)]...) :
                seg.pind_good_inS
        elseif detection_type==4
            r.res_tonal.pind_good_inS
        else
            nothing
        end
        t!==nothing && length(t)==nrows && return Float64.(t)
    catch err
        @debug "exact detection times unavailable" detection_type exception=err
    end
    return nothing
end

# Chapter table from the sync JSON's video_timeline.files (T5): where each chapter starts on the global video
# timeline and how many frames it adds to the concatenated overlay video. `nothing` when it does not describe `vidlist`.
function _rc2_timeline(info,vidlist)
    files=try info["video_timeline"]["files"] catch; return nothing end
    files isa AbstractVector && length(files)==length(vidlist) && !isempty(files) || return nothing
    rows=NamedTuple[]; cum=0
    for (f,v) in zip(files,vidlist)
        keys_ok=all(k->get(f,k,nothing) isa Real && isfinite(f[k]),("global_start_s","first_pts_s","n_frames"))
        keys_ok && get(f,"path","")!="" && abspath(string(f["path"]))==abspath(string(v)) || return nothing
        push!(rows,(global_start_s=Float64(f["global_start_s"]),audio_start_s=Float64(something(get(f,"audio_start_s",0.0),0.0)),
            first_pts_s=Float64(f["first_pts_s"]),n_frames=Int(f["n_frames"]),cum_before=cum))
        cum+=Int(f["n_frames"])
    end
    issorted(r.global_start_s for r in rows) || return nothing
    return rows
end

# (scale, offset_s, mode, timeline): t_video = scale*t_recorder + offset. mode is :affine or :constant.
# timeline is the chapter table (affine only); the constant delay keeps the original flat frame mapping.
function _rc2_sync_from_info(info,delays,vidlist=String[])
    constant=(scale=1.0,offset=Float64(delays),mode=:constant,timeline=nothing)
    if info!==nothing && get(info,"status","")=="valid" && all(k->get(info,k,nothing) isa Real && isfinite(info[k]),("scale","offset_s"))
        @info "Affine sync: scale=$(info["scale"]) offset=$(info["offset_s"]) s, drift=$(get(info,"drift_ppm",missing)) ppm, rms=$(get(info,"rms_residual_s",missing)) s"
        timeline=_rc2_timeline(info,vidlist)
        if timeline===nothing && length(vidlist)>1
            @warn "Sync JSON has no usable video_timeline for these chapters; using the flat frame mapping (frames may be off by about 2 after chapter 1)"
        end
        return (scale=Float64(info["scale"]),offset=Float64(info["offset_s"]),mode=:affine,timeline=timeline)
    end
    info===nothing || @warn "Affine sync not valid (status=$(get(info,"status","missing"))); using the constant delay $(delays) s"
    length(vidlist)>1 && @warn "Constant delay: flat frame mapping over $(length(vidlist)) chapters; chapter audio and video durations differ by tens of ms, so frames may be off by about 2 after chapter 1"
    return constant
end

# Frames (0-based Int, index in the concatenated overlay video) for recorder times `t_rec` (s on the concatenated
# audio). Affine with a chapter table: global video time -> (chapter, local time) like prep_match_frame
# (local = t - global_start + audio_start), frame = frames of earlier chapters + floor((local - first_pts)*fps),
# clamped into the chapter. Otherwise the flat floor(t_video*fps). Returns (frames, n_clamped).
function _rc2_frames(model,t_rec,fps)
    tv=model.scale .* t_rec .+ model.offset
    model.timeline===nothing && return (floor.(Int,tv .* fps),0)
    files=model.timeline; starts=[f.global_start_s for f in files]; clamped=0
    frames=map(tv) do t
        i=searchsortedlast(starts,t)
        i==0 && return floor(Int,t*fps)   # before the first chapter: flat
        f=files[i]
        k=floor(Int,(t-f.global_start_s+f.audio_start_s-f.first_pts_s)*fps)
        (k<0 || k>=f.n_frames) && (clamped+=1)
        f.cum_before+clamp(k,0,f.n_frames-1)
    end
    return (frames,clamped)
end

# A cached sync_affine.json is reused only when it is valid and its recorded inputs (path, size, mtime of every
# video and recorder file) equal the current files.
function _rc2_sync_cache_ok(cached,vidlist,audlist)
    cached isa AbstractDict && get(cached,"status","")=="valid" || return false
    inputs=get(cached,"inputs",nothing); inputs isa AbstractDict || return false
    for (key,list) in (("videos",vidlist),("recorders",audlist))
        recs=get(inputs,key,nothing)
        recs isa AbstractVector && length(recs)==length(list) || return false
        for (r,p) in zip(recs,list)
            current=_clicker_call(:_fsync_file_record,p)
            r isa AbstractDict && all(k->get(r,k,nothing)==current[k],("path","size_bytes","mtime_utc")) || return false
        end
    end
    return true
end

function _rc2_sync_model(sync,delays,vidlist,audlist,res_dir2;refit=false)
    sync in (:affine,:constant) || throw(ArgumentError("sync must be :affine or :constant"))
    sync===:constant && return _rc2_sync_from_info(nothing,delays,vidlist)
    out_path=joinpath(res_dir2,"sync_affine.json")
    info=nothing
    try
        _ensure_clicker_module()
        if !refit && isfile(out_path)
            cached=try _clicker_call(:read_json,out_path) catch; nothing end
            if _rc2_sync_cache_ok(cached,vidlist,audlist)
                @info "Reusing cached $out_path (inputs unchanged; refit=true forces a new fit)"
                info=cached
            end
        end
        info===nothing && (info=_clicker_call(:fit_folder_sync_affine,vidlist,audlist;out_path=out_path))
    catch err
        @warn "fit_folder_sync_affine failed; using the constant delay" exception=(err,catch_backtrace())
    end
    return _rc2_sync_from_info(info,delays,vidlist)
end

# recursively run all folder if given an array of folders
function run_contiguous_folders_2(folname::AbstractVector; flag_return_results_array=false, kwargs...)
    results = Vector{Any}(undef, length(folname))
    for (ind, fol) in enumerate(folname)
        try
            if flag_return_results_array
                results[ind] = run_contiguous_folders_2(fol; kwargs...)
            else
                run_contiguous_folders_2(fol; kwargs...)
            end
        catch err
            @error "run_contiguous_folders_2 failed for.........................\n$fol" exception=(err, catch_backtrace())
        end
    end
    return results
end

"""
    run_contiguous_folders_2(folname; sync=:affine, pixel_filter=true, kwargs...)

Copy of `run_contiguous_folders` with frames from exact detection times through the sync model.
`sync=:affine` fits `t_video = scale·t_recorder + offset` for the whole folder
(`fit_folder_sync_affine`, written to `sync_affine.json`) and falls back to the constant delay with a
warning unless its status is `"valid"`; `sync=:constant` is the original constant delay. With
`pixel_filter=true` a NaN-free `detection_pixels_valid.csv` (plus a status CSV) is written and the overlay
reads it; the overlay never sees the raw CSV, and `pixel_filter=false` throws if any row has a non-finite
frame/px/py (one NaN would abort `overlay_annotations_on_video` with an `InexactError`).

Frame convention (differs from `run_contiguous_folders`, which draws 1-2 frames late): the written `frame` is a
0-based integer index into the concatenated overlay video, no `+1`. With `t_video = scale·(t_in_file +
cumulative_audio_duration) + offset` (`video_time_from_recorder`), `sync=:affine` maps `t_video` to a chapter and
local time through the sync JSON's `video_timeline.files` like `prep_match_frame` (`local = t_video - global_start +
audio_start`) and writes `frames of earlier chapters + floor((local - first_pts)·fps)`, so chapter audio/video
duration differences do not accumulate. The constant delay (`sync=:constant`, or the fallback, with a warning for
several chapters) keeps the flat `floor((t_in_file + cumulative_audio_duration + delay)·fps)`. `t_in_file` is the exact
detection time (`pind_good_inS`), not a rounded frame. A valid `sync_affine.json` in the results folder whose recorded
inputs (path, size, mtime) still match is reused instead of the ~13 min refit; `sync_refit=true` forces a refit.
The overlay's temporary frames go to `res_dir/<folder>/tmp_overlay` (not `/tmp`; keyword `tmpdir`), removed afterwards.
Returns, with `flag_return=true`, the original's named tuple plus `sync_model` and `overlay_pixels`.
"""
function run_contiguous_folders_2(folname; res_dir="", overlay_radius=32, detection_types = 1:2,
    flag_overlayvideo=true, flag_overlayimages=false, flag_return = false,
    flag_spectro=false, spectro_height=200, sync=:affine, pixel_filter=true, tmpdir=nothing, sync_refit=false, kwargs...)
    sync in (:affine,:constant) || throw(ArgumentError("sync must be :affine or :constant"))
    isdefined(@__MODULE__,:write_valid_detection_pixels) ||
        error("run_contiguous_folders_2 needs write_valid_detection_pixels (src/detection_pixels_valid.jl)")
    @info "Processing folder: $folname ............."

    delays, conf, output_vidname, vidlist, audlist, res_dir2 = run_analysis_split_vidau(folname; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true, kwargs...)
    results = Vector{Any}(undef, length(audlist))
    Threads.@threads for i in eachindex(audlist)
        try
            results[i] = process_detections(audlist[i], vidlist[1]; res_dir=res_dir2, flag_return=true)
        catch err
            @error "process_detections failed for $(audlist[i])" exception=(err, catch_backtrace())
            results[i] = nothing
        end
    end

    sync_model = _rc2_sync_model(sync, delays, vidlist, audlist, res_dir2; refit=sync_refit)
    n_clamped = 0
    fps = get_fps(vidlist[1])
    detection_pixels = joinpath(res_dir2, "detection_pixels.csv")
    open( detection_pixels, "w") do io
        writedlm(io, ["frame" "px" "py" "radius" "r" "g" "b" "a" "shape"], ',')
        for detection_type in detection_types
            cum_duration = 0.0
            for (ind, res) in enumerate(results)
                if !isnothing(res)
                    rounded_frames = res[detection_type][1]
                    times = _rc2_detection_times(res, detection_type, length(rounded_frames), fps)
                    # exact time when known; otherwise undo the original's 1-origin rounded frame
                    times === nothing && (times = (rounded_frames .- 1) ./ fps)
                    frames, nc = _rc2_frames(sync_model, cum_duration .+ times, fps)   # 0-based, index in the concatenated video
                    n_clamped += nc
                    writedlm(io, [frames res[detection_type][2]], ',')
                end
                cum_duration += get_duration(audlist[ind])
            end
        end
    end

    detection_pixels_df = CSV.read(detection_pixels, DataFrame)
    detection_pixels_df.type = detection_pixels_df.shape
    detection_pixels_df[:, ["r", "g", "b"]] = detection_pixels_df[:, ["r", "g", "b"]] .* 255.0 # convert to 0-255 range
    CSV.write(splitext(detection_pixels) |> x-> x[1]*"__color255_webui"*x[2], detection_pixels_df)

    n_clamped > 0 && @warn "$n_clamped detection rows fell outside their chapter's frames (audio is slightly longer than video) and were clamped to the first/last frame"
    # The overlay never reads the raw CSV: a single NaN row aborts it with InexactError.
    info = get_media_info(vidlist[1])
    video_stream = first(filter(s -> get(s, "codec_type", "") == "video", info["streams"]))
    valid = write_valid_detection_pixels(detection_pixels; imsize=(Int(video_stream["height"]), Int(video_stream["width"])))
    if !pixel_filter && valid.n_nonfinite > 0
        error("pixel_filter=false but $(valid.n_nonfinite) of $(valid.n_input) detection rows have a non-finite frame/px/py ($detection_pixels); the overlay would abort. Use pixel_filter=true (see $(valid.status_path)).")
    end
    overlay_pixels = valid.out_path
    @info "Valid detection pixels: $(valid.n_written)/$(valid.n_input) rows ($(valid.n_nonfinite) non-finite dropped, $(valid.n_outside_image) outside the image kept)"

    # Only process valid (non-nothing) results for audacity file combination
    valid_idx = findall(!isnothing, results)
    all_durations = get_duration.(audlist)
    cum_durations_before = cumsum([0.0; all_durations[1:end-1]])  # cumulative duration before each file

    dfs_impulse = [CSV.read(joinpath(res_dir2, results[i].res.res_impulsetrain.outfname * ".txt"), DataFrame; header=false) for i in valid_idx]
    dfs_tonal = [CSV.read(joinpath(res_dir2, results[i].res.res_tonalsegment.outfname * ".txt"), DataFrame; header=false) for i in valid_idx]
    dfs_impulsetrain = [CSV.read(joinpath(res_dir2, results[i].res.res_impulsetrain.outfname * "_train-only.txt"), DataFrame; header=false) for i in valid_idx]

    cum_index_impulse = 0; cum_index_tonal = 0; cum_index_impulsetrain = 0;
    for (j, i) in enumerate(valid_idx)
        t_offset = cum_durations_before[i]
        dfs_impulse[j][:, 1:2] .+= t_offset
        dfs_tonal[j][:, 1:2] .+= t_offset
        dfs_impulsetrain[j][:, 1:2] .+= t_offset

        dfs_impulse[j][:, 3] .+= cum_index_impulse
        dfs_tonal[j][:, 3] .+= cum_index_tonal
        dfs_impulsetrain[j][:, 3] .+= cum_index_impulsetrain

        cum_index_impulse += nrow(dfs_impulse[j])
        cum_index_tonal += nrow(dfs_tonal[j])
        cum_index_impulsetrain += nrow(dfs_impulsetrain[j])
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
    overlay_tmpdir = isnothing(tmpdir) ? joinpath(res_dir2, "tmp_overlay") : tmpdir
    if flag_overlayvideo
        GC.gc()
        try
            # If flag_spectro, launch spectrogram video creation in parallel with overlay
            spectro_task = nothing
            if flag_spectro
                info = get_media_info(vidlist[1])
                video_stream = first(filter(s -> get(s, "codec_type", "") == "video", info["streams"]))
                vid_width = get(video_stream, "width", 1280)
                vid_fps_val = something(get_fps(vidlist[1]), 30.0)
                spectro_out = splitext(output_vidname)[1] * "_spectro.mp4"
                spectro_task = @async create_spectro_video(audlist, spectro_out;
                    delays=max(0.0, Float64(delays)), vid_width=vid_width,
                    vid_fps=vid_fps_val, spectro_height=spectro_height)
            end

            out_vid_path = overlay_annotations_on_video(overlay_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaid.mkv";
                mode=:VideoIO, radius=overlay_radius, tmpdir=overlay_tmpdir,
                kwargs...)

            # Wait for spectrogram (if running) and collect its path
            spectro_vidpath = nothing
            if !isnothing(spectro_task)
                try
                    spectro_vidpath = fetch(spectro_task)
                catch err
                    @error "Spectrogram video creation failed" exception=(err, catch_backtrace())
                end
            end

            combine_vidau(out_vid_path, audlist; vidau_syncdiff=delays, MERGE_VID_AU_DYNAMIC_NORM=true, rx_vect=rx_vect,
                spectro_vidpath=spectro_vidpath, kwargs...)
        catch err
            @error "Failed to overlay boxes on video($output_vidname)" exception=(err, catch_backtrace())
        finally
            isnothing(tmpdir) && rm(overlay_tmpdir; recursive=true, force=true)   # only the directory this function chose
        end
    end
    flag_overlayimages &&
    (overlay_boxes_on_video_imageonly(overlay_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaidIMG";
        radius=(overlay_radius isa Number ? overlay_radius : OVERLAY_RADIUS));
    pic2vid(splitext(output_vidname)[1]*"_overlaidIMG", splitext(output_vidname)[1]*"_overlaidIMG.mp4"; auto_mode=true)
    )

    if flag_return
        return (;results, detection_pixels_df, delays, conf, output_vidname, vidlist, audlist, res_dir2, sync_model, overlay_pixels)
    end
end
