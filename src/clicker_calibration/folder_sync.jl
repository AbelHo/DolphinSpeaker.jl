# Drift-aware audio-video synchronization for one recording folder (T5).
# Reuses prep_wav_info, prep_boundaries, prep_video_info and prep_synchronize from
# prepare.jl; no sync algorithm lives here. Including this file performs no work.
#
# Convention:  t_video_global = scale * t_recorder + offset_s
#   t_recorder     : seconds on the concatenated `audlist` (sum of samples/fs, no BEXT gaps)
#   t_video_global : seconds on the concatenated chapters, defined by prep_synchronize as the
#                    sum of declared camera audio-track durations (codec padding excluded),
#                    sampled at the 2 kHz coarse decode rate for chapter start times.
# The existing constant-delay pipeline is the special case scale = 1, offset_s = delays
# (frame = (t + cum_duration + delays) * fps).
const FOLDER_SYNC_SCHEMA = "clicker-folder-sync-v1"
const FOLDER_SYNC_SYNC_KEYS = (:ffmpeg, :sync_fs, :anchor_window_s, :sync_min_score,
    :sync_min_ambiguity, :sync_residual_limit_s)
const FOLDER_SYNC_RECORDER_CHANNEL = 1   # must equal prep_synchronize's physical channel

"""Map recorder time (s on the concatenated `audlist`) to global video time (s on the
concatenated chapters): `scale*t + offset_s`. Broadcasts over arrays. `sync` is the
Dict returned by `fit_folder_sync_affine` (or its parsed JSON)."""
function video_time_from_recorder(t_rec_s, sync)
    scale = get(sync, "scale", nothing); offset = get(sync, "offset_s", nothing)
    (scale isa Real && offset isa Real && isfinite(scale) && isfinite(offset)) ||
        throw(ArgumentError("sync needs finite numeric \"scale\" and \"offset_s\" (status=$(get(sync, "status", missing)))"))
    t_rec_s isa AbstractArray ? scale .* t_rec_s .+ offset : scale * t_rec_s + offset
end

_fsync_file_record(p) = Dict{String,Any}("path"=>abspath(p),
    "size_bytes"=>isfile(p) ? filesize(p) : nothing,
    "mtime_utc"=>isfile(p) ? string(Dates.unix2datetime(mtime(p))) * "Z" : nothing)

function _fsync_recorder_timeline(infos, boundaries)
    Dict{String,Any}(
        "definition"=>"concatenation of audlist in the given order; each file contributes samples/fs seconds; BEXT time_reference is diagnostic only",
        "sync_physical_channel"=>FOLDER_SYNC_RECORDER_CHANNEL,
        "fs"=>infos[1]["fs"],
        "total_duration_s"=>sum(q["samples"] for q in infos)/infos[1]["fs"],
        "files"=>[Dict{String,Any}("path"=>q["path"],"start_s"=>q["start_sample"]/q["fs"],
            "start_sample"=>q["start_sample"],"samples"=>q["samples"],"fs"=>q["fs"],
            "channels"=>q["channels"],"duration_s"=>q["duration_s"]) for q in infos],
        "boundaries"=>boundaries)
end

function _fsync_video_timeline(videos)
    files = Dict{String,Any}[]
    for v in videos
        pts = v["pts"]
        push!(files, Dict{String,Any}("path"=>v["path"],"global_start_s"=>get(v,"global_start_s",nothing),
            "audio_duration_s"=>v["audio_duration_s"],"decoded_audio_duration_s"=>get(v,"decoded_audio_duration_s",nothing),
            "video_duration_s"=>v["duration_s"],"audio_start_s"=>v["audio_start_s"],"video_start_s"=>v["video_start_s"],
            "video_minus_audio_duration_s"=>v["duration_s"]-v["audio_duration_s"],
            "n_frames"=>length(pts),"first_pts_s"=>first(pts),"last_pts_s"=>last(pts),
            "width"=>v["width"],"height"=>v["height"]))
    end
    Dict{String,Any}(
        "definition"=>"concatenation of chapters (in the given order); global time = sum of declared camera audio-track durations, codec padding excluded; each chapter keeps its original local PTS (local = global - global_start_s + audio_start_s, see prep_match_frame)",
        "camera_rate_hz"=>48000,
        "total_duration_s"=>isempty(videos) ? 0.0 : sum(get(v,"decoded_audio_duration_s",v["audio_duration_s"]) for v in videos),
        "files"=>files)
end

"""Assemble the folder-sync Dict from an already-computed `prep_synchronize` /
`prep_fit_sync` result. Pure: no media access."""
function _fsync_assemble(infos, videos, boundaries, sync; vidlist=[v["path"] for v in videos],
        audlist=[q["path"] for q in infos], settings=Dict{String,Any}())
    prep_status = string(get(sync, "status", "missing"))
    has_fit = get(sync, "offset_s", nothing) isa Real && get(sync, "scale", nothing) isa Real
    status = prep_status == "valid" && has_fit ? "valid" : "insufficient_anchors"
    anchors = get(sync, "anchors", Any[])
    inliers = [a for a in anchors if get(a, "fit_inlier", false)]
    out = Dict{String,Any}(
        "schema"=>FOLDER_SYNC_SCHEMA, "status"=>status, "prep_status"=>prep_status,
        "convention"=>"t_video_global = scale * t_recorder + offset_s",
        "scale"=>get(sync,"scale",nothing), "offset_s"=>get(sync,"offset_s",nothing),
        "drift_ppm"=>get(sync,"drift_ppm",nothing),
        "rms_residual_s"=>get(sync,"rms_residual_s",nothing), "max_residual_s"=>get(sync,"max_residual_s",nothing),
        "n_anchors"=>length(inliers), "n_anchors_candidates"=>length(anchors),
        "anchor_min_time_s"=>get(sync,"anchor_min_time_s",nothing),
        "anchor_max_time_s"=>get(sync,"anchor_max_time_s",nothing),
        "clock_scale_uncertainty"=>get(sync,"clock_scale_uncertainty",nothing),
        "template_scale_bias_bound_s"=>get(sync,"template_scale_bias_bound_s",nothing),
        "recorder_sync_physical_channel"=>FOLDER_SYNC_RECORDER_CHANNEL,
        "video_timeline"=>_fsync_video_timeline(videos),
        "recorder_timeline"=>_fsync_recorder_timeline(infos, boundaries),
        "inputs"=>Dict{String,Any}("videos"=>[_fsync_file_record(p) for p in vidlist],
            "recorders"=>[_fsync_file_record(p) for p in audlist],"settings"=>settings),
        "segment_diagnostics"=>get(sync,"segment_diagnostics",nothing),
        "refinement_history"=>get(sync,"refinement_history",nothing),
        "anchors"=>anchors, "coarse_anchors"=>get(sync,"coarse_anchors",nothing),
        "method"=>get(sync,"method","prep_synchronize"),
        "created_utc"=>string(Dates.now(Dates.UTC)) * "Z")
    out
end

function _fsync_failed(vidlist, audlist, exception, bt, settings)
    log_failure("fit_folder_sync_affine", exception, bt)
    Dict{String,Any}("schema"=>FOLDER_SYNC_SCHEMA, "status"=>"failed",
        "error"=>sprint(showerror, exception),
        "convention"=>"t_video_global = scale * t_recorder + offset_s",
        "scale"=>nothing, "offset_s"=>nothing, "drift_ppm"=>nothing, "rms_residual_s"=>nothing,
        "max_residual_s"=>nothing, "n_anchors"=>0, "anchor_min_time_s"=>nothing, "anchor_max_time_s"=>nothing,
        "video_timeline"=>nothing, "recorder_timeline"=>nothing,
        "inputs"=>Dict{String,Any}("videos"=>[abspath(string(p)) for p in vidlist],
            "recorders"=>[abspath(string(p)) for p in audlist],"settings"=>settings),
        "created_utc"=>string(Dates.now(Dates.UTC)) * "Z")
end

"""
    fit_folder_sync_affine(vidlist, audlist; out_path=nothing, ffprobe="/usr/bin/ffprobe", kwargs...)

Fit the robust affine map `t_video_global = scale*t_recorder + offset_s` (drift = (scale-1)*1e6 ppm)
between the concatenated recorder files `audlist` (physical channel 1) and the concatenated camera
chapters `vidlist`, by waveform correlation (`prep_synchronize`). `kwargs` accepted:
$(join(string.(FOLDER_SYNC_SYNC_KEYS), ", ")). Never throws: insufficient/inconsistent anchors give
status "insufficient_anchors" (original reason in "prep_status"); internal failures give status
"failed" plus "error" and an entry in docs/error.log. Writes JSON to `out_path` when given.
"""
function fit_folder_sync_affine(vidlist::AbstractVector{<:AbstractString}, audlist::AbstractVector{<:AbstractString};
        out_path=nothing, ffprobe="/usr/bin/ffprobe", kwargs...)
    settings = Dict{String,Any}(string(k)=>v for (k,v) in kwargs); settings["ffprobe"] = ffprobe
    local result
    try
        unknown = [k for k in keys(kwargs) if !(k in FOLDER_SYNC_SYNC_KEYS)]
        isempty(unknown) || throw(ArgumentError("unknown keywords $(unknown); allowed: $(FOLDER_SYNC_SYNC_KEYS)"))
        isempty(vidlist) && throw(ArgumentError("empty vidlist"))
        isempty(audlist) && throw(ArgumentError("empty audlist"))
        for p in vcat(collect(vidlist), collect(audlist))
            isfile(p) || throw(ArgumentError("missing input file: $p"))
        end
        infos = prep_wav_info.(audlist)
        length(unique(q["fs"] for q in infos)) == 1 || throw(ErrorException("recorder files differ in sample rate: $(unique(q["fs"] for q in infos))"))
        all(q["channels"] >= FOLDER_SYNC_RECORDER_CHANNEL for q in infos) || throw(ErrorException("recorder WAV without channel $(FOLDER_SYNC_RECORDER_CHANNEL)"))
        boundaries = prep_boundaries(infos)     # sets start_sample on each info
        videos = [prep_video_info(p; ffprobe=ffprobe) for p in vidlist]
        syncargs = (k=>v for (k,v) in kwargs)
        sync = prep_synchronize(infos, videos; syncargs...)
        result = _fsync_assemble(infos, videos, boundaries, sync; vidlist=vidlist, audlist=audlist, settings=settings)
    catch exception
        result = _fsync_failed(vidlist, audlist, exception, catch_backtrace(), settings)
    end
    if out_path !== nothing
        try
            write_json(string(out_path), result)
            result["out_path"] = abspath(string(out_path))
        catch exception
            log_failure("fit_folder_sync_affine write $(out_path)", exception, catch_backtrace())
            result["status"] == "failed" || (result["write_error"] = sprint(showerror, exception))
        end
    end
    result
end
