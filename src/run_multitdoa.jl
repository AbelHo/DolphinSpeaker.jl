# run_contiguous_folders with one webui overlay CSV per TDOA method (Impulse1-6), for comparing methods in the viewer
# plan: reviews/plans/multitdoa-webui-csvs_2026-10-22.md, background: improvements/tdoa_localization_2026-10-21/README.md
# usage: include("run_example.jl"); set_device__ophk_acoustic_D3(); include("run_multitdoa.jl")
#        run_contiguous_folders_multitdoa(folder; res_dir=..., closure_max=1.0)
# existing functions are untouched: process_detections gives Impulse1-4, Impulse5/6 are computed here
using CSV, DataFrames, DelimitedFiles, Statistics, Printf, FFMPEG

# Impulse1-4 are the pipeline's methods, Impulse5/6 the new ones (get_tdoa_gcc / get_tdoa_srp + tdoa2dir_ls);
# "g" = gated twin: clicks flagged as probably wrong are dropped (ungated files draw them with flag_shape)
const MULTITDOA_VARIANTS = ["Impulse1_minmax", "Impulse2_xcorrMaxEnergy", "Impulse3_envelopeFiltered", "Impulse4_xcorrMaxPeak",
    "Impulse5_gccPhatGuided", "Impulse6_srpPhat", "Impulse5g_gccPhatGuided-gated", "Impulse6g_srpPhat-gated"]
const MULTITDOA_COLUMNS = ["frame", "px", "py", "radius", "r", "g", "b", "a", "shape"]

# overlay rows [px py radius r g b a shape] of one method; the shape of the flagged clicks becomes flag_shape
# (in the webui CSV `type` is a copy of `shape`)
function impulse_variant_rows(px, colour_cols, flagged; flag_shape="triangle")
    rows = Matrix{Any}(undef, size(px, 1), size(px, 2) + size(colour_cols, 2))
    rows[:, 1:size(px, 2)] = px
    rows[:, size(px, 2)+1:end] = colour_cols
    rows[flagged, end] .= flag_shape
    return rows
end

# clicks that are probably wrong: GCC closure |τ12+τ23+τ31| above closure_max samples, or TDOAs no plane wave can
# produce (u²+v² > 1, see tdoa2uv). Returns (flagged, infeasible); closure=nothing skips the closure test
function flag_clicks(tdoas, closure, rx, fs; closure_max=1.0, ref_channel=ref_channel)
    uv = tdoa2uv(tdoas, rx, fs; ref_channel=ref_channel)
    infeasible = vec(.!(hypot.(uv[:, 1], uv[:, 2]) .<= 1 + 1e-9))   # also catches NaN
    flagged = isnothing(closure) ? infeasible : infeasible .| .!(closure .<= closure_max)
    return flagged, infeasible
end

# process_detections (Impulse1-4, unchanged) plus Impulse5 = GCC-PHAT guided by the envelope, Impulse6 = SRP-PHAT,
# both with the closed-form solver; returns process_detections' fields, then `variants`:
# Vector of (; key, frames, rows, flagged) with rows = [px py radius r g b a shape]
function process_detections_multitdoa(aufname, vidfname; res_dir=nothing, closure_max=1.0, flag_shape="triangle", kwargs...)
    out = process_detections(aufname, vidfname; res_dir=res_dir, flag_return=true, flag_color=true, kwargs...)
    res = out.res
    d, pk, fs = res.res_impulse.data_filt, res.res_impulsetrain.pind_good, res.fs
    rx = rx_vect[:, get_relevant_channels(rx_vect)]
    size(rx, 2) == 3 || error("closure = 3 x pair rms assumes 3 hydrophones, got $(size(rx, 2))")

    frames, p1 = out.pixel_related_impulsive
    colour = p1[:, 3:end]   # radius r g b a shape: the same click set in every method, so the colours are shared
    @info "Impulse5/6: $(length(pk)) clicks of $(basename(aufname))"
    T5, q = get_tdoa_gcc(d, pk; window=window_impulsive, ref_channel=ref_channel, fs=fs, rx_vect=rx, weighting=:phat, guide=:envelope, return_quality=true)
    closure = 3 .* q.pair_rms
    T6 = get_tdoa_srp(d, pk; window=window_impulsive, ref_channel=ref_channel, fs=fs, rx_vect=rx)
    flag5, _ = flag_clicks(T5, closure, rx, fs; closure_max=closure_max)
    flag6, _ = flag_clicks(T6, closure, rx, fs; closure_max=closure_max)   # SRP is feasible by construction; closure from the GCC pass
    px5 = angle2px(tdoa2dir_ls(T5, rx, fs), fov_angle)
    px6 = angle2px(tdoa2dir_ls(T6, rx, fs), fov_angle)

    none_flagged = falses(length(frames))
    px_of(p) = p[:, 1:2]
    ungated = [
        (; key=MULTITDOA_VARIANTS[1], frames=frames, rows=impulse_variant_rows(px_of(p1), colour, none_flagged), flagged=none_flagged),
        (; key=MULTITDOA_VARIANTS[2], frames=out.pixel_related_impulsive2[1], rows=impulse_variant_rows(px_of(out.pixel_related_impulsive2[2]), colour, none_flagged), flagged=none_flagged),
        (; key=MULTITDOA_VARIANTS[3], frames=out.pixel_related_impulsive3[1], rows=impulse_variant_rows(px_of(out.pixel_related_impulsive3[2]), colour, none_flagged), flagged=none_flagged),
        (; key=MULTITDOA_VARIANTS[4], frames=out.pixel_related_impulsive4[1], rows=impulse_variant_rows(px_of(out.pixel_related_impulsive4[2]), colour, none_flagged), flagged=none_flagged),
        (; key=MULTITDOA_VARIANTS[5], frames=frames, rows=impulse_variant_rows(px5, colour, flag5; flag_shape=flag_shape), flagged=flag5),
        (; key=MULTITDOA_VARIANTS[6], frames=frames, rows=impulse_variant_rows(px6, colour, flag6; flag_shape=flag_shape), flagged=flag6),
    ]
    gate(v, key) = (; key=key, frames=v.frames[.!v.flagged], rows=v.rows[.!v.flagged, :], flagged=falses(count(.!v.flagged)))
    variants = [ungated..., gate(ungated[5], MULTITDOA_VARIANTS[7]), gate(ungated[6], MULTITDOA_VARIANTS[8])]

    # the detection-filtered data is only needed above; drop it like process_detections does when flag_return=false
    res = (; Base.structdiff(res, NamedTuple{(:res_impulse,)})..., res_impulse=Base.structdiff(res.res_impulse, NamedTuple{(:data_filt,)}))
    return merge(Base.structdiff(out, NamedTuple{(:res,)}), (; res, variants))
end

# frame, px, py, ... rows of one detection type over all audio files, shifted by the cumulative duration and the A/V sync delay
function _webui_rows(get_frames_rows, results, audlist, vidlist, delays)
    fps = get_fps(vidlist[1])
    cum_duration = 0.0
    parts = Matrix{Any}[]
    for (ind, res) in enumerate(results)
        if !isnothing(res)
            frames, rows = get_frames_rows(res)
            push!(parts, [frames .+ ((cum_duration + delays) * fps) rows])
        end
        cum_duration += get_duration(audlist[ind])
    end
    return isempty(parts) ? Matrix{Any}(undef, 0, length(MULTITDOA_COLUMNS)) : reduce(vcat, parts)
end

# same file format as detection_pixels__color255_webui.csv: r, g, b on 0-255 and type = shape
function _write_webui_csv(fname, mat)
    df = DataFrame([Symbol(c) => (c == "shape" ? String.(mat[:, k]) : Float64.(mat[:, k])) for (k, c) in enumerate(MULTITDOA_COLUMNS)])
    df.type = df.shape
    df[:, ["r", "g", "b"]] = df[:, ["r", "g", "b"]] .* 255.0
    CSV.write(fname, df)
    return df
end

# run_contiguous_folders without the video rendering: besides the standard detection_pixels.csv and
# detection_pixels__color255_webui.csv it writes detection_pixels__color255_webui__<variant>.csv for every entry of
# MULTITDOA_VARIANTS (that method's impulse rows + the tonal rows, so each is a drop-in for the standard file)
# and multitdoa_summary.csv (clicks, flagged, inside the video frame per variant)
function run_contiguous_folders_multitdoa(folname; res_dir="", detection_types=1:2, closure_max=1.0, flag_shape="triangle", max_parallel_files=1, kwargs...)
    @info "Processing folder: $folname ............. (impulse methods: $(join(MULTITDOA_VARIANTS, ", ")))"
    delays, conf, output_vidname, vidlist, audlist, res_dir2 = run_analysis_split_vidau(folname; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true, kwargs...)
    results = Vector{Any}(undef, length(audlist))
    # a 2 GB WAV needs roughly 15-20 GB while it is processed, so only max_parallel_files run at once
    # (the kernel OOM-killed the run with 3 in parallel while other programs held ~40 GB)
    gate = Base.Semaphore(max_parallel_files)
    Threads.@threads for i in eachindex(audlist)
        Base.acquire(gate)
        try
            results[i] = process_detections_multitdoa(audlist[i], vidlist[1]; res_dir=res_dir2, closure_max=closure_max, flag_shape=flag_shape)
        catch err
            @error "process_detections_multitdoa failed for $(audlist[i])" exception=(err, catch_backtrace())
            results[i] = nothing
        finally
            GC.gc()
            Base.release(gate)
        end
    end

    # standard files, written exactly as run_contiguous_folders does
    detection_pixels = joinpath(res_dir2, "detection_pixels.csv")
    open(detection_pixels, "w") do io
        writedlm(io, ["frame" "px" "py" "radius" "r" "g" "b" "a" "shape"], ',')
        for detection_type in detection_types
            cum_duration = 0.0
            for (ind, res) in enumerate(results)
                if !isnothing(res)
                    writedlm(io, [res[detection_type][1] .+ ((cum_duration + delays) * get_fps(vidlist[1])) res[detection_type][2]], ',')
                end
                cum_duration += get_duration(audlist[ind])
            end
        end
    end
    detection_pixels_df = CSV.read(detection_pixels, DataFrame)
    detection_pixels_df.type = detection_pixels_df.shape
    detection_pixels_df[:, ["r", "g", "b"]] = detection_pixels_df[:, ["r", "g", "b"]] .* 255.0
    CSV.write(splitext(detection_pixels) |> x -> x[1] * "__color255_webui" * x[2], detection_pixels_df)

    # one webui file per method
    summary = DataFrame(variant=String[], clicks=Int[], flagged=Int[], inside_frame=Int[], inside_frame_pct=Float64[])
    for key in MULTITDOA_VARIANTS
        mats = [dt == 1 ? _webui_rows(res -> let v = res.variants[findfirst(v -> v.key == key, res.variants)]; (v.frames, v.rows) end, results, audlist, vidlist, delays) :
                _webui_rows(res -> res[dt], results, audlist, vidlist, delays) for dt in detection_types]
        df = _write_webui_csv(joinpath(res_dir2, "detection_pixels__color255_webui__$(key).csv"), reduce(vcat, mats))
        imp = df[1:size(mats[1], 1), :]   # impulse rows come first (detection_types starts with 1)
        inside = count(r -> 0 <= r.px <= imsize[2] && 0 <= r.py <= imsize[1], eachrow(imp))
        push!(summary, (key, nrow(imp), count(==(flag_shape), imp.shape), inside, 100 * inside / max(nrow(imp), 1)))
    end
    CSV.write(joinpath(res_dir2, "multitdoa_summary.csv"), summary)
    show(stdout, MIME("text/plain"), summary); println()

    # combined audacity labels, as run_contiguous_folders
    valid_idx = findall(!isnothing, results)
    all_durations = get_duration.(audlist)
    cum_durations_before = cumsum([0.0; all_durations[1:end-1]])
    dfs_impulse = [CSV.read(joinpath(res_dir2, results[i].res.res_impulsetrain.outfname * ".txt"), DataFrame; header=false) for i in valid_idx]
    dfs_tonal = [CSV.read(joinpath(res_dir2, results[i].res.res_tonalsegment.outfname * ".txt"), DataFrame; header=false) for i in valid_idx]
    dfs_impulsetrain = [CSV.read(joinpath(res_dir2, results[i].res.res_impulsetrain.outfname * "_train-only.txt"), DataFrame; header=false) for i in valid_idx]
    cum_index_impulse = 0; cum_index_tonal = 0; cum_index_impulsetrain = 0
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

    # the combined video is written by an @async ffmpeg started in run_analysis_split_vidau: wait until it is readable
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
    @info "Combined video ready: $output_vidname"
    return results
end
