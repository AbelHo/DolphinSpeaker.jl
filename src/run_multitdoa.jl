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

# ---- memory-adaptive scheduling of the per-file work -------------------------------------------------------------
# Probes use Julia / the libuv bundled with it, so they work on Linux, macOS and Windows. Only the exact peak
# measurement (kernel VmHWM, resettable) is Linux-specific; elsewhere Sys.maxrss gives an upper bound (see run_memory_adaptive)
function _proc_bytes(file, key)
    for line in eachline(file)
        startswith(line, key) && return 1024 * parse(Int, split(line)[2])
    end
    error("$key not found in $file")
end
function rss_bytes()                # resident memory of this process (libuv: /proc, task_info or the Windows working set)
    r = Ref{Csize_t}(0)
    ccall(:uv_resident_set_memory, Cint, (Ptr{Csize_t},), r) == 0 || return Int(Sys.maxrss())
    return Int(r[])
end
function memory_limit_bytes()       # container (cgroup) limit, or typemax(Int) if there is none
    try
        lim = ccall(:uv_get_constrained_memory, UInt64, ())
        return 0 < lim < Sys.total_memory() ? Int(lim) : typemax(Int)
    catch
        return typemax(Int)
    end
end
function available_memory_bytes()   # what the system can still give to programs (libuv: MemAvailable on Linux)
    free = Int(Sys.free_memory())
    lim = memory_limit_bytes()
    return lim == typemax(Int) ? free : min(free, lim - rss_bytes())
end
function peak_rss_bytes()           # peak resident memory: Linux VmHWM (resettable), else Sys.maxrss (since process start)
    try
        Sys.islinux() && return _proc_bytes("/proc/self/status", "VmHWM:")
    catch
    end
    return Int(Sys.maxrss())
end
function reset_peak_rss!()          # VmHWM := current RSS; false if not Linux or not allowed
    Sys.islinux() || return false
    try
        write("/proc/self/clear_refs", "5")
        return true
    catch
        return false
    end
end
function release_memory!()          # full GC, then hand freed heap back to the OS (glibc)
    GC.gc()
    if Sys.islinux()
        try
            ccall(:malloc_trim, Cint, (Cint,), 0)
        catch
        end
    end
    return nothing
end
default_memory_margin_bytes() = max(4 * 2^30, Sys.total_memory() ÷ 10)
_spawn_item(f, k) = Threads.@spawn f(k)   # k is an argument: a task spawned from a loop variable would see its later values

const MEMORY_PROBES = (; available=available_memory_bytes, rss=rss_bytes, peak=peak_rss_bytes, reset_peak=reset_peak_rss!, release=release_memory!)

"""
    run_memory_adaptive(f, n; weights=ones(n), max_parallel=nothing, adaptive=true,
        margin_bytes=default_memory_margin_bytes(), safety=1.25, poll_s=1.0, probes=MEMORY_PROBES)

Call `f(i)` for `i in 1:n` (it stores its own result) with as many items in flight as the memory allows.
`weights[i]` is the relative size of item i (here the audio duration); memory use is assumed proportional to it.

- `adaptive=true`: the largest item runs alone first and its peak resident memory is measured (Linux: kernel VmHWM,
  reset before it starts; elsewhere Sys.maxrss, which can only overestimate it, so the result stays safe but may be
  too cautious). Per unit of weight that, times `safety`, is the estimate. Another item then starts only if
  `sum(estimates of the running items) + estimate(new) <= capacity`, with
  `capacity = (rss - rss_idle) + available - margin_bytes`: what this process already holds beyond its idle size plus
  what the system can still give, less the margin. Capacity does not change while a running item ramps up its memory
  (rss rises exactly as available falls), and it follows other programs coming and going. One item always runs.
- `adaptive=false`: no measurement, only the cap.
- `max_parallel`: upper bound on items in flight; `nothing` or `Inf` = no bound.
`probes` can be replaced to simulate memory (see test/test_multitdoa.jl). Returns `(; per_weight, max_concurrent, order)`.
"""
function run_memory_adaptive(f, n::Integer; weights=ones(n), max_parallel=nothing, adaptive=true,
    margin_bytes=default_memory_margin_bytes(), safety=1.25, poll_s=1.0, probes=MEMORY_PROBES)

    cap = max_parallel isa Integer ? Int(max_parallel) : (isnothing(max_parallel) || isinf(max_parallel)) ? typemax(Int) : floor(Int, max_parallel)
    cap = max(cap, 1)
    gb(x) = x / 2^30
    order = sortperm(weights; rev=true)   # the largest first: it sets the estimate, the others are no larger
    queue = collect(order)
    per_weight = 0.0
    rss_idle = 0
    max_concurrent = 0

    if adaptive && n > 0
        k1 = popfirst!(queue)
        probes.reset_peak()
        rss_before = probes.rss()
        f(k1)
        per_weight = safety * max(probes.peak() - rss_before, 0) / weights[k1]
        probes.release()
        rss_idle = probes.rss()
        max_concurrent = 1
        @info @sprintf("memory-adaptive: item %d (weight %.0f) peaked at %.1f GB -> estimate %.3f GB per unit (safety x%.2f); margin %.1f GB, available now %.1f GB",
            k1, weights[k1], gb(max(probes.peak() - rss_before, 0)), gb(per_weight), safety, gb(margin_bytes), gb(probes.available()))
    end

    est(i) = per_weight * weights[i]
    capacity() = (probes.rss() - rss_idle) + probes.available() - margin_bytes
    running = Dict{Int,Task}()
    blocked = 0
    while !isempty(queue) || !isempty(running)
        for (i, t) in collect(running)
            istaskdone(t) || continue
            delete!(running, i)
            try
                fetch(t)
            catch err
                @error "item $i failed" exception=(err, catch_backtrace())
            end
            probes.release()
        end
        while !isempty(queue)
            k = first(queue)
            fits = isempty(running) || (length(running) < cap &&
                (!adaptive || sum(est, keys(running); init=0.0) + est(k) <= capacity()))
            if !fits
                if blocked != k
                    blocked = k
                    @info @sprintf("memory-adaptive: item %d waits (needs %.1f GB more, capacity %.1f GB, %d running)",
                        k, gb(est(k)), gb(adaptive ? capacity() : 0), length(running))
                end
                break
            end
            popfirst!(queue)
            running[k] = _spawn_item(f, k)
            max_concurrent = max(max_concurrent, length(running))
            @info @sprintf("memory-adaptive: started item %d (estimate %.1f GB), %d running", k, gb(est(k)), length(running))
        end
        (isempty(queue) && isempty(running)) || sleep(poll_s)
    end
    return (; per_weight, max_concurrent, order)
end

# run_contiguous_folders without the video rendering: besides the standard detection_pixels.csv and
# detection_pixels__color255_webui.csv it writes detection_pixels__color255_webui__<variant>.csv for every entry of
# MULTITDOA_VARIANTS (that method's impulse rows + the tonal rows, so each is a drop-in for the standard file)
# and multitdoa_summary.csv (clicks, flagged, inside the video frame per variant)
function run_contiguous_folders_multitdoa(folname; res_dir="", detection_types=1:2, closure_max=1.0, flag_shape="triangle",
    max_parallel_files=nothing, adaptive_memory=true, memory_margin_gb=nothing, memory_safety=1.25, kwargs...)
    @info "Processing folder: $folname ............. (impulse methods: $(join(MULTITDOA_VARIANTS, ", ")))"
    delays, conf, output_vidname, vidlist, audlist, res_dir2 = run_analysis_split_vidau(folname; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true, kwargs...)
    results = Vector{Any}(undef, length(audlist))
    # files are processed in parallel as far as the RAM allows (run_memory_adaptive): the longest one runs alone first
    # and its measured peak memory sets how many of the others may start, leaving memory_margin_gb free
    # (all 3 files at once were OOM-killed while other programs held ~40 GB). max_parallel_files=nothing/Inf: no cap
    process_one(i) = try
        results[i] = process_detections_multitdoa(audlist[i], vidlist[1]; res_dir=res_dir2, closure_max=closure_max, flag_shape=flag_shape)
    catch err
        @error "process_detections_multitdoa failed for $(audlist[i])" exception=(err, catch_backtrace())
        results[i] = nothing
    end
    margin = isnothing(memory_margin_gb) ? default_memory_margin_bytes() : round(Int, memory_margin_gb * 2^30)
    run_memory_adaptive(process_one, length(audlist); weights=Float64.(get_duration.(audlist)), max_parallel=max_parallel_files,
        adaptive=adaptive_memory, margin_bytes=margin, safety=memory_safety)

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
