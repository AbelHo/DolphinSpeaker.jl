# Compare TDOA methods on real dolphin clicks (no ground truth) and on semi-synthetic clicks (known truth).
# plan: reviews/plans/trace-the-function-run-contiguous-folder-eager-puppy.md
# usage: julia --project=. -t 16 src/compare_tdoa.jl [wav files...]
#    or: include("compare_tdoa.jl"); compare_tdoa([aufname]; res_dir=...)
include("run_example.jl")
set_device__ophk_acoustic_D3()
using Statistics, Printf, Random, Dates, CSV, DataFrames, Plots, FFTW

COMPARE_TDOA_FOLDER = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
COMPARE_TDOA_RESULTS = "/media/spin/anas2/data/calf/upload/results"

# each method: f(data, peaks, fs) -> (tdoa N×m, closure per click or nothing); windows as used in the pipeline
function tdoa_methods(rx)
    W = window_impulsive
    plain(f) = (d, p, fs) -> (f(d, p, fs), nothing)
    gcc_m(; win=W, kw...) = (d, p, fs) -> begin
        t, q = get_tdoa_gcc(d, p; window=win, ref_channel=ref_channel, fs=fs, rx_vect=rx, return_quality=true, kw...)
        (t, 3 .* q.pair_rms)   # 3 sensors: closure |τ12+τ23+τ31| = 3 x rms pair misfit
    end
    return [
        (name="minmax [Impulse, overlay]", group="existing", f=plain((d, p, fs) -> get_tdoa_minmax(d, p; window=W, ref_channel=ref_channel))),
        (name="xcorr MaxEnergyRef [Impulse2]", group="existing", f=plain((d, p, fs) -> get_tdoa_raw_MaxEnergyRefChannel(d, p; window=-130:130, ref_channel=ref_channel))),
        (name="envelope_filtered [Impulse3]", group="existing", f=plain((d, p, fs) -> get_tdoa_envelope_filtered(d, p; window=W, ref_channel=ref_channel))),
        (name="xcorr MaxPeakRef [Impulse4]", group="existing", f=plain((d, p, fs) -> get_tdoa_raw_MaxPeakRefChannel(d, p; window=W, ref_channel=ref_channel))),
        (name="xcorr MaxEnergyRef resample x10", group="existing", f=plain((d, p, fs) -> get_tdoa_raw_MaxEnergyRefChannel_resample(d, p; window=W, ref_channel=ref_channel))),
        (name="minmax2 (vote fixed)", group="fixed copy", f=plain((d, p, fs) -> get_tdoa_minmax2(d, p; window=W, ref_channel=ref_channel))),
        (name="envelope_filtered2 (fs fixed)", group="fixed copy", f=plain((d, p, fs) -> get_tdoa_envelope_filtered2(d, p; window=W, ref_channel=ref_channel, fs=fs))),
        (name="gcc none (ablation)", group="new", f=gcc_m(weighting=:none)),
        (name="gcc phat", group="new", f=gcc_m(weighting=:phat)),
        (name="gcc phat +guide", group="new", f=gcc_m(weighting=:phat, guide=:envelope)),
        (name="gcc phat 20kHz+", group="new", f=gcc_m(weighting=:phat, band=[20_000, Inf])),
        (name="gcc scot", group="new", f=gcc_m(weighting=:scot)),
        (name="gcc scot +guide", group="new", f=gcc_m(weighting=:scot, guide=:envelope)),
        (name="gcc ml", group="new", f=gcc_m(weighting=:ml)),
        (name="gcc ml +guide", group="new", f=gcc_m(weighting=:ml, guide=:envelope)),
        (name="srp phat", group="new", f=plain((d, p, fs) -> get_tdoa_srp(d, p; window=W, ref_channel=ref_channel, fs=fs, rx_vect=rx))),
        # tighter window: keeps every direct arrival (|τ| <= 22.4 samples) but fewer late echoes / overlapping clicks
        (name="gcc phat [window -30:30]", group="new", f=gcc_m(weighting=:phat, win=-30:30)),
        (name="srp phat [window -30:30]", group="new", f=plain((d, p, fs) -> get_tdoa_srp(d, p; window=-30:30, ref_channel=ref_channel, fs=fs, rx_vect=rx))),
    ]
end

# TDOAs (samples, relative to ref_channel) of a plane wave from (az, el) in degrees, same model as cost_tdoa2ang
function model_tdoa(az, el, rx, fs; c=default_c)
    t = vec(-[sind(az) * cosd(el) sind(el) cosd(az) * cosd(el)] * rx ./ c .* fs)
    return t .- t[ref_channel]
end

median_or_nan(x) = isempty(x) ? NaN : median(x)
dirvec_deg(az, el) = [sind(az) * cosd(el), sind(el), cosd(az) * cosd(el)]
angle_between(a1, e1, a2, e2) = acosd(clamp(dot(dirvec_deg(a1, e1), dirvec_deg(a2, e2)), -1, 1))

# az/el (deg) from both solvers and the unclamped direction cosines; NelderMead runs in threaded chunks
function solve_directions(T, rx, fs)
    uv = tdoa2uv(T, rx, fs)
    ls = rad2deg.(tdoa2dir_ls(T, rx, fs))
    nm = zeros(size(T, 1), 2)
    if !isempty(T)
        chunks = collect(Iterators.partition(axes(T, 1), cld(size(T, 1), 4Threads.nthreads())))
        Threads.@threads for ch in chunks
            nm[ch, :] = rad2deg.(tdoa2dir(T[ch, :], rx, fs))
        end
    end
    return (; uv, ls, nm)
end

# consecutive clicks of the same train: |Δaz|, |Δel| and the angle between their directions (deg)
function train_steps(ang, idx_trains)
    pairs = [(a, b) for tr in idx_trains for (a, b) in zip(tr[1:end-1], tr[2:end])]
    Δ = [abs.(ang[b, :] .- ang[a, :]) for (a, b) in pairs]
    return Δ, [angle_between(ang[a, 1], ang[a, 2], ang[b, 1], ang[b, 2]) for (a, b) in pairs], pairs
end

# x delayed by δ samples (FFT phase shift), zero-padded by `pad` on both sides
function fracdelay(x, δ; pad=64)
    xp = [zeros(pad); x; zeros(pad)]
    N = length(xp)
    return irfft(rfft(xp) .* exp.(-im * 2π .* rfftfreq(N) .* δ), N)
end

# semi-synthetic benchmark: real high-SNR clicks (ref channel, tapered) delayed by the exact TDOAs of known
# directions inside the camera FOV, added to real click-free 3-channel noise at a given window SNR
function synthetic_benchmark(d, fs, pk, ppeak, all_det, rx, methods; snrs=[0, 5, 10, 15, 20, 30],
    ntrials=300, ntemplates=40, seglen=2000, gap=2000, half=64, seed=1)

    rng = MersenneTwister(seed)
    W = window_impulsive
    n, m = size(d)
    top = pk[sortperm(ppeak; rev=true)][1:min(ntemplates, end)]
    templates = [d[p-half:p+half, ref_channel] .* tukey(2half + 1, 0.5) for p in top]
    det = sort(all_det)
    clickfree(q) = (i = searchsortedfirst(det, q - gap); i > length(det) || det[i] > q + seglen + gap)
    cand = filter(clickfree, rand(rng, gap+1:n-seglen-gap, 50ntrials))
    length(cand) < ntrials && error("only $(length(cand)) click-free noise segments found")
    noise_pos = cand[1:ntrials]
    En = mean(sum(abs2, d[q+seglen÷2 .+ W, ch]) for q in noise_pos for ch in 1:m)   # noise energy in one window
    dirs = [(rand(rng) * 70 - 35, rand(rng) * 50 - 25) for _ in 1:ntrials]
    tmpl = rand(rng, eachindex(templates), ntrials)

    rows = DataFrame()
    for snr in snrs
        data = zeros(seglen * ntrials, m)
        peaks = zeros(Int, ntrials)
        Ttrue = zeros(ntrials, m)
        for t in 1:ntrials
            base = (t - 1) * seglen
            data[base+1:base+seglen, :] = d[noise_pos[t]:noise_pos[t]+seglen-1, :]
            s = templates[tmpl[t]]
            s = s .* sqrt(10^(snr / 10) * En / sum(abs2, s[half+1 .+ W]))
            Ttrue[t, :] = model_tdoa(dirs[t]..., rx, fs)
            c0 = base + seglen ÷ 2
            peaks[t] = c0
            for ch in 1:m
                data[c0-half-64:c0+half+64, ch] .+= fracdelay(s, Ttrue[t, ch])
            end
        end
        for meth in methods
            T̂ = try
                Float64.(meth.f(data, peaks, fs)[1])
            catch err
                @error "synthetic: $(meth.name) failed" exception = err
                continue
            end
            T̂ .-= T̂[:, ref_channel]
            others = setdiff(1:m, ref_channel)
            terr = vec(T̂[:, others] .- Ttrue[:, others])
            ls = rad2deg.(tdoa2dir_ls(T̂, rx, fs))
            aerr = [angle_between(ls[t, 1], ls[t, 2], dirs[t]...) for t in 1:ntrials]
            push!(rows, (method=meth.name, group=meth.group, snr_db=snr,
                tdoa_rmse=sqrt(mean(abs2, terr)), tdoa_med_abs=median(abs.(terr)),
                angle_err_median=median(aerr), angle_err_p90=quantile(aerr, 0.9), outlier_pct=100 * mean(aerr .> 5)))
        end
    end
    return rows, En
end

# metrics on the strongest detections (by ppeak, independent of any TDOA method), which are probably true clicks;
# false detections and multi-arrival clicks make every method look inconsistent on the full set
function strength_table(per_click; fractions=[1.0, 0.5, 0.25, 0.1])
    rows = DataFrame()
    ppk = per_click.ppeak[per_click.method.==per_click.method[1]]
    for f in fractions
        idx = findall(ppk .>= quantile(ppk, 1 - f))
        for g in groupby(per_click, :method; sort=false)   # rows of each method are in click order
            steps = [angle_between(g.az_ls[a], g.el_ls[a], g.az_ls[b], g.el_ls[b])
                     for (a, b) in zip(idx[1:end-1], idx[2:end]) if g.train[a] == g.train[b] && g.train[a] > 0]
            rng = MersenneTwister(1)
            spread = median(angle_between(g.az_ls[a], g.el_ls[a], g.az_ls[b], g.el_ls[b])
                            for (a, b) in zip(rand(rng, idx, 20_000), rand(rng, idx, 20_000)))
            push!(rows, (top_pct=100f, method=String(g.method[1]), n=length(idx), feasible_pct=100 * mean(g.feasible[idx]),
                jitter_deg=median_or_nan(steps), spread_deg=spread, consistency_ratio=median_or_nan(steps) / spread))
        end
    end
    return rows
end

function print_strength_table(io, tab)
    println(io, "\nProbable true clicks: the strongest detections by ppeak (this selection does not depend on any TDOA method).")
    println(io, "  ratio = jitter / spread as above; integer-TDOA methods can show jitter 0 on strong clicks because consecutive")
    println(io, "  integer TDOAs repeat: that is quantization, not accuracy (see the semi-synthetic table).")
    tops = sort(unique(tab.top_pct); rev=true)
    @printf(io, "%-32s", "method")
    foreach(t -> @printf(io, " %9s", t == 100 ? "ratio all" : @sprintf("top%.0f%%", t)), tops)
    @printf(io, " %12s %14s\n", "jitter°@min", "feasible%@min")
    for g in groupby(tab, :method; sort=false)
        @printf(io, "%-32s", g.method[1])
        foreach(t -> @printf(io, " %9.2f", g.consistency_ratio[g.top_pct.==t][1]), tops)
        last = g[g.top_pct.==minimum(tops), :]
        @printf(io, " %12.2f %14.1f\n", last.jitter_deg[1], last.feasible_pct[1])
    end
end

# video frame offset of this audio file in the combined video, as run_contiguous_folders computes it:
# frame = round(t * fps) + 1 + (duration of the earlier audio files in the folder + A/V sync delay) * fps
# the sync delay comes from sync_delay_s, else sync_dir/sync_delay.csv (written by run_contiguous_folders), else is computed
function webui_frame_offset(aufname; sync_delay_s=nothing, sync_dir=nothing, res_dir=nothing)
    vidlist, audlist = split_vid_au(dirname(aufname); readdir_func=readdir_all)
    fps = get_fps(vidlist[1])
    i = findfirst(a -> basename(a) == basename(aufname), audlist)
    isnothing(i) && error("$aufname not among the audio files of its folder")
    cum = sum(get_duration.(audlist[1:i-1]); init=0.0)
    if isnothing(sync_delay_s) && !isnothing(sync_dir)
        sync_delay_s = CSV.read(joinpath(sync_dir, "sync_delay.csv"), DataFrame).delay_s[1]
    end
    if isnothing(sync_delay_s)
        @info "computing the A/V sync delay of $(basename(dirname(aufname)))"
        sync_delay_s, conf = find_vid_vs_audio_syncdiff_timesegment(vidlist[1], audlist[1]; flag_return_conf=true)
        isnothing(res_dir) || write(joinpath(res_dir, "sync_delay.csv"), "foldername,delay_s,confidence\n$(basename(dirname(aufname))),$(sync_delay_s),$(conf)\n")
    end
    @info @sprintf("webui frames: fps %.3f, earlier audio %.1f s, sync delay %.3f s", fps, cum, sync_delay_s)
    return fps, (cum + sync_delay_s) * fps
end

# one detection_pixels__color255_webui__<method>.csv per method, same columns and click colouring as the pipeline's
# detection_pixels__color255_webui.csv (impulse rows), so each method can be overlaid on the combined video
function write_webui_csvs(res_dir, tag, d, fs, pk, t_s, per_click, aufname; solver=:nm, radius=OVERLAY_RADIUS, webui_kwargs...)
    fps, offset = webui_frame_offset(aufname; res_dir=res_dir, webui_kwargs...)
    frames = round.(Int, t_s .* fps) .+ 1 .+ offset
    clips = extract_clips(d, pk, 100)   # colour and alpha per click, as in process_detections
    coloring = sig2rgb.(clips; fs=fs, rgb_bands=[[1000, 30000], [30000, 60000], [60000, 96000]])
    rgb = 255 .* hcat((coloring .|> x -> x[1][1:3, argmax(x[1][end, :])])...)'
    pow = @pipe extrema.(clips; dims=1) .|> reduce.(-, _) .|> minimum .|> log10(-_)
    alphas = stretch(pow, 0.4, 0.8)
    az, el = solver == :nm ? (:az_nm, :el_nm) : (:az_ls, :el_ls)
    for g in groupby(per_click, :method; sort=false)
        px = angle2pxd(hcat(g[!, az], g[!, el]), fov_angle, imsize)
        name = replace(String(g.method[1]), r"[^A-Za-z0-9]+" => "-") |> x -> strip(x, '-')
        CSV.write(joinpath(res_dir, "detection_pixels__color255_webui__$(tag)__$(name).csv"),
            DataFrame(frame=frames, px=px[:, 1], py=px[:, 2], radius=radius, r=rgb[:, 1], g=rgb[:, 2], b=rgb[:, 3], a=alphas,
                shape="circle", type="circle"))
    end
end

# export_webui=true also writes one webui overlay CSV per method (see write_webui_csvs; webui_kwargs: sync_delay_s, sync_dir)
function compare_tdoa_file(aufname; res_dir, edge_margin=600, run_synthetic=true, synthetic_kwargs=(;),
    export_webui=false, webui_solver=:nm, webui_kwargs=(;))
    tag = splitext(basename(aufname))[1]
    rx = rx_vect[:, get_relevant_channels(rx_vect)]
    methods = tdoa_methods(rx)
    W = window_impulsive

    @info "reading $aufname"
    data, fs, _, _, _ = readAudio(aufname)
    res_impulse = detect_impulse((aufname, data, fs), nothing; ref_channel=ref_channel, return_datafilt=true)
    res_train = detect_impulsetrain(res_impulse, nothing)
    data = nothing
    GC.gc()
    d = res_impulse.data_filt
    n = size(d, 1)

    keep = findall(p -> edge_margin < p <= n - edge_margin, res_train.pind_good)
    pk = res_train.pind_good[keep]
    ppeak = res_train.ppeak[keep]
    t_s = res_train.pind_good_inS[keep]
    idx_of = Dict(p => i for (i, p) in enumerate(pk))
    idx_trains = [[idx_of[p] for p in tr if haskey(idx_of, p)] for tr in res_train.train_impulse_list]
    train_id = zeros(Int, length(pk))
    for (k, tr) in enumerate(idx_trains), i in tr
        train_id[i] = k
    end
    @info @sprintf("%s: fs %.0f Hz, threshold %.3f, %d clicks in %d trains (%d kept away from the file edges)",
        tag, fs, res_impulse.threshold, length(res_train.pind_good), length(res_train.train_impulse_list), length(pk))

    # weaker click of each consecutive in-train pair -> ppeak quartile
    qedges = quantile(ppeak, [0.25, 0.5, 0.75])
    quart(x) = 1 + count(x .> qedges)
    # random click pairs across the file: how far apart directions are when they need not agree
    rng = MersenneTwister(1)
    rand_pairs = collect(zip(rand(rng, eachindex(pk), 20_000), rand(rng, eachindex(pk), 20_000)))

    per_click = DataFrame()
    summary = DataFrame()
    by_quartile = DataFrame()
    ls_angles = Dict{String,Matrix{Float64}}()
    closures = Dict{String,Vector{Float64}}()
    for meth in methods
        @info "real data: $(meth.name)"
        local T, clo
        try
            meth.f(d, pk[1:min(20, end)], fs)   # compile before timing
            elapsed = @elapsed (T, clo) = meth.f(d, pk, fs)
            T = Float64.(T)
            T .-= T[:, ref_channel]
            sol = solve_directions(T, rx, fs)
            feasible = vec(hypot.(sol.uv[:, 1], sol.uv[:, 2]) .<= 1)
            infov = (abs.(sol.ls[:, 1]) .<= fov_angle[1] / 2) .& (abs.(sol.ls[:, 2]) .<= fov_angle[2] / 2)
            nm_differs = vec(maximum(abs.(sol.nm .- sol.ls); dims=2) .> 0.1)
            steps, Δgc, step_pairs = train_steps(sol.ls, idx_trains)
            Δaz = getindex.(steps, 1)
            Δel = getindex.(steps, 2)
            spread = median(angle_between(sol.ls[a, 1], sol.ls[a, 2], sol.ls[b, 1], sol.ls[b, 2]) for (a, b) in rand_pairs)
            ls_angles[meth.name] = sol.ls
            isnothing(clo) || (closures[meth.name] = clo)

            append!(per_click, DataFrame(file=tag, click=eachindex(pk), pind=pk, t_s=t_s, ppeak=ppeak, train=train_id,
                method=meth.name, group=meth.group,
                tdoa1=T[:, 1], tdoa2=T[:, 2], tdoa3=T[:, 3], u=sol.uv[:, 1], v=sol.uv[:, 2], feasible=feasible,
                az_ls=sol.ls[:, 1], el_ls=sol.ls[:, 2], az_nm=sol.nm[:, 1], el_nm=sol.nm[:, 2],
                closure=isnothing(clo) ? fill(NaN, length(pk)) : clo))
            push!(summary, (method=meth.name, group=meth.group, n=length(pk),
                feasible_pct=100 * mean(feasible), in_fov_pct=100 * mean(infov),
                jitter_deg=median_or_nan(Δgc), spread_deg=spread, consistency_ratio=median_or_nan(Δgc) / spread,
                jitter_az=median_or_nan(Δaz), jitter_el=median_or_nan(Δel),
                closure_median=isnothing(clo) ? NaN : median(clo), closure_p90=isnothing(clo) ? NaN : quantile(clo, 0.9),
                closure_over3_pct=isnothing(clo) ? NaN : 100 * mean(clo .> 3),
                nm_differs_pct=100 * mean(nm_differs), ms_per_click=1e3 * elapsed / length(pk)))
            weak = [quart(min(ppeak[a], ppeak[b])) for (a, b) in step_pairs]
            for q in 1:4
                sel = weak .== q
                clicks_q = quart.(ppeak) .== q
                push!(by_quartile, (method=meth.name, ppeak_quartile=q, n_pairs=count(sel),
                    jitter_deg=median_or_nan(Δgc[sel]), jitter_az=median_or_nan(Δaz[sel]), jitter_el=median_or_nan(Δel[sel]),
                    feasible_pct=100 * mean(feasible[clicks_q]), closure_median=isnothing(clo) ? NaN : median(clo[clicks_q])))
            end
        catch err
            @error "real data: $(meth.name) failed" exception = (err, catch_backtrace())
        end
    end

    CSV.write(joinpath(res_dir, "per_click_$tag.csv"), per_click)
    CSV.write(joinpath(res_dir, "summary_$tag.csv"), summary)
    CSV.write(joinpath(res_dir, "by_ppeak_quartile_$tag.csv"), by_quartile)
    strength = strength_table(per_click)
    CSV.write(joinpath(res_dir, "by_strength_$tag.csv"), strength)
    export_webui && write_webui_csvs(res_dir, tag, d, fs, pk, t_s, per_click, aufname; solver=webui_solver, webui_kwargs...)

    # pairwise agreement: median angle (deg) between the LS directions of two methods
    names_ok = [meth.name for meth in methods if haskey(ls_angles, meth.name)]
    agree = [median(angle_between(ls_angles[a][i, 1], ls_angles[a][i, 2], ls_angles[b][i, 1], ls_angles[b][i, 2]) for i in eachindex(pk))
             for a in names_ok, b in names_ok]
    agreement = DataFrame(agree, names_ok)
    insertcols!(agreement, 1, :method => names_ok)

    # real-click window SNR, for reading the synthetic results at realistic levels
    syn, En = run_synthetic ?
        synthetic_benchmark(d, fs, pk, ppeak, res_impulse.pind_good, rx, methods; synthetic_kwargs...) : (DataFrame(), NaN)
    snr_real = run_synthetic ? [10log10(max(sum(abs2, d[p .+ W, ref_channel]) / En - 1, 1e-3)) for p in pk] : Float64[]

    CSV.write(joinpath(res_dir, "agreement_$tag.csv"), agreement)
    run_synthetic && CSV.write(joinpath(res_dir, "synthetic_$tag.csv"), syn)

    report = sprint() do io
        println(io, "TDOA method comparison: $tag")
        @printf(io, "fs %.0f Hz | %d clicks in trains (%d analysed) | window %s | ref_channel %d\n",
            fs, length(res_train.pind_good), length(pk), string(W), ref_channel)
        if run_synthetic
            @printf(io, "real-click window SNR (dB): p10 %.1f | median %.1f | p90 %.1f\n", quantile(snr_real, [0.1, 0.5, 0.9])...)
        end
        println(io, "\nReal clicks (no ground truth):")
        println(io, "  jitter = median angle between consecutive clicks of a train; spread = the same for random click pairs;")
        println(io, "  ratio = jitter / spread (lower is better; 1 = consecutive clicks agree no better than random ones;")
        println(io, "  a method that barely moves gets a low jitter, the ratio does not reward that);")
        println(io, "  feasible = TDOAs a plane wave can produce; closure = |τ12+τ23+τ31| before reconciling (gcc only, samples);")
        println(io, "  NM≠LS = NelderMead az/el differ from the closed-form solver (mirror or infeasible).")
        @printf(io, "%-32s %-10s %9s %7s %8s %8s %6s %15s %8s %7s %8s\n", "method", "group", "feasible%", "inFOV%",
            "jitter°", "spread°", "ratio", "closure med/p90", "clo>3%", "NM≠LS%", "ms/click")
        fmtnan(x, f) = isnan(x) ? "-" : Printf.format(Printf.Format(f), x)
        for r in eachrow(summary)
            @printf(io, "%-32s %-10s %9.1f %7.1f %8.1f %8.1f %6.2f %7s/%-7s %8s %7.1f %8.3f\n", r.method, r.group, r.feasible_pct,
                r.in_fov_pct, r.jitter_deg, r.spread_deg, r.consistency_ratio, fmtnan(r.closure_median, "%.2f"),
                fmtnan(r.closure_p90, "%.2f"), fmtnan(r.closure_over3_pct, "%.1f"), r.nm_differs_pct, r.ms_per_click)
        end
        print_strength_table(io, strength)
        if run_synthetic
            println(io, "\nSemi-synthetic (known truth): median angle error (deg) by window SNR; outliers = error > 5°")
            snrs = sort(unique(syn.snr_db))
            @printf(io, "%-32s", "method")
            foreach(s -> @printf(io, " %7s", "$(s)dB"), snrs)
            @printf(io, " %12s\n", "outl%@10dB")
            for g in groupby(syn, :method; sort=false)
                @printf(io, "%-32s", g.method[1])
                foreach(s -> @printf(io, " %7.2f", g.angle_err_median[g.snr_db.==s][1]), snrs)
                o = g.outlier_pct[g.snr_db.==10]
                @printf(io, " %12s\n", isempty(o) ? "-" : @sprintf("%.1f", o[1]))
            end
        end
    end
    print(report)
    write(joinpath(res_dir, "report_$tag.txt"), report)

    plot_comparison(res_dir, tag, per_click, summary, syn, closures)
    return (; summary, by_quartile, agreement, syn, per_click)
end

#~ plots (static PNGs; colours from the dataviz reference palette)
viz_ink = (surface="#fcfcfb", primary="#0b0b0b", secondary="#52514e", muted="#898781", grid="#e1e0d9", axis="#c3c2b7")
group_colors = Dict("new" => "#2a78d6", "existing" => "#eb6834", "fixed copy" => "#1baf7a")
seq_blue = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plot_style() = (background_color=viz_ink.surface, foreground_color_axis=viz_ink.axis, foreground_color_border=viz_ink.axis,
    gridcolor=viz_ink.grid, gridalpha=1, tickfontcolor=viz_ink.secondary, guidefontcolor=viz_ink.secondary,
    titlefontcolor=viz_ink.primary, titlefontsize=9, tickfontsize=7, guidefontsize=8, legendfontsize=8,
    foreground_color_legend=nothing, background_color_legend=viz_ink.surface)

function small_multiples(xs, ys, titles, ylab, fname; ylims=:auto, kind=:scatter, xlab="time (s)", bins=range(0, 30; length=61))
    nc = 4
    nr = cld(length(titles), nc)
    plt = plot(; layout=(nr, nc), size=(1800, 330nr), legend=false, left_margin=6Plots.mm, bottom_margin=4Plots.mm, plot_style()...)
    for k in eachindex(titles)
        bottom = k > length(titles) - nc   # lowest panel of its column
        if kind == :scatter
            scatter!(plt[k], xs[k], ys[k]; markersize=1.2, markerstrokewidth=0, color=group_colors["new"], alpha=0.35,
                title=titles[k], ylims=ylims, xlabel=bottom ? xlab : "", ylabel=(k - 1) % nc == 0 ? ylab : "")
        else
            histogram!(plt[k], ys[k]; bins=bins, color=group_colors["new"], linecolor=viz_ink.surface, yformatter=:plain,
                title=titles[k], xlabel=bottom ? xlab : "", ylabel=(k - 1) % nc == 0 ? ylab : "")
        end
    end
    for k in length(titles)+1:nr*nc   # unused grid cells
        plot!(plt[k]; framestyle=:none, grid=false, xticks=false, yticks=false)
    end
    savefig(plt, fname)
end

function bars_by_group(summary, col, title, fmt)
    pos = nrow(summary):-1:1
    # horizontal bars: set both limits explicitly (with only xlims, Plots also stretches the category axis)
    plt = plot(; title=title, legend=:outertop, legend_columns=3, yticks=(collect(pos), summary.method),
        xlims=(0, 1.18maximum(filter(!isnan, summary[!, col]))), ylims=(0.4, nrow(summary) + 0.6), plot_style()...)
    for g in ["new", "existing", "fixed copy"]
        idx = findall(summary.group .== g)
        isempty(idx) || bar!(plt, collect(pos)[idx], summary[idx, col]; orientation=:h, bar_width=0.7, color=group_colors[g],
            linecolor=viz_ink.surface, label=g)
    end
    annotate!(plt, [(v, p, text(" " * fmt(v), 7, viz_ink.secondary, :left)) for (v, p) in zip(summary[!, col], pos)])
    return plt
end

function plot_comparison(res_dir, tag, per_click, summary, syn, closures)
    gr()
    names = summary.method
    sub = [per_click[per_click.method.==nm, :] for nm in names]
    small_multiples([s.t_s for s in sub], [s.az_ls for s in sub], names, "azimuth (°)", joinpath(res_dir, "az_vs_time_$tag.png"); ylims=(-90, 90))
    small_multiples([s.t_s for s in sub], [s.el_ls for s in sub], names, "elevation (°)", joinpath(res_dir, "el_vs_time_$tag.png"); ylims=(-90, 90))

    p1 = bars_by_group(summary, :jitter_deg, "Within-train jitter: median angle between consecutive clicks (°)", v -> @sprintf("%.1f", v))
    p2 = bars_by_group(summary, :consistency_ratio, "Jitter / spread of random click pairs (lower is better, 1 = none)", v -> @sprintf("%.2f", v))
    p3 = bars_by_group(summary, :feasible_pct, "Feasible TDOAs (%), higher is better", v -> @sprintf("%.1f", v))
    savefig(plot(p1, p2, p3; layout=(1, 3), size=(2300, 700), left_margin=8Plots.mm, plot_style()...), joinpath(res_dir, "metrics_$tag.png"))

    cnames = [nm for nm in names if haskey(closures, nm)]
    ctitles = [@sprintf("%s: %.0f%% of clicks > 3 samples", nm, 100 * mean(closures[nm] .> 3)) for nm in cnames]
    isempty(cnames) || small_multiples(fill(nothing, length(cnames)), [closures[nm] for nm in cnames], ctitles, "clicks",
        joinpath(res_dir, "closure_$tag.png"); kind=:hist, xlab="closure |τ12+τ23+τ31| (samples, > 30 not shown)")

    if !isempty(syn)
        snrs = sort(unique(syn.snr_db))
        mnames = unique(syn.method)
        M = [syn.angle_err_median[(syn.method.==nm).&(syn.snr_db.==s)][1] for nm in mnames, s in snrs]
        L = log10.(clamp.(M, 0.01, 100))
        # every cell carries its value, so no colorbar; colour is on a log scale
        hm = Plots.heatmap(string.(snrs) .* " dB", mnames, L; color=Plots.cgrad(seq_blue), clims=(log10(0.03), log10(60)), yflip=true,
            colorbar=false, yticks=(0.5:1:length(mnames)-0.5, mnames),
            title="Semi-synthetic clicks: median angle error (°) by window SNR, lower is better (colour on a log scale)",
            size=(1100, 700), left_margin=8Plots.mm, plot_style()...)
        annotate!(hm, [(j - 0.5, i - 0.5, text(@sprintf("%.2f", M[i, j]), 7, L[i, j] > 0.3 ? "#ffffff" : viz_ink.primary))
                       for i in axes(M, 1) for j in axes(M, 2)])
        savefig(hm, joinpath(res_dir, "synthetic_$tag.png"))
    end
end

# rebuild the figures from a finished run's CSVs (no recomputation)
function replot_comparison(res_dir, tag)
    per_click = CSV.read(joinpath(res_dir, "per_click_$tag.csv"), DataFrame)
    summary = CSV.read(joinpath(res_dir, "summary_$tag.csv"), DataFrame)
    synf = joinpath(res_dir, "synthetic_$tag.csv")
    syn = isfile(synf) ? CSV.read(synf, DataFrame) : DataFrame()
    closures = Dict(String(g.method[1]) => Vector{Float64}(g.closure) for g in groupby(per_click, :method) if !all(isnan, g.closure))
    plot_comparison(res_dir, tag, per_click, summary, syn, closures)
end

function compare_tdoa(aufnames=[joinpath(COMPARE_TDOA_FOLDER, "251006_001_0001.WAV")];
    res_dir=joinpath(COMPARE_TDOA_RESULTS, "compare_tdoa_" * Dates.format(now(), "yyyymmdd_HHMMSS")), kwargs...)
    mkpath(res_dir)
    @info "results -> $res_dir ($(Threads.nthreads()) threads)"
    return [compare_tdoa_file(f; res_dir=res_dir, kwargs...) for f in aufnames]
end

if abspath(PROGRAM_FILE) == @__FILE__
    haskey(ENV, "GKSwstype") || (ENV["GKSwstype"] = "100")   # headless GR when run as a script
    webui = "--webui" in ARGS   # also write per-method webui overlay CSVs
    files = filter(!=("--webui"), ARGS)
    compare_tdoa(isempty(files) ? [joinpath(COMPARE_TDOA_FOLDER, "251006_001_0001.WAV")] : files; export_webui=webui)
end
