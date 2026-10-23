# End-to-end benchmark on simulated reverberant pool recordings (src/pool_sim.jl): the pipeline's detector
# (detect_impulse + detect_impulsetrain), SIC candidates, tracker and selection rules, scored against the truth table.
# Also computes the real-data proxies on the simulated output, to check whether they rank the rules like the truth.
# needs: results/data/click_features_New_device.csv (round 1; picks the click templates)
# usage: julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/bench_sequences.jl calibrate
#        julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/bench_sequences.jl run [n_scenarios] [A0]
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "pool_sim.jl"))
include(joinpath(@__DIR__, "bench_common.jl"))
using Printf, Random
mode = length(ARGS) >= 1 ? ARGS[1] : "run"
out_dir = joinpath(@__DIR__, "..", "results"); mkpath(joinpath(out_dir, "data"))
fs = BENCH_FS; rx = BENCH_RX
const DUR = 4.0                                   # seconds per scenario
const REAL_RATIO = (1.09, 1.79, 7.0)              # real in-train ppeak / threshold, p10 / p50 / p90 (New_device)

templates, _ = bench_material(joinpath(out_dir, "data", "click_features_New_device.csv"))
noise = read_filtered_block(BENCH_WAV, round(Int, 53.5fs):round(Int, 97.5fs); fs=fs)     # 44 s without detections
# detection threshold of the real file the noise comes from: its auto-threshold (9 x the 75% envelope quantile over the
# whole file) is 1.8x the noise-only value, because clicks and reverberation raise the quantile
thr_noise = impulsive_autothreshold_median_ratio * quantile(abs.(hilbert(noise[:, 3])), 0.75)
thr = 490.376461598132
@info @sprintf("noise %.1f s; detection threshold %.1f (file 0001; noise-only auto-threshold would be %.1f)", size(noise, 1) / fs, thr, thr_noise)

function simulate_and_detect(s, A0)
    rng = Xoshiro(10_000 + s)
    pool, array_pos = random_pool_setup(rng)
    nsrc = rand(rng, 1:2)
    sim = simulate_scenario(rng, templates, noise; fs=fs, dur=DUR, nsrc=nsrc, pool=pool, array_pos=array_pos, A0=A0, rx=rx)
    res = detect_impulse(("sim$s", sim.data, fs), nothing; threshold=thr, return_datafilt=true)
    rt = detect_impulsetrain(res, nothing)
    return (; sim, res, rt, pool, array_pos, nsrc)
end

if mode == "calibrate"
    for A0 in (6e3, 9e3, 1.2e4, 1.8e4, 2.7e4)
        r = Float64[]; fol = Float64[]; nd = 0
        for s in 1:6
            d = simulate_and_detect(s, A0)
            append!(r, d.rt.ppeak ./ thr); nd += length(d.rt.pind_good)
            append!(fol, follower_flags(d.rt.pind_good, fs).follower)
        end
        @printf("A0 %7.0f: in-train detections %5d, ppeak/threshold p10/p50/p90 %.2f/%.2f/%.2f (real %.2f/%.2f/%.2f), followers %.0f%% (real 51%%)\n",
            A0, nd, quantile(r, [0.1, 0.5, 0.9])..., REAL_RATIO..., 100mean(fol))
    end
    exit()
end

n_scen = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 30
A0 = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e4
const WINDOWS = Dict("sic" => window_impulsive, "sicw" => -110:50)
const TRACKERS = [   # (name, kwargs) - the first is round 2's default
    ("default", (;)), ("gate 4", (; gate_deg=4.0)), ("gate 9", (; gate_deg=9.0)), ("birth 0.1", (; birth_min_score=0.1)),
    ("birth 0.3", (; birth_min_score=0.3)), ("delete 0.5 s", (; t_delete=0.5)), ("q 10", (; sigma_q_deg=10.0)), ("q 40", (; sigma_q_deg=40.0))]
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))

det_rows = DataFrame(scenario=Int[], window=String[], det=Int[], t=Float64[], gap_prev=Float64[], follower=Bool[], matched=Bool[],
    src=Int[], click=Int[], order=Int[], true_az=Float64[], true_el=Float64[], nfound=Int[], err_rank1=Float64[], err_best=Float64[])
rule_rows = DataFrame(scenario=Int[], window=String[], tracker=String[], rule=String[], det=Int[], err=Float64[], echo_label=Bool[], az=Float64[], el=Float64[])
sim_stats = DataFrame(scenario=Int[], nsrc=Int[], L=Float64[], W=Float64[], D=Float64[], depth=Float64[], clicks=Int[], detections=Int[],
    in_train=Int[], matched_pct=Float64[], direct_pct=Float64[], followers_pct=Float64[], ratio_p50=Float64[])
for s in 1:n_scen
    d = simulate_and_detect(s, A0)
    P = d.rt.pind_good; n = length(P)
    n < 20 && continue
    arr = d.sim.arrivals; clk = d.sim.clicks
    # truth: the strongest arrival within ±12 samples of each detection
    order_ = zeros(Int, n); srcm = zeros(Int, n); clkm = zeros(Int, n); tz = fill(NaN, n); te = fill(NaN, n); matched = falses(n)
    srt = sortperm(arr.sample); samples = arr.sample[srt]
    for i in 1:n
        lo = searchsortedfirst(samples, P[i] - 12); hi = searchsortedlast(samples, P[i] + 12)
        lo > hi && continue
        j = srt[lo:hi][argmax(abs.(arr.amp[srt[lo:hi]]))]
        matched[i] = true; order_[i] = arr.order[j]; srcm[i] = arr.src[j]; clkm[i] = arr.click[j]
        row = findfirst((clk.src .== srcm[i]) .& (clk.click .== clkm[i]))
        tz[i] = clk.direct_az[row]; te[i] = clk.direct_el[row]
    end
    ff = follower_flags(P, fs); tsec = (P .- 1) ./ fs
    push!(sim_stats, (s, d.nsrc, d.pool.L, d.pool.W, d.pool.D, d.array_pos[3], nrow(clk), length(d.res.pind_good), n, 100mean(matched),
        100mean(order_[matched] .== 0), 100mean(ff.follower), median(d.rt.ppeak ./ thr)))
    for (wname, W) in WINDOWS
        c = get_srp_candidates_sic(d.res.data_filt, P; window=W, fs=fs, rx_vect=rx, min_rel=0.7)
        dirs = [[(isnan(c.uv[i, k, 1]) ? (NaN, NaN) : rad2deg.(uv2azel(c.uv[i, k, 1], c.uv[i, k, 2]))) for k in 1:3] for i in 1:n]
        errk(i, k) = isnan(dirs[i][k][1]) || !matched[i] ? NaN : angd(vec3(dirs[i][k]...), vec3(tz[i], te[i]))
        for i in 1:n
            push!(det_rows, (s, wname, i, tsec[i], ff.gap_prev[i], ff.follower[i], matched[i], srcm[i], clkm[i], order_[i], tz[i], te[i],
                c.nfound[i], errk(i, 1), minimum(k -> (e = errk(i, k); isnan(e) ? Inf : e), 1:3)))
        end
        earliest = [(ks = findall(k -> !isnan(c.tarr[i, k]), 1:3); isempty(ks) ? 1 : ks[argmin(c.tarr[i, ks])]) for i in 1:n]
        for (tname, tkw) in TRACKERS
            T = track_candidates(tsec, c.uv, c.rscore; tkw...)
            conf = Dict(zip(T.tracks.id, T.tracks.confirmed)); hits = Dict(zip(T.tracks.id, T.tracks.hits))
            St = select_candidates(T.track, c.tarr, c.score, conf, hits)
            longest = [(ks = [k for k in 1:3 if T.track[i, k] > 0 && conf[T.track[i, k]]]; isempty(ks) ? 1 : ks[argmax([hits[T.track[i, k]] for k in ks])]) for i in 1:n]
            sels = Dict("rank 1" => ones(Int, n), "earliest" => earliest, "track lag" => St.sel, "longest track" => longest)
            tname == "default" || (sels = Dict("track lag" => St.sel, "longest track" => longest))   # tracker variants only affect these
            for (rname, sel) in sels
                az = [dirs[i][sel[i]][1] for i in 1:n]; el = [dirs[i][sel[i]][2] for i in 1:n]
                lab = rname == "track lag" ? St.status .!= "direct" : falses(n)
                for (rn, a, e) in [(rname, az, el), ("$rname + followers inherit", copy(az), copy(el))]
                    if endswith(rn, "inherit")
                        for i in 2:n; ff.follower[i] && (a[i] = a[i-1]; e[i] = e[i-1]); end
                    end
                    for i in 1:n
                        push!(rule_rows, (s, wname, tname, rn, i, matched[i] ? angd(vec3(a[i], e[i]), vec3(tz[i], te[i])) : NaN, lab[i], a[i], e[i]))
                    end
                end
            end
        end
    end
    @info @sprintf("scenario %d/%d: %d sources, pool %.0fx%.0fx%.1f m, %d clicks, %d in-train detections (%.0f%% direct-path, %.0f%% followers)",
        s, n_scen, d.nsrc, d.pool.L, d.pool.W, d.pool.D, nrow(clk), n, last(sim_stats).direct_pct, last(sim_stats).followers_pct)
end
CSV.write(joinpath(out_dir, "data", "bench_sequences_detections.csv"), det_rows)
CSV.write(joinpath(out_dir, "data", "bench_sequences_rules.csv"), rule_rows)
CSV.write(joinpath(out_dir, "bench_sequences_scenarios.csv"), sim_stats)

# ---- summary: accuracy per rule vs the proxies ---------------------------------------------------------------------
dr = det_rows
summ = DataFrame(window=String[], tracker=String[], rule=String[], acc_all=Float64[], acc_direct_det=Float64[], acc_echo_det=Float64[],
    med_err=Float64[], proxy_burst=Float64[], proxy_2_10=Float64[], proxy_10_30=Float64[], in_frame=Float64[])
for g in groupby(rule_rows, [:window, :tracker, :rule])
    w, tn, rn = g.window[1], g.tracker[1], g.rule[1]
    info = dr[dr.window .== w, :]
    key = Dict((r.scenario, r.det) => r for r in eachrow(info))
    m = [key[(r.scenario, r.det)] for r in eachrow(g)]
    ok = [x.matched for x in m]
    err = g.err[ok]; ord = [x.order for x in m][ok]
    # proxies: consecutive isolated detections (both gap >= 1.6 ms) in the same scenario
    function proxy(lo, hi)
        a = Float64[]
        for i in 2:nrow(g)
            (g.scenario[i] == g.scenario[i-1] && !m[i].follower && !m[i-1].follower && lo <= m[i].gap_prev < hi) || continue
            push!(a, angd(vec3(g.az[i], g.el[i]), vec3(g.az[i-1], g.el[i-1])))
        end
        return isempty(a) ? NaN : median(a)
    end
    push!(summ, (w, tn, rn, 100mean(err .< 3), 100mean(err[ord .== 0] .< 3), 100mean(err[ord .> 0] .< 3), median(err),
        proxy(1.6, 2.2), proxy(2.2, 10.0), proxy(10.0, 30.0), 100mean((abs.(g.az) .<= 40) .& (abs.(g.el) .<= 27.23))))
end
for w in keys(WINDOWS)   # oracle: the best of the candidates
    info = dr[(dr.window .== w) .& dr.matched, :]
    push!(summ, (w, "-", "oracle (best candidate)", 100mean(info.err_best .< 3), 100mean(info.err_best[info.order .== 0] .< 3),
        100mean(info.err_best[info.order .> 0] .< 3), median(info.err_best), NaN, NaN, NaN, NaN))
end
sort!(summ, [:window, :tracker, order(:acc_all; rev=true)])
CSV.write(joinpath(out_dir, "bench_sequences_summary.csv"), summ)

# echo label of the track rule: precision / recall against detections that are image arrivals
io = IOBuffer()
@printf(io, "%d scenarios x %.0f s, A0 %.0f, threshold %.1f; %d in-train detections, %.1f%% matched to an arrival, %.1f%% of them direct-path\n",
    nrow(sim_stats), DUR, A0, thr, sum(sim_stats.in_train), mean(sim_stats.matched_pct), mean(sim_stats.direct_pct))
@printf(io, "simulated vs real: followers %.0f%% (real 51%%); median ppeak/threshold %.2f (real %.2f); SIC >= 2 candidates %.0f%% (real 78%%)\n\n",
    mean(sim_stats.followers_pct), median(sim_stats.ratio_p50), REAL_RATIO[2], 100mean(dr.nfound[dr.window .== "sic"] .>= 2))
println(io, "selected direction within 3 deg of the emitting click's direct path (%), and the real-data proxies on the simulation")
@printf(io, "%-5s %-13s %-36s %7s %7s %7s %7s | %7s %7s %7s %6s\n", "win", "tracker", "rule", "all", "direct", "echo", "med err", "burst", "2-10ms", "10-30ms", "frame")
for r in eachrow(summ)
    (r.tracker in ("default", "-") || r.rule in ("track lag", "track lag + followers inherit")) || continue
    @printf(io, "%-5s %-13s %-36s %7.1f %7.1f %7.1f %7.1f | %7.1f %7.1f %7.1f %6.1f\n", r.window, r.tracker, r.rule, r.acc_all, r.acc_direct_det,
        r.acc_echo_det, r.med_err, r.proxy_burst, r.proxy_2_10, r.proxy_10_30, r.in_frame)
end
for w in keys(WINDOWS)
    g = rule_rows[(rule_rows.window .== w) .& (rule_rows.tracker .== "default") .& (rule_rows.rule .== "track lag"), :]
    info = dr[dr.window .== w, :]; key = Dict((r.scenario, r.det) => r for r in eachrow(info))
    m = [key[(r.scenario, r.det)] for r in eachrow(g)]; ok = [x.matched for x in m]
    truth_echo = [x.order > 0 for x in m][ok]; lab = g.echo_label[ok]
    @printf(io, "\n%s track rule echo label: precision %.1f%%, recall %.1f%% (image detections: %.1f%%)", w, 100mean(truth_echo[lab]), 100mean(lab[truth_echo]), 100mean(truth_echo))
end
# do the proxies rank the default-tracker rules like the truth? (Spearman over rules, per window)
spearman(a, b) = (ra = invperm(sortperm(a)); rb = invperm(sortperm(b)); cor(Float64.(ra), Float64.(rb)))
println(io, "\n\nrank correlation over the rules (default tracker) between accuracy and each proxy (proxy: lower change / higher frame share = better):")
for w in keys(WINDOWS)
    g = summ[(summ.window .== w) .& (summ.tracker .== "default"), :]
    @printf(io, "  %-5s burst %+.2f, 2-10 ms %+.2f, 10-30 ms %+.2f, in frame %+.2f (n rules %d)\n", w, spearman(g.acc_all, -g.proxy_burst),
        spearman(g.acc_all, -g.proxy_2_10), spearman(g.acc_all, -g.proxy_10_30), spearman(g.acc_all, g.in_frame), nrow(g))
end
body = String(take!(io))
open(joinpath(out_dir, "bench_sequences.md"), "w") do f; println(f, "# End-to-end benchmark on simulated pool recordings\n\n```\n", body, "```"); end
print(body)
