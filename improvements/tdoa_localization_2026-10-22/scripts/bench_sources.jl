# Round-4 end-to-end benchmark on the calibrated (realistic) pool simulation: per-detection direction with no
# inheritance / inherit-all followers / smart inheritance (separate_sources), and separation of n dolphins
# (n given, n estimated, k-means baseline), scored against the truth. Scenarios with 2 dolphins (main) and 3 (stress).
# needs: results/data/click_features_New_device.csv (round 1; click templates)
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/bench_sources.jl [n_scen_2] [n_scen_3]
#   realism via environment variables: A0 (6000), SCATTER (0.7), REVERB (0.1)
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "pool_sim.jl")); include(joinpath(src_dir, "source_separation.jl"))
include(joinpath(@__DIR__, "bench_common.jl"))
using Printf, Random
n2 = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 30
n3 = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 15
A0 = parse(Float64, get(ENV, "A0", "6000")); SCATTER = parse(Float64, get(ENV, "SCATTER", "0.7")); REVERB = parse(Float64, get(ENV, "REVERB", "0.1"))
out_dir = joinpath(@__DIR__, "..", "results"); mkpath(joinpath(out_dir, "data"))
fs = BENCH_FS; rx = BENCH_RX; DUR = 4.0; W = -110:50
templates, _ = bench_material(joinpath(out_dir, "data", "click_features_New_device.csv"))
noise = read_filtered_block(BENCH_WAV, round(Int, 53.5fs):round(Int, 97.5fs); fs=fs)
thr = 490.376461598132
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))

det = DataFrame(set=Int[], scenario=Int[], i=Int[], t=Float64[], follower=Bool[], matched=Bool[], src=Int[], order=Int[],
    other_follower=Bool[], method=String[], err=Float64[], role=String[])
sep = DataFrame(set=Int[], scenario=Int[], method=String[], precision=Float64[], recall=Float64[], f1=Float64[], ari=Float64[],
    n_true=Int[], n_found=Int[], unassigned_pct=Float64[])
proxy = DataFrame(set=Int[], scenario=Int[], rule=String[], burst=Float64[], p2_10=Float64[], p10_30=Float64[], acc=Float64[])
stats = DataFrame(set=Int[], scenario=Int[], detections=Int[], followers_pct=Float64[], ge2_pct=Float64[], closure_le1_pct=Float64[],
    direct_det_pct=Float64[], echo_det_pct=Float64[], diffuse_det_pct=Float64[])
# per-detection inputs of the separation (longest-track directions) and truth, for tune_sources.jl
inputs = DataFrame(set=Int[], scenario=Int[], t=Float64[], az=Float64[], el=Float64[], follower=Bool[], matched=Bool[], src=Int[],
    order=Int[], true_az=Float64[], true_el=Float64[], other_follower=Bool[])
for (set, nscen) in ((2, n2), (3, n3)), s in 1:nscen
    rng = Xoshiro(30_000 + 1000set + s)
    pool, array_pos = random_pool_setup(rng)
    tpl = templates[randperm(rng, length(templates))[1:10]]
    sim = simulate_scenario(rng, tpl, noise; fs=fs, dur=DUR, nsrc=set, pool=pool, array_pos=array_pos, A0=A0, rx=rx,
        hyd_scatter=SCATTER, reverb_level=REVERB, offaxis_lowpass=true)
    res = detect_impulse(("sim$s", sim.data, fs), nothing; threshold=thr, return_datafilt=true)
    rt = detect_impulsetrain(res, nothing)
    P = rt.pind_good; n = length(P)
    n < 30 && continue
    tsec = (P .- 1) ./ fs
    # truth: the strongest arrival in [p-35, p+12] (scattering taps delay the envelope peak by up to 30 samples)
    arr = sim.arrivals; clk = sim.clicks; srt = sortperm(arr.sample); smp = arr.sample[srt]
    matched = falses(n); src = zeros(Int, n); ord = fill(-2, n); tz = fill(NaN, n); te = fill(NaN, n)
    for i in 1:n
        lo = searchsortedfirst(smp, P[i] - 35); hi = searchsortedlast(smp, P[i] + 12)
        lo > hi && continue
        j = srt[lo:hi][argmax(abs.(arr.amp[srt[lo:hi]]))]
        matched[i] = true; src[i] = arr.src[j]; ord[i] = arr.order[j]
        r = findfirst((clk.src .== src[i]) .& (clk.click .== arr.click[j]))
        tz[i] = clk.direct_az[r]; te[i] = clk.direct_el[r]
    end
    ff = follower_flags(P, fs)
    # followers that are the other dolphin's click: matched to a direct path of a source other than the previous detection's
    other = [i > 1 && ff.follower[i] && ord[i] == 0 && src[i] != src[i-1] for i in 1:n]
    c = get_srp_candidates_sic(res.data_filt, P; window=W, fs=fs, rx_vect=rx, min_rel=0.7)
    T = track_candidates(tsec, c.uv, c.rscore)
    conf = Dict(zip(T.tracks.id, T.tracks.confirmed)); hits = Dict(zip(T.tracks.id, T.tracks.hits))
    St = select_candidates(T.track, c.tarr, c.score, conf, hits)
    longest = [(ks = [k for k in 1:3 if T.track[i, k] > 0 && conf[T.track[i, k]]]; isempty(ks) ? 1 : ks[argmax([hits[T.track[i, k]] for k in ks])]) for i in 1:n]
    dir(k, i) = isnan(c.uv[i, k, 1]) ? (NaN, NaN) : rad2deg.(uv2azel(c.uv[i, k, 1], c.uv[i, k, 2]))
    push!(stats, (set, s, n, 100mean(ff.follower), 100mean(c.nfound .>= 2), 100mean(c.closure .<= 1), 100mean(ord[matched] .== 0),
        100mean(ord[matched] .> 0), 100mean(ord[matched] .== -1)))
    err(az, el) = [matched[i] ? angd(vec3(az[i], el[i]), vec3(tz[i], te[i])) : NaN for i in 1:n]
    for (rname, sel) in (("rank 1", ones(Int, n)), ("earliest", [(ks = findall(k -> !isnan(c.tarr[i, k]), 1:3); ks[argmin(c.tarr[i, ks])]) for i in 1:n]),
        ("track lag", St.sel), ("longest track", longest))
        az = [dir(sel[i], i)[1] for i in 1:n]; el = [dir(sel[i], i)[2] for i in 1:n]
        variants = [("$rname", az, el, fill("-", n))]
        if rname in ("rank 1", "longest track")
            ai, ei = copy(az), copy(el)
            for i in 2:n; ff.follower[i] && (ai[i] = ai[i-1]; ei[i] = ei[i-1]); end
            push!(variants, ("$rname + inherit all", ai, ei, fill("-", n)))
            for nn in (set, nothing)
                ss = separate_sources(tsec, az, el; follower=ff.follower, n=nn)
                nn === nothing || push!(variants, ("$rname + smart inherit", ss.az, ss.el, ss.role))
                if rname == "longest track"   # source separation on the longest-track directions
                    sc = pairwise_scores(src, ss.label, tsec; max_dt=0.2)
                    nf = isnothing(nn) ? count(ss.sources.confirmed) : length(unique(filter(>(0), ss.label)))
                    push!(sep, (set, s, isnothing(nn) ? "tracker, n estimated" : "tracker, n given", sc.precision, sc.recall, sc.f1,
                        adjusted_rand(src, ss.label), length(unique(src[matched])), nf, 100mean(ss.label .== 0)))
                end
            end
            if rname == "longest track"
                for i in 1:n; push!(inputs, (set, s, tsec[i], az[i], el[i], ff.follower[i], matched[i], src[i], ord[i], tz[i], te[i], other[i])); end
                km = kmeans_sources(tsec, az, el, set; follower=ff.follower)
                sc = pairwise_scores(src, km, tsec; max_dt=0.2)
                push!(sep, (set, s, "k-means baseline (k = n)", sc.precision, sc.recall, sc.f1, adjusted_rand(src, km), length(unique(src[matched])), set, 0.0))
            end
        end
        for (vname, a, e, role) in variants
            ev = err(a, e)
            for i in 1:n
                push!(det, (set, s, i, tsec[i], ff.follower[i], matched[i], src[i], ord[i], other[i], vname, ev[i], role[i]))
            end
        end
        # proxies on the simulation (consecutive isolated detections), per rule
        iso = .!ff.follower
        pr(lo, hi) = (a = [angd(vec3(az[i], el[i]), vec3(az[i-1], el[i-1])) for i in 2:n if iso[i] && iso[i-1] && lo <= ff.gap_prev[i] < hi]; isempty(a) ? NaN : median(a))
        push!(proxy, (set, s, rname, pr(1.6, 2.2), pr(2.2, 10.0), pr(10.0, 30.0), 100mean(filter(!isnan, err(az, el)) .< 3)))
    end
    @info @sprintf("set %d scenario %d: %d detections, %.0f%% followers, %.0f%% with >= 2 candidates, %.0f%% direct-path detections",
        set, s, n, last(stats).followers_pct, last(stats).ge2_pct, last(stats).direct_det_pct)
end
CSV.write(joinpath(out_dir, "data", "bench_sources_detections.csv"), det)
CSV.write(joinpath(out_dir, "bench_sources_separation.csv"), sep)
CSV.write(joinpath(out_dir, "bench_sources_proxies.csv"), proxy)
CSV.write(joinpath(out_dir, "bench_sources_scenarios.csv"), stats)
CSV.write(joinpath(out_dir, "data", "bench_sources_inputs.csv"), inputs)

io = IOBuffer()
@printf(io, "realistic simulation: A0 %.0f, hydrophone scatter %.2f, diffuse reverb %.2f, off-axis low-pass on; window %s; SIC min_rel 0.7\n", A0, SCATTER, REVERB, string(W))
for set in (2, 3)
    st = stats[stats.set .== set, :]
    nrow(st) == 0 && continue
    @printf(io, "\n=== %d dolphins: %d scenarios x %.0f s, %d detections; followers %.0f%% (real 51%%), >= 2 SIC candidates %.0f%% (real 79%%), closure <= 1 %.0f%% (real 35%%)\n",
        set, nrow(st), DUR, sum(st.detections), mean(st.followers_pct), mean(st.ge2_pct), mean(st.closure_le1_pct))
    @printf(io, "    matched detections: direct path %.0f%%, specular echo %.0f%%, diffuse %.0f%%\n", mean(st.direct_det_pct), mean(st.echo_det_pct), mean(st.diffuse_det_pct))
    println(io, "  direction within 3 deg of the emitting click's direct path (%):")
    @printf(io, "  %-36s %7s %8s %7s %10s %8s\n", "method", "all", "direct", "echo", "other-dol.", "median°")
    d = det[(det.set .== set) .& det.matched, :]
    for g in sort(collect(groupby(d, :method)); by=g -> -mean(g.err .< 3))
        @printf(io, "  %-36s %7.1f %8.1f %7.1f %10.1f %8.1f\n", g.method[1], 100mean(g.err .< 3), 100mean(g.err[g.order .== 0] .< 3),
            100mean(g.err[g.order .> 0] .< 3), 100mean(g.err[g.other_follower] .< 3), median(g.err))
    end
    sm = d[d.method .== "longest track + smart inherit", :]
    @printf(io, "  smart inherit roles: echo label precision %.1f%% / recall %.1f%% vs specular+diffuse echo detections; 'own' on other-dolphin followers %.1f%%\n",
        100mean(sm.order[sm.role .== "echo"] .!= 0), 100mean(sm.role[sm.order .!= 0] .== "echo"), 100mean(sm.role[sm.other_follower] .== "own"))
    println(io, "  source separation (pairs of detections <= 0.2 s apart):")
    @printf(io, "  %-28s %9s %7s %6s %6s %8s %8s %11s\n", "method", "precision", "recall", "F1", "ARI", "n true", "n found", "unassigned%")
    for g in groupby(sep[sep.set .== set, :], :method)
        @printf(io, "  %-28s %9.3f %7.3f %6.3f %6.3f %8.1f %8.1f %11.1f\n", g.method[1], mean(g.precision), mean(g.recall), mean(g.f1), mean(g.ari),
            mean(g.n_true), mean(g.n_found), mean(g.unassigned_pct))
    end
    px = combine(groupby(proxy[proxy.set .== set, :], :rule), :acc => mean => :acc, :burst => (x -> median(filter(!isnan, x))) => :burst,
        :p2_10 => (x -> median(filter(!isnan, x))) => :p2_10, :p10_30 => (x -> median(filter(!isnan, x))) => :p10_30)
    println(io, "  rules: accuracy vs the real-data proxies computed on the simulation (median over scenarios):")
    for r in eachrow(sort(px, :acc; rev=true))
        @printf(io, "  %-16s acc %5.1f%% | burst %5.1f°  2-10 ms %5.1f°  10-30 ms %5.1f°\n", r.rule, r.acc, r.burst, r.p2_10, r.p10_30)
    end
end
body = String(take!(io))
open(joinpath(out_dir, "bench_sources.md"), "w") do f; println(f, "# Round-4 benchmark: realistic simulation, inheritance, n dolphins\n\n```\n", body, "```"); end
print(body)
