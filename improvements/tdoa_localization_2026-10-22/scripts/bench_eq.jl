# Round-5 benchmark: hydrophone equalization on the calibrated realistic simulation (2 dolphins). Per scenario the
# equalizers are estimated from the scenario's own strong isolated clicks (4 s of data): "aligned" (estimate_hydrophone_eq,
# self-calibration by click directions), "blind" (estimate_hydrophone_eq_blind; pure delays unknown = 0) and
# "blind + true delays" (the pure delays a calibration recording would give, here from the truth); compared with no
# equalization and an oracle built from the true responses: within-window statistics, direction accuracy (longest
# track + smart inheritance) and dolphin separation (SEPARATION_PARAMS_D3, n = 2).
# needs: results/data/click_features_New_device.csv (round 1; click templates)
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/bench_eq.jl [n_scenarios]
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "pool_sim.jl"))
include(joinpath(src_dir, "source_separation.jl")); include(joinpath(src_dir, "hydrophone_eq.jl"))
include(joinpath(@__DIR__, "bench_common.jl"))
using Printf, Random
n_scen = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
A0, SCATTER, REVERB = 6000.0, 0.7, 0.1      # calibrated realism (round 4)
out_dir = joinpath(@__DIR__, "..", "results"); mkpath(joinpath(out_dir, "data"))
fs = BENCH_FS; rx = BENCH_RX; DUR = 4.0; W = -110:50
templates, _ = bench_material(joinpath(out_dir, "data", "click_features_New_device.csv"))
noise = read_filtered_block(BENCH_WAV, round(Int, 53.5fs):round(Int, 97.5fs); fs=fs)
thr = 490.376461598132
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))

# oracle equalizer from the simulator's true responses (no linear-phase removal: the full inverse)
function oracle_eq(responses, ref; N=256, reg=0.1)
    f = rfftfreq(N, fs)
    H = [[1 + sum(a * cis(-2π * fk * d / fs) for (d, a) in responses[k]; init=0.0im) for fk in f] for k in eachindex(responses)]
    G = hcat([H[k] ./ H[ref] for k in eachindex(H)]...)
    E = conj.(G) ./ (abs2.(G) .+ reg) .* (1 + reg); E[:, ref] .= 1
    return (; N, fs, f, G, coh=ones(length(f), length(H)), E, ref, band=default_gcc_band(nothing, fs), n_clicks=0, delay_removed=zeros(length(H)))
end
# the pure delays (samples) separating an estimated relative response from the truth (what calibration data would give)
function true_delays(eq, responses, ref)
    f = eq.f; ib = (f .>= eq.band[1]) .& (f .<= eq.band[2]); d = zeros(length(responses))
    for k in eachindex(responses)
        k == ref && continue
        Gt = [(1 + sum(a * cis(-2π * fk * dd / fs) for (dd, a) in responses[k]; init=0.0im)) /
              (1 + sum(a * cis(-2π * fk * dd / fs) for (dd, a) in responses[ref]; init=0.0im)) for fk in f]
        φ = unwrap(angle.(Gt[ib] ./ eq.G[ib, k])); X = [ones(count(ib)) f[ib]]
        d[k] = -(X \ φ)[2] * fs / 2π
    end
    return d
end
# phase error (rad) of an estimated relative response vs the truth, in the band; linear phase removed from both
function phase_error(eq, responses, ref)
    f = eq.f; ib = (f .>= eq.band[1]) .& (f .<= eq.band[2])
    errs = Float64[]
    for k in eachindex(responses)
        k == ref && continue
        Gt = [(1 + sum(a * cis(-2π * fk * d / fs) for (d, a) in responses[k]; init=0.0im)) /
              (1 + sum(a * cis(-2π * fk * d / fs) for (d, a) in responses[ref]; init=0.0im)) for fk in f]
        φ = unwrap(angle.(Gt[ib] ./ eq.G[ib, k]))
        X = [ones(count(ib)) f[ib]]; ab = X \ φ
        push!(errs, median(abs.(φ .- X * ab)))
    end
    return mean(errs)
end

rows = DataFrame(scenario=Int[], version=String[], detections=Int[], ge2=Float64[], closure_le1=Float64[], ncc_med=Float64[],
    acc_all=Float64[], acc_direct=Float64[], med_err=Float64[], sep_f1=Float64[], phase_err=Float64[], removed_delay=Float64[])
for s in 1:n_scen
    rng = Xoshiro(50_000 + s)
    pool, array_pos = random_pool_setup(rng)
    tpl = templates[randperm(rng, length(templates))[1:10]]
    sim = simulate_scenario(rng, tpl, noise; fs=fs, dur=DUR, nsrc=2, pool=pool, array_pos=array_pos, A0=A0, rx=rx,
        hyd_scatter=SCATTER, reverb_level=REVERB, offaxis_lowpass=true)
    res = detect_impulse(("sim$s", sim.data, fs), nothing; threshold=thr, return_datafilt=true)
    rt = detect_impulsetrain(res, nothing)
    P = rt.pind_good; n = length(P)
    n < 30 && continue
    tsec = (P .- 1) ./ fs; ff = follower_flags(P, fs)
    arr = sim.arrivals; clk = sim.clicks; srt = sortperm(arr.sample); smp = arr.sample[srt]
    matched = falses(n); src = zeros(Int, n); ord = fill(-2, n); tz = fill(NaN, n); te = fill(NaN, n)
    for i in 1:n
        lo = searchsortedfirst(smp, P[i] - 35); hi = searchsortedlast(smp, P[i] + 12)
        lo > hi && continue
        j = srt[lo:hi][argmax(abs.(arr.amp[srt[lo:hi]]))]
        matched[i] = true; src[i] = arr.src[j]; ord[i] = arr.order[j]
        r = findfirst((clk.src .== src[i]) .& (clk.click .== arr.click[j])); tz[i] = clk.direct_az[r]; te[i] = clk.direct_el[r]
    end
    # equalizer from the strong isolated clicks of this scenario (upper half of the amplitudes)
    iso = findall(.!ff.follower)
    strong = P[iso[rt.ppeak[iso] .>= median(rt.ppeak[iso])]]
    eq_est = estimate_hydrophone_eq(res.data_filt, strong; fs=fs, rx_vect=rx, ref=ref_channel)
    eq_bl = estimate_hydrophone_eq_blind(res.data_filt, strong; fs=fs, rx_vect=rx, ref=ref_channel)
    eq_orc = oracle_eq(sim.responses, ref_channel)
    eq_blc = with_delays(eq_bl, true_delays(eq_bl, sim.responses, ref_channel))
    for (ver, d, pe, rd) in (("none", res.data_filt, NaN, NaN),
            ("aligned", apply_hydrophone_eq(res.data_filt, eq_est), phase_error(eq_est, sim.responses, ref_channel), maximum(abs, eq_est.delay_removed)),
            ("blind", apply_hydrophone_eq(res.data_filt, eq_bl), phase_error(eq_bl, sim.responses, ref_channel), NaN),
            ("blind + true delays", apply_hydrophone_eq(res.data_filt, eq_blc), phase_error(eq_bl, sim.responses, ref_channel), NaN),
            ("oracle", apply_hydrophone_eq(res.data_filt, eq_orc), 0.0, 0.0))
        c = get_srp_candidates_sic(d, P; window=W, fs=fs, rx_vect=rx, min_rel=0.7)
        T = track_candidates(tsec, c.uv, c.rscore)
        conf = Dict(zip(T.tracks.id, T.tracks.confirmed)); hits = Dict(zip(T.tracks.id, T.tracks.hits))
        sel = [(ks = [k for k in 1:3 if T.track[i, k] > 0 && conf[T.track[i, k]]]; isempty(ks) ? 1 : ks[argmax([hits[T.track[i, k]] for k in ks])]) for i in 1:n]
        az = [rad2deg(uv2azel(c.uv[i, sel[i], 1], c.uv[i, sel[i], 2])[1]) for i in 1:n]
        el = [rad2deg(uv2azel(c.uv[i, sel[i], 1], c.uv[i, sel[i], 2])[2]) for i in 1:n]
        ss = separate_sources(tsec, az, el; follower=ff.follower, n=2, SEPARATION_PARAMS_D3...)
        err = [matched[i] ? angd(vec3(ss.az[i], ss.el[i]), vec3(tz[i], te[i])) : NaN for i in 1:n]
        m = matched
        push!(rows, (s, ver, n, 100mean(c.nfound .>= 2), 100mean(c.closure .<= 1), median(c.ncc_min), 100mean(err[m] .< 3),
            100mean(err[m .& (ord .== 0)] .< 3), median(err[m]), pairwise_scores(src, ss.label, tsec; max_dt=0.2).f1, pe, rd))
    end
    @info @sprintf("scenario %d: %d detections; acc none / aligned / blind / blind+delays / oracle %.0f / %.0f / %.0f / %.0f / %.0f%%", s, n,
        rows.acc_all[end-4:end]...)
end
CSV.write(joinpath(out_dir, "bench_eq_scenarios.csv"), rows)

io = IOBuffer()
@printf(io, "%d scenarios x %.0f s, 2 dolphins, calibrated realism (A0 %.0f, scatter %.1f, reverb %.1f); SIC %s, longest track + smart inherit, separation n = 2 (tuned)\n\n",
    length(unique(rows.scenario)), DUR, A0, SCATTER, REVERB, string(W))
@printf(io, "%-19s | %7s %9s %7s | %8s %10s %8s | %7s | %9s %9s\n", "version", ">=2 cand", "closure<=1", "ncc", "acc3° all", "acc3° dir.", "med err", "sep F1",
    "phase err", "|delay|")
for v in ("none", "aligned", "blind", "blind + true delays", "oracle")
    g = rows[rows.version .== v, :]
    @printf(io, "%-19s | %7.1f%% %8.1f%% %7.2f | %7.1f%% %9.1f%% %7.1f° | %7.3f | %8.2f %8.1f\n", v, mean(g.ge2), mean(g.closure_le1), median(g.ncc_med),
        mean(g.acc_all), mean(g.acc_direct), median(g.med_err), mean(g.sep_f1), (x = filter(!isnan, g.phase_err); isempty(x) ? NaN : mean(x)), (x = filter(!isnan, g.removed_delay); isempty(x) ? NaN : mean(x)))
end
# per-scenario gains over no equalization
gn = rows[rows.version .== "none", :]; go = rows[rows.version .== "oracle", :]
for v in ("aligned", "blind", "blind + true delays")
    ge = rows[rows.version .== v, :]
    gain = go.acc_all .- gn.acc_all
    @printf(io, "\n%-19s - none: accuracy %+.1f points (median; better in %d of %d), separation F1 %+.3f (median); share of the oracle's gain recovered %.0f%% (scenarios where the oracle gains >= 5 points)",
        v, median(ge.acc_all .- gn.acc_all), count(ge.acc_all .> gn.acc_all), nrow(ge), median(ge.sep_f1 .- gn.sep_f1),
        100mean(((ge.acc_all .- gn.acc_all) ./ gain)[gain .>= 5]))
end
println(io)
body = String(take!(io))
open(joinpath(out_dir, "bench_eq.md"), "w") do f; println(f, "# Round-5 benchmark: hydrophone equalization (simulation)\n\n```\n", body, "```"); end
print(body)
