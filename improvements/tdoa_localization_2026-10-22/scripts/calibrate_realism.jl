# Calibrate the simulator's realism options (src/pool_sim.jl: hyd_scatter, reverb_level, offaxis_lowpass) and source
# level against the real within-window statistics of New_device (SIC -50:50 feature table).
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/calibrate_realism.jl [n_scenarios] [A0 ...]
#   grid via environment variables (defaults = the first sweep): SCATTERS="0,0.3,0.6" REVERBS="0,0.1,0.25" TAG="" (output suffix)
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "pool_sim.jl"))
include(joinpath(@__DIR__, "bench_common.jl"))
using Printf, Random
n_scen = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8
A0s = length(ARGS) >= 2 ? parse.(Float64, ARGS[2:end]) : [9000.0]
out_dir = joinpath(@__DIR__, "..", "results")
fs = BENCH_FS; rx = BENCH_RX; DUR = 3.0
templates, _ = bench_material(joinpath(out_dir, "data", "click_features_New_device.csv"))
noise = read_filtered_block(BENCH_WAV, round(Int, 53.5fs):round(Int, 97.5fs); fs=fs)
thr = 490.376461598132   # file 0001's detection threshold (see bench_sequences.jl)

# real targets (SIC -50:50, New_device): name => (value, scale for the distance)
const TARGETS = [("ge2 candidates %", 78.0, 10.0), ("closure<=1 %", 35.1, 10.0), ("ncc_min median", 0.21, 0.1),
    ("efrac_1 median", 0.30, 0.1), ("env_peaks>=2 %", 95.1, 10.0), ("followers %", 51.4, 10.0), ("ppeak/thr median", 1.79, 0.5)]

function stats_for(realism, A0)
    st = Dict{String,Vector{Float64}}(k => Float64[] for k in ("nf", "clo", "ncc", "efr", "env", "fol", "rat", "dir"))
    for s in 1:n_scen
        rng = Xoshiro(20_000 + s)
        pool, array_pos = random_pool_setup(rng)
        tpl = templates[randperm(rng, length(templates))[1:10]]
        sim = simulate_scenario(rng, tpl, noise; fs=fs, dur=DUR, nsrc=2, pool=pool, array_pos=array_pos, A0=A0, rx=rx, realism...)
        res = detect_impulse(("sim$s", sim.data, fs), nothing; threshold=thr, return_datafilt=true)
        rt = detect_impulsetrain(res, nothing)
        P = rt.pind_good
        length(P) < 20 && continue
        c = get_srp_candidates_sic(res.data_filt, P; window=window_impulsive, fs=fs, rx_vect=rx, min_rel=0.7)
        ea = envelope_arrivals(res.data_filt, P; window=window_impulsive)
        append!(st["nf"], c.nfound .>= 2); append!(st["clo"], c.closure); append!(st["ncc"], c.ncc_min); append!(st["efr"], c.efrac[:, 1])
        append!(st["env"], vec(maximum(ea; dims=2)) .>= 2); append!(st["fol"], follower_flags(P, fs).follower); append!(st["rat"], rt.ppeak ./ thr)
        # truth: share of detections whose strongest arrival within ±12 samples is a direct path
        arr = sim.arrivals; srt = sortperm(arr.sample); smp = arr.sample[srt]
        for p in P
            lo = searchsortedfirst(smp, p - 12); hi = searchsortedlast(smp, p + 12)
            lo > hi && continue
            push!(st["dir"], arr.order[srt[lo:hi][argmax(abs.(arr.amp[srt[lo:hi]]))]] == 0)
        end
    end
    v = [100mean(st["nf"]), 100mean(st["clo"] .<= 1), median(st["ncc"]), median(st["efr"]), 100mean(st["env"]), 100mean(st["fol"]), median(st["rat"])]
    return v, 100mean(st["dir"]), length(st["nf"])
end

rows = DataFrame(A0=Float64[], hyd_scatter=Float64[], reverb_level=Float64[], ge2=Float64[], closure_le1=Float64[], ncc_med=Float64[],
    efrac_med=Float64[], env_ge2=Float64[], followers=Float64[], ratio_med=Float64[], direct_det=Float64[], n=Int[], distance=Float64[])
scatters = parse.(Float64, split(get(ENV, "SCATTERS", "0,0.3,0.6"), ","))
reverbs = parse.(Float64, split(get(ENV, "REVERBS", "0,0.1,0.25"), ","))
tag = get(ENV, "TAG", "")
for A0 in A0s, hs in scatters, rl in reverbs
    v, dirpct, n = stats_for((; hyd_scatter=hs, reverb_level=rl, offaxis_lowpass=true), A0)
    dist = sqrt(sum(((v[i] - TARGETS[i][2]) / TARGETS[i][3])^2 for i in eachindex(v)))
    push!(rows, (A0, hs, rl, v..., dirpct, n, dist))
    @info @sprintf("A0 %.0f scatter %.1f reverb %.2f: ge2 %.0f%%, closure<=1 %.0f%%, ncc %.2f, efrac %.2f, env>=2 %.0f%%, followers %.0f%%, ratio %.2f, direct %.0f%% (n %d) -> distance %.2f",
        A0, hs, rl, v..., dirpct, n, dist)
end
sort!(rows, :distance)
CSV.write(joinpath(out_dir, "calibrate_realism$tag.csv"), rows)
io = IOBuffer()
@printf(io, "%d scenarios x %.0f s per setting (2 sources, off-axis low-pass on); real targets: %s\n\n", n_scen, DUR,
    join([@sprintf("%s %.2f", t[1], t[2]) for t in TARGETS], ", "))
@printf(io, "%7s %7s %7s | %6s %8s %6s %6s %7s %6s %6s | %7s %6s %8s\n", "A0", "scatter", "reverb", "ge2%", "clo<=1%", "ncc", "efrac", "env>=2%", "fol%", "ratio", "direct%", "n", "distance")
for r in eachrow(rows)
    @printf(io, "%7.0f %7.1f %7.2f | %6.1f %8.1f %6.2f %6.2f %7.1f %6.1f %6.2f | %7.1f %6d %8.2f\n", r.A0, r.hyd_scatter, r.reverb_level, r.ge2,
        r.closure_le1, r.ncc_med, r.efrac_med, r.env_ge2, r.followers, r.ratio_med, r.direct_det, r.n, r.distance)
end
body = String(take!(io))
open(joinpath(out_dir, "calibrate_realism$tag.md"), "w") do f; println(f, "# Simulator realism calibration $tag\n\n```\n", body, "```"); end
print(body)
