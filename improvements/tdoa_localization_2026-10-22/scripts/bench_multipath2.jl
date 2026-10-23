# Round-2 multipath benchmark: which candidate generator keeps the direct path and the echoes while adding few false
# candidates? Plain local maxima (round 1, with the pair-support filter), GCC-kernel CLEAN and time-domain successive
# interference cancellation (SIC), on the same trials; a third of the trials have no echo (false-candidate rate).
# needs: results/data/click_features_New_device.csv (round 1)
# usage: julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/bench_multipath2.jl [n_trials]
src = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src, "multiarrival.jl"))
include(joinpath(@__DIR__, "bench_common.jl"))
using Printf
n_trials = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3000
out_dir = joinpath(@__DIR__, "..", "results")
templates, noises = bench_material(joinpath(out_dir, "data", "click_features_New_device.csv"))
fs = BENCH_FS; rx = BENCH_RX; W = window_impulsive
const W_WIDE = -110:50   # sic_wide: same as sic with a window reaching 110 samples (570 us) before the detection

# one row per (trial, method, candidate)
rows = DataFrame(trial=Int[], nech=Int[], delay=Int[], gain=Float64[], kind=String[], method=String[], k=Int[],
    err_direct=Float64[], err_echo=Float64[], rel=Float64[], support_rel=Float64[], tarr=Float64[])
lk = ReentrantLock()
Threads.@threads for trial in 1:n_trials
    rng = Xoshiro(1000 + trial)
    s = templates[rand(rng, eachindex(templates))]; nz = noises[rand(rng, eachindex(noises))]
    nech = rand(rng, 0:2)
    arr, info = bench_arrivals(rng, nech)
    data = bench_window(s, nz, arr); pk = bench_peak(data)
    truth = bench_uv(arr[1][1], arr[1][2]); echoes = [bench_uv(a[1], a[2]) for a in arr[2:end]]
    plain = get_srp_candidates(data, [pk]; window=W, fs=fs, rx_vect=rx, min_rel=0.0)
    sup = srp_support(data, [pk], plain.uv; window=W, fs=fs, rx_vect=rx)
    clean = get_srp_candidates_clean(data, [pk]; window=W, fs=fs, rx_vect=rx, threaded=false)
    sic = get_srp_candidates_sic(data, [pk]; window=W, fs=fs, rx_vect=rx, threaded=false)
    sicw = get_srp_candidates_sic(data, [pk]; window=W_WIDE, fs=fs, rx_vect=rx, threaded=false)   # looks further back
    local_rows = []
    for (nm, c, rel, supp) in [("plain", plain, plain.score, sup), ("clean", clean, clean.rscore, clean.pmin), ("sic", sic, sic.rscore, sic.pmin),
            ("sic_wide", sicw, sicw.rscore, sicw.pmin)]
        for k in 1:3
            isnan(c.uv[1, k, 1]) && continue
            u = Tuple(c.uv[1, k, :])
            push!(local_rows, (trial, nech, info.delay, info.gain, info.kind, nm, k, _angle_uv(u, truth),
                isempty(echoes) ? Inf : minimum(e -> _angle_uv(u, e), echoes), rel[1, k] / rel[1, 1], supp[1, k] / supp[1, 1], c.tarr[1, k]))
        end
    end
    lock(lk) do; foreach(r -> push!(rows, r), local_rows); end
end
CSV.write(joinpath(out_dir, "bench_multipath2_candidates.csv"), rows)

# accepted candidates for a method and thresholds: rank 1 always, others need rel >= r (and support_rel >= sp)
accept(df, r, sp) = df[(df.k .== 1) .| ((df.rel .>= r) .& (df.support_rel .>= sp)), :]
function stats(meth, r, sp)
    d = accept(rows[rows.method .== meth, :], r, sp)
    g = groupby(d, :trial)
    found = combine(g, :err_direct => (e -> any(e .< 2)) => :found, :nech => first => :nech,
        [:k, :err_direct, :err_echo] => ((k, ed, ee) -> count((k .> 1) .& (ed .>= 3) .& (ee .>= 3))) => :false_extra,
        [:k, :err_echo] => ((k, ee) -> any((k .> 1) .& (ee .< 3))) => :echo_found,
        [:tarr, :err_direct] => ((t, e) -> e[argmin(t)] < 2) => :earliest_ok,
        [:k, :err_direct] => ((k, e) -> e[k .== 1][1] < 2) => :rank1_ok)
    e0 = found.nech .== 0; e1 = .!e0
    return (; recall=100mean(found.found), recall_echo=100mean(found.found[e1]), rank1=100mean(found.rank1_ok[e1]),
        earliest=100mean(found.earliest_ok[e1]), false_per_click_noecho=mean(found.false_extra[e0]), false_per_click_echo=mean(found.false_extra[e1]),
        echo_found=100mean(found.echo_found[e1]), extra_noecho=100mean(found.false_extra[e0] .> 0))
end

io = IOBuffer()
@printf(io, "trials %d (a third without echo), templates %d real clean clicks, real noise, SNR %.0f dB, window %s (sic_wide: %s)\n", n_trials, length(templates), BENCH_SNR_DB, string(W), string(W_WIDE))
println(io, "recall = direct path within 2 deg among the accepted candidates (echo trials); echo found = an accepted extra candidate within 3 deg of an echo;")
println(io, "false = accepted extra candidates within 3 deg of no true arrival (per click); rank1 / earliest = selection by rank 1 / by earliest beam arrival (echo trials)\n")
@printf(io, "%-7s %-6s %-6s | %7s %7s %8s | %10s %12s %10s\n", "method", "rel>=", "supp>=", "recall", "rank1", "earliest", "echo found", "false/click", "false/click")
@printf(io, "%-7s %-6s %-6s | %7s %7s %8s | %10s %12s %10s\n", "", "", "", "(echo)", "", "", "", "(no echo)", "(echo)")
sweep = DataFrame(method=String[], rel_min=Float64[], support_min=Float64[], recall=Float64[], rank1=Float64[], earliest=Float64[],
    echo_found=Float64[], false_noecho=Float64[], false_echo=Float64[])
for (meth, sps) in [("plain", [-Inf, 0.3, 0.5]), ("clean", [-Inf]), ("sic", [-Inf]), ("sic_wide", [-Inf])], sp in sps, r in [0.0, 0.3, 0.5, 0.7, 0.9, 2.0]
    st = stats(meth, r, sp)
    push!(sweep, (meth, r, sp, st.recall_echo, st.rank1, st.earliest, st.echo_found, st.false_per_click_noecho, st.false_per_click_echo))
    @printf(io, "%-7s %-6.1f %-6s | %7.1f %7.1f %8.1f | %10.1f %12.2f %10.2f\n", meth, r, sp == -Inf ? "-" : string(sp), st.recall_echo, st.rank1, st.earliest,
        st.echo_found, st.false_per_click_noecho, st.false_per_click_echo)
end
CSV.write(joinpath(out_dir, "bench_multipath2_sweep.csv"), sweep)
body = String(take!(io))
open(joinpath(out_dir, "bench_multipath2.md"), "w") do f; println(f, "# Round-2 multipath benchmark: candidate generators\n\n```\n", body, "```"); end
print(body)
