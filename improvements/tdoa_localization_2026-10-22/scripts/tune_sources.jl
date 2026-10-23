# Tune separate_sources (src/source_separation.jl) on the 2-dolphin scenarios of the realistic simulation and check the
# chosen parameters on the 3-dolphin scenarios (not used for tuning). Inputs from bench_sources.jl
# (results/data/bench_sources_inputs.csv: longest-track directions and the truth per detection).
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/tune_sources.jl
using CSV, DataFrames, Statistics, Printf, Random
include(joinpath(@__DIR__, "..", "..", "..", "src", "source_separation.jl"))
res = joinpath(@__DIR__, "..", "results")
inp = CSV.read(joinpath(res, "data", "bench_sources_inputs.csv"), DataFrame)
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))
scen = [g for g in groupby(inp, [:set, :scenario])]

function evaluate(set, kw, nmode)
    f1 = Float64[]; pr = Float64[]; rc = Float64[]; acc = Float64[]; un = Float64[]; nf = Float64[]
    for g in scen
        g.set[1] == set || continue
        n = nmode == :given ? set : nothing
        s = separate_sources(g.t, g.az, g.el; follower=g.follower, n=n, kw...)
        sc = pairwise_scores(g.src, s.label, g.t; max_dt=0.2)
        push!(f1, sc.f1); push!(pr, sc.precision); push!(rc, sc.recall); push!(un, 100mean(s.label .== 0))
        push!(nf, nmode == :given ? length(unique(filter(>(0), s.label))) : count(s.sources.confirmed))
        m = g.matched
        push!(acc, 100mean(angd.(vec3.(s.az[m], s.el[m]), vec3.(g.true_az[m], g.true_el[m])) .< 3))
    end
    return (; f1=mean(f1), precision=mean(pr), recall=mean(rc), acc=mean(acc), unassigned=mean(un), n_found=mean(nf))
end
function kmeans_eval(set)
    f1 = Float64[]
    for g in scen
        g.set[1] == set || continue
        push!(f1, pairwise_scores(g.src, kmeans_sources(g.t, g.az, g.el, set; follower=g.follower), g.t; max_dt=0.2).f1)
    end
    return mean(f1)
end

grid = [(; sigma_meas_deg=sm, gate_deg=ga, force_gate_deg=fg, t_drop=td, follow_gate_deg=fo, merge_deg=md)
        for sm in (2.0, 5.0, 8.0, 12.0), ga in (8.0, 15.0, 25.0), fg in (12.0, 25.0, 40.0), td in (1.0, 3.0), fo in (5.0, 12.0), md in (4.0, 8.0)]
rows = Vector{Any}(undef, length(grid))
Threads.@threads for k in eachindex(grid)
    rows[k] = (grid[k], evaluate(2, grid[k], :given), evaluate(2, grid[k], :estimated))
end
tab = DataFrame(sigma_meas=Float64[], gate=Float64[], force_gate=Float64[], t_drop=Float64[], follow_gate=Float64[], merge=Float64[],
    f1_given=Float64[], f1_est=Float64[], recall_given=Float64[], precision_given=Float64[], acc_given=Float64[], unassigned_given=Float64[], n_est=Float64[])
for (kw, a, b) in rows
    push!(tab, (kw..., a.f1, b.f1, a.recall, a.precision, a.acc, a.unassigned, b.n_found))
end
sort!(tab, :f1_given; rev=true)
CSV.write(joinpath(res, "tune_sources.csv"), tab)
best = first(eachrow(tab))
bkw = (; sigma_meas_deg=best.sigma_meas, gate_deg=best.gate, force_gate_deg=best.force_gate, t_drop=best.t_drop, follow_gate_deg=best.follow_gate, merge_deg=best.merge)
best_est = first(eachrow(sort(tab, :f1_est; rev=true)))
ekw = (; sigma_meas_deg=best_est.sigma_meas, gate_deg=best_est.gate, force_gate_deg=best_est.force_gate, t_drop=best_est.t_drop, follow_gate_deg=best_est.follow_gate, merge_deg=best_est.merge)

io = IOBuffer()
@printf(io, "grid of %d parameter sets, tuned on %d two-dolphin scenarios, checked on %d three-dolphin scenarios (not used for tuning)\n",
    length(grid), count(g -> g.set[1] == 2, scen), count(g -> g.set[1] == 3, scen))
@printf(io, "best for n given: %s\nbest for n estimated: %s\n\n", string(bkw), string(ekw))
@printf(io, "%-34s %6s %9s %7s %6s %11s %8s\n", "", "F1", "precision", "recall", "acc3°", "unassigned%", "n found")
for (set, nm) in ((2, "2 dolphins (tuning set)"), (3, "3 dolphins (held out)"))
    println(io, nm)
    for (lab, kw, mode) in (("defaults, n given", (;), :given), ("defaults, n estimated", (;), :estimated),
            ("tuned, n given", bkw, :given), ("tuned (for estimated), n estimated", ekw, :estimated))
        e = evaluate(set, kw, mode)
        @printf(io, "  %-32s %6.3f %9.3f %7.3f %6.1f %11.1f %8.1f\n", lab, e.f1, e.precision, e.recall, e.acc, e.unassigned, e.n_found)
    end
    @printf(io, "  %-32s %6.3f\n", "k-means baseline (k = n)", kmeans_eval(set))
end
println(io, "\ntop 10 parameter sets (n given, 2 dolphins):")
show(io, MIME("text/plain"), first(tab, 10)); println(io)
body = String(take!(io))
open(joinpath(res, "tune_sources.md"), "w") do f; println(f, "# Tuning separate_sources on the realistic simulation\n\n```\n", body, "```"); end
print(body)
