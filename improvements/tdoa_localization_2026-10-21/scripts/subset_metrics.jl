# Real-data metrics on subsets of probable true clicks (false detections inflate inconsistency for every method)
using CSV, DataFrames, Statistics, Printf, Random
pc = CSV.read(ARGS[1], DataFrame)
dirv(a, e) = (sind(a) * cosd(e), sind(e), cosd(a) * cosd(e))
gc(a1, e1, a2, e2) = (d1 = dirv(a1, e1); d2 = dirv(a2, e2); acosd(clamp(sum(d1 .* d2), -1, 1)))
methods = unique(pc.method)
ref = pc[pc.method .== "gcc phat", :]                    # per-click quantities shared by all methods
ppk = ref.ppeak
gate_closure = ref.closure .< 1.0
function metrics(g, sel)
    idx = findall(sel)
    same = [(g.train[idx[i]] == g.train[idx[i+1]] && g.train[idx[i]] > 0) for i in 1:length(idx)-1]
    jit = [gc(g.az_ls[idx[i]], g.el_ls[idx[i]], g.az_ls[idx[i+1]], g.el_ls[idx[i+1]]) for i in findall(same)]
    rng = MersenneTwister(1)
    a = rand(rng, idx, 20000); b = rand(rng, idx, 20000)
    spread = median(gc(g.az_ls[i], g.el_ls[i], g.az_ls[j], g.el_ls[j]) for (i, j) in zip(a, b))
    return 100 * mean(g.feasible[idx]), median(jit), spread, median(jit) / spread
end
subsets = [("all clicks", trues(nrow(ref))),
           ("top 50% ppeak", ppk .>= quantile(ppk, 0.5)),
           ("top 25% ppeak", ppk .>= quantile(ppk, 0.75)),
           ("top 10% ppeak", ppk .>= quantile(ppk, 0.9)),
           ("gcc-phat closure < 1 smp", gate_closure)]
for (lab, sel) in subsets
    @printf("\n== %s: %d clicks (%.0f%%) ==\n", lab, count(sel), 100 * mean(sel))
    @printf("%-32s %9s %8s %8s %6s\n", "method", "feasible%", "jitter°", "spread°", "ratio")
    rows = []
    for m in methods
        g = pc[pc.method .== m, :]
        push!(rows, (m, metrics(g, sel)...))
    end
    for r in sort(rows; by=r -> r[5])
        @printf("%-32s %9.1f %8.1f %8.1f %6.2f\n", r...)
    end
end
println("\nshare of clicks passing the gcc-phat closure < 1 gate, by ppeak quartile (Q1 weakest):")
q = quantile(ppk, [0.25, 0.5, 0.75])
for k in 1:4
    sel = k == 1 ? ppk .<= q[1] : k == 4 ? ppk .> q[3] : (ppk .> q[k-1]) .& (ppk .<= q[k])
    @printf("  Q%d: %.1f%%\n", k, 100 * mean(gate_closure[sel]))
end
