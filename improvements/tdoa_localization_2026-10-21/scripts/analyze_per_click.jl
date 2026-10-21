using CSV, DataFrames, Statistics, Printf, Random
pc = CSV.read(ARGS[1], DataFrame)
dirv(a, e) = (sind(a) * cosd(e), sind(e), cosd(a) * cosd(e))
gc(a1, e1, a2, e2) = (d1 = dirv(a1, e1); d2 = dirv(a2, e2); acosd(clamp(sum(d1 .* d2), -1, 1)))
Random.seed!(1)
@printf("%-32s %8s %8s %7s %8s %10s %10s\n", "method", "jit°", "spread°", "ratio", "ident%", "clo>3 %", "az=±90 %")
for g in groupby(pc, :method; sort=false)
    n = nrow(g)
    same = [(g.train[i] == g.train[i+1] && g.train[i] > 0) for i in 1:n-1]
    idx = findall(same)
    jit = [gc(g.az_ls[i], g.el_ls[i], g.az_ls[i+1], g.el_ls[i+1]) for i in idx]
    ident = mean([g.tdoa1[i] == g.tdoa1[i+1] && g.tdoa2[i] == g.tdoa2[i+1] && g.tdoa3[i] == g.tdoa3[i+1] for i in idx])
    r1 = rand(1:n, 20000); r2 = rand(1:n, 20000)
    spread = median(gc(g.az_ls[a], g.el_ls[a], g.az_ls[b], g.el_ls[b]) for (a, b) in zip(r1, r2))
    clo = all(isnan, g.closure) ? NaN : 100 * mean(g.closure .> 3)
    edge = 100 * mean(abs.(g.az_ls) .> 89.9)
    @printf("%-32s %8.2f %8.2f %7.3f %8.1f %10s %10.1f\n", g.method[1], median(jit), spread, median(jit) / spread, 100ident,
        isnan(clo) ? "-" : @sprintf("%.1f", clo), edge)
end
g = pc[pc.method .== "gcc phat", :]
q = quantile(g.ppeak, [0.25, 0.5, 0.75])
println("\ngcc phat closure by ppeak quartile (median, % > 3 samples):")
for k in 1:4
    sel = (k == 1 ? g.ppeak .<= q[1] : k == 4 ? g.ppeak .> q[3] : (g.ppeak .> q[k-1]) .& (g.ppeak .<= q[k]))
    @printf("  Q%d: median %.2f, >3: %.1f%%\n", k, median(g.closure[sel]), 100 * mean(g.closure[sel] .> 3))
end
println("gcc phat closure quantiles (10,25,50,75,90): ", round.(quantile(g.closure, [0.1, 0.25, 0.5, 0.75, 0.9]); digits=2))
