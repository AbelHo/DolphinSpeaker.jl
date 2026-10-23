# Detection timing on a multitdoa result folder: interval and all-pair lag histograms, amplitude ratios,
# direction agreement by lag, follower geometry, flag rates.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/diag_detection_timing.jl <res_folder> [out_dir]
include(joinpath(@__DIR__, "common.jl"))
res_folder = ARGS[1]
out_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "results")
df, info = load_clicks(res_folder)
N = nrow(df); fs = info.fs
V = vec3.(df.az6, df.el6); V5 = vec3.(df.az5, df.el5)
same(g) = g > 1 && df.file[g] == df.file[g-1]
io = IOBuffer()
@printf(io, "in-train clicks: %d (files: %s)\n", N, join(basename.(info.aufnames), ", "))
@printf(io, "flagged (closure>1 or infeasible): %.1f%%\n", 100mean(df.flag))

edges = [0.4, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 5, 7.5, 10, 15, 20, 30, 50, 75, 100, 200, 500, 1e9]
idi = filter(isfinite, df.gap_prev)
println(io, "\n1. consecutive interval (ms), n=$(length(idi))")
for k in 1:length(edges)-1
    c = count(x -> edges[k] <= x < edges[k+1], idi)
    @printf(io, "  %7.2f-%-8.2f %7d  %5.1f%%\n", edges[k], min(edges[k+1], 9999), c, 100c / length(idi))
end

bw = 0.1; nb = 150; hist = zeros(Int, nb)
for g in 1:N
    j = g + 1
    while j <= N && df.file[j] == df.file[g] && (df.sample[j] - df.sample[g]) / fs * 1000 < nb * bw
        hist[floor(Int, (df.sample[j] - df.sample[g]) / fs * 1000 / bw) + 1] += 1; j += 1
    end
end
println(io, "\n2. all-pair lag histogram, counts per 0.1 ms (row = 1 ms)")
for r in 0:nb÷10-1; @printf(io, "  %4.1f-%4.1f ms:%s\n", r, r + 1, join(lpad.(hist[10r+1:10r+10], 7))); end
CSV.write(joinpath(out_dir, "lag_histogram.csv"), DataFrame(lag_ms_lo=(0:nb-1) .* bw, count=hist))

lags = [(0.4,0.75),(0.75,1.5),(1.5,1.8),(1.8,2.2),(2.2,3.5),(3.5,4.5),(4.5,10),(10,30),(30,1e9)]
println(io, "\n3. consecutive pairs by lag: amplitude ratio later/earlier, direction change (SRP, GCC; both unflagged SRP)")
for (lo, hi) in lags
    gs = [g for g in 2:N if same(g) && lo <= df.gap_prev[g] < hi]
    r = df.ppeak[gs] ./ df.ppeak[gs .- 1]
    a6 = [angd(V[g], V[g-1]) for g in gs]; a5 = [angd(V5[g], V5[g-1]) for g in gs]
    au = [angd(V[g], V[g-1]) for g in gs if !df.flag[g] && !df.flag[g-1]]
    @printf(io, "  %5.2f-%-7.2f n=%6d  weaker %4.1f%% ratio med %.2f p10 %.2f | SRP med %5.1f <2deg %4.1f%% | GCC med %5.1f | unflagged n=%5d med %5.1f <2deg %4.1f%%\n",
        lo, min(hi, 9999), length(gs), 100mean(r .< 1), median(r), quantile(r, 0.1), median(a6), 100mean(a6 .< 2), median(a5),
        length(au), isempty(au) ? NaN : median(au), isempty(au) ? NaN : 100mean(au .< 2))
end
Random.seed!(1)
rp = [angd(V[rand(1:N)], V[rand(1:N)]) for _ in 1:100000]
@printf(io, "  random pairs SRP: median %.1f deg, <5deg %.1f%%\n", median(rp), 100mean(rp .< 5))

println(io, "\n4. follower (0.4-1.6 ms after previous) vs predecessor, both unflagged")
fg = [g for g in 2:N if same(g) && 0.4 <= df.gap_prev[g] < 1.6 && !df.flag[g] && !df.flag[g-1]]
dz = df.az6[fg] .- df.az6[fg .- 1]; de = df.el6[fg] .- df.el6[fg .- 1]
rr = rand(1:N, length(fg)); dzr = df.az6[fg] .- df.az6[rr]
@printf(io, "  n=%d |daz|<3: %.1f%% (random %.1f%%); same az & |del|>5: %.1f%%\n", length(fg), 100mean(abs.(dz) .< 3), 100mean(abs.(dzr) .< 3), 100mean((abs.(dz) .< 3) .& (abs.(de) .> 5)))
@printf(io, "  del p10/p25/p50/p75/p90: %.1f/%.1f/%.1f/%.1f/%.1f\n", quantile(de, [.1,.25,.5,.75,.9])...)

println(io, "\n5. populations and flag rates")
iso = df.gap_prev .>= 1.6
for (nm, s) in [("isolated (gap>=1.6ms)", iso), ("follower", .!iso)]
    @printf(io, "  %-22s n=%6d (%.1f%%) flagged %.1f%% median ppeak %.0f\n", nm, count(s), 100mean(s), 100mean(df.flag[s]), median(df.ppeak[s]))
end
q = quantile(df.ppeak, [0.25, 0.5, 0.75])
for (lo, hi, nm) in [(0, q[1], "Q1 weakest"), (q[1], q[2], "Q2"), (q[2], q[3], "Q3"), (q[3], Inf, "Q4 strongest")]
    s = (df.ppeak .>= lo) .& (df.ppeak .< hi); @printf(io, "  amplitude %-12s flagged %.1f%%\n", nm, 100mean(df.flag[s]))
end
clean = iso .& .!df.flag
@printf(io, "  clean (isolated + unflagged): %d (%.1f%%)\n", count(clean), 100mean(clean))
println(io, "  consecutive clean clicks, direction change by lag:")
for (lo, hi) in [(1.6,2.2),(2.2,5),(5,10),(10,30),(30,100),(100,1e9)]
    a = Float64[]
    for f in unique(df.file)
        ix = findall((df.file .== f) .& clean)
        for m in 1:length(ix)-1
            dt = (df.sample[ix[m+1]] - df.sample[ix[m]]) / fs * 1000
            lo <= dt < hi && push!(a, angd(V[ix[m]], V[ix[m+1]]))
        end
    end
    @printf(io, "    %6.1f-%-8.1f n=%6d median %5.1f deg <2deg %5.1f%%\n", lo, min(hi, 9999), length(a), median(a), 100mean(a .< 2))
end
write_report(joinpath(out_dir, "diag_detection_timing.md"), "Detection timing diagnostics", String(take!(io));
    notes="Source: `$res_folder`. Directions: Impulse6 (SRP-PHAT) and Impulse5 (GCC-PHAT guided) pixels inverted to az/el; flag = Impulse5 triangle.")
