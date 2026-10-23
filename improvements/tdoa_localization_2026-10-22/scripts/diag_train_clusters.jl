# Direction clusters of clean clicks (isolated + unflagged) per click train, and whether the second cluster is
# time-locked to the first (echo) or has its own rhythm (other animal / movement).
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/diag_train_clusters.jl <res_folder> [out_dir]
include(joinpath(@__DIR__, "common.jl"))
res_folder = ARGS[1]
out_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "results")
df, info = load_clicks(res_folder)
N = nrow(df); fs = info.fs
clean = (df.gap_prev .>= 1.6) .& .!df.flag
V = vec3.(df.az6, df.el6); tms = df.t .* 1000
io = IOBuffer()
@printf(io, "clean clicks: %d of %d\n\n2D histogram of clean SRP directions (%% of clean, 10 deg bins; rows el, cols az)\n%8s", count(clean), N, "el\\az")
for a in -90:10:80; @printf(io, "%5d", a); end; println(io)
for e in 80:-10:-90
    @printf(io, "%8d", e)
    for a in -90:10:80
        s = 100count(i -> clean[i] && a <= df.az6[i] < a + 10 && e <= df.el6[i] < e + 10, 1:N) / count(clean)
        print(io, s < 0.3 ? "    ." : @sprintf("%5.1f", s))
    end
    println(io)
end

# greedy clustering (8 deg radius) of the clean clicks of each train with >= 20 of them
ncl = Int[]; top = Float64[]; sec = Float64[]
locked = Float64[]; null = Float64[]; ici1 = Float64[]; ici2 = Float64[]; dazs = Float64[]; dels = Float64[]
Random.seed!(3)
for grp in groupby(DataFrame(g=1:N, file=df.file, train=df.train), [:file, :train])
    all_idx = grp.g; idx = filter(i -> clean[i], all_idx)
    length(idx) < 20 && continue
    rem = copy(idx); sizes = Int[]; centres = Int[]
    while !isempty(rem) && length(sizes) < 6
        s = rem[argmax([count(j -> angd(V[i], V[j]) < 8, rem) for i in rem])]
        mem = filter(j -> angd(V[s], V[j]) < 8, rem)
        push!(sizes, length(mem)); push!(centres, s); rem = setdiff(rem, mem)
    end
    fr = sizes ./ length(idx)
    push!(ncl, count(fr .>= 0.10)); push!(top, fr[1]); push!(sec, length(fr) > 1 ? fr[2] : 0.0)
    (length(fr) > 1 && fr[2] >= 0.10) || continue
    c1, c2 = centres[1], centres[2]
    m1 = filter(j -> angd(V[c1], V[j]) < 8, all_idx)            # all detections near cluster 1, flagged too
    m2 = filter(j -> clean[j] && angd(V[c2], V[j]) < 8, all_idx)
    push!(dazs, df.az6[c2] - df.az6[c1]); push!(dels, df.el6[c2] - df.el6[c1])
    t1 = sort(tms[m1]); span = tms[all_idx[end]] - tms[all_idx[1]]
    for j in m2
        p = searchsortedlast(t1, tms[j] - 1e-9); push!(locked, p == 0 ? Inf : tms[j] - t1[p])
        tr_ = tms[all_idx[1]] + rand() * span; p = searchsortedlast(t1, tr_); push!(null, p == 0 ? Inf : tr_ - t1[p])
    end
    append!(ici1, diff(sort(tms[filter(j -> clean[j], m1)]))); append!(ici2, diff(sort(tms[m2])))
end
println(io, "\ntrains analysed (>= 20 clean clicks): ", length(ncl))
for c in 0:6; @printf(io, "  %d clusters >= 10%%: %d trains\n", c, count(==(c), ncl)); end
@printf(io, "  share of largest cluster: median %.2f; second: median %.2f\n", median(top), median(sec))
@printf(io, "\nsecond cluster (n trains=%d): daz p25/p50/p75 %.1f/%.1f/%.1f, del %.1f/%.1f/%.1f, |daz|<5: %.0f%%\n", length(dazs),
    quantile(dazs, [.25,.5,.75])..., quantile(dels, [.25,.5,.75])..., 100mean(abs.(dazs) .< 5))
println(io, "lag from the latest cluster-1 detection to each clean cluster-2 click (random times in the same train):")
for (lo, hi) in [(0,1.6),(1.6,2.5),(2.5,5),(5,10),(10,30),(30,Inf)]
    @printf(io, "  %5.1f-%-5.1f ms: %5.1f%%  (random %5.1f%%)\n", lo, hi, 100mean((locked .>= lo) .& (locked .< hi)), 100mean((null .>= lo) .& (null .< hi)))
end
@printf(io, "ICI (ms) cluster1 p10/p25/p50/p75/p90: %.2f/%.2f/%.2f/%.2f/%.2f\n", quantile(ici1, [.1,.25,.5,.75,.9])...)
@printf(io, "ICI (ms) cluster2 p10/p25/p50/p75/p90: %.2f/%.2f/%.2f/%.2f/%.2f\n", quantile(ici2, [.1,.25,.5,.75,.9])...)
write_report(joinpath(out_dir, "diag_train_clusters.md"), "Direction clusters per click train", String(take!(io));
    notes="Source: `$res_folder`. Clean = gap to previous detection >= 1.6 ms and Impulse5 not flagged; directions from Impulse6 (SRP).")
