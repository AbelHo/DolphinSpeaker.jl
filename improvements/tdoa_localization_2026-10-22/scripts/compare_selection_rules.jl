# Offline comparison of per-click selection rules on a round-2 feature table (no recomputation).
# Proxy for "same source, same path": consecutive isolated clicks should keep their direction (burst pulses 1.6-2.2 ms
# apart come from one animal within milliseconds). Also the share of points inside the video frame.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/compare_selection_rules.jl sic|sicw [folder name]
using CSV, DataFrames, Statistics, Printf
tag = length(ARGS) >= 1 ? ARGS[1] : "sicw"
res = joinpath(@__DIR__, "..", "results")
folder = length(ARGS) >= 2 ? ARGS[2] : "New_device"; sfx = folder == "New_device" ? "" : "_$folder"
f = CSV.read(joinpath(res, "data", "click_features_$(tag)_$folder.csv"), DataFrame)
N = nrow(f); K = 3
W, H = 2740, 1520; fovh, fovv = 80.0, 54.46086642417594
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))
col(nm, k) = f[!, "$(nm)_$k"]
cand(i) = [k for k in 1:K if !isnan(col("az", k)[i])]

# rules: candidate index per click
rules = Dict{String,Vector{Int}}()
rules["rank 1 (SRP maximum)"] = ones(Int, N)
rules["earliest arrival"] = f.selected_earliest
rules["track lag rule"] = f.selected_track
for x in (0.85, 0.95)
    rules[@sprintf("earliest, rscore >= %.2f x best", x)] = [(ks = filter(k -> col("rscore", k)[i] >= x * f.rscore_1[i], cand(i)); ks[argmin(col("tarr", k)[i] for k in ks)]) for i in 1:N]
end
for dmax in (20, 40)
    rules["earliest, arriving <= $dmax smp before the detection"] = [(ks = filter(k -> col("tarr", k)[i] >= -dmax, cand(i)); isempty(ks) ? 1 : ks[argmin(col("tarr", k)[i] for k in ks)]) for i in 1:N]
end
rules["largest energy removed (efrac)"] = [(ks = cand(i); ks[argmax(col("efrac", k)[i] for k in ks)]) for i in 1:N]

iso = .!f.follower
function consistency(sel, lo, hi)
    a = Float64[]
    for i in 2:N
        (f.file[i] == f.file[i-1] && iso[i] && iso[i-1] && lo <= f.gap_prev_ms[i] < hi) || continue
        push!(a, angd(vec3(col("az", sel[i])[i], col("el", sel[i])[i]), vec3(col("az", sel[i-1])[i-1], col("el", sel[i-1])[i-1])))
    end
    return median(a), 100mean(a .< 2), length(a)
end
inside(sel) = 100mean(i -> (az = col("az", sel[i])[i]; el = col("el", sel[i])[i]; abs(az) <= fovh / 2 && abs(el) <= fovv / 2), 1:N)

io = IOBuffer()
@printf(io, "feature table: click_features_%s_%s.csv (%d clicks); consecutive isolated clicks (no detection 0.4-1.6 ms before either)\n", tag, folder, N)
@printf(io, "%-48s %9s | %-16s | %-16s | %-16s | %6s\n", "rule", "≠ rank 1", "1.6-2.2 ms", "2.2-10 ms", "10-30 ms", "in")
@printf(io, "%-48s %9s | %-16s | %-16s | %-16s | %6s\n", "", "", "median / <2°", "median / <2°", "median / <2°", "frame")
out = DataFrame(rule=String[], changed_pct=Float64[], med_burst=Float64[], lt2_burst=Float64[], med_2_10=Float64[], lt2_2_10=Float64[],
    med_10_30=Float64[], lt2_10_30=Float64[], inside_pct=Float64[])
for nm in sort(collect(keys(rules)); by=x -> x == "rank 1 (SRP maximum)" ? "" : x)
    sel = rules[nm]
    c1 = consistency(sel, 1.6, 2.2); c2 = consistency(sel, 2.2, 10.0); c3 = consistency(sel, 10.0, 30.0)
    ch = 100mean(sel .!= 1); ins = inside(sel)
    push!(out, (nm, ch, c1[1], c1[2], c2[1], c2[2], c3[1], c3[2], ins))
    @printf(io, "%-48s %8.1f%% | %5.1f° / %5.1f%% | %5.1f° / %5.1f%% | %5.1f° / %5.1f%% | %5.1f%%\n", nm, ch, c1[1], c1[2], c2[1], c2[2], c3[1], c3[2], ins)
end
@printf(io, "pairs: %d (1.6-2.2 ms), %d (2.2-10 ms), %d (10-30 ms)\n", consistency(ones(Int, N), 1.6, 2.2)[3], consistency(ones(Int, N), 2.2, 10.0)[3], consistency(ones(Int, N), 10.0, 30.0)[3])
CSV.write(joinpath(res, "selection_rules_$tag$sfx.csv"), out)
body = String(take!(io))
open(joinpath(res, "selection_rules_$tag$sfx.md"), "w") do g; println(g, "# Selection rules on real data ($tag)\n\n```\n", body, "```"); end
print(body)
