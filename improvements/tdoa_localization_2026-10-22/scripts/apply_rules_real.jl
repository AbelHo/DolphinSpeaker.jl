# Apply a tracker setting and selection rule (as benchmarked in bench_sequences.jl) to the real New_device candidates
# of a round-2 feature table, without recomputing the candidates: reruns the tracker per file, computes every rule's
# real-data proxies, and writes the viewer CSV of the chosen rule.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/apply_rules_real.jl <sic|sicw> "<tracker>" "<rule>" [res_folder]
#   tracker: default | gate 4 | gate 9 | birth 0.1 | birth 0.3 | delete 0.5 s | q 10 | q 40
#   rule:    rank 1 | earliest | track lag | longest track, optionally + " + followers inherit"
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl"))
using Printf
tag = ARGS[1]; tracker = ARGS[2]; chosen = ARGS[3]
res_folder = length(ARGS) >= 4 ? ARGS[4] : "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device"
const TRACKERS = Dict("default" => (;), "gate 4" => (; gate_deg=4.0), "gate 9" => (; gate_deg=9.0), "birth 0.1" => (; birth_min_score=0.1),
    "birth 0.3" => (; birth_min_score=0.3), "delete 0.5 s" => (; t_delete=0.5), "q 10" => (; sigma_q_deg=10.0), "q 40" => (; sigma_q_deg=40.0))
res = joinpath(@__DIR__, "..", "results")
folder = basename(rstrip(res_folder, '/')); sfx = folder == "New_device" ? "" : "_$folder"   # output names for other recordings get the folder suffix
f = CSV.read(joinpath(res, "data", "click_features_$(tag)_$folder.csv"), DataFrame)
N = nrow(f); K = 3
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))

sel = Dict{String,Vector{Int}}(); status = fill("", N)
for r in ("rank 1", "earliest", "track lag", "longest track"); sel[r] = ones(Int, N); end
next_id = 1
for fl in unique(f.file)
    ix = findall(f.file .== fl); n = length(ix)
    uv = fill(NaN, n, K, 2); rs = fill(NaN, n, K); sc = fill(NaN, n, K); ta = fill(NaN, n, K)
    for k in 1:K
        az = f[ix, "az_$k"]; el = f[ix, "el_$k"]
        uv[:, k, 1] = sind.(az) .* cosd.(el); uv[:, k, 2] = sind.(el)
        rs[:, k] = f[ix, "rscore_$k"]; sc[:, k] = f[ix, "score_$k"]; ta[:, k] = f[ix, "tarr_$k"]
    end
    T = track_candidates(f.t_s[ix], uv, rs; first_id=next_id, TRACKERS[tracker]...)
    isempty(T.tracks.id) || (global next_id = maximum(T.tracks.id) + 1)
    conf = Dict(zip(T.tracks.id, T.tracks.confirmed)); hits = Dict(zip(T.tracks.id, T.tracks.hits))
    St = select_candidates(T.track, ta, sc, conf, hits)
    sel["earliest"][ix] = [(ks = findall(k -> !isnan(ta[i, k]), 1:K); isempty(ks) ? 1 : ks[argmin(ta[i, ks])]) for i in 1:n]
    sel["track lag"][ix] = St.sel; status[ix] = St.status
    sel["longest track"][ix] = [(ks = [k for k in 1:K if T.track[i, k] > 0 && conf[T.track[i, k]]]; isempty(ks) ? 1 : ks[argmax([hits[T.track[i, k]] for k in ks])]) for i in 1:n]
end
dirs(s) = ([f[i, "az_$(s[i])"] for i in 1:N], [f[i, "el_$(s[i])"] for i in 1:N])
function inherit(az, el)
    a, e = copy(az), copy(el)
    for i in 2:N; (f.follower[i] && f.file[i] == f.file[i-1]) && (a[i] = a[i-1]; e[i] = e[i-1]); end
    return a, e
end
rules = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
for (r, s) in sel; rules[r] = dirs(s); rules["$r + followers inherit"] = inherit(dirs(s)...); end

iso = .!f.follower
function proxy(az, el, lo, hi)
    a = [angd(vec3(az[i], el[i]), vec3(az[i-1], el[i-1])) for i in 2:N if f.file[i] == f.file[i-1] && iso[i] && iso[i-1] && lo <= f.gap_prev_ms[i] < hi]
    return median(a)
end
inframe(az, el) = 100mean((abs.(az) .<= 40) .& (abs.(el) .<= 27.23))
io = IOBuffer()
@printf(io, "real %s candidates (%s), tracker \"%s\"; proxies: median direction change of consecutive isolated clicks; share inside the frame\n", folder, tag, tracker)
@printf(io, "%-36s %8s %8s %8s %7s\n", "rule", "burst", "2-10 ms", "10-30 ms", "frame")
out = DataFrame(rule=String[], burst=Float64[], p2_10=Float64[], p10_30=Float64[], in_frame=Float64[])
for r in sort(collect(keys(rules)))
    az, el = rules[r]
    push!(out, (r, proxy(az, el, 1.6, 2.2), proxy(az, el, 2.2, 10.0), proxy(az, el, 10.0, 30.0), inframe(az, el)))
    @printf(io, "%-36s %7.1f° %7.1f° %7.1f° %6.1f%%\n", r, out.burst[end], out.p2_10[end], out.p10_30[end], out.in_frame[end])
end
slug = replace(tracker, " " => "", "." => "") * "_" * replace(chosen, " + followers inherit" => "-inherit", " " => "")
CSV.write(joinpath(res, "real_rules_$(tag)_$(slug)$sfx.csv"), out)
body = String(take!(io))
open(joinpath(res, "real_rules_$(tag)_$(slug)$sfx.md"), "w") do g; println(g, "# Selection rules on real data: $tag, tracker \"$tracker\"\n\n```\n", body, "```"); end
print(body)

# viewer CSV of the chosen rule: click colours / frames from Impulse6; triangle = follower or (track rule) not "direct"
w6 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse6_srpPhat.csv"), DataFrame)
@assert all(f.frame .== w6.frame[1:N])
az, el = rules[chosen]
W, H = 2740, 1520; fovh, fovv = 80.0, 54.46086642417594
px = az ./ (fovh / 2) .* (W / 2) .+ W / 2; py = el ./ (-fovv / 2) .* (H / 2) .+ H / 2
tri = startswith(chosen, "track lag") ? (status .!= "direct") .| f.follower : f.follower
imp = DataFrame(frame=w6.frame[1:N], px=px, py=py, radius=w6.radius[1:N], r=w6.r[1:N], g=w6.g[1:N], b=w6.b[1:N], a=w6.a[1:N],
    shape=[t ? "triangle" : "circle" for t in tri])
imp.type = imp.shape
# per-click directions of the chosen rule (for export_audacity.jl)
CSV.write(joinpath(res, "data", "real_rules_$(tag)_$(slug)_$folder.csv"), DataFrame(file=f.file, t_s=f.t_s, az=az, el=el, triangle=tri))
fname = joinpath(res_folder, "detection_pixels__color255_webui__Impulse10$(tag)_$(slug).csv")
CSV.write(fname, vcat(imp, w6[N+1:end, names(imp)]))
@printf("wrote %s (%.1f%% inside the frame, %.1f%% triangles)\n", fname, 100mean((0 .<= px .<= W) .& (0 .<= py .<= H)), 100mean(tri))
