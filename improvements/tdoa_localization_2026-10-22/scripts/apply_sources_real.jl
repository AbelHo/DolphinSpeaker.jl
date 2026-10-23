# Separate the dolphins in the real New_device recording: longest-track directions from the round-2 SIC candidates
# (window -110:50), then separate_sources with n given (default 2: "Two-way_2males") and n estimated. Writes viewer
# CSVs coloured by dolphin (echo followers as triangles at their click's direction, unassigned clicks grey) and
# plausibility statistics (no ground truth: consistency and ICI regularity within each dolphin vs the unseparated stream).
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/apply_sources_real.jl [n] [params] [res_folder]
#   params: "tuned" (SEPARATION_PARAMS_D3, default) or "defaults" (separate_sources' own defaults)
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "source_separation.jl"))
using Printf
n_given = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
params = length(ARGS) >= 2 ? ARGS[2] : "tuned"
sep_kw = params == "tuned" ? SEPARATION_PARAMS_D3 : (;)
res_folder = length(ARGS) >= 3 ? ARGS[3] : "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device"
res = joinpath(@__DIR__, "..", "results")
folder = basename(rstrip(res_folder, '/')); sfx = folder == "New_device" ? "" : "_$folder"   # output names for other recordings get the folder suffix
f = CSV.read(joinpath(res, "data", "click_features_sicw_$folder.csv"), DataFrame)
N = nrow(f); K = 3
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))

# longest-track direction per detection (tracker rerun per file on the saved candidates)
az = fill(NaN, N); el = fill(NaN, N)
out = Dict(:given => (label=zeros(Int, N), role=fill("", N), az=zeros(N), el=zeros(N)),
           :auto => (label=zeros(Int, N), role=fill("", N), az=zeros(N), el=zeros(N)))
nsrc_auto = Dict{String,Int}(); next_auto = 0
for fl in unique(f.file)
    ix = findall(f.file .== fl); n = length(ix)
    uv = fill(NaN, n, K, 2); rs = fill(NaN, n, K)
    for k in 1:K
        a = f[ix, "az_$k"]; e = f[ix, "el_$k"]
        uv[:, k, 1] = sind.(a) .* cosd.(e); uv[:, k, 2] = sind.(e); rs[:, k] = f[ix, "rscore_$k"]
    end
    T = track_candidates(f.t_s[ix], uv, rs)
    conf = Dict(zip(T.tracks.id, T.tracks.confirmed)); hits = Dict(zip(T.tracks.id, T.tracks.hits))
    sel = [(ks = [k for k in 1:K if T.track[i, k] > 0 && conf[T.track[i, k]]]; isempty(ks) ? 1 : ks[argmax([hits[T.track[i, k]] for k in ks])]) for i in 1:n]
    az[ix] = [f[ix[i], "az_$(sel[i])"] for i in 1:n]; el[ix] = [f[ix[i], "el_$(sel[i])"] for i in 1:n]
    for (key, nn) in ((:given, n_given), (:auto, nothing))
        s = separate_sources(f.t_s[ix], az[ix], el[ix]; follower=f.follower[ix], n=nn, sep_kw...)
        lab = s.label
        if key == :auto   # labels unique across files
            lab = [l > 0 ? l + next_auto : 0 for l in lab]
            nsrc_auto[fl] = count(s.sources.confirmed)
        end
        out[key].label[ix] = lab; out[key].role[ix] = s.role; out[key].az[ix] = s.az; out[key].el[ix] = s.el
    end
    global next_auto = maximum(out[:auto].label)
end
det = DataFrame(file=f.file, click=f.click, t_s=f.t_s, frame=f.frame, follower=f.follower, az_own=az, el_own=el,
    label_n=out[:given].label, role_n=out[:given].role, az_out=out[:given].az, el_out=out[:given].el,
    label_auto=out[:auto].label, role_auto=out[:auto].role)
CSV.write(joinpath(res, "data", "real_sources_sicw_$(params)_$folder.csv"), det)

# ---- plausibility statistics -----------------------------------------------------------------------------------------
io = IOBuffer()
@printf(io, "%s, %d detections; directions: longest-track rule on SIC -110:50 candidates; separate_sources parameters: %s %s\n", folder, N, params, string(sep_kw))
for (key, nm) in ((:given, "n = $n_given"), (:auto, "n estimated"))
    o = out[key]
    @printf(io, "\n=== %s: unassigned %.1f%%; roles: %s\n", nm, 100mean(o.label .== 0),
        join([@sprintf("%s %.1f%%", r, 100mean(o.role .== r)) for r in ("leader", "own", "repeat", "echo", "unassigned")], ", "))
    if key == :auto
        @printf(io, "  confirmed sources per file: %s\n", join(["$(fl[end-7:end-4]) $(nsrc_auto[fl])" for fl in sort(collect(keys(nsrc_auto)))], ", "))
        # simultaneously active sources: distinct labels with >= 3 genuine detections per 1-s bin (bins with detections)
        act = Int[]
        for fl in unique(f.file)
            ix = findall((f.file .== fl) .& (o.label .> 0) .& (o.role .!= "echo"))
            isempty(ix) && continue
            bins = floor.(Int, f.t_s[ix])
            for b in unique(bins)
                labs = o.label[ix[bins .== b]]
                push!(act, count(l -> count(==(l), labs) >= 3, unique(labs)))
            end
        end
        @printf(io, "  simultaneously active sources per 1 s with clicks: %s\n", join([@sprintf("%d: %.0f%%", k, 100mean(act .== k)) for k in 0:4], ", "))
    end
    # per label: detections, median direction, ICI of genuine clicks
    labs = sort(filter(>(0), unique(o.label)))
    big = sort(labs; by=l -> -count(==(l), o.label))[1:min(length(labs), 6)]
    @printf(io, "  %-6s %9s %8s %8s %10s %10s %10s\n", "label", "detections", "az med", "el med", "ICI p25", "ICI med", "ICI p75")
    for l in big
        ix = findall(o.label .== l); g = filter(i -> o.role[i] != "echo", ix)
        ici = Float64[]
        for i in 2:length(g); f.file[g[i]] == f.file[g[i-1]] && push!(ici, 1000 * (f.t_s[g[i]] - f.t_s[g[i-1]])); end
        ici = filter(<(200), ici)
        @printf(io, "  %-6d %9d %7.1f° %7.1f° %8.1f ms %8.1f ms %8.1f ms\n", l, length(ix), median(o.az[ix]), median(o.el[ix]),
            quantile(ici, [0.25, 0.5, 0.75])...)
    end
    # within-label vs unseparated: direction change and ICI regularity of consecutive genuine clicks (2-30 ms apart)
    gen = findall((o.label .> 0) .& (o.role .!= "echo"))
    function chain_stats(groups)
        dch = Float64[]; cv = Float64[]
        for g in groups
            for i in 2:length(g)
                dt = 1000 * (f.t_s[g[i]] - f.t_s[g[i-1]])
                (f.file[g[i]] == f.file[g[i-1]] && 2 <= dt <= 30) || continue
                push!(dch, angd(vec3(o.az[g[i]], o.el[g[i]]), vec3(o.az[g[i-1]], o.el[g[i-1]])))
                if i >= 3
                    dt0 = 1000 * (f.t_s[g[i-1]] - f.t_s[g[i-2]])
                    2 <= dt0 <= 30 && push!(cv, abs(dt - dt0) / ((dt + dt0) / 2))
                end
            end
        end
        return median(dch), 100mean(dch .< 2), median(cv)
    end
    w = chain_stats([filter(i -> o.label[i] == l, gen) for l in unique(o.label[gen])])
    a = chain_stats([gen])
    @printf(io, "  consecutive genuine clicks 2-30 ms apart: within a dolphin median change %.1f° (<2°: %.1f%%), ICI change %.0f%% | unseparated %.1f° (<2°: %.1f%%), ICI change %.0f%%\n",
        w[1], w[2], 100w[3], a[1], a[2], 100a[3])
    @printf(io, "  inside the frame: %.1f%% (own directions %.1f%%)\n", 100mean((abs.(o.az) .<= 40) .& (abs.(o.el) .<= 27.23)),
        100mean((abs.(az) .<= 40) .& (abs.(el) .<= 27.23)))
end
body = String(take!(io))
open(joinpath(res, "real_sources_sicw_$params$sfx.md"), "w") do g; println(g, "# Dolphin separation on $folder ($params parameters)\n\n```\n", body, "```"); end
print(body)

# ---- viewer CSVs: colour = dolphin (n given: reference palette slots; estimated: track palette), grey = unassigned ----
w6 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse6_srpPhat.csv"), DataFrame)
@assert all(f.frame .== w6.frame[1:N])
W, H = 2740, 1520; fovh, fovv = 80.0, 54.46086642417594
hex(h) = (parse(Int, h[2:3]; base=16), parse(Int, h[4:5]; base=16), parse(Int, h[6:7]; base=16))
slots = hex.(["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"])
for (key, nm) in ((:given, "n$(n_given)"), (:auto, "auto"))
    o = out[key]
    px = o.az ./ (fovh / 2) .* (W / 2) .+ W / 2; py = o.el ./ (-fovv / 2) .* (H / 2) .+ H / 2
    col = [o.label[i] == 0 ? (154, 154, 154) : key == :given ? slots[mod1(o.label[i], 8)] : round.(Int, 255 .* TRACK_PALETTE[mod1(o.label[i], 8)]) for i in 1:N]
    imp = DataFrame(frame=w6.frame[1:N], px=px, py=py, radius=w6.radius[1:N], r=Float64.(first.(col)), g=Float64.(getindex.(col, 2)),
        b=Float64.(last.(col)), a=[o.label[i] == 0 ? 0.5w6.a[i] : w6.a[i] for i in 1:N], shape=[r == "echo" ? "triangle" : "circle" for r in o.role])
    imp.type = imp.shape
    fname = joinpath(res_folder, "detection_pixels__color255_webui__Impulse11sicw_sources_$(nm)_$params.csv")
    CSV.write(fname, vcat(imp, w6[N+1:end, names(imp)]))
    @printf("wrote %s\n", fname)
end
