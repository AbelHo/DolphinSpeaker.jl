# Figures for the round-4 report: simulator realism, direction accuracy per method, dolphin separation, and the
# real-data separation (direction vs time coloured by dolphin). Palette: dataviz reference slots, light mode.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/report_figures_round4.jl
using CSV, DataFrames, Statistics, Plots, Printf
gr()
root = joinpath(@__DIR__, "..")
res = joinpath(root, "results"); fig = joinpath(root, "reports", "figures"); mkpath(fig)
const SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
const INK, INK2, SURF, GRIDC, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#d9d8d4", "#a3a29e"
base(; kw...) = plot(; background_color=SURF, foreground_color_text=INK, foreground_color_guide=INK2, foreground_color_axis=INK2,
    foreground_color_border=INK2, gridcolor=GRIDC, gridalpha=0.6, framestyle=:axes, grid=:y, titlefontsize=11, guidefontsize=10,
    tickfontsize=9, legendfontsize=9, legend_background_color=SURF, legend_foreground_color=SURF, size=(860, 460),
    left_margin=5Plots.mm, bottom_margin=6Plots.mm, kw...)
function grouped!(p, cats, series; width=0.8)   # series: [(label, values, colour)]
    k = length(series); w = width / k
    for (j, (lab, v, col)) in enumerate(series)
        bar!(p, (1:length(cats)) .+ (j - (k + 1) / 2) * w, v; bar_width=w * 0.92, color=col, linewidth=0, label=lab)
    end
    plot!(p; xticks=(1:length(cats), cats))
end

# 1. realism: statistics of the simulated windows vs the real recording
cal = vcat(CSV.read(joinpath(res, "calibrate_realism.csv"), DataFrame), CSV.read(joinpath(res, "calibrate_realism_target.csv"), DataFrame))
row(a0, hs, rl) = cal[(cal.A0 .== a0) .& (cal.hyd_scatter .== hs) .& (cal.reverb_level .== rl), :][1, :]
dry = row(9000, 0.0, 0.0); calib = row(6000, 0.7, 0.1)
cats = [">= 2 SIC candidates", "closure <= 1 sample", "followers"]
p = base(; ylabel="share of detections (%)", ylims=(0, 110), legend=:topright,
    title="Hydrophone scattering makes the simulated windows as messy as the real ones")
grouped!(p, cats, [("simulation without scattering", [dry.ge2, dry.closure_le1, dry.followers], GREY),
    ("calibrated simulation (scatter 0.7)", [calib.ge2, calib.closure_le1, calib.followers], SLOT[1]),
    ("real New_device", [78.0, 35.1, 51.4], SLOT[2])])
savefig(p, joinpath(fig, "round4_realism.png"))

# 2. direction accuracy per method (realistic simulation)
det = CSV.read(joinpath(res, "data", "bench_sources_detections.csv"), DataFrame)
methods = ["rank 1", "longest track", "longest track + inherit all", "longest track + smart inherit"]
acc(set, m) = (d = det[(det.set .== set) .& det.matched .& (det.method .== m), :]; 100mean(d.err .< 3))
p = base(; ylabel="direction within 3° of the direct path (%)", ylims=(0, 70), legend=:topleft,
    title="Realistic simulation: direction within 3° of the direct path, by method")
grouped!(p, ["rank 1", "longest track", "+ inherit all", "+ smart inherit"],
    [("2 dolphins", [acc(2, m) for m in methods], SLOT[1]), ("3 dolphins", [acc(3, m) for m in methods], SLOT[3])])
savefig(p, joinpath(fig, "round4_accuracy.png"))

# 3. dolphin separation: pairwise F1, defaults vs tuned vs k-means
txt = read(joinpath(res, "tune_sources.md"), String)
function f1s(set)
    sec = split(txt, set == 2 ? "2 dolphins (tuning set)" : "3 dolphins (held out)")[2]
    get(lab) = parse(Float64, split(strip(split(split(sec, lab)[2], "\n")[1]))[1])
    return [get("defaults, n given"), get("tuned, n given"), get("tuned (for estimated), n estimated"), get("k-means baseline (k = n)")]
end
p = base(; ylabel="pairwise F1 (pairs <= 0.2 s apart)", ylims=(0.5, 1.0), legend=:topright,
    title="Separating dolphins: tuning on 2 dolphins also helps the held-out 3-dolphin case")
grouped!(p, ["defaults, n given", "tuned, n given", "tuned, n estimated", "k-means (k = n)"],
    [("2 dolphins (tuning set)", f1s(2), SLOT[1]), ("3 dolphins (held out)", f1s(3), SLOT[3])])
savefig(p, joinpath(fig, "round4_separation.png"))

# 4. real data: busiest 3 s of file 0001, own vs output direction coloured by dolphin (n = 2, tuned)
r = CSV.read(joinpath(res, "data", "real_sources_sicw_tuned_New_device.csv"), DataFrame)
f1 = r[r.file .== r.file[1], :]
t = f1.t_s; t0 = t[argmax([count(x -> s <= x < s + 3, t) for s in t[1:50:end]]) * 50 - 49]
w = (t0 .<= f1.t_s .< t0 + 3)
p = base(; xlabel="time (s)", ylabel="azimuth (°)", ylims=(-60, 60), legend=:outerbottom, legend_column=3, size=(860, 520),
    title=@sprintf("Real data, %s, %.0f–%.0f s: 2 dolphins (triangles: echoes)", f1.file[1], t0, t0 + 3))
for (lab, col, nm) in ((1, SLOT[1], "dolphin 1"), (2, SLOT[2], "dolphin 2"), (0, GREY, "unassigned"))
    s = w .& (f1.label_n .== lab) .& (f1.role_n .!= "echo"); e = w .& (f1.label_n .== lab) .& (f1.role_n .== "echo")
    scatter!(p, f1.t_s[s], f1.az_out[s]; ms=2.5, msw=0, color=col, label=nm, alpha=0.8)
    scatter!(p, f1.t_s[e], f1.az_out[e]; ms=2.5, msw=0, color=col, marker=:utriangle, label="", alpha=0.5)
end
savefig(p, joinpath(fig, "round4_real_sources.png"))
println("round-4 figures written to ", fig)
