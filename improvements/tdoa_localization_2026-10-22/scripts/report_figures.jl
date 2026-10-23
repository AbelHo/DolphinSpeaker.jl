# Figures for the round-1 report (reports/): from results/lag_histogram.csv, results/bench_multipath.csv and
# results/data/click_features_New_device.csv. Palette: dataviz reference categorical slots, light mode, fixed order.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/report_figures.jl
using CSV, DataFrames, Statistics, Plots
gr()
root = joinpath(@__DIR__, "..")
res = joinpath(root, "results"); fig = joinpath(root, "reports", "figures"); mkpath(fig)
const SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]      # blue, orange, aqua, yellow
const INK, INK2, SURF, GRIDC = "#0b0b0b", "#52514e", "#fcfcfb", "#d9d8d4"
base(; kw...) = plot(; background_color=SURF, foreground_color_text=INK, foreground_color_guide=INK2, foreground_color_axis=INK2,
    foreground_color_border=INK2, gridcolor=GRIDC, gridalpha=0.6, framestyle=:axes, grid=:y, titlefontsize=11, guidefontsize=10,
    tickfontsize=9, legendfontsize=9, legend_background_color=SURF, legend_foreground_color=SURF, size=(860, 440),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm, kw...)

# 1. all-pair lag histogram
lh = CSV.read(joinpath(res, "lag_histogram.csv"), DataFrame)
p = base(; xlabel="time between two detections (ms)", ylabel="pairs per 0.1 ms", legend=false, yformatter=y -> string(round(Int, y)),
    title="Pairs pile up at multiples of 2 ms (burst pulses) and just above 0.4 ms")
bar!(p, lh.lag_ms_lo .+ 0.05, lh.count; bar_width=0.08, color=SLOT[1], linewidth=0)
savefig(p, joinpath(fig, "lag_histogram.png"))

# 2. direction change between consecutive detections, by gap
f = CSV.read(joinpath(res, "data", "click_features_New_device.csv"), DataFrame)
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))
V = vec3.(f.az_1, f.el_1)
classes = [(0.4, 0.75), (0.75, 1.5), (1.5, 1.8), (1.8, 2.2), (2.2, 3.5), (3.5, 4.5), (4.5, 10), (10, 30), (30, 100)]
lab = ["0.4–0.75", "0.75–1.5", "1.5–1.8", "1.8–2.2", "2.2–3.5", "3.5–4.5", "4.5–10", "10–30", "30–100"]
all_med = Float64[]; clean_med = Float64[]
for (lo, hi) in classes
    g = [i for i in 2:nrow(f) if f.file[i] == f.file[i-1] && lo <= f.gap_prev_ms[i] < hi]
    push!(all_med, median([angd(V[i], V[i-1]) for i in g]))
    gu = [i for i in g if !f.impulse5_flagged[i] && !f.impulse5_flagged[i-1]]
    push!(clean_med, median([angd(V[i], V[i-1]) for i in gu]))
end
x = 1:length(classes)
p = base(; xlabel="time since the previous detection (ms)", ylabel="median direction change (°)", xticks=(x, lab), ylims=(0, 60),
    legend=:top, title="Detections 0.4–1.5 ms after another point elsewhere; 1.8–30 ms apart they agree")
hline!(p, [45.5]; color=INK2, ls=:dash, lw=1, label="random pairs (45.5°)")
plot!(p, x, all_med; color=SLOT[1], lw=2, marker=:circle, ms=6, msw=0, label="all consecutive pairs")
plot!(p, x, clean_med; color=SLOT[2], lw=2, marker=:circle, ms=6, msw=0, label="both pass the closure gate")
savefig(p, joinpath(fig, "direction_change_by_gap.png"))

# 3. closure flag rate by amplitude quartile and by isolation
q = quantile(f.ppeak, [0.25, 0.5, 0.75])
grp = [("Q1 weakest", f.ppeak .< q[1]), ("Q2", q[1] .<= f.ppeak .< q[2]), ("Q3", q[2] .<= f.ppeak .< q[3]), ("Q4 strongest", f.ppeak .>= q[3]),
    ("isolated", .!f.follower), ("follower", f.follower)]
p = base(; ylabel="flagged by closure / feasibility (%)", ylims=(0, 100), legend=false, xticks=(1:6, first.(grp)),
    title="Weak clicks and followers fail the closure gate most often")
vals = [100mean(f.impulse5_flagged[s]) for (_, s) in grp]
bar!(p, 1:6, vals; color=[fill(SLOT[1], 4); fill(SLOT[2], 2)], linewidth=0, bar_width=0.6)
annotate!(p, [(i, vals[i] + 4, text(string(round(vals[i]; digits=1), "%"), 9, INK)) for i in 1:6])
savefig(p, joinpath(fig, "flag_rate.png"))

# 4. benchmark: direct path within 2 deg by echo gain
b = CSV.read(joinpath(res, "bench_multipath.csv"), DataFrame)
gains = [0.5, 0.8, 1.0, 1.25]
ser = [("SRP single estimate (Impulse6)", r -> r.err_srp < 2), ("GCC guided (Impulse5)", r -> r.err_gcc < 2),
    ("any of the 3 SRP candidates", r -> r.recall3), ("earliest beam arrival", r -> r.earliest_ok)]
p = base(; xlabel="echo amplitude relative to the direct path", ylabel="direct path found within 2° (%)", ylims=(0, 100),
    xticks=gains, legend=:bottomleft, title="Benchmark: once the echo is as strong as the direct path, single estimates fail")
for (k, (nm, fn)) in enumerate(ser)
    y = [100mean(fn.(eachrow(b[b.gain .== g, :]))) for g in gains]
    plot!(p, gains, y; color=SLOT[k], lw=2, marker=:circle, ms=6, msw=0, label=nm)
end
savefig(p, joinpath(fig, "bench_by_gain.png"))
println("figures written to ", fig)
