# Figure for the round-5 report: hydrophone equalization on the realistic simulation (results/bench_eq_scenarios.csv).
# The real-data figures come from estimate_eq_real.jl and echo_structure_real.jl.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/report_figures_round5.jl
using CSV, DataFrames, Statistics, Plots
gr()
root = joinpath(@__DIR__, "..")
res = joinpath(root, "results"); fig = joinpath(root, "reports", "figures")
const SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
const INK, INK2, SURF, GRIDC, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#d9d8d4", "#a3a29e"
b = CSV.read(joinpath(res, "bench_eq_scenarios.csv"), DataFrame)
vers = [("none", "no equalization", GREY), ("aligned", "estimated (pure delays unknown)", SLOT[1]),
    ("blind + true delays", "blind estimate + true pure delays", SLOT[3]), ("oracle", "oracle (true responses)", SLOT[2])]
metrics = [("acc_all", "within 3° of the direct path"), ("closure_le1", "GCC closure <= 1 sample"), ("sep_f1", "dolphin separation F1 (x100)")]
p = plot(; background_color=SURF, foreground_color_text=INK, foreground_color_guide=INK2, foreground_color_axis=INK2, foreground_color_border=INK2,
    gridcolor=GRIDC, gridalpha=0.6, framestyle=:axes, grid=:y, ylims=(0, 110), ylabel="% (mean over 15 scenarios)", size=(900, 470),
    title="Simulation: equalizing restores coherence, but directions need the pure delays", titlefontsize=11, guidefontsize=10, tickfontsize=9,
    legend=:outerbottom, legend_column=2, legendfontsize=9, legend_background_color=SURF, legend_foreground_color=SURF, left_margin=5Plots.mm)
k = length(vers); w = 0.8 / k
for (j, (v, lab, col)) in enumerate(vers)
    g = b[b.version .== v, :]
    vals = [m == "sep_f1" ? 100mean(g[!, m]) : mean(g[!, m]) for (m, _) in metrics]
    bar!(p, (1:length(metrics)) .+ (j - (k + 1) / 2) * w, vals; bar_width=w * 0.92, color=col, linewidth=0, label=lab)
end
plot!(p; xticks=(1:length(metrics), last.(metrics)))
savefig(p, joinpath(fig, "round5_eq_sim.png"))
println("round-5 figure written to ", fig)
