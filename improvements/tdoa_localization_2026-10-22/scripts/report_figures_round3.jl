# Figures for the round-3 report: simulated accuracy per rule, and the real-data proxies per rule.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/report_figures_round3.jl
using CSV, DataFrames, Statistics, Plots
gr()
root = joinpath(@__DIR__, "..")
res = joinpath(root, "results"); fig = joinpath(root, "reports", "figures"); mkpath(fig)
const SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
const INK, INK2, SURF, GRIDC = "#0b0b0b", "#52514e", "#fcfcfb", "#d9d8d4"
base(; kw...) = plot(; background_color=SURF, foreground_color_text=INK, foreground_color_guide=INK2, foreground_color_axis=INK2,
    foreground_color_border=INK2, gridcolor=GRIDC, gridalpha=0.6, framestyle=:axes, grid=:y, titlefontsize=11, guidefontsize=10,
    tickfontsize=9, legendfontsize=9, legend_background_color=SURF, legend_foreground_color=SURF, size=(860, 440),
    left_margin=5Plots.mm, bottom_margin=8Plots.mm, kw...)
rules = ["rank 1", "earliest", "track lag", "longest track"]
x = 1:length(rules)

# 1. simulation: accuracy per rule, with and without followers inheriting (window -110:50, default tracker)
s = CSV.read(joinpath(res, "bench_sequences_summary.csv"), DataFrame)
acc(r) = s.acc_all[(s.window .== "sicw") .& (s.tracker .== "default") .& (s.rule .== r)][1]
oracle = s.acc_all[(s.window .== "sicw") .& (s.rule .== "oracle (best candidate)")][1]
p = base(; ylabel="direct path within 3° (%)", ylims=(0, 100), xticks=(x, rules), legend=:outerbottom, legend_column=3, size=(860, 500),
    title="Simulation (55 pools, SIC −110:50): followers inheriting the previous direction adds ~12 points")
bar!(p, x .- 0.18, acc.(rules); bar_width=0.34, color=SLOT[1], linewidth=0, label="per detection")
bar!(p, x .+ 0.18, [acc("$r + followers inherit") for r in rules]; bar_width=0.34, color=SLOT[2], linewidth=0, label="followers inherit")
hline!(p, [oracle]; color=INK2, ls=:dash, lw=1, label="oracle: best candidate")
savefig(p, joinpath(fig, "round3_sim_rules.png"))

# 2. real data: proxies per rule (window -110:50, default tracker)
r = CSV.read(joinpath(res, "real_rules_sicw_default_tracklag-inherit.csv"), DataFrame)
get(rule, col) = r[r.rule .== rule, col][1]
p = base(; ylabel="median direction change (°)", ylims=(0, 40), xlims=(0.6, 4.4), xticks=(x, rules), legend=:topright, right_margin=5Plots.mm,
    title="Real data (SIC −110:50): consecutive isolated clicks change least with the longest-track rule")
for (k, (col, nm)) in enumerate([(:burst, "1.6–2.2 ms (burst pulses)"), (:p2_10, "2.2–10 ms"), (:p10_30, "10–30 ms")])
    plot!(p, x, [get(ru, col) for ru in rules]; color=SLOT[k], lw=2, marker=:circle, ms=6, msw=0, label=nm)
end
savefig(p, joinpath(fig, "round3_real_rules.png"))
println("round-3 figures written to ", fig)
