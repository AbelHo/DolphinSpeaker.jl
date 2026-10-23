# Figures for the round-2 report (reports/round2_sic_selection.md): benchmark v2 trade-off and selection, and the
# real-data consecutive-click agreement / direction-vs-time views. Palette: dataviz reference slots, light mode.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/report_figures_round2.jl
using CSV, DataFrames, Statistics, Plots, Printf
gr()
root = joinpath(@__DIR__, "..")
res = joinpath(root, "results"); fig = joinpath(root, "reports", "figures"); mkpath(fig)
const SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]      # blue, orange, aqua, yellow
const INK, INK2, SURF, GRIDC = "#0b0b0b", "#52514e", "#fcfcfb", "#d9d8d4"
base(; kw...) = plot(; background_color=SURF, foreground_color_text=INK, foreground_color_guide=INK2, foreground_color_axis=INK2,
    foreground_color_border=INK2, gridcolor=GRIDC, gridalpha=0.6, framestyle=:axes, grid=:y, titlefontsize=11, guidefontsize=10,
    tickfontsize=9, legendfontsize=9, legend_background_color=SURF, legend_foreground_color=SURF, size=(860, 440),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm, kw...)

# A. recall of the direct path vs false candidates, per generator, over the acceptance threshold
sw = CSV.read(joinpath(res, "bench_multipath2_sweep.csv"), DataFrame)
p = base(; xlabel="false candidates per echo-free click", ylabel="direct path among the candidates (%), echo trials",
    ylims=(50, 90), legend=:bottomright, grid=:both, title="Successive cancellation keeps more direct paths for the same false-candidate rate")
for (k, (m, nm)) in enumerate([("plain", "round 1: local maxima"), ("clean", "CLEAN on the GCC curves"), ("sic", "SIC, window −50:50"), ("sic_wide", "SIC, window −110:50")])
    d = sort(sw[(sw.method .== m) .& isinf.(sw.support_min) .& (sw.rel_min .< 2), :], :rel_min)
    plot!(p, d.false_noecho, d.recall; color=SLOT[k], lw=2, marker=:circle, ms=5, msw=0, label=nm)
end
for m in ("sic", "sic_wide")   # mark the chosen threshold 0.7
    d = sw[(sw.method .== m) .& (sw.rel_min .== 0.7), :]
    annotate!(p, d.false_noecho[1] + 0.05, d.recall[1] + 1.2, text("0.7", 8, INK2, :left))
end
savefig(p, joinpath(fig, "round2_tradeoff.png"))

# B. one direction per click: share within 2 deg of the direct path by echo amplitude
cand = CSV.read(joinpath(res, "bench_multipath2_candidates.csv"), DataFrame)
function by_gain(meth, rule; thr=0.7)
    d = cand[(cand.method .== meth) .& ((cand.k .== 1) .| (cand.rel .>= thr)), :]
    g = combine(groupby(d, :trial), :gain => first => :gain, :nech => first => :nech,
        [:k, :err_direct] => ((k, e) -> e[k .== 1][1] < 2) => :rank1, [:tarr, :err_direct] => ((t, e) -> e[argmin(t)] < 2) => :earliest)
    return [100mean(g[(g.nech .== 0) .| (g.gain .== 0), rule]); [100mean(g[(g.nech .> 0) .& (g.gain .== x), rule]) for x in [0.5, 0.8, 1.0, 1.25]]]
end
labs = ["no echo", "0.5", "0.8", "1.0", "1.25"]; x = 1:5
p = base(; xlabel="echo amplitude relative to the direct path", ylabel="selected direction = direct path within 2° (%)",
    xticks=(x, labs), ylims=(0, 100), legend=:bottomleft, title="Earliest SIC arrival vs the single SRP estimate")
for (k, (nm, meth, rule)) in enumerate([("single SRP estimate (Impulse6)", "sic", :rank1), ("SIC, earliest arrival, window −50:50", "sic", :earliest),
        ("SIC, earliest arrival, window −110:50", "sic_wide", :earliest)])
    plot!(p, x, by_gain(meth, rule); color=SLOT[k], lw=2, marker=:circle, ms=6, msw=0, label=nm)
end
savefig(p, joinpath(fig, "round2_selection_by_gain.png"))

# C/D. real data (only when the round-2 runs exist)
fs_ = Dict(t => joinpath(res, "data", "click_features_$(t)_New_device.csv") for t in ("sic", "sicw"))
if all(isfile, values(fs_))
    F = Dict(t => CSV.read(f, DataFrame) for (t, f) in fs_)
    f6 = CSV.read(joinpath(res, "data", "click_features_New_device.csv"), DataFrame)   # round 1: az_1/el_1 = Impulse6
    vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
    angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))
    classes = [(1.6, 2.2), (2.2, 10), (10, 30)]; lab = ["1.6–2.2 (burst pulses)", "2.2–10", "10–30"]   # isolated: gap >= 1.6 ms
    function med_change(az, el, file, gap, iso)
        [median([angd(vec3(az[i], el[i]), vec3(az[i-1], el[i-1])) for i in 2:length(az)
            if file[i] == file[i-1] && iso[i] && iso[i-1] && lo <= gap[i] < hi]) for (lo, hi) in classes]
    end
    iso = .!f6.follower
    fw = F["sicw"]
    trk(c) = [c[i, "az_$(c.selected_track[i])"] for i in 1:nrow(c)], [c[i, "el_$(c.selected_track[i])"] for i in 1:nrow(c)]
    azt, elt = trk(fw)
    ser = [("single SRP estimate (Impulse6)", f6.az_1, f6.el_1), ("SIC −110:50, earliest arrival per click", fw.az_sel, fw.el_sel),
        ("SIC −110:50, track lag rule", azt, elt)]
    p = base(; xlabel="time between consecutive isolated clicks (ms)", ylabel="median direction change (°)", xticks=(1:3, lab),
        legend=:topleft, ylims=(0, 50), title="Real data: per-click earliest arrival is noisier, the track rule steadier")
    for (k, (nm, az, el)) in enumerate(ser)
        plot!(p, 1:3, med_change(az, el, f6.file, f6.gap_prev_ms, iso); color=SLOT[k], lw=2, marker=:circle, ms=6, msw=0, label=nm)
    end
    savefig(p, joinpath(fig, "round2_real_agreement.png"))

    # direction vs time in the busiest 3 s of file 0001: Impulse6 vs the wide-window selection
    t = f6.t_s[f6.file .== f6.file[1]]
    t0 = t[argmax([count(x -> s <= x < s + 3, t) for s in t[1:50:end]]) * 50 - 49]
    sel = (f6.file .== f6.file[1]) .& (t0 .<= f6.t_s .< t0 + 3)
    pa = base(; ylabel="azimuth (°)", ylims=(-90, 90), legend=false, title=@sprintf("Single SRP estimate (Impulse6), %s, %.0f–%.0f s", f6.file[1], t0, t0 + 3))
    scatter!(pa, f6.t_s[sel], f6.az_1[sel]; ms=2, msw=0, color=SLOT[1], alpha=0.6)
    pb = base(; xlabel="time (s)", ylabel="azimuth (°)", ylims=(-90, 90), legend=false, title="Selected direction (SIC −110:50, track lag rule)")
    scatter!(pb, fw.t_s[sel], azt[sel]; ms=2, msw=0, color=SLOT[3], alpha=0.6)
    savefig(plot(pa, pb; layout=(2, 1), size=(860, 640)), joinpath(fig, "round2_az_vs_time.png"))
end
println("round-2 figures written to ", fig)
