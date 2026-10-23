# Effect of the hydrophone equalizers on the real New_device clicks (before a full rerun): within-window statistics of
# the strongest isolated clicks of each file, raw vs blind vs aligned equalization (the pooled estimates of
# estimate_eq_real.jl), and the shift of their SRP directions (a pure-delay error would shift all directions alike).
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/eq_effect_real.jl [per_file]
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "hydrophone_eq.jl"))
using Printf, JLD2
per_file = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1500
res = joinpath(@__DIR__, "..", "results")
wavdir = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
fs = 192000; rx = rx_vect[:, get_relevant_channels(rx_vect)]
eqs = load(joinpath(res, "data", "hydrophone_eq_New_device.jld2"))
feat = CSV.read(joinpath(res, "data", "click_features_sicw_New_device.csv"), DataFrame)
half = 400
# all isolated clicks (not only the strong ones used for the estimate), sampled evenly
g = feat[.!feat.follower, :]
g = g[round.(Int, range(1, nrow(g); length=min(3per_file, nrow(g)))), :]
parts = Vector{Matrix{Float64}}(undef, nrow(g))
Threads.@threads for i in 1:nrow(g)
    parts[i] = read_filtered_block(joinpath(wavdir, g.file[i]), (g.sample[i] - half):(g.sample[i] + half); fs=fs)
end
data = reduce(vcat, parts); pk = [(i - 1) * (2half + 1) + half + 1 for i in 1:nrow(g)]
vec3(u, v) = (u, v, sqrt(max(1 - u^2 - v^2, 0.0))); angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))
io = IOBuffer()
@printf(io, "%d isolated clicks sampled evenly over New_device (all amplitudes), SIC -50:50, min_rel 0.7\n\n", nrow(g))
@printf(io, "%-9s %9s %11s %9s %11s %9s %10s\n", "version", ">=2 cand", "closure<=1", "ncc med", "SRP score", "efrac_1", "shift°")
c0 = nothing
for (nm, d) in (("raw", data), ("blind", apply_hydrophone_eq(data, eqs["blind"])), ("aligned", apply_hydrophone_eq(data, eqs["aligned"])))
    c = get_srp_candidates_sic(d, pk; window=window_impulsive, fs=fs, rx_vect=rx, min_rel=0.7)
    isnothing(c0) && (global c0 = c)
    shift = median([angd(vec3(c.uv[i, 1, :]...), vec3(c0.uv[i, 1, :]...)) for i in eachindex(pk)])
    @printf(io, "%-9s %8.1f%% %10.1f%% %9.2f %11.3f %9.2f %10.2f\n", nm, 100mean(c.nfound .>= 2), 100mean(c.closure .<= 1), median(c.ncc_min),
        median(c.rscore[:, 1]), median(c.efrac[:, 1]), shift)
end
body = String(take!(io))
open(joinpath(res, "eq_effect_real.md"), "w") do f; println(f, "# Effect of hydrophone equalization on real clicks\n\n```\n", body, "```"); end
print(body)
