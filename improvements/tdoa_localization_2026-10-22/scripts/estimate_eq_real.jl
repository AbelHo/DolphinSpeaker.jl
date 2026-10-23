# Estimate the hydrophone equalizers on the real New_device recording, per file and pooled, from the strongest isolated
# clicks (upper quarter of the amplitudes, no followers; up to `per_file` per file), and check that the files agree (a
# fixed mount/housing response should give the same estimate everywhere). Saves the pooled equalizers for run_multiarrival2.
# needs: results/data/click_features_sicw_New_device.csv (round 2)
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/estimate_eq_real.jl [per_file]
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl")); include(joinpath(src_dir, "hydrophone_eq.jl"))
using Printf, JLD2, Plots
gr()
per_file = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4000
res = joinpath(@__DIR__, "..", "results"); fig = joinpath(@__DIR__, "..", "reports", "figures")
wavdir = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
fs = 192000; rx = rx_vect[:, get_relevant_channels(rx_vect)]
feat = CSV.read(joinpath(res, "data", "click_features_sicw_New_device.csv"), DataFrame)
half = 300
snips = Dict{String,Matrix{Float64}}(); peaks = Dict{String,Vector{Int}}()
for fl in unique(feat.file)
    g = feat[(feat.file .== fl) .& .!feat.follower, :]
    g = g[g.ppeak .>= quantile(g.ppeak, 0.75), :]
    g = g[round.(Int, range(1, nrow(g); length=min(per_file, nrow(g)))), :]
    parts = Vector{Matrix{Float64}}(undef, nrow(g))
    Threads.@threads for i in 1:nrow(g)
        parts[i] = read_filtered_block(joinpath(wavdir, fl), (g.sample[i] - half):(g.sample[i] + half); fs=fs)
    end
    snips[fl] = reduce(vcat, parts); peaks[fl] = [(i - 1) * (2half + 1) + half + 1 for i in 1:nrow(g)]
    @info "$fl: $(nrow(g)) strong isolated clicks"
end
files = sort(collect(keys(snips)))
pooled = reduce(vcat, [snips[f] for f in files])
offs = cumsum([0; [size(snips[f], 1) for f in files[1:end-1]]])
ppool = reduce(vcat, [peaks[f] .+ o for (f, o) in zip(files, offs)])

est = Dict{Tuple{String,String},Any}()
for (nm, fn) in (("blind", estimate_hydrophone_eq_blind), ("aligned", estimate_hydrophone_eq))
    for f in files; est[(nm, f)] = fn(snips[f], peaks[f]; fs=fs, rx_vect=rx, ref=ref_channel); end
    est[(nm, "pooled")] = fn(pooled, ppool; fs=fs, rx_vect=rx, ref=ref_channel)
end

# agreement between files: dispersive phase (after removing a linear trend) and magnitude, coherence-weighted band
function compare(a, b, k)
    ib = (a.f .>= a.band[1]) .& (a.f .<= a.band[2]) .& (a.coh[:, k] .>= 0.3) .& (b.coh[:, k] .>= 0.3)
    count(ib) < 4 && return (NaN, NaN)
    r = a.G[ib, k] ./ b.G[ib, k]
    φ = unwrap(angle.(r)); X = [ones(count(ib)) a.f[ib]]; ab = X \ φ
    return median(abs.(φ .- X * ab)), median(abs.(log.(abs.(r) ./ median(abs.(r)))))
end
io = IOBuffer()
@printf(io, "New_device: strongest isolated clicks per file: %s (snippets ±%d samples, filtered as the pipeline)\n", join(["$(f[end-7:end-4]) $(length(peaks[f]))" for f in files], ", "), half)
for nm in ("blind", "aligned"), k in 1:3
    k == ref_channel && continue
    @printf(io, "\n%-8s hydrophone %d vs reference %d:\n", nm, k, ref_channel)
    for f in [files; "pooled"]
        e = est[(nm, f)]
        ib = (e.f .>= e.band[1]) .& (e.f .<= e.band[2])
        @printf(io, "  %-22s coherence median %.2f, bins corrected %3.0f%%, |G| range p10-p90 %.2f-%.2f%s\n", f, median(e.coh[ib, k]),
            100mean(e.E[ib, k] .!= 1), quantile(abs.(e.G[ib, k]), [0.1, 0.9])..., nm == "aligned" ? @sprintf(", pure delay removed %.2f samples", e.delay_removed[k]) : "")
    end
    for i in eachindex(files), j in i+1:length(files)
        ph, mg = compare(est[(nm, files[i])], est[(nm, files[j])], k)
        @printf(io, "  %s vs %s: dispersive phase differs by %.2f rad (median), magnitude by %.0f%%\n", files[i][end-7:end-4], files[j][end-7:end-4], ph, 100(exp(mg) - 1))
    end
end
jldsave(joinpath(res, "data", "hydrophone_eq_New_device.jld2"); blind=est[("blind", "pooled")], aligned=est[("aligned", "pooled")],
    per_file=Dict(string(k) => v for (k, v) in est))
body = String(take!(io))
open(joinpath(res, "eq_real_New_device.md"), "w") do g; println(g, "# Hydrophone equalizers estimated on New_device\n\n```\n", body, "```"); end
print(body)

# figure: dispersive phase and magnitude of the blind estimate per file, hydrophones 1 and 2 relative to the reference
const SLOT = ["#2a78d6", "#eb6834", "#1baf7a"]
plt = []
for k in (1, 2)
    pa = plot(; ylabel="phase (rad), trend removed", title="hydrophone $k vs reference $(ref_channel)", legend=k == 1 ? :bottomleft : false,
        background_color="#fcfcfb", gridalpha=0.4, titlefontsize=10, guidefontsize=9, tickfontsize=8, legendfontsize=8)
    pm = plot(; ylabel="|G| (median 1)", xlabel="frequency (kHz)", legend=false, background_color="#fcfcfb", gridalpha=0.4, guidefontsize=9, tickfontsize=8)
    for (i, f) in enumerate(files)
        e = est[("blind", f)]; ib = (e.f .>= e.band[1]) .& (e.f .<= e.band[2])
        φ = unwrap(angle.(e.G[ib, k])); X = [ones(count(ib)) e.f[ib]]; φ .-= X * (X \ φ)
        plot!(pa, e.f[ib] ./ 1000, φ; color=SLOT[i], lw=1.5, label="file $(f[end-7:end-4])")
        plot!(pm, e.f[ib] ./ 1000, abs.(e.G[ib, k]) ./ median(abs.(e.G[ib, k])); color=SLOT[i], lw=1.5)
    end
    push!(plt, pa, pm)
end
savefig(plot(plt[1], plt[3], plt[2], plt[4]; layout=(2, 2), size=(980, 620), left_margin=4Plots.mm, bottom_margin=4Plots.mm),
    joinpath(fig, "round5_eq_real_files.png"))
