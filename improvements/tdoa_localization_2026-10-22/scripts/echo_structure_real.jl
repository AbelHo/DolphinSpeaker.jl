# What makes the real windows incoherent, if not a fixed hydrophone response? Look for in-window echoes whose delay depends
# on the arrival direction (e.g. a reflecting plate or wall behind the array at distance d: delay 2 d cos θ / c for a
# source θ off the camera axis). For strong isolated clicks of New_device, the real cepstrum and the envelope
# autocorrelation of each hydrophone (lags 3-100 samples after the click) are averaged in bins of the SIC direction.
# usage: julia --project=. -t 24 improvements/tdoa_localization_2026-10-22/scripts/echo_structure_real.jl [n_clicks]
src_dir = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src_dir, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src_dir, "multiarrival.jl"))
using Printf, Plots
gr()
n_clicks = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6000
res = joinpath(@__DIR__, "..", "results"); fig = joinpath(@__DIR__, "..", "reports", "figures")
wavdir = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
fs = 192000; rx = rx_vect[:, get_relevant_channels(rx_vect)]
feat = CSV.read(joinpath(res, "data", "click_features_sicw_New_device.csv"), DataFrame)
g = feat[.!feat.follower, :]
g = g[g.ppeak .>= quantile(g.ppeak, 0.5), :]                      # the stronger half
g = g[round.(Int, range(1, nrow(g); length=min(n_clicks, nrow(g)))), :]
half = 300
snip = Vector{Matrix{Float64}}(undef, nrow(g))
Threads.@threads for i in 1:nrow(g)
    snip[i] = read_filtered_block(joinpath(wavdir, g.file[i]), (g.sample[i] - half):(g.sample[i] + half); fs=fs)
end
# direction of each click: SRP rank 1 on the pipeline window
c = get_srp_candidates_sic(reduce(vcat, snip), [(i - 1) * (2half + 1) + half + 1 for i in eachindex(snip)];
    window=window_impulsive, fs=fs, rx_vect=rx, min_rel=0.7, K=1)
u, v = c.uv[:, 1, 1], c.uv[:, 1, 2]
θ = acosd.(clamp.(sqrt.(max.(1 .- u .^ 2 .- v .^ 2, 0.0)), -1, 1))   # angle from the camera axis (front hemisphere)
el = asind.(clamp.(v, -1, 1))
maxlag = 100; N = 512
w = 30   # click window: the click starts ~10 samples before the detected peak; echoes up to maxlag after
function cepstrum(x)
    X = rfft([x .* tukey(length(x), 0.1); zeros(N - length(x))])
    irfft(log.(abs.(X) .+ 1e-6 * maximum(abs.(X))), N)[1:maxlag+1]
end
function env_acorr(x)
    e = abs.(hilbert(x)); e .-= mean(e)
    a = [sum(e[1:end-l] .* e[1+l:end]) for l in 0:maxlag]
    a ./ a[1]
end
bins = [(0, 10), (10, 20), (20, 30), (30, 45), (45, 60), (60, 90)]
io = IOBuffer()
@printf(io, "%d strong isolated clicks of New_device; cepstrum / envelope autocorrelation of each hydrophone, window -%d:+%d samples around the detected peak\n", nrow(g), w, maxlag + 20)
qmin = 8   # below ~8 samples (42 µs) the cepstrum shows the click's own spectral shape (peaks at 3-6 samples in every bin)
@printf(io, "peaks (lag in samples and µs) of the average real cepstrum (quefrency %d-%d), per off-axis bin and hydrophone:\n", qmin, maxlag)
cep_img = zeros(maxlag + 1, length(bins)); ac_img = zeros(maxlag + 1, length(bins))
for (bi, (lo, hi)) in enumerate(bins)
    ix = findall((θ .>= lo) .& (θ .< hi))
    length(ix) < 30 && continue
    for ch in 1:3
        C = mean(cepstrum(snip[i][half+1-w:half+1+maxlag+20, ch]) for i in ix)
        A = mean(env_acorr(snip[i][half+1-w:half+1+maxlag+20, ch]) for i in ix)
        ch == ref_channel && (cep_img[:, bi] = C; ac_img[:, bi] = A)
        q = qmin+1:maxlag+1   # lags qmin..maxlag
        top = q[sortperm(abs.(C[q]); rev=true)[1:3]] .- 1
        @printf(io, "  θ %2d-%2d° (n %4d) hydrophone %d: cepstral peaks at %s samples (%s µs), values %s\n", lo, hi, length(ix), ch,
            join(top, ", "), join(round.(Int, top ./ fs .* 1e6), ", "), join(round.(C[top .+ 1]; digits=3), ", "))
    end
end
# elevation dependence (horizontal surfaces: surface / bottom reflections depend on elevation and range, not on the array)
@printf(io, "\nreference hydrophone, by elevation:\n")
for (lo, hi) in [(-90, -20), (-20, -5), (-5, 5), (5, 20), (20, 90)]
    ix = findall((el .>= lo) .& (el .< hi))
    length(ix) < 30 && continue
    C = mean(cepstrum(snip[i][half+1-w:half+1+maxlag+20, ref_channel]) for i in ix)
    q = qmin+1:maxlag+1; top = q[sortperm(abs.(C[q]); rev=true)[1:3]] .- 1
    @printf(io, "  el %3d to %3d° (n %4d): cepstral peaks at %s samples, values %s\n", lo, hi, length(ix), join(top, ", "), join(round.(C[top .+ 1]; digits=3), ", "))
end
body = String(take!(io))
open(joinpath(res, "echo_structure_real.md"), "w") do f; println(f, "# In-window echo structure of real clicks (New_device)\n\n```\n", body, "```"); end
print(body)

labels = ["$(b[1])-$(b[2])°" for b in bins]
p1 = Plots.heatmap(1:length(bins), 0:maxlag, cep_img; xticks=(1:length(bins), labels), xlabel="angle from the camera axis", ylabel="quefrency (samples)",
    title="average real cepstrum, reference hydrophone", color=:balance, clims=(-maximum(abs, cep_img[qmin+1:end, :]), maximum(abs, cep_img[qmin+1:end, :])),
    ylims=(qmin, maxlag), titlefontsize=10, guidefontsize=9, tickfontsize=8)
p2 = Plots.heatmap(1:length(bins), 0:maxlag, ac_img; xticks=(1:length(bins), labels), xlabel="angle from the camera axis", ylabel="lag (samples)",
    title="average envelope autocorrelation", color=:viridis, ylims=(qmin, maxlag), clims=(-0.2, 0.6), titlefontsize=10, guidefontsize=9, tickfontsize=8)
savefig(plot(p1, p2; layout=(1, 2), size=(1000, 460), left_margin=4Plots.mm, bottom_margin=5Plots.mm), joinpath(fig, "round5_echo_structure_real.png"))
