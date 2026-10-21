# Why are pairwise TDOAs inconsistent on real clicks? Per-channel SNR (all 5 channels), pairwise waveform
# similarity, and example clicks.
include(joinpath(@__DIR__, "..", "..", "..", "src", "compare_tdoa.jl"))
using Random
out = ARGS[1]
aufname = joinpath(COMPARE_TDOA_FOLDER, "251006_001_0001.WAV")
data, fs, _, _, _ = readAudio(aufname)
res_impulse = detect_impulse((aufname, data, fs), nothing; ref_channel=ref_channel, return_datafilt=true)
res_train = detect_impulsetrain(res_impulse, nothing)
all5 = filter_simple(Float64.(data), impulsive_band_pass; fs=fs)   # same band as the detector, all 5 channels
data = nothing; GC.gc()
d = res_impulse.data_filt
rel = get_relevant_channels(rx_vect)
@printf("relevant channels (rx_vect columns without NaN): %s; ref_channel %d of them = physical channel %d\n", rel, ref_channel, rel[ref_channel])
@printf("check: data_filt equals filtered physical channels %s: %s\n", rel, maximum(abs, d[1:10^6, :] .- all5[1:10^6, rel]) < 1e-6 * maximum(abs, d[1:10^6, :]))

W = window_impulsive
pk = res_train.pind_good
n = size(d, 1)
Random.seed!(1)
det = sort(res_impulse.pind_good)
clickfree(q) = (i = searchsortedfirst(det, q - 4000); i > length(det) || det[i] > q + 4000)
noise_pos = filter(clickfree, rand(5000:n-5000, 20000))[1:2000]
En = [mean(sum(abs2, all5[q .+ W, ch]) for q in noise_pos) for ch in 1:5]
snr(p, ch) = 10log10(max(sum(abs2, all5[p .+ W, ch]) / En[ch] - 1, 1e-3))
S = [snr(p, ch) for p in pk, ch in 1:5]
println("\nwindow SNR (dB) at the detected clicks, per physical channel (p10 / median / p90):")
for ch in 1:5
    @printf("  ch%d%s: %6.1f %6.1f %6.1f\n", ch, ch in rel ? " (used)" : "       ", quantile(S[:, ch], [0.1, 0.5, 0.9])...)
end

# pairwise similarity: normalized cross-correlation peak within the physical lag limit
rx = rx_vect[:, rel]
τmax = tau_max_samples(rx, fs)
pairs = [(1, 2), (1, 3), (2, 3)]
function ncc(x, y, maxlag)
    x = x .- mean(x); y = y .- mean(y)
    s = xcorr(x, y; padmode=:none)
    c = length(y)
    lags = -round(Int, maxlag):round(Int, maxlag)
    return maximum(abs, s[c .+ lags]) / sqrt(sum(abs2, x) * sum(abs2, y))
end
sub = pk[1:10:end]
C = [ncc(d[p .+ W, j], d[p .+ W, k], τmax[j, k]) for p in sub, (j, k) in pairs]
println("\npeak normalized cross-correlation within ±τmax between used channels (p10 / median / p90; 1 = identical waveforms):")
for (c, (j, k)) in enumerate(pairs)
    @printf("  ch%d-ch%d: %.2f %.2f %.2f\n", rel[j], rel[k], quantile(C[:, c], [0.1, 0.5, 0.9])...)
end

# envelope peak position of each used channel relative to the window centre (the ref channel's peak)
offs = [argmax(abs.(hilbert(d[p .+ W, ch]))) - (1 - first(W)) for p in sub, ch in 1:3]
println("\nenvelope-peak offset from the ref-channel peak (samples), share within the physical limit ±22.4:")
for ch in 1:3
    @printf("  ch%d: median |offset| %.0f, within ±23: %.1f%%\n", rel[ch], median(abs.(offs[:, ch])), 100 * mean(abs.(offs[:, ch]) .<= 23))
end

# closure of gcc phat on the same clicks
_, q = get_tdoa_gcc(d, sub; window=W, fs=fs, rx_vect=rx, return_quality=true)
clo = 3 .* q.pair_rms
minC = vec(minimum(C; dims=2))
println("\ngcc phat closure (samples) vs the weakest pair similarity:")
for (lo, hi) in [(0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
    sel = lo .<= minC .< hi
    any(sel) && @printf("  min ncc %.1f-%.1f: %5d clicks, median closure %.2f\n", lo, hi, count(sel), median(clo[sel]))
end

# example clicks: 3 consistent (closure < 0.3) and 3 inconsistent (closure > 5), strongest of each kind
gr()
ord = sortperm(res_train.ppeak[1:10:end]; rev=true)
good = [i for i in ord if clo[i] < 0.3][1:3]
bad = [i for i in ord if clo[i] > 5][1:3]
chcol = ["#2a78d6", "#eb6834", "#1baf7a"]
plt = plot(; layout=(2, 3), size=(1600, 700), plot_style()..., left_margin=6Plots.mm, bottom_margin=5Plots.mm)
for (k, i) in enumerate([good; bad])
    p = sub[i]
    for ch in 1:3
        plot!(plt[k], collect(W) ./ fs .* 1e6, d[p .+ W, ch]; color=chcol[ch], linewidth=1.5, label=k == 1 ? "ch$(rel[ch])" : "",
            title=@sprintf("%s: t=%.3f s, closure %.1f smp, min ncc %.2f", k <= 3 ? "consistent" : "inconsistent", (p - 1) / fs, clo[i], minC[i]),
            xlabel=k > 3 ? "time from ref-channel peak (µs)" : "")
    end
end
savefig(plt, joinpath(out, "example_clicks.png"))
println("\nsaved ", joinpath(out, "example_clicks.png"))
