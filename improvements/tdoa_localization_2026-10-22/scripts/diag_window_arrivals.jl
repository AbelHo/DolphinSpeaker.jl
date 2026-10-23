# Arrivals inside the TDOA window from raw WAV snippets (first file of the folder), and the earliest-peak test:
# does a window on the earliest envelope peak give a better direction than the window on the detection?
# usage: julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/diag_window_arrivals.jl <res_folder> [out_dir]
include(joinpath(@__DIR__, "..", "..", "..", "src", "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(@__DIR__, "common.jl"))
using WAV, DSP, Peaks
res_folder = ARGS[1]
out_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "results")
df, info = load_clicks(res_folder)
fs = info.fs; wavf = info.aufnames[1]
d1 = df[df.file .== 1, :]; N1 = nrow(d1)
rx = rx_vect[:, get_relevant_channels(rx_vect)]
hp = digitalfilter(Highpass(impulsive_band_pass[1] / fs * 2), Butterworth(butterworth_size))
pre, post = 400, 600; c0 = pre + 1
snippet(p) = (x = wavread(wavf; subrange=(p - pre):(p + post), format="native")[1];
    reduce(hcat, [filtfilt(hp, Float64.(x[:, ch])) for ch in get_relevant_channels(rx_vect)]))
function arrivals(y)            # envelope peaks >= 30% of the channel max in -150:+60, 8 samples apart
    map(1:3) do k
        w = (c0 - 150):(c0 + 60); e = abs.(hilbert(y[:, k]))[w]; pk = findmaxima(e, 8)
        isempty(pk.heights) && return (Int[], Float64[])
        keep = pk.heights .>= 0.3maximum(pk.heights)
        (pk.indices[keep] .+ (first(w) - c0), pk.heights[keep] ./ maximum(pk.heights))
    end
end
io = IOBuffer()

# A. every 15th click: arrival counts and flag rate by window content
sel = 1:15:N1
A = Vector{Any}(undef, length(sel))
Threads.@threads for m in eachindex(sel)
    i = sel[m]; ar = arrivals(snippet(d1.sample[i]))
    any(a -> isempty(a[1]), ar) && (A[m] = nothing; continue)
    a3 = ar[3][1]
    A[m] = (; n50=count(abs.(a3) .<= 50), nany=maximum(k -> count(abs.(ar[k][1]) .<= 50), 1:3), pre=any(a3 .< -8),
        lag=-minimum(a3), rel=ar[3][2][argmin(a3)], flag=d1.flag[i], iso=d1.gap_prev[i] >= 1.6)
end
A = filter(!isnothing, A)
n50 = getfield.(A, :n50); nany = getfield.(A, :nany); prec = getfield.(A, :pre); flg = getfield.(A, :flag); iso = getfield.(A, :iso)
@printf(io, "A. sampled clicks of %s: %d (every 15th)\n", basename(wavf), length(A))
println(io, "envelope peaks (>= 30% of channel max) inside -50:50, ref channel:")
for c in 0:4; @printf(io, "  %s: %5.1f%%\n", c < 4 ? string(c) : "4+", 100mean(c < 4 ? n50 .== c : n50 .>= 4)); end
@printf(io, ">= 2 peaks on some channel: %.1f%%\n", 100mean(nany .>= 2))
@printf(io, "preceded by an earlier peak (8-150 samples): %.1f%%; lag p25/p50/p75 %.0f/%.0f/%.0f samples; height med %.2f\n", 100mean(prec),
    quantile(getfield.(A[prec], :lag), [.25,.5,.75])..., median(getfield.(A[prec], :rel)))
for (nm, s) in [("1 peak on every channel", nany .== 1), (">= 2 peaks on some channel", nany .>= 2), ("not preceded", .!prec), ("preceded", prec)]
    @printf(io, "  %-28s n=%5d flagged %5.1f%%\n", nm, count(s), 100mean(flg[s]))
end
for (nm, s) in [("isolated", iso), ("follower", .!iso)]
    @printf(io, "  %-28s n=%5d >=2 peaks %5.1f%% preceded %5.1f%% flagged %5.1f%%\n", nm, count(s), 100mean(nany[s] .>= 2), 100mean(prec[s]), 100mean(flg[s]))
end

# B. three blocks of 1000 consecutive clicks: SRP/GCC on the detection window vs a window on the earliest peak
blocks = [10001:11000, 30001:31000, 50001:51000]
idxs = reduce(vcat, collect.(blocks))
R = Vector{Any}(undef, length(idxs))
Threads.@threads for m in eachindex(idxs)
    i = idxs[m]; y = snippet(d1.sample[i]); ar = arrivals(y)
    t1 = minimum(a -> isempty(a[1]) ? 0 : minimum(a[1]), ar)
    _, qd = get_tdoa_srp(y, [c0]; window=-50:50, ref_channel=3, fs=fs, rx_vect=rx, return_quality=true)
    _, gd = get_tdoa_gcc(y, [c0]; window=-50:50, ref_channel=3, fs=fs, rx_vect=rx, guide=:envelope, return_quality=true)
    _, qf = get_tdoa_srp(y, [c0 + t1]; window=-15:40, ref_channel=3, fs=fs, rx_vect=rx, return_quality=true)
    _, gf = get_tdoa_gcc(y, [c0 + t1]; window=-15:40, ref_channel=3, fs=fs, rx_vect=rx, guide=:envelope, return_quality=true)
    R[m] = (; uvd=qd.uv[1, :], uvf=qf.uv[1, :], clod=3gd.pair_rms[1], clof=3gf.pair_rms[1], iso=d1.gap_prev[i] >= 1.6, t1)
end
uv2v(uv) = (u = uv[1]; v = clamp(uv[2], -1, 1); (u, v, sqrt(max(1 - u^2 - v^2, 0.0))))
Vd = [uv2v(r.uvd) for r in R]; Vf = [uv2v(r.uvf) for r in R]
@printf(io, "\nB. %d consecutive clicks (3 blocks); earliest peak > 20 samples before detection: %.1f%%\n", length(R), 100mean([r.t1 < -20 for r in R]))
dd = [angd(Vd[k], Vf[k]) for k in eachindex(R)]
@printf(io, "  earliest-peak vs detection direction: median %.1f deg; closure<=1: detection %.0f%%, earliest %.0f%%\n",
    median(dd), 100mean([r.clod <= 1 for r in R]), 100mean([r.clof <= 1 for r in R]))
for (nm, VV) in [("detection window -50:50", Vd), ("earliest-peak window -15:40", Vf)]
    a = Float64[]
    for b in blocks
        ks = [m for m in eachindex(idxs) if idxs[m] in b && R[m].iso]
        for q in 1:length(ks)-1
            dt = (d1.sample[idxs[ks[q+1]]] - d1.sample[idxs[ks[q]]]) / fs * 1000
            1.6 <= dt < 30 && push!(a, angd(VV[ks[q]], VV[ks[q+1]]))
        end
    end
    @printf(io, "  consecutive isolated clicks, %-28s n=%4d median %5.1f deg <2deg %5.1f%%\n", nm, length(a), median(a), 100mean(a .< 2))
end
write_report(joinpath(out_dir, "diag_window_arrivals.md"), "Arrivals inside the TDOA window", String(take!(io));
    notes="Source: `$res_folder`, audio `$wavf`. Snippets are high-passed like the pipeline (Butterworth, `impulsive_band_pass`). Ref channel = 3rd hydrophone column.")
