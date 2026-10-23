# Shared generator for the semi-synthetic multipath benchmarks (round 2): clean real clicks of 251006_001_0001 as
# templates, real noise from quiet stretches, plane-wave arrivals with fractional TDOAs for the D3 array.
# Needs: include run_example.jl + set_device__ophk_acoustic_D3() + multiarrival.jl, and the round-1 feature table.
using Random, CSV, DataFrames, Statistics
const BENCH_FS = 192000
const BENCH_RX = rx_vect[:, get_relevant_channels(rx_vect)]
const BENCH_WAV = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device/251006_001_0001.WAV"
const BENCH_C0 = 500                     # sample of the direct arrival in a 1000-sample window
const BENCH_SNR_DB = 19.0                # median real-click window SNR (2026-10-21 report)

# 300 templates (ref channel, -40:60 around the click, peak 1) and 300 noise stretches of 1000 samples
function bench_material(feature_csv; n=300, seed=11)
    feat = CSV.read(feature_csv, DataFrame)
    f1 = feat[feat.file .== basename(BENCH_WAV), :]
    clean = f1[(.!f1.follower) .& (f1.closure_phat .<= 0.3) .& (f1.nfound .== 1) .& (f1.score_1 .>= quantile(f1.score_1, 0.75)), :]
    Random.seed!(seed)
    tsel = clean.sample[randperm(nrow(clean))[1:min(n, nrow(clean))]]
    blk(p, a, b) = read_filtered_block(BENCH_WAV, (p + a):(p + b); fs=BENCH_FS)
    templates = [(x = blk(p, -400, 400)[401-40:401+60, 3]; x ./ maximum(abs, x)) for p in tsel]
    allp = sort(f1.sample); noises = Matrix{Float64}[]
    while length(noises) < n
        p = rand(allp[1]:allp[end]); i = searchsortedfirst(allp, p)
        minimum(abs.(allp[max(i - 1, 1):min(i, length(allp))] .- p)) > 0.02BENCH_FS && push!(noises, blk(p, 0, 999))
    end
    return templates, noises
end

bench_tdoa(az, el) = (t = vec(-[sind(az)*cosd(el) sind(el) cosd(az)*cosd(el)] * BENCH_RX ./ default_c .* BENCH_FS); t .- t[ref_channel])
bench_uv(az, el) = (sind(az) * cosd(el), sind(el))

# window with the given arrivals (az, el, delay in samples, signed gain) of template s, plus scaled noise
function bench_window(s, nz, arrivals; snr_db=BENCH_SNR_DB)
    data = zeros(1000, 3)
    for (az, el, d, g) in arrivals, ch in 1:3
        seg = zeros(1000); seg[BENCH_C0-40+d:BENCH_C0+60+d] .= g .* s
        data[:, ch] .+= _frac_delay(seg, bench_tdoa(az, el)[ch])
    end
    return data .+ nz ./ sqrt(mean(abs2, nz)) .* 10^(-snr_db / 20)
end

# detection point as in the pipeline: the strongest ref-channel envelope peak
bench_peak(data) = (e3 = abs.(hilbert(data[:, 3])); BENCH_C0 - 20 + argmax(e3[BENCH_C0-20:BENCH_C0+120]) - 1)

# random direct direction in the field of view plus nech echoes (delay 10-58 samples, gain 0.5-1.25, random polarity),
# from a vertical image (same azimuth) or a random azimuth
function bench_arrivals(rng, nech; delays=[10, 20, 40, 58], gains=[0.5, 0.8, 1.0, 1.25])
    az0, el0 = 70rand(rng) - 35, 44rand(rng) - 22
    delay = rand(rng, delays); gain = rand(rng, gains); kind = rand(rng, ["vertical image", "random azimuth"])
    arr = [(az0, el0, 0, 1.0)]
    for e in 1:nech
        d = e == 1 ? delay : delay + rand(rng, 5:30)
        dir = kind == "vertical image" ? (az0, clamp(el0 + rand(rng, [-1, 1]) * (15 + 35rand(rng)), -70, 70)) :
              (clamp(az0 + rand(rng, [-1, 1]) * (15 + 45rand(rng)), -80, 80), el0 + 10randn(rng))
        push!(arr, (dir..., d, (e == 1 ? gain : 0.7gain) * rand(rng, [-1, 1])))
    end
    return arr, (; delay=nech == 0 ? 0 : delay, gain=nech == 0 ? 0.0 : gain, kind=nech == 0 ? "none" : kind)
end
