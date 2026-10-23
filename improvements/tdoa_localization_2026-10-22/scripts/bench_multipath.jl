# Semi-synthetic multipath benchmark: a clean real click arrives from a known direction, plus 1-2 delayed copies
# (echoes) from other directions inside the TDOA window, plus real noise. How often does each method give the
# direct path, and is the direct path among the SRP candidates?
# needs: results/data/click_features_New_device.csv (run_multiarrival_newdevice.jl)
# usage: julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/bench_multipath.jl [n_trials]
src = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(src, "run_example.jl")); set_device__ophk_acoustic_D3()
include(joinpath(src, "multiarrival.jl"))
using Random, Plots
n_trials = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4000
out_dir = joinpath(@__DIR__, "..", "results")
feat = CSV.read(joinpath(out_dir, "data", "click_features_New_device.csv"), DataFrame)
wavdir = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device"
fs = 192000; rx = rx_vect[:, get_relevant_channels(rx_vect)]
wavf = joinpath(wavdir, "251006_001_0001.WAV")
f1 = feat[feat.file .== "251006_001_0001.WAV", :]

# templates: clean, strong, single-candidate clicks; the ref-channel waveform around the click
clean = f1[(.!f1.follower) .& (f1.closure_phat .<= 0.3) .& (f1.nfound .== 1) .& (f1.score_1 .>= quantile(f1.score_1, 0.75)), :]
Random.seed!(11)
tsel = clean.sample[randperm(nrow(clean))[1:min(300, nrow(clean))]]
blk(p, a, b) = read_filtered_block(wavf, (p + a):(p + b); fs=fs)
templates = [(x = blk(p, -400, 400)[401-40:401+60, 3]; x ./ maximum(abs, x)) for p in tsel]   # 101 samples, peak 1
# noise: 1000-sample stretches at least 20 ms away from any detection
allp = sort(f1.sample); noises = Matrix{Float64}[]
while length(noises) < 300
    p = rand(allp[1]:allp[end])
    i = searchsortedfirst(allp, p)
    near = minimum(abs.(allp[max(i-1, 1):min(i, length(allp))] .- p))
    near > 0.02fs && push!(noises, blk(p, 0, 999))
end
# real window SNR: click peak vs noise rms, from the templates' snippets
snr_db = 19.0                                              # median real-click window SNR (2026-10-21 report)

arr_t(az, el) = (t = vec(-[sind(az)*cosd(el) sind(el) cosd(az)*cosd(el)] * rx ./ default_c .* fs); t .- t[ref_channel])
uvd(az, el) = (sind(az) * cosd(el), sind(el))
c0 = 500
rows = DataFrame(delay=Int[], gain=Float64[], nech=Int[], kind=String[], err_srp=Float64[], err_gcc=Float64[],
    recall1=Bool[], recall3=Bool[], rank=Int[], earliest_ok=Bool[], nfound=Int[])
lk = ReentrantLock()
Threads.@threads for trial in 1:n_trials
    rng = Xoshiro(trial)
    s = templates[rand(rng, eachindex(templates))]; nz = noises[rand(rng, eachindex(noises))]
    az0, el0 = 70rand(rng) - 35, 44rand(rng) - 22
    nech = rand(rng, 1:2); delay = rand(rng, [10, 20, 40, 58]); gain = rand(rng, [0.5, 0.8, 1.0, 1.25])
    kind = rand(rng, ["vertical image", "random azimuth"])
    arrivals = [(az0, el0, 0, 1.0)]
    for e in 1:nech
        d = e == 1 ? delay : delay + rand(rng, 5:30)
        dir = kind == "vertical image" ? (az0, clamp(el0 + rand(rng, [-1, 1]) * (15 + 35rand(rng)), -70, 70)) : (clamp(az0 + rand(rng, [-1, 1]) * (15 + 45rand(rng)), -80, 80), el0 + 10randn(rng))
        push!(arrivals, (dir..., d, (e == 1 ? gain : 0.7gain) * rand(rng, [-1, 1])))
    end
    data = zeros(1000, 3)
    for (az, el, d, g) in arrivals, ch in 1:3
        seg = zeros(1000); seg[c0-40+d:c0+60+d] .= g .* s
        data[:, ch] .+= _frac_delay(seg, arr_t(az, el)[ch])
    end
    data .+= nz ./ sqrt(mean(abs2, nz)) .* 10^(-snr_db / 20)
    # the detector sits on the strongest envelope peak of the ref channel
    e3 = abs.(hilbert(data[:, 3])); pk = c0 - 20 + argmax(e3[c0-20:c0+120]) - 1
    _, q = get_tdoa_srp(data, [pk]; window=window_impulsive, fs=fs, rx_vect=rx, return_quality=true)
    Tg = get_tdoa_gcc(data, [pk]; window=window_impulsive, fs=fs, rx_vect=rx, guide=:envelope)
    ug = tdoa2uv(Tg, rx, fs)
    c = get_srp_candidates(data, [pk]; window=window_impulsive, fs=fs, rx_vect=rx)
    truth = uvd(az0, el0)
    errs = [isnan(c.uv[1, k, 1]) ? Inf : _angle_uv(Tuple(c.uv[1, k, :]), truth) for k in 1:3]
    found = findfirst(<(2.0), errs)
    ks = findall(isfinite, errs)
    early = ks[argmin(c.tarr[1, ks])]
    gerr = hypot(ug[1], ug[2]) > 1 ? 90.0 : _angle_uv((ug[1], ug[2]), truth)
    lock(lk) do
        push!(rows, (delay, gain, nech, kind, _angle_uv(Tuple(q.uv[1, :]), truth), gerr, errs[1] < 2, !isnothing(found),
            something(found, 0), errs[early] < 2, c.nfound[1]))
    end
end
CSV.write(joinpath(out_dir, "bench_multipath.csv"), rows)

io = IOBuffer()
@printf(io, "trials %d, templates %d clean clicks of 251006_001_0001, noise from quiet stretches, SNR %.0f dB\n", nrow(rows), length(templates), snr_db)
println(io, "share of trials with the direct path within 2 deg (%):")
@printf(io, "%-24s %6s %6s %6s %8s %8s %8s %7s\n", "group", "n", "SRP", "GCC", "cand#1", "any of 3", "earliest", "found>=2")
grp(nm, s) = @printf(io, "%-24s %6d %6.1f %6.1f %6.1f %8.1f %8.1f %7.1f\n", nm, count(s), 100mean(rows.err_srp[s] .< 2), 100mean(rows.err_gcc[s] .< 2),
    100mean(rows.recall1[s]), 100mean(rows.recall3[s]), 100mean(rows.earliest_ok[s]), 100mean(rows.nfound[s] .>= 2))
grp("all", trues(nrow(rows)))
for d in [10, 20, 40, 58]; grp("delay $d smp ($(round(Int, d / fs * 1e6)) us)", rows.delay .== d); end
for g in [0.5, 0.8, 1.0, 1.25]; grp("echo gain $g", rows.gain .== g); end
for k in ["vertical image", "random azimuth"]; grp(k, rows.kind .== k); end
for n in 1:2; grp("$n echo(es)", rows.nech .== n); end
write_md = (fname, title, body) -> (open(fname, "w") do f; println(f, "# ", title, "\n\n```\n", body, "```"); end; print(body))
write_md(joinpath(out_dir, "bench_multipath.md"), "Semi-synthetic multipath benchmark", String(take!(io)))

# figure: share within 2 deg vs echo delay, one line per method (Okabe-Ito order)
ds = [10, 20, 40, 58]; dus = ds ./ fs .* 1e6
series = [("SRP (Impulse6)", :err_srp), ("GCC guided (Impulse5)", :err_gcc), ("candidate 1", :recall1), ("any of 3 candidates", :recall3), ("earliest beam", :earliest_ok)]
cols = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"]
p = plot(; xlabel="echo delay (µs)", ylabel="direct path within 2° (%)", ylims=(0, 100), legend=:bottomright, size=(800, 480), grid=:y, framestyle=:box)
for ((nm, col), cc) in zip(series, cols)
    y = [col in (:err_srp, :err_gcc) ? 100mean(rows[rows.delay .== d, col] .< 2) : 100mean(rows[rows.delay .== d, col]) for d in ds]
    plot!(p, dus, y; label=nm, color=cc, lw=2, marker=:circle, ms=5, msw=0)
end
savefig(p, joinpath(out_dir, "bench_multipath.png"))
