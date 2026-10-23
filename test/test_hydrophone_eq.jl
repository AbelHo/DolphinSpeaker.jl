# Unit tests for src/hydrophone_eq.jl (synthetic clicks; no recordings needed).
# run: julia --project=. -t auto test/test_hydrophone_eq.jl
# details: improvements/tdoa_localization_2026-10-22/reports/round5_hydrophone_eq.md
using Test, DSP, FFTW, Random, Statistics, LinearAlgebra
@isdefined(fname2dt_zoomf6) || (fname2dt_zoomf6(_) = nothing)   # only referenced by set_device; not needed here
include(joinpath(@__DIR__, "..", "src", "localization.jl"))
set_device__ophk_acoustic_D3()
include(joinpath(@__DIR__, "..", "src", "hydrophone_eq.jl"))

const RX = rx_vect[:, get_relevant_channels(rx_vect)]
const FS = 192_000.0
true_tdoa(az, el) = (sv = [sind(az)*cosd(el) sind(el) cosd(az)*cosd(el)]; dt = vec(-sv * RX ./ default_c .* FS); dt .- dt[ref_channel])
fracdelay(x, d) = (N = length(x); irfft(rfft(x) .* exp.(-im * 2π .* rfftfreq(N) .* d), N))
click(n, c0; f0=60e3, σ=6e-6) = [exp(-((i - c0) / (σ * FS))^2) * cos(2π * f0 * (i - c0) / FS) for i in 1:n]
withresp(x, taps) = (y = copy(x); for (d, a) in taps; y[d+1:end] .+= a .* x[1:end-d]; end; y)
angle_between(a, b) = acosd(clamp(dot([sind(a[1])*cosd(a[2]), sind(a[2]), cosd(a[1])*cosd(a[2])], [sind(b[1])*cosd(b[2]), sind(b[2]), cosd(b[1])*cosd(b[2])]), -1, 1))

# clicks from random directions in front, one every 2000 samples, each hydrophone filtered by its own response
function synth(rng, nclicks, resp; snr_db=40)   # click peak 1 vs noise: the correlation metric's ceiling is about 0.98
    sp = 2000; n = sp * (nclicks + 1)
    data = zeros(n, 3); peaks = Int[]; dirs = Tuple{Float64,Float64}[]
    for i in 1:nclicks
        az, el = 80rand(rng) - 40, 50rand(rng) - 25
        p = i * sp; push!(peaks, p); push!(dirs, (az, el))
        seg = click(sp, sp ÷ 2; f0=40e3 + 30e3rand(rng))
        T = true_tdoa(az, el)
        for ch in 1:3
            data[p-sp÷2+1:p+sp÷2, ch] .+= withresp(fracdelay(seg, T[ch]), resp[ch])
        end
    end
    return data .+ 10^(-snr_db / 20) .* randn(rng, n, 3), peaks, dirs
end
# correlation between hydrophone k and the reference of each click, aligned by the true delays
function aligned_ncc(data, peaks, dirs, k)
    map(eachindex(peaks)) do i
        T = true_tdoa(dirs[i]...)
        x = fracdelay(data[peaks[i]-60:peaks[i]+80, k], -(T[k] - T[3])); y = data[peaks[i]-60:peaks[i]+80, 3]
        abs(sum(x .* y)) / sqrt(sum(abs2, x) * sum(abs2, y))
    end
end
direction_errors(data, peaks, dirs) = map(eachindex(peaks)) do i
    _, q = get_tdoa_srp(data, [peaks[i]]; window=-50:50, fs=FS, rx_vect=RX, return_quality=true)
    u, v = q.uv[1, :]; el = asind(v); angle_between((asind(u / cosd(el)), el), dirs[i])
end

@testset "hydrophone equalization" begin
    resp = [[(7, 0.6), (19, -0.4)], [(5, -0.5), (23, 0.3)], Tuple{Int,Float64}[]]   # minimum phase; reference undistorted
    data, peaks, dirs = synth(Xoshiro(1), 400, resp)
    eq = estimate_hydrophone_eq(data, peaks; fs=FS, rx_vect=RX, ref=3, iters=3)
    deq = apply_hydrophone_eq(data, eq)

    @testset "recovers the responses" begin
        f = eq.f; ib = (f .>= eq.band[1]) .& (f .<= eq.band[2])
        for k in 1:2
            Gtrue = [1 + sum(a * cis(-2π * fk * d / FS) for (d, a) in resp[k]) for fk in f]
            Gest = eq.G[:, k] .* median(abs.(Gtrue[ib]))                 # the estimate is normalized to median magnitude 1
            @test median(abs.(angle.(Gest[ib] ./ Gtrue[ib]))) < 0.15      # phase within about 9°
            @test median(abs.(abs.(Gest[ib]) ./ abs.(Gtrue[ib]) .- 1)) < 0.15
            @test abs(eq.delay_removed[k]) < 0.5                            # minimum-phase responses: no pure delay
        end
        @test all(eq.E[:, 3] .== 1) && all(eq.G[:, 3] .== 1)
    end

    @testset "restores coherence and directions" begin
        # (at 30 dB even the oracle equalizer, built from the true responses, only reaches 0.85: the window is mostly noise)
        for k in 1:2
            before, after = median(aligned_ncc(data, peaks, dirs, k)), median(aligned_ncc(deq, peaks, dirs, k))
            @test before < 0.9
            @test after > 0.95
        end
        e0, e1 = median(direction_errors(data, peaks, dirs)), median(direction_errors(deq, peaks, dirs))
        @test e1 < e0
        @test e1 < 0.5
        @test deq[:, 3] == data[:, 3]                                     # the reference hydrophone is untouched
    end

    @testset "identity without distortion" begin
        d0, p0, _ = synth(Xoshiro(2), 200, [Tuple{Int,Float64}[] for _ in 1:3])
        eq0 = estimate_hydrophone_eq(d0, p0; fs=FS, rx_vect=RX, ref=3, iters=2)
        ib = (eq0.f .>= eq0.band[1]) .& (eq0.f .<= eq0.band[2])
        @test median(abs.(angle.(eq0.G[ib, 1:2]))) < 0.05
        dd = apply_hydrophone_eq(d0, eq0)
        @test norm(dd .- d0) / norm(d0) < 0.1
    end
end

@testset "blind estimate under strong distortion" begin
    # a scattered copy stronger than the direct arrival on hydrophone 1, and a distorted reference
    resp = [[(7, 1.3), (19, -0.4)], [(5, -0.5), (23, 0.3)], [(11, 0.6)]]
    data, peaks, dirs = synth(Xoshiro(3), 600, resp)
    f(eq) = eq.f; inband(eq) = (eq.f .>= eq.band[1]) .& (eq.f .<= eq.band[2])
    H(k, fk) = 1 + sum(a * cis(-2π * fk * d / FS) for (d, a) in resp[k]; init=0.0im)
    eqb = estimate_hydrophone_eq_blind(data, peaks; fs=FS, rx_vect=RX, ref=3)
    ib = inband(eqb)
    for k in 1:2
        Gt = [H(k, fk) / H(3, fk) for fk in eqb.f]
        r = eqb.G[ib, k] ./ Gt[ib]                                    # equal up to a pure delay and a scale
        φ = unwrap(angle.(r)); X = [ones(count(ib)) eqb.f[ib]]; ab = X \ φ
        @test median(abs.(φ .- X * ab)) < 0.2                          # dispersive phase recovered
        @test abs(rem(ab[1], 2π, RoundNearest)) < 0.3                    # constant phase recovered (no polarity error)
        @test median(abs.(log.(abs.(r) ./ median(abs.(r))))) < 0.2       # relative magnitude recovered
    end
    # with the true pure delays (from calibration data) the directions are restored
    dtrue = zeros(3)
    for k in 1:2
        Gt = [H(k, fk) / H(3, fk) for fk in eqb.f]
        φ = unwrap(angle.(Gt[ib] ./ eqb.G[ib, k])); X = [ones(count(ib)) eqb.f[ib]]
        dtrue[k] = -(X \ φ)[2] * FS / 2π
    end
    e_raw = median(direction_errors(data, peaks, dirs))
    e_cal = median(direction_errors(apply_hydrophone_eq(data, with_delays(eqb, dtrue)), peaks, dirs))
    @test e_cal < 0.5
    @test e_cal < e_raw
end
