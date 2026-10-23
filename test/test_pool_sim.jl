# Unit tests for src/pool_sim.jl (image-source pool simulation; no recordings needed).
# run: julia --project=. -t auto test/test_pool_sim.jl
# details: improvements/tdoa_localization_2026-10-22/reports/round3_pool_simulation.md
using Test, DSP, FFTW, Random, Statistics, LinearAlgebra, DataFrames
@isdefined(fname2dt_zoomf6) || (fname2dt_zoomf6(_) = nothing)   # only referenced by set_device; not needed here
include(joinpath(@__DIR__, "..", "src", "localization.jl"))
set_device__ophk_acoustic_D3()
include(joinpath(@__DIR__, "..", "src", "pool_sim.jl"))

const RX = rx_vect[:, get_relevant_channels(rx_vect)]
const FS = 192_000

@testset "pool_sim" begin
    @testset "image counts" begin
        pool = Pool(L=20, W=10, D=4)
        @test length(pool_images([5, 5, 2], [1, 0, 0], pool; order=0)) == 1
        @test length(pool_images([5, 5, 2], [1, 0, 0], pool; order=1)) == 7
        @test length(pool_images([5, 5, 2], [1, 0, 0], pool; order=2)) == 25   # 1 + 6 + 18 distinct images
        imgs = pool_images([5, 5, 2], [1, 0, 0], pool; order=2)
        @test length(unique(round.(im.pos; digits=6) for im in imgs)) == length(imgs)
    end

    @testset "surface echo: closed form" begin
        pool = Pool(L=40, W=40, D=50, refl_bottom=0, refl_wall=0)   # only the surface reflects
        za, zs, R = 1.5, 2.5, 8.0
        array_pos = [2.0, 20.0, za]; src = [2.0 + R, 20.0, zs]
        arr = click_arrivals(src, [-1.0, 0, 0], array_pos, pool; order=1, rel_min=0)
        direct = only(filter(a -> a.order == 0, arr)); surf = only(filter(a -> a.order == 1 && a.gain != 0, arr))
        rd, rs = hypot(R, zs - za), hypot(R, zs + za)
        @test direct.delay ≈ rd / default_c
        @test surf.delay ≈ rs / default_c
        @test perceived_azel(direct.dir)[2] ≈ asind((za - zs) / rd)     # source deeper than the array: below
        @test perceived_azel(surf.dir)[2] ≈ asind((za + zs) / rs)       # image above the surface
        @test abs(perceived_azel(direct.dir)[1]) < 1e-9 && abs(perceived_azel(surf.dir)[1]) < 1e-9   # straight ahead
        @test surf.gain < 0                                              # polarity inverted
        @test direct.gain ≈ beam_gain(acosd(R / rd)) / rd   # beam along -x; the ray to the array is tilted by the depth difference
    end

    @testset "fold" begin
        for x in -30:0.37:30
            p, s = _fold(x, 1.0, 4.0)
            @test 1.0 <= p <= 4.0
            @test s in (-1.0, 1.0)
        end
        @test _fold(5.0, 1.0, 4.0) == (3.0, -1.0)
    end

    @testset "simulated clicks localize to the truth" begin
        tmpl = [exp(-((i - 41) / 1.5)^2) * cos(2π * 60e3 * (i - 41) / FS) for i in 1:101]
        tmpl ./= maximum(abs, tmpl)
        pool = Pool(L=20, W=10, D=5, refl_surface=0, refl_bottom=0, refl_wall=0)
        rng = Xoshiro(3)
        sim = simulate_scenario(rng, [tmpl], zeros(200_000, 3); dur=0.8, nsrc=1, pool=pool, array_pos=[1.5, 5.0, 2.0],
            A0=1000.0, rx=RX, order=0)
        @test size(sim.data, 2) == 5 && all(sim.data[:, 1:2] .== 0)
        @test nrow(sim.arrivals) == nrow(sim.clicks) > 5
        @test all(1 .<= sim.arrivals.sample .<= size(sim.data, 1))
        for i in 1:min(nrow(sim.arrivals), 20)
            a = sim.arrivals[i, :]
            p = round(Int, a.sample)
            _, q = get_tdoa_srp(sim.data[:, 3:5], [p]; window=-50:50, fs=FS, rx_vect=RX, return_quality=true)
            u, v = q.uv[1, :]
            el = asind(v); az = asind(u / cosd(el))
            @test hypot(az - a.az, el - a.el) < 0.5
            @test abs(argmax(abs.(sim.data[p-20:p+20, 5])) - 21) <= 1   # peak on the reference hydrophone at the truth sample
        end
        @test sim.clicks.direct_az ≈ sim.arrivals.az
    end

    # ---- round 4: realism options ----
    tm = [exp(-((i - 41) / 2.0)^2) * cos(2π * 50e3 * (i - 41) / FS) for i in 1:101]
    tm2 = [exp(-((i - 41) / 3.0)^2) * cos(2π * 35e3 * (i - 41) / FS) for i in 1:101]
    @testset "defaults reproduce round 3" begin
        noise = 0.01 .* randn(Xoshiro(5), 400_000, 3)
        rng = Xoshiro(42); pool, ap = random_pool_setup(rng)
        sim = simulate_scenario(rng, [tm, tm2], noise; dur=1.5, nsrc=2, pool=pool, array_pos=ap, A0=100.0, rx=RX)
        # fingerprint of the round-3 simulator (commit b752612 + random reflections of 3ecf536) for this seed
        @test size(sim.data, 1) == 290000 && nrow(sim.arrivals) == 20362
        @test sum(abs, sim.data) ≈ 1.4406981242e+05 rtol = 1e-9
        @test sum(sim.arrivals.sample) ≈ 3.4249072927e+09 rtol = 1e-9
    end

    @testset "realism options" begin
        dry = Pool(L=20, W=10, D=5, refl_surface=0, refl_bottom=0, refl_wall=0)
        ncc(sim) = begin   # correlation between hydrophones 1 and 3 of the first click, aligned by the true fractional delays
            a = sim.arrivals[1, :]; p = round(Int, a.sample)
            u, v = sind(a.az) * cosd(a.el), sind(a.el)
            τ = -(RX' * [u, v, sqrt(max(1 - u^2 - v^2, 0.0))]) ./ default_c .* FS
            x = _sim_frac_delay(sim.data[p-60:p+80, 3], -(τ[1] - τ[3])); y = sim.data[p-60:p+80, 5]
            abs(sum(x .* y)) / sqrt(sum(abs2, x) * sum(abs2, y))
        end
        s0 = simulate_scenario(Xoshiro(7), [tm], zeros(300_000, 3); dur=0.5, nsrc=1, pool=dry, array_pos=[1.5, 5.0, 2.0], A0=1000.0, rx=RX, order=0)
        s1 = simulate_scenario(Xoshiro(7), [tm], zeros(300_000, 3); dur=0.5, nsrc=1, pool=dry, array_pos=[1.5, 5.0, 2.0], A0=1000.0, rx=RX, order=0,
            hyd_scatter=0.6)
        @test ncc(s0) > 0.97
        @test ncc(s1) < ncc(s0) - 0.05                              # hydrophone-specific scattering decorrelates the channels
        s2 = simulate_scenario(Xoshiro(7), [tm], zeros(300_000, 3); dur=0.5, nsrc=1, pool=dry, array_pos=[1.5, 5.0, 2.0], A0=1000.0, rx=RX, order=0,
            reverb_level=0.1, reverb_n=40, reverb_tmax=0.03)
        diffuse = s2.arrivals[s2.arrivals.order .== -1, :]; direct = s2.arrivals[s2.arrivals.order .== 0, :]
        @test nrow(diffuse) == 40 * nrow(direct)
        for c in 1:nrow(direct)   # diffuse arrivals of each click come after its direct arrival, within reverb_tmax (+ array delays)
            d = diffuse[diffuse.click .== direct.click[c], :]
            @test all(direct.sample[c] - 25 .<= d.sample .<= direct.sample[c] + 0.03FS + 25)
        end
        x = randn(Xoshiro(1), 400)
        lp = _sim_lowpass(x, 25e3, FS)
        P(v) = abs2.(rfft(v)); f = rfftfreq(400, FS)
        @test sum(P(lp)[f .> 50e3]) < 0.01 * sum(P(x)[f .> 50e3])  # off-axis low-pass removes the high band
    end
end
