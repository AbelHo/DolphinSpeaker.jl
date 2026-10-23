# Unit tests for src/multiarrival.jl (synthetic data only, no recordings needed).
# run: julia --project=. -t auto test/test_multiarrival.jl
# details: improvements/tdoa_localization_2026-10-22/README.md
using Test, DSP, FFTW, Random, Statistics, LinearAlgebra, Printf, DataFrames
@isdefined(fname2dt_zoomf6) || (fname2dt_zoomf6(_) = nothing)   # only referenced by set_device; not needed here
include(joinpath(@__DIR__, "..", "src", "localization.jl"))
set_device__ophk_acoustic_D3()
include(joinpath(@__DIR__, "..", "src", "multiarrival.jl"))

const RX = rx_vect[:, get_relevant_channels(rx_vect)]
const FS = 192_000.0
const M = size(RX, 2)

true_tdoa(az, el) = (sv = [sind(az)*cosd(el) sind(el) cosd(az)*cosd(el)]; dt = vec(-sv * RX ./ default_c .* FS); dt .- dt[ref_channel])
uv_of(az, el) = (sind(az) * cosd(el), sind(el))
fracdelay(x, d) = (N = length(x); irfft(rfft(x) .* exp.(-im * 2π .* rfftfreq(N) .* d), N))
# broadband like a dolphin click (narrower bands raise the SRP sidelobes)
click(n, c0; f0=60e3, σ=8e-6) = [exp(-((i - c0) / (σ * FS))^2) * cos(2π * f0 * (i - c0) / FS) for i in 1:n]

# one window of 3-channel data: plane-wave arrivals (az, el, delay in samples, gain) of the same click, plus noise
function synth_window(arrivals; n=1000, snr_db=30, seed=1)
    Random.seed!(seed)
    data = zeros(n, M)
    for (az, el, delay, gain) in arrivals
        seg = gain .* click(n, n ÷ 2 + delay)
        T = true_tdoa(az, el)
        for ch in 1:M; data[:, ch] .+= fracdelay(seg, T[ch]); end
    end
    return data .+ 10^(-snr_db / 20) .* randn(n, M)
end
angle_uv(a, b) = _angle_uv(a, b)

@testset "multiarrival" begin
    @testset "one plane wave: candidate 1 = get_tdoa_srp" begin
        for (seed, (az, el)) in enumerate([(10.0, -5.0), (-25.0, 12.0), (33.0, 20.0)])
            d = synth_window([(az, el, 0, 1.0)]; seed=seed)
            c = get_srp_candidates(d, [500]; window=-50:50, fs=FS, rx_vect=RX)
            _, q = get_tdoa_srp(d, [500]; window=-50:50, fs=FS, rx_vect=RX, return_quality=true)
            @test c.uv[1, 1, :] ≈ q.uv[1, :]
            @test c.score[1, 1] ≈ q.srp[1]
            @test angle_uv(Tuple(c.uv[1, 1, :]), uv_of(az, el)) < 0.5
            @test c.nfound[1] == 1
            @test c.closure[1] < 0.3
            @test c.ncc_min[1] > 0.9
            @test abs(c.tarr[1, 1]) <= 2
        end
    end

    @testset "two plane waves 30 samples apart" begin
        d1, d2 = (5.0, -10.0), (5.0, 30.0)   # same bearing, vertical offset: an image-like pair
        d = synth_window([(d1..., 0, 1.0), (d2..., 30, 0.8)])
        c = get_srp_candidates(d, [500]; window=-50:50, fs=FS, rx_vect=RX)
        @test c.nfound[1] >= 2
        found(dir) = findfirst(k -> !isnan(c.uv[1, k, 1]) && angle_uv(Tuple(c.uv[1, k, :]), uv_of(dir...)) < 1, 1:3)
        k1, k2 = found(d1), found(d2)
        @test !isnothing(k1) && !isnothing(k2)
        @test c.tarr[1, k1] < c.tarr[1, k2]                 # the direct arrival is earlier in its beam
        @test abs(c.tarr[1, k2] - c.tarr[1, k1] - 30) <= 3
    end

    @testset "envelope_arrivals" begin
        d = synth_window([(0.0, 0.0, 0, 1.0), (0.0, 0.0, 30, 0.8)])
        @test all(envelope_arrivals(d, [500]; window=-50:50) .== 2)
        @test all(envelope_arrivals(synth_window([(0.0, 0.0, 0, 1.0)]), [500]; window=-50:50) .== 1)
    end

    @testset "follower_flags" begin
        s = [1000, 1100, 1500, 2000, 2400]   # fs 1000 Hz: gaps 100, 400, 500, 400 ms
        f = follower_flags(s, 1000; echo_window_ms=450)
        @test f.gap_prev ≈ [Inf, 100, 400, 500, 400]
        @test f.gap_next ≈ [100, 400, 500, 400, Inf]
        @test f.follower == [false, true, true, false, true]
    end

    @testset "tracker: two sources, their images and clutter" begin
        Random.seed!(7)
        # source A clicks every 10 ms, B every 13 ms, for 2 s; each click also shows a vertical image; plus clutter
        t = Float64[]; truth = Vector{Vector{Symbol}}(); uv = Vector{Vector{NTuple{2,Float64}}}(); sc = Vector{Vector{Float64}}()
        dirA(t) = (-20 + 5t, -10.0); dirB(t) = (25 - 4t, 5.0)
        noisy((az, el)) = uv_of(az + 0.8randn(), el + 0.8randn())
        events = sort([[(ti, :A) for ti in 0:0.010:2]; [(ti, :B) for ti in 0.003:0.013:2]])
        for (ti, s) in events
            dir = s == :A ? dirA(ti) : dirB(ti)
            img = (dir[1], dir[2] + 35)
            cl = (rand() * 120 - 60, rand() * 80 - 40)
            push!(t, ti); push!(truth, [s, Symbol(s, :img), :clutter])
            push!(uv, [noisy(dir), noisy(img), uv_of(cl...)]); push!(sc, [0.5, 0.35, 0.25])
        end
        n = length(t); U = fill(NaN, n, 3, 2); S = fill(NaN, n, 3)
        for i in 1:n, k in 1:3; U[i, k, :] .= uv[i][k]; S[i, k] = sc[i][k]; end
        T = track_candidates(t, U, S)
        conf = Dict(zip(T.tracks.id, T.tracks.confirmed))
        # majority track per true path, then the share of that path's candidates it holds
        for path in (:A, :B, :Aimg, :Bimg)
            ids = [T.track[i, k] for i in 1:n for k in 1:3 if truth[i][k] == path]
            main = argmax(id -> count(==(id), ids), unique(filter(>(0), ids)))
            @test conf[main]
            @test count(==(main), ids) / length(ids) >= 0.95
            # and it holds (almost) nothing else
            other = [truth[i][k] for i in 1:n for k in 1:3 if T.track[i, k] == main]
            @test count(==(path), other) / length(other) >= 0.95
        end
        clutter_ids = unique(filter(>(0), [T.track[i, 3] for i in 1:n]))
        @test count(id -> conf[id], clutter_ids) <= 0.02 * length(clutter_ids)

        # a gap longer than t_delete ends a track
        tt = [0.0, 0.01, 0.02, 0.03, 0.5, 0.51]
        U2 = fill(NaN, 6, 1, 2); U2[:, 1, 1] .= 0.1; U2[:, 1, 2] .= 0.0
        T2 = track_candidates(tt, U2, fill(0.5, 6, 1); t_delete=0.25)
        @test T2.track[1:4, 1] == fill(T2.track[1, 1], 4)
        @test T2.track[5, 1] != T2.track[1, 1]
    end

    @testset "track_colours" begin
        Random.seed!(3)
        tr = DataFrame(t_start=rand(200) .* 10)
        tr.t_end = tr.t_start .+ rand(200) .* 0.5
        col = track_colours(tr)
        for i in 1:200, j in i+1:200
            overlap = tr.t_start[j] <= tr.t_end[i] && tr.t_end[j] >= tr.t_start[i]
            if overlap && count(q -> tr.t_start[q] <= tr.t_end[i] && tr.t_end[q] >= tr.t_start[i], 1:200) <= length(TRACK_PALETTE)
                @test col[i] != col[j]
            end
        end
        @test all(1 .<= col .<= length(TRACK_PALETTE))
    end

    # ---- round 2 ----
    @testset "srp_support: true direction vs sidelobe" begin
        d = synth_window([(10.0, -5.0, 0, 1.0)])
        c = get_srp_candidates(d, [500]; window=-50:50, fs=FS, rx_vect=RX, min_rel=0.0)   # keeps sidelobe candidates
        sp = srp_support(d, [500], c.uv; window=-50:50, fs=FS, rx_vect=RX)
        @test c.nfound[1] >= 2
        @test sp[1, 1] > 0.7
        @test all(sp[1, 2:c.nfound[1]] .< 0.3)
    end

    @testset "CLEAN candidates: single wave" begin
        d = synth_window([(10.0, -5.0, 0, 1.0)])
        c = get_srp_candidates_clean(d, [500]; window=-50:50, fs=FS, rx_vect=RX, min_rel=0.35)
        _, q = get_tdoa_srp(d, [500]; window=-50:50, fs=FS, rx_vect=RX, return_quality=true)
        @test c.uv[1, 1, :] ≈ q.uv[1, :]
        @test c.nfound[1] == 1
    end

    @testset "SIC candidates" begin
        # single wave: candidate 1 = get_tdoa_srp, nothing else above 0.5 of it
        for (seed, (az, el)) in enumerate([(10.0, -5.0), (-30.0, 20.0)])
            d = synth_window([(az, el, 0, 1.0)]; seed=seed)
            c = get_srp_candidates_sic(d, [500]; window=-50:50, fs=FS, rx_vect=RX, min_rel=0.5)
            _, q = get_tdoa_srp(d, [500]; window=-50:50, fs=FS, rx_vect=RX, return_quality=true)
            @test c.uv[1, 1, :] ≈ q.uv[1, :]
            @test c.nfound[1] == 1
            @test c.efrac[1, 1] > 0.8
        end
        # direct + image 30 samples later: both found, timed, and the earliest is the direct path
        d1, d2 = (5.0, -10.0), (5.0, 30.0)
        d = synth_window([(d1..., 0, 1.0), (d2..., 30, 0.8)])
        c = get_srp_candidates_sic(d, [500]; window=-50:50, fs=FS, rx_vect=RX, min_rel=0.5)
        found(dir) = findfirst(k -> !isnan(c.uv[1, k, 1]) && angle_uv(Tuple(c.uv[1, k, :]), uv_of(dir...)) < 1, 1:3)
        k1, k2 = found(d1), found(d2)
        @test !isnothing(k1) && !isnothing(k2)
        @test abs(c.tarr[1, k1]) <= 2 && abs(c.tarr[1, k2] - 30) <= 2
        # equally strong arrivals 10 samples apart: the earliest beam arrival is still the direct path
        d = synth_window([(d1..., 0, 1.0), (d2..., 10, 1.0)])
        c = get_srp_candidates_sic(d, [500]; window=-50:50, fs=FS, rx_vect=RX, min_rel=0.5)
        ks = findall(k -> !isnan(c.uv[1, k, 1]), 1:3)
        @test length(ks) >= 2
        @test angle_uv(Tuple(c.uv[1, ks[argmin(c.tarr[1, ks])], :]), uv_of(d1...)) < 1
    end

    @testset "track_lags and select_candidates" begin
        # 12 clicks; tracks 1 (direct, tarr 0) and 2 (image, tarr 20) share clicks 1-10, click 11 has only track 2,
        # click 12 only an unconfirmed track 3
        track = zeros(Int, 12, 3); tarr = fill(NaN, 12, 3); score = fill(NaN, 12, 3)
        for i in 1:10
            track[i, 1] = 2; tarr[i, 1] = 20; score[i, 1] = 0.6     # the image is candidate 1 (stronger)
            track[i, 2] = 1; tarr[i, 2] = 0; score[i, 2] = 0.5
        end
        track[11, 1] = 2; tarr[11, 1] = 21; score[11, 1] = 0.6
        track[12, 1] = 3; tarr[12, 1] = 0; score[12, 1] = 0.4
        conf = Dict(1 => true, 2 => true, 3 => false); hits = Dict(1 => 10, 2 => 11, 3 => 1)
        lg = track_lags(track, tarr, conf)
        @test lg[1] == (0.0, 10) && lg[2] == (20.0, 10)
        s = select_candidates(track, tarr, score, conf, hits; lags=lg)
        @test s.sel[1:10] == fill(2, 10) && all(s.status[1:10] .== "direct")
        @test s.status[11] == "echo"
        @test s.sel[12] == 1 && s.status[12] == "untracked"
    end
end
