# Unit tests for the TDOA functions in src/localization.jl (synthetic data only, no recordings needed).
# run: julia --project=. -t auto test/test_tdoa.jl
# details: improvements/tdoa_localization_2026-10-21/README.md
using Test, DSP, FFTW, Random, Statistics, LinearAlgebra, Printf
@isdefined(fname2dt_zoomf6) || (fname2dt_zoomf6(_) = nothing)   # only referenced by set_device; not needed here
include(joinpath(@__DIR__, "..", "src", "localization.jl"))
@isdefined(energy) || (energy(x) = sum(abs2, x; dims=1))         # SignalAnalysis.energy, loaded elsewhere in the pipeline
set_device__ophk_acoustic_D3()

const W = window_impulsive   # D3 click window (-50:50)
const RX = rx_vect[:, get_relevant_channels(rx_vect)]
const FS = 192_000.0
const M = size(RX, 2)

# TDOAs (samples, relative to ref_channel) of a plane wave from (az, el) in degrees, same model as cost_tdoa2ang
true_tdoa(az, el) = (sv = [sind(az)*cosd(el) sind(el) cosd(az)*cosd(el)]; dt = vec(-sv * RX ./ default_c .* FS); dt .- dt[ref_channel])
dirvec(a) = [sin(a[1])*cos(a[2]), sin(a[2]), cos(a[1])*cos(a[2])]
fracdelay(x, d) = (N = length(x); irfft(rfft(x) .* exp.(-im * 2π .* rfftfreq(N) .* d), N))
click(n, c0; f0=60e3, σ=15e-6) = [exp(-((i - c0) / (σ * FS))^2) * cos(2π * f0 * (i - c0) / FS) for i in 1:n]

# 3-channel data with one Gabor click per direction, delayed by the exact fractional TDOAs, plus white noise
function synth(dirs, snr_db; spacing=2000, seed=2)
	Random.seed!(seed)
	n = spacing * (length(dirs) + 1)
	data = zeros(n, M)
	peaks = Int[]
	Ttrue = zeros(length(dirs), M)
	for (i, (az, el)) in enumerate(dirs)
		p = i * spacing
		push!(peaks, p)
		Ttrue[i, :] = true_tdoa(az, el)
		seg = click(spacing, spacing ÷ 2)
		for ch in 1:M
			data[p-spacing÷2+1:p+spacing÷2, ch] .+= fracdelay(seg, Ttrue[i, ch])
		end
	end
	data .+= 10^(-snr_db / 20) .* randn(n, M)   # click peak amplitude is 1
	return data, peaks, Ttrue
end

@testset "TDOA localization" begin
	@testset "pairs_to_delays" begin
		t_true = [0.3, -1.7, 1.4]; t_true .-= mean(t_true)
		pairs = [(j, k) for j in 1:M for k in j+1:M]
		t, rms = pairs_to_delays([t_true[j] - t_true[k] for (j, k) in pairs], pairs, M)
		@test maximum(abs, t .- t_true) < 1e-12
		@test rms < 1e-12
		_, rms = pairs_to_delays([t_true[j] - t_true[k] for (j, k) in pairs] .+ [0.3, 0, 0], pairs, M)
		@test 3rms ≈ 0.3 atol = 1e-12   # 3 sensors: closure = 3 x rms misfit
	end

	@testset "tdoa2dir_ls vs NelderMead tdoa2dir" begin
		Random.seed!(1)
		dirs = [(rand() * 120 - 60, rand() * 100 - 50) for _ in 1:200]
		T = reduce(vcat, [true_tdoa(a, e)' for (a, e) in dirs])
		a_true = reduce(vcat, [[a e] for (a, e) in dirs])
		@test maximum(abs, rad2deg.(tdoa2dir_ls(T, RX, FS)) .- a_true) < 1e-6
		# NelderMead can return the mirrored back-hemisphere az/el (same TDOAs): compare direction cosines (u, v)
		uvgap(A, B) = maximum(i -> maximum(abs, dirvec(A[i, :])[1:2] .- dirvec(B[i, :])[1:2]), axes(A, 1))
		infov = vec(abs.(a_true[:, 1]) .<= 40 .&& abs.(a_true[:, 2]) .<= 27)
		for TT in (T, round.(T))   # exact and integer TDOAs
			ls, nm = tdoa2dir_ls(TT, RX, FS), tdoa2dir(TT, RX, FS)
			@test uvgap(ls, nm) < 1e-3
			@test maximum(abs, rad2deg.(ls[infov, :] .- nm[infov, :])) < 0.1   # identical az/el inside the camera FOV
		end
		_, res = tdoa2dir_ls(T, RX, FS; return_residual=true)
		@test maximum(res) < 1e-12
	end

	@testset "fractional delays of synthetic clicks (40 dB)" begin
		Random.seed!(3)
		dirs = [(rand() * 70 - 35, rand() * 50 - 25) for _ in 1:100]
		data, peaks, Ttrue = synth(dirs, 40)
		methods = [
			"gcc none" => (d, p) -> get_tdoa_gcc(d, p; window=W, fs=FS, rx_vect=RX, weighting=:none),
			"gcc phat" => (d, p) -> get_tdoa_gcc(d, p; window=W, fs=FS, rx_vect=RX, weighting=:phat),
			"gcc phat guided" => (d, p) -> get_tdoa_gcc(d, p; window=W, fs=FS, rx_vect=RX, weighting=:phat, guide=:envelope),
			"gcc scot" => (d, p) -> get_tdoa_gcc(d, p; window=W, fs=FS, rx_vect=RX, weighting=:scot),
			"gcc ml" => (d, p) -> get_tdoa_gcc(d, p; window=W, fs=FS, rx_vect=RX, weighting=:ml),
			"gcc phat refpairs" => (d, p) -> get_tdoa_gcc(d, p; window=W, fs=FS, rx_vect=RX, all_pairs=false),
			"srp phat" => (d, p) -> get_tdoa_srp(d, p; window=W, fs=FS, rx_vect=RX),
		]
		for (name, f) in methods
			T̂ = f(data, peaks)
			T̂ = T̂ .- T̂[:, ref_channel]
			rmse = sqrt(mean(abs2, T̂ .- Ttrue))
			@info @sprintf("%-18s tdoa rmse %.3f samples", name, rmse)
			@test rmse < 0.1
		end
		# same sign convention as the existing integer method (positive = later than ref_channel)
		R = Float64.(get_tdoa_raw(data, peaks; window=W, ref_channel=ref_channel))
		R .-= R[:, ref_channel]
		others = setdiff(1:M, ref_channel)
		big = abs.(Ttrue[:, others]) .> 2
		@test all(sign.(R[:, others][big]) .== sign.(Ttrue[:, others][big]))
	end

	@testset "fixed copies vs originals" begin
		Random.seed!(4)
		data = randn(5000, M)
		peaks = collect(200:300:4800)
		o = get_tdoa_minmax(data, peaks; window=W, ref_channel=ref_channel)
		n2 = get_tdoa_minmax2(data, peaks; window=W, ref_channel=ref_channel)
		votes_max = [count(ch -> maximum(data[W .+ p, ch]) > -minimum(data[W .+ p, ch]), 1:M) > M / 2 for p in peaks]
		@test o[votes_max, :] == n2[votes_max, :]
		@test o[.!votes_max, :] != n2[.!votes_max, :]
		Random.seed!(3)
		dirs = [(rand() * 70 - 35, rand() * 50 - 25) for _ in 1:20]
		data_c, peaks_c, _ = synth(dirs, 30)
		@test get_tdoa_envelope_filtered(data_c, peaks_c; window=W, ref_channel=ref_channel) ==
			get_tdoa_envelope_filtered2(data_c, peaks_c; window=W, ref_channel=ref_channel, fs=500_000)
	end

	@testset "clicks at the file edges" begin
		data = randn(5000, M)
		edge = [10, 4990]
		@test size(get_tdoa_minmax2(data, edge; window=W)) == (2, M)
		@test size(get_tdoa_envelope_filtered2(data, edge; window=W, fs=FS)) == (2, M)
		@test size(get_tdoa_gcc(data, edge; window=W, fs=FS, rx_vect=RX)) == (2, M)
		@test size(get_tdoa_srp(data, edge; window=W, fs=FS, rx_vect=RX)) == (2, M)
	end
end
