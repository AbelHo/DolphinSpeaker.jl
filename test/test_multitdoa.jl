# Unit tests for the pure helpers of src/run_multitdoa.jl (synthetic data only, no recordings needed).
# run: julia --project=. -t auto test/test_multitdoa.jl
# details: reviews/plans/multitdoa-webui-csvs_2026-10-22.md
using Test, DSP, FFTW, Random, Statistics, LinearAlgebra, Printf
@isdefined(fname2dt_zoomf6) || (fname2dt_zoomf6(_) = nothing)   # only referenced by set_device; not needed here
include(joinpath(@__DIR__, "..", "src", "localization.jl"))
set_device__ophk_acoustic_D3()
include(joinpath(@__DIR__, "..", "src", "run_multitdoa.jl"))

const RX = rx_vect[:, get_relevant_channels(rx_vect)]
const FS = 192_000.0

# TDOAs (samples, relative to ref_channel) of a plane wave from (az, el) in degrees, same model as cost_tdoa2ang
true_tdoa(az, el) = (sv = [sind(az)*cosd(el) sind(el) cosd(az)*cosd(el)]; dt = vec(-sv * RX ./ default_c .* FS); dt .- dt[ref_channel])

@testset "multi-TDOA overlay helpers" begin
	@testset "impulse_variant_rows" begin
		n = 5
		px = rand(n, 2)
		colour = Any[fill(20, n) rand(n) rand(n) rand(n) rand(n) fill("circle", n)]
		flagged = Bool[0, 1, 0, 0, 1]
		rows = impulse_variant_rows(px, colour, flagged)
		@test size(rows) == (n, 8)
		@test rows[:, 1:2] == px
		@test rows[:, 3:7] == colour[:, 1:5]
		@test rows[:, 8] == ["circle", "triangle", "circle", "circle", "triangle"]
		@test impulse_variant_rows(px, colour, flagged; flag_shape="cross")[2, 8] == "cross"
		@test impulse_variant_rows(px, colour, falses(n)) == [px colour]
		@test colour[2, 6] == "circle"   # the input is not modified
	end

	@testset "flag_clicks" begin
		T = reduce(vcat, [true_tdoa(az, el)' for (az, el) in [(10, 5), (-30, -10), (0, 0), (35, 20)]])
		# consistent and inside the plane-wave range: nothing flagged
		flagged, infeasible = flag_clicks(T, zeros(4), RX, FS)
		@test !any(flagged) && !any(infeasible)
		# closure above the limit is flagged, up to the limit is not
		flagged, _ = flag_clicks(T, [0.2, 3.0, 1.0, 1.5], RX, FS; closure_max=1.0)
		@test flagged == Bool[0, 1, 0, 1]
		flagged, _ = flag_clicks(T, [0.2, 3.0, 1.0, 1.5], RX, FS; closure_max=2.0)
		@test flagged == Bool[0, 1, 0, 0]
		# TDOAs larger than the array allows (|τ| > 22.4 samples on the long side) are infeasible whatever the closure
		Tbad = [T; 45.0 -45.0 0.0]
		flagged, infeasible = flag_clicks(Tbad, zeros(5), RX, FS)
		@test infeasible == Bool[0, 0, 0, 0, 1]
		@test flagged == infeasible
		# no closure: only feasibility; NaN is flagged
		flagged, _ = flag_clicks([T; NaN NaN NaN], nothing, RX, FS)
		@test flagged == Bool[0, 0, 0, 0, 1]
		@test flag_clicks(T, [0.0, NaN, 0.0, 0.0], RX, FS)[1] == Bool[0, 1, 0, 0]
	end
end
