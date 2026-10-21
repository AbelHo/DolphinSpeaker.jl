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

	# simulated machine: `total` units of RAM, `base` held by the idle process, `other` by other programs;
	# a running item holds its weight. Capacity = total - other - base - margin, whatever is running
	function simulate(weights; total=20.0, base=1.0, other=3.0, margin=1.0, kwargs...)
		lk = ReentrantLock()
		held, peak, running, maxrunning, maxheld = Ref(0.0), Ref(base), Ref(0), Ref(0), Ref(0.0)
		done = Int[]
		rss() = base + held[]
		probes = (; available=() -> total - other - rss(), rss=rss, peak=() -> peak[],
			reset_peak=() -> (peak[] = rss(); true), release=() -> nothing)
		function f(i)
			lock(lk) do
				held[] += weights[i]; running[] += 1
				peak[] = max(peak[], rss()); maxrunning[] = max(maxrunning[], running[]); maxheld[] = max(maxheld[], held[])
			end
			sleep(0.05)
			lock(lk) do
				held[] -= weights[i]; running[] -= 1; push!(done, i)
			end
		end
		out = run_memory_adaptive(f, length(weights); weights=weights, margin_bytes=margin, safety=1.0, poll_s=0.01, probes=probes, kwargs...)
		return (; out..., maxrunning=maxrunning[], maxheld=maxheld[], done=sort(done))
	end

	@testset "run_memory_adaptive" begin
		w = [5.0, 4, 4, 3, 3, 2]
		# capacity 15: the largest runs alone first, then 4+4+3+3 = 14 together, and the last (2) would exceed 15 and waits
		r = simulate(w)
		@test r.done == 1:6
		@test r.per_weight ≈ 1.0
		@test r.maxheld == 14
		@test r.maxrunning == 4
		# other programs hold more RAM: capacity 4 -> one item at a time, all still processed
		r = simulate(w; other=14.0)
		@test r.done == 1:6
		@test r.maxrunning == 1
		# a cap on the items in flight (Inf = no cap)
		@test simulate(w; max_parallel=2).maxrunning == 2
		@test simulate(w; max_parallel=Inf, other=0.0).maxrunning >= 4
		@test simulate(w; max_parallel=1).maxrunning == 1
		# no measurement: everything starts at once
		r = simulate(w; adaptive=false)
		@test r.done == 1:6 && r.maxrunning == 6
		@test simulate(w; adaptive=false, max_parallel=3).maxrunning == 3
		# one item, and no items
		@test simulate([2.0]).done == [1]
		@test run_memory_adaptive(identity, 0).max_concurrent == 0
	end
end
