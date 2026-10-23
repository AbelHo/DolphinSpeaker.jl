# Unit tests for src/source_separation.jl (synthetic detection streams; no recordings needed).
# run: julia --project=. -t auto test/test_source_separation.jl
# details: improvements/tdoa_localization_2026-10-22/reports/round4_realism_sources.md
using Test, Random, Statistics, LinearAlgebra, DataFrames, Clustering
include(joinpath(@__DIR__, "..", "src", "source_separation.jl"))

# k sources clicking at their own ICI, directions spread in azimuth and drifting slowly; optional echo followers
# (a vertical image 0.6-1.2 ms after each click) and clutter detections from random directions
function synth_detections(rng, k; dur=3.0, noise_deg=1.0, clutter=0.03, echoes=false)
    ici = [10e-3, 13e-3, 17e-3, 23e-3, 29e-3][1:k]
    base = [(-40 + 80 * (s - 1) / max(k - 1, 1), -8.0 + 16 * (s - 1) / max(k - 1, 1)) for s in 1:k]
    rows = Tuple{Float64,Float64,Float64,Int,Bool}[]   # t, az, el, source (0 = clutter), echo
    for s in 1:k
        t = rand(rng) * ici[s]
        while t < dur
            az = base[s][1] + 4.0 * (-1)^s * t + noise_deg * randn(rng); el = base[s][2] + noise_deg * randn(rng)
            push!(rows, (t, az, el, s, false))
            echoes && push!(rows, (t + 6e-4 + 6e-4rand(rng), az + 0.5randn(rng), el + 35 + randn(rng), s, true))
            t += ici[s] * (1 + 0.05randn(rng))
        end
    end
    for _ in 1:round(Int, clutter * length(rows))
        push!(rows, (dur * rand(rng), 120rand(rng) - 60, 80rand(rng) - 40, 0, false))
    end
    sort!(rows; by=first)
    return (t=first.(rows), az=[r[2] for r in rows], el=[r[3] for r in rows], src=[r[4] for r in rows], echo=[r[5] for r in rows])
end

@testset "source separation" begin
    @testset "evaluation helpers" begin
        t = collect(0:0.01:1); lab = repeat([1, 2], 51)[1:length(t)]
        @test pairwise_scores(lab, lab, t).f1 ≈ 1
        @test pairwise_scores(lab, 3 .- lab, t).f1 ≈ 1                  # permutation-invariant
        @test pairwise_scores(lab, ones(Int, length(t)), t).precision < 0.6
        @test adjusted_rand(lab, lab) ≈ 1
        @test adjusted_rand(lab, 3 .- lab) ≈ 1
    end

    for k in (2, 3), n in (k, nothing)
        @testset "$k sources, n = $(something(n, "estimated"))" begin
            d = synth_detections(Xoshiro(10k), k)
            s = separate_sources(d.t, d.az, d.el; n=n)
            real = d.src .> 0
            sc = pairwise_scores(d.src, s.label, d.t; max_dt=0.2)
            @test sc.precision >= 0.95
            @test sc.recall >= 0.9
            @test adjusted_rand(d.src, s.label) >= 0.85
            isnothing(n) || @test all(0 .<= s.label .<= n)
            # most clutter stays unassigned or in a source of its own
            @test count(s.label[.!real] .== 0) >= 0.5 * count(.!real)
            if isnothing(n)
                @test count(s.sources.confirmed) <= k + 2                 # the estimated number of sources is close to k
            end
        end
    end

    @testset "followers: echoes inherit, the other dolphin keeps its own direction" begin
        d = synth_detections(Xoshiro(5), 2; echoes=true, clutter=0.0)
        s = separate_sources(d.t, d.az, d.el; n=2)
        ech = findall(d.echo)
        @test count(s.role[ech] .== "echo") >= 0.9 * length(ech)
        # an echo takes the direction of its leader: the source's own click just before it
        ok = [abs(s.el[i] - d.el[i]) > 20 for i in ech]                    # moved off the image (35° above)
        @test count(ok) >= 0.9 * length(ech)
        # clicks of the other source that follow a click within 1.6 ms keep their own direction and source
        fol = [i for i in 2:length(d.t) if !d.echo[i] && (d.t[i] - d.t[i-1]) * 1000 < 1.6]
        @test length(fol) >= 5
        own_ok = [s.role[i] in ("own", "repeat") && s.az[i] == d.az[i] for i in fol]
        @test count(own_ok) >= 0.8 * length(fol)
        sc = pairwise_scores(d.src, s.label, d.t; max_dt=0.2)
        @test sc.precision >= 0.95 && sc.recall >= 0.9
    end

    @testset "k-means baseline on separated sources" begin
        # the baseline separates the event leaders; its followers inherit naively (the other dolphin's clicks within
        # 1.6 ms of a click get the wrong label), which is what separate_sources' follower handling improves on
        d = synth_detections(Xoshiro(3), 2; clutter=0.0)
        lab = kmeans_sources(d.t, d.az, d.el, 2)
        lead = [i == 1 || (d.t[i] - d.t[i-1]) * 1000 >= 1.6 for i in eachindex(d.t)]
        @test adjusted_rand(d.src[lead], lab[lead]) >= 0.9
        s = separate_sources(d.t, d.az, d.el; n=2)
        @test adjusted_rand(d.src, s.label) > adjusted_rand(d.src, lab)
    end
end
