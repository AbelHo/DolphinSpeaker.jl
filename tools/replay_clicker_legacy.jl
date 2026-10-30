#= Replay the original D3 numerical path on authenticated, already detected clicks.

Usage: julia --project=. tools/replay_clicker_legacy.jl OUTPUT_ROOT
Source functions are taken verbatim as Julia syntax trees from the pre-calibration
Git snapshot. This does not rerun the detector or change click/frame identities.
=#
using CSV, DataFrames, JSON, SHA, LinearAlgebra, Optim, Test
include(joinpath(@__DIR__,"../src/clicker_calibration/ClickerCalibration.jl"))
const C = ClickerCalibration
const REPOSITORY = normpath(joinpath(@__DIR__,".."))
const SOURCE_COMMIT = "d22056c39b949ee5bead6ddd3b9644e4261ab267"
const SOURCE_FILES = ("src/localization.jl","src/dsp.jl","src/config.jl")

function definition_name(node)
    node isa Expr || return nothing
    node.head in (:function, :(=)) || return nothing
    head=node.args[1]
    while head isa Expr && head.head in (:where, :(::))
        head=head.args[1]
    end
    head isa Expr && head.head == :call && head.args[1] isa Symbol ? head.args[1] : nothing
end

const SOURCES=Dict(path=>read(`git -C $REPOSITORY show $SOURCE_COMMIT:$path`,String) for path in SOURCE_FILES)
const FROZEN=Module(:FrozenOriginalD3)
Core.eval(FROZEN, :(using Optim, LinearAlgebra))
Core.eval(FROZEN, quote
    default_tdoa2dir_solver=NelderMead
    default_window=-5000:10000
    default_c=1540.0
    ref_channel=3
    fname2dt_zoomf6(path)=nothing # unused timestamp callback; no PCM or time mapping
end)
const REQUIRED=Dict(
    "src/dsp.jl"=>Set([:extrema_and_indices]),
    "src/localization.jl"=>Set([:cost_tdoa2ang,:mod_sign,:tdoa2dir,:get_tdoa_minmax,:angle2px,:angle2pxd]),
    "src/config.jl"=>Set([:set_device__ophk_acoustic_D3_clicker]))
const DEFINITION_COUNTS=Dict{String,Any}()
for path in ("src/dsp.jl","src/localization.jl","src/config.jl")
    counts=Dict{String,Int}()
    for node in Meta.parseall(SOURCES[path]).args
        name=definition_name(node)
        name in REQUIRED[path] || continue
        # The unused old positional minmax convenience overload references
        # get_tdoa_max; retaining its definition is harmless and faithful.
        Core.eval(FROZEN,node)
        counts[string(name)]=get(counts,string(name),0)+1
    end
    Set(Symbol.(keys(counts))) == REQUIRED[path] || error("Missing frozen definitions in $path")
    DEFINITION_COUNTS[path]=counts
end
Core.eval(FROZEN, :(set_device__ophk_acoustic_D3_clicker()))

function replay_snippet(snip,fs,width,height)
    size(snip)==(1001,3) || error("Exact original window requires complete1001x3 snippet")
    # Original get_tdoa_minmax returns window indices, not relative delays.
    raw=FROZEN.get_tdoa_minmax(snip,[301];window=FROZEN.window_impulsive,ref_channel=3)
    angles,residual=FROZEN.tdoa2dir(raw,FROZEN.rx_vect[:,3:5],fs;return_residual=true)
    pixels=FROZEN.angle2px(angles)
    independent=optimize(x->FROZEN.cost_tdoa2ang(raw[1,:],x,FROZEN.rx_vect[:,3:5],fs),zeros(2),NelderMead())
    @assert isapprox(FROZEN.mod_sign.(Optim.minimizer(independent),2pi),vec(angles);atol=1e-12,rtol=0)
    delay=Float64.(vec(raw).-raw[1,3])
    geometry=permutedims(FROZEN.rx_vect[:,3:5])
    baseline=geometry[1:2,1:2].-permutedims(geometry[3,1:2])
    uv=-FROZEN.default_c .* (baseline \ (delay[1:2]./fs))
    feasible=sum(abs2,uv)<=1
    valid=all(isfinite,pixels)
    inside=valid && 0<=pixels[1,1]<width && 0<=pixels[1,2]<height
    converged=Optim.converged(independent)
    statuses=String[]
    valid || push!(statuses,"nonfinite_legacy_output")
    feasible || push!(statuses,"legacy_physically_inconsistent_tdoa")
    inside || push!(statuses,"outside_image")
    converged || push!(statuses,"legacy_solver_not_converged")
    (old_x=pixels[1,1],old_y=pixels[1,2],old_valid=valid,old_in_frame=inside,
     old_status=isempty(statuses) ? "ok" : join(statuses,";"),
     old_d1_samples=delay[1],old_d2_samples=delay[2],old_d3_samples=delay[3],
     old_azimuth_rad=angles[1,1],old_elevation_rad=angles[1,2],
     old_residual_samples2=residual[1],old_solver_converged=converged,
     old_physical_feasible=feasible)
end

function controls()
    @testset "Frozen legacy polarity, solver and declared image width" begin
        @test FROZEN.imsize==(1520,2740)
        @test FROZEN.fov_angle[1:2]==[80.0,54.46086642417594]
        @test FROZEN.impulsive_band_pass==[1000,Inf]
        @test FROZEN.window_impulsive==(-300:700)
        s=zeros(1001,3)
        for (j,at) in enumerate((305,299,301));s[at,j]=1;end
        p=replay_snippet(s,192000,2704,1520)
        @test (p.old_d1_samples,p.old_d2_samples,p.old_d3_samples)==(4.,-2.,0.)
        @test p.old_physical_feasible
        @test p.old_residual_samples2<1e-6
        @test FROZEN.angle2px(zeros(1,2))==[1370.0 760.0]
        # Preserve, rather than repair, the legacy positive-extremum behavior.
        for (j,at) in enumerate((400,450,500));s[at,j]=-10;end
        p=replay_snippet(s,192000,2704,1520)
        @test (p.old_d1_samples,p.old_d2_samples)==(4.,-2.)
        bad=zeros(1001,3);bad[100,1]=1;bad[900,2]=1;bad[301,3]=1
        p=replay_snippet(bad,192000,2704,1520)
        @test !p.old_physical_feasible
        @test p.old_valid # old optimizer emits a finite direction for this case
        @test occursin("physically_inconsistent",p.old_status)
    end
end

function main(output_root)
    controls()
    mkpath(output_root)
    snapshots=joinpath(output_root,"legacy_source_snapshot");mkpath(snapshots)
    for (path,source) in SOURCES
        write(joinpath(snapshots,basename(path)),source)
    end
    for device in ("A","B")
        original=joinpath(REPOSITORY,"calibration_runs/clicker_20260209/20261030_v2",device)
        manifest=C.read_json(joinpath(original,"annotations.complete.json"))
        @assert manifest["frames_sha256"]==C.sha256_file(joinpath(original,"frames.csv"))
        @assert manifest["annotations_sha256"]==C.sha256_file(joinpath(original,"annotations.csv"))
        selection=C.read_json(joinpath(original,"human_review_selection.json"))
        ids=String.(selection["accepted_frame_ids"])
        events=CSV.read(joinpath(original,"events.csv"),DataFrame)
        byframe=Dict(String(r.frame_id)=>NamedTuple(r) for r in eachrow(events))
        hashes=Dict(name=>C.sha256_file(joinpath(original,name)) for name in
                    ("events.csv","frames.csv","annotations.csv","annotations.complete.json","human_review_selection.json","prepare_manifest.json","event_snippets.f64"))
        rows=Vector{Any}(undef,length(ids))
        # Authentication helpers retain mutable verification caches: read them
        # serially, then parallelize only independent immutable computations.
        snippets=[C.read_event_snippet(original,byframe[id];verify=true) for id in ids]
        # Each event independent; use available native Julia threads.
        Threads.@threads for i in eachindex(ids)
            event=byframe[ids[i]]
            snip=snippets[i]
            p=replay_snippet(snip,Float64(event.fs),event.width,event.height)
            # Preparation's original minmax diagnostic must match frozen source.
            @assert [p.old_d1_samples,p.old_d2_samples,p.old_d3_samples] ==
                    [event.minmax_d1,event.minmax_d2,event.minmax_d3]
            rows[i]=merge((event_id=String(event.event_id),frame_id=ids[i],),p)
        end
        out=joinpath(output_root,device);mkpath(out)
        table=DataFrame(convert(Vector{typeof(first(rows))},rows))
        C._prep_atomic_csv(joinpath(out,"legacy_predictions.csv"),table)
        for (name,hash) in hashes
            @assert C.sha256_file(joinpath(original,name))==hash
        end
        C.write_json(joinpath(out,"legacy_manifest.json"),Dict(
            "status"=>"complete","device"=>device,"events"=>length(ids),"input_hashes"=>hashes,
            "legacy_source_commit"=>SOURCE_COMMIT,"source_sha256"=>Dict(k=>bytes2hex(sha256(v)) for (k,v) in SOURCES),
            "definition_counts"=>DEFINITION_COUNTS,"replay_source_sha256"=>C.sha256_file(@__FILE__),
            "csv_sha256"=>C.sha256_file(joinpath(out,"legacy_predictions.csv")),
            "method"=>"Exact frozen get_tdoa_minmax -> zero-initialized Optim.NelderMead tdoa2dir -> angle2px",
            "declared_imsize_height_width"=>[1520,2740],"actual_imsize_height_width"=>[1520,2704],
            "fov_degrees"=>FROZEN.fov_angle[1:2],"sound_speed_m_s"=>1540,
            "reference_physical_channel"=>5,"window"=>[-300,700],
            "conditioning"=>"Exact authenticated whole-file HP1kHz fourth-order filtfilt prepared snippet; no refiltering",
            "old_finite"=>count(table.old_valid),"old_in_frame"=>count(table.old_in_frame),
            "old_physically_feasible"=>count(table.old_physical_feasible),
            "old_solver_converged"=>count(table.old_solver_converged),
            "limitations"=>["Same prepared detected click population and sync mapping, not a rerun of original detection/train selection.",
                "Original2740 width intentionally retained, no rescaling to actual2704 image.",
                "Legacy finite outputs for impossible delays retained and explicitly labelled.",
                "Optim runtime is current installed version; original code and defaults are frozen, historical Julia environment was not reconstructed."]))
        println("READY ",device," ",joinpath(out,"legacy_predictions.csv")," count=",nrow(table))
    end
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(abspath(only(ARGS)))
end
