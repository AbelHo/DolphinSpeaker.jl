# Orchestration controls only: stage doubles record calls, never claim scientific
# correctness, render media, fit a model, spawn agents, or contact Telegram.
module FullWorkflowTests
using Test,JSON,SHA,Dates
const SOURCE_DIR=joinpath(@__DIR__,"..","..","src","clicker_calibration")
include(joinpath(SOURCE_DIR,"support.jl"))
include(joinpath(SOURCE_DIR,"workflow.jl"))
include(joinpath(SOURCE_DIR,"full_workflow.jl"))

function fixture(root;ready=Dict("A"=>true,"B"=>true),fail=nothing,prepare_complete=true,pause_ready=false)
    trace=NamedTuple[];clock=Ref(0.0)
    record(stage,args...;kwargs...)=begin
        push!(trace,(;stage,args,kwargs=(;kwargs...)))
        stage==fail && error("injected orchestration failure: $stage")
        Dict("status"=>"complete")
    end
    operations=(
        prepare=(source,destination;kwargs...)->begin
            record(:prepare,source,destination;kwargs...);mkpath(destination)
            Dict("status"=>prepare_complete ? "complete" : "bounded_or_unmapped")
        end,
        sync=(args...;kwargs...)->record(:sync,args...;kwargs...),
        annotate=(args...;kwargs...)->begin
            record(:annotate,args...;kwargs...);Dict("status"=>"launched")
        end,
        annotation_current=dir->get(ready,basename(dir),false),
        validate_annotation=(args...;kwargs...)->record(:validate,args...;kwargs...),
        gallery=(args...;kwargs...)->record(get(kwargs,:fitted,false) ? :comparison_gallery : :annotation_gallery,args...;kwargs...),
        video=(args...;kwargs...)->record(kwargs[:annotation_only] ? :annotation_video : :comparison_video,args...;kwargs...),
        package=(original,destination;kwargs...)->begin
            final=kwargs[:final];record(final ? :final_package : :initial_package,original,destination;kwargs...)
            Dict("status"=>final ? "complete" : "review_available","destination"=>destination)
        end,
        notify=(delivery;kwargs...)->record(kwargs[:phase]==:annotation ? :initial_notify : :final_notify,delivery;kwargs...),
        fit_input=(args...;kwargs...)->record(:adapter,args...;kwargs...),
        policy=(args...;kwargs...)->record(:native_policy,args...;kwargs...),
        fit=(args...;kwargs...)->record(:fit,args...;kwargs...),
        review=(args...;kwargs...)->record(:review,args...;kwargs...),
        refine=(args...;kwargs...)->record(:refine,args...;kwargs...),
        project=(args...;kwargs...)->record(:project,args...;kwargs...),
        diagnostics=(args...;kwargs...)->record(:diagnostics,args...;kwargs...),
        clock=()->clock[],pause=seconds->begin
            clock[]+=seconds
            pause_ready && foreach(name->ready[name]=true,keys(ready))
        end)
    common=(output_dir=joinpath(root,"work"),results_dir=joinpath(root,"delivery"),
        error_log=joinpath(root,"expected_errors.jsonl"),operations=operations)
    return (;trace,clock,operations,common,input=joinpath(root,"source"))
end

@testset "One-call order delivers original review before fitting" begin
    mktempdir() do root
        f=fixture(root)
        result=run_full_clicker_calibration(f.input;f.common...,devices=(:A,:B),notify=:none,
            annotation_fit_policy=Dict("A"=>:user_confirmed_uncertain_coordinates),
            observation_policy=Dict(:B=>:terra_estimated_underwater_coordinates),
            fit_input_root=joinpath(root,"fit_inputs"),prepare_options=(ratio=12.0,),
            fit_options=(nstarts=6,),video_options=(threads=0,encoder="auto",decoder="auto"))
        @test result["status"]=="complete"
        @test all(d["status"]=="complete" for d in values(result["devices"]))
        stages=getproperty.(f.trace,:stage)
        @test count(==(:prepare),stages)==2
        @test count(==(:annotate),stages)==2
        @test maximum(findall(==(:initial_notify),stages))<minimum(findall(==(:fit),stages))
        @test maximum(findall(==(:annotate),stages))<minimum(findall(==(:annotation_gallery),stages))
        @test count(==(:adapter),stages)==1 && count(==(:native_policy),stages)==1
        for name in ("A","B")
            original=joinpath(root,"work",name)
            fit=name=="A" ? joinpath(root,"fit_inputs",name) : original
            science=[row.stage for row in f.trace if row.stage in (:fit,:review,:refine) && first(row.args)==fit]
            @test science==[:fit,:review,:refine]
            projection=only(filter(row->row.stage==:project && first(row.args)==original,f.trace))
            @test projection.args[2]==joinpath(fit,"refinement_v1")
            final=only(filter(row->row.stage==:final_package && first(row.args)==original,f.trace))
            @test final.kwargs.final && final.kwargs.model_dir==joinpath(fit,"refinement_v1")
        end
        videos=filter(row->row.stage in (:annotation_video,:comparison_video),f.trace)
        @test length(videos)==4
        @test all(row.kwargs.threads==0 && row.kwargs.encoder=="auto" && row.kwargs.decoder=="auto" for row in videos)
        @test all(row.kwargs.resume for row in videos)
        @test findfirst(==(:project),stages)<findfirst(==(:comparison_gallery),stages)<findfirst(==(:comparison_video),stages)<findfirst(==(:final_package),stages)<findfirst(==(:final_notify),stages)
        @test read_json(joinpath(root,"work","full_workflow_manifest.json"))["status"]=="complete"
    end
end

@testset "Awaiting annotation, timeout, and preparation failure remain explicit" begin
    mktempdir() do root
        f=fixture(root;ready=Dict("A"=>true,"B"=>false))
        result=run_full_clicker_calibration(f.input;f.common...,wait_annotations=false,notify=:none)
        @test result["status"]=="awaiting_annotation"
        @test result["pending_devices"]==["B"]
        @test result["devices"]["A"]["status"]=="complete"
        @test result["devices"]["B"]["status"]=="awaiting_annotation"
        @test count(row->row.stage==:fit,f.trace)==1
        @test f.clock[]==0
    end
    mktempdir() do root
        f=fixture(root;ready=Dict("A"=>false),pause_ready=true)
        result=run_full_clicker_calibration(f.input;f.common...,devices=(:A,),notify=:none)
        @test result["status"]=="complete" && f.clock[]==5
        @test count(row->row.stage==:prepare,f.trace)==1
    end
    mktempdir() do root
        f=fixture(root;ready=Dict("A"=>false))
        result=run_full_clicker_calibration(f.input;f.common...,devices=(:A,),notify=:none,annotation_timeout_s=5)
        @test result["status"]=="awaiting_annotation" && f.clock[]==5
        @test !any(row->row.stage==:fit,f.trace)
    end
    mktempdir() do root
        f=fixture(root;prepare_complete=false)
        result=run_full_clicker_calibration(f.input;f.common...,devices=(:A,),notify=:none)
        @test result["status"]=="needs_review"
        @test result["devices"]["A"]["status"]=="needs_preparation_review"
        @test getproperty.(f.trace,:stage)==[:prepare]
    end
end

@testset "Frozen delivery bypasses every scientific mutation stage" begin
    mktempdir() do root
        f=fixture(root);original=joinpath(root,"existing","A");model=joinpath(root,"old_fit","refinement_v1")
        mkpath(original);mkpath(model)
        for name in ("calibration.json","predictions.csv","report.json");write(joinpath(model,name),"frozen fixture $name\n");end
        write(joinpath(original,"annotations.csv"),"unchanged original Terra fixture\n")
        paths=vcat([joinpath(model,n) for n in readdir(model)],[joinpath(original,"annotations.csv")])
        before=Dict(path=>sha256_file(path) for path in paths)
        result=run_full_clicker_calibration(f.input;f.common...,devices=(:A,),notify=:none,delivery_only=true,
            original_dirs=Dict(:A=>original),model_dirs=Dict("A"=>model))
        @test result["status"]=="complete"
        @test !any(row->row.stage in (:prepare,:sync,:annotate,:adapter,:native_policy,:fit,:review,:refine),f.trace)
        @test only(filter(row->row.stage==:project,f.trace)).args==(original,model)
        @test Dict(path=>sha256_file(path) for path in paths)==before
        @test count(row->row.stage==:initial_package,f.trace)==1
        @test count(row->row.stage==:final_package,f.trace)==1
    end
end

@testset "Failures stop dependent stages and protect old-source science" begin
    for failure in (:initial_notify,:fit,:project,:comparison_video)
        mktempdir() do root
            f=fixture(root;fail=failure)
            @test_throws ErrorException run_full_clicker_calibration(f.input;f.common...,devices=(:A,),notify=:telegram)
            state=read_json(joinpath(root,"work","full_workflow_manifest.json"))
            @test state["status"]=="failed" && state["failed_device"]=="A"
            @test isfile(joinpath(root,"expected_errors.jsonl"))
            @test !any(row->row.stage==:final_package,f.trace)
            failure==:initial_notify && @test !any(row->row.stage==:fit,f.trace)
        end
    end
    mktempdir() do root
        f=fixture(root);original=joinpath(root,"work","A");mkpath(original)
        report=joinpath(original,"fit_report.json")
        write_json(report,Dict("hashes"=>Dict("physical_fit_source_sha256"=>"older source")))
        before=sha256_file(report)
        @test_throws ErrorException run_full_clicker_calibration(f.input;f.common...,devices=(:A,),notify=:none)
        @test isempty(f.trace) && sha256_file(report)==before
    end
    mktempdir() do root
        f=fixture(root)
        @test_throws ArgumentError run_full_clicker_calibration(f.input;f.common...,devices=(:A,),
            annotation_fit_policy=:user_confirmed_uncertain_coordinates,observation_policy=:terra_estimated_underwater_coordinates)
        @test_throws ArgumentError run_full_clicker_calibration(f.input;f.common...,devices=(:A,),model_dirs="existing")
        @test_throws ErrorException run_full_clicker_calibration(f.input;output_dir=joinpath(f.input,"A","output"),results_dir=f.common.results_dir,operations=f.operations)
        @test_throws ArgumentError _full_operations((unknown=identity,))
    end
end

@testset "Dataset-root results are allowed; recordings stay protected" begin
    mktempdir() do root
        f=fixture(root)
        result=run_full_clicker_calibration(f.input;output_dir=f.common.output_dir,
            results_dir=joinpath(f.input,"results"),devices=(:A,),notify=:none,
            error_log=f.common.error_log,operations=f.operations)
        @test result["status"]=="complete"
        @test result["devices"]["A"]["destination"]==joinpath(f.input,"results","A")
        @test all(row.args[2]==joinpath(f.input,"results","A","videos") for row in f.trace if row.stage in (:annotation_video,:comparison_video))
        @test_throws ErrorException run_full_clicker_calibration(f.input;output_dir=f.common.output_dir,
            results_dir=joinpath(f.input,"B","results"),devices=(:A,),operations=f.operations)
    end
end

@testset "Real report source-key contracts are recognized without rerunning science" begin
    mktempdir() do root
        fit=joinpath(root,"fit");model=joinpath(fit,"refinement_v1");mkpath(model)
        first=Dict("physical_fit_source_sha256"=>sha256_file(joinpath(SOURCE_DIR,"physical_fit.jl")),
            "support_source_sha256"=>sha256_file(joinpath(SOURCE_DIR,"support.jl")))
        refined=Dict("source/"*name=>sha256_file(joinpath(SOURCE_DIR,name)) for name in
            ("physical_fit.jl","prepare.jl","srp_review.jl","method_refinement.jl","support.jl"))
        write_json(joinpath(fit,"fit_report.json"),Dict("hashes"=>first))
        write_json(joinpath(model,"report.json"),Dict("hashes"=>refined))
        @test _full_existing_science(fit,model)
        refined["source/prepare.jl"]="stale"
        write_json(joinpath(model,"report.json"),Dict("hashes"=>refined))
        @test_throws ErrorException _full_existing_science(fit,model)
    end
end

@testset "Mode-aware video resume restores pointers and still invokes exporter" begin
    mktempdir() do root
        original=joinpath(root,"original");model=joinpath(root,"model");output=joinpath(root,"results","A","videos")
        mkpath(original);mkpath(model);mkpath(output)
        write(joinpath(original,"annotations.csv"),"original truth fixture")
        write(joinpath(model,"calibration.json"),"frozen model fixture")
        video_paths=Dict();manifests=Dict()
        for annotation_only in (true,false)
            kind=annotation_only ? "annotation_only" : "calibrated_comparison"
            video=joinpath(output,kind*".mp4");audio=joinpath(output,kind*".wav")
            write(video,"media pointer fixture $kind");write(audio,"audio pointer fixture $kind")
            hashes=Dict(joinpath(original,"annotations.csv")=>sha256_file(joinpath(original,"annotations.csv")))
            annotation_only || (hashes[joinpath(model,"calibration.json")]=sha256_file(joinpath(model,"calibration.json")))
            manifest=Dict("schema_version"=>"clicker-video-export/1.0","status"=>"complete","mode"=>kind,
                "settings"=>Dict("annotation_only"=>annotation_only),"inputs_sha256"=>hashes,
                "video_path"=>video,"reference_audio_path"=>audio)
            name=annotation_only ? "review_export_manifest.json" : "comparison_export_manifest.json"
            path=joinpath(output,name);write_json(path,manifest);manifests[annotation_only]=manifest;video_paths[annotation_only]=path
        end
        pointer=joinpath(output,"export_manifest.json")
        cp(video_paths[false],pointer)
        commands=Vector{String}[]
        runner=(command;operation,capture)->begin
            push!(commands,collect(command.exec))
            desired="--annotation-only" in command.exec
            # This double checks orchestration only; the real exporter performs
            # source, settings, media and checksum checks after pointer restore.
            @test read(pointer)==read(video_paths[desired])
            Vector{UInt8}(codeunits("{\"resource_selection\":{}}\n"*JSON.json(manifests[desired],2)))
        end
        review=_full_export_video(original,output;annotation_only=true,python="python3",runner=runner)
        @test review==manifests[true]
        final=_full_export_video(original,output;model_dir=model,python="python3",runner=runner)
        @test final==manifests[false]
        @test length(commands)==2 # Restoring a pointer never skips actual exporter checks.
        @test all(args[findfirst(==("--threads"),args)+1]=="0" for args in commands)
        @test all(args[findfirst(==("--encoder"),args)+1]=="auto" for args in commands)
        @test all(args[findfirst(==("--decoder"),args)+1]=="auto" for args in commands)
        before=read(pointer)
        stale=deepcopy(manifests[true]);stale["inputs_sha256"][joinpath(original,"annotations.csv")]="stale"
        write_json(video_paths[true],stale)
        @test !_full_restore_video_manifest(original,output;annotation_only=true)
        @test read(pointer)==before
        wrong=deepcopy(manifests[false]);wrong["video_path"]=joinpath(root,"outside.mp4");write(wrong["video_path"],"outside fixture")
        @test !_full_video_matches(wrong,original,output;model_dir=model)
        # Captured stdout can contain multiple JSON records: the atomic file
        # supplies the return value, while wrong-mode output still fails closed.
        @test_throws ErrorException _full_export_video(original,output;annotation_only=true,resume=false,
            runner=(args...;kwargs...)->UInt8[])
    end
end
end # module
