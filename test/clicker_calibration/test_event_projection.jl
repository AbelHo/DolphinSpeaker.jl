# Bounded synthetic replay controls; no optimization, real input, or waveform job.
module EventProjectionTests
using Test,CSV,DataFrames,JSON,SHA,Dates,Statistics,LinearAlgebra,Random,Optim,DSP,FFTW,WAV,Peaks
const SOURCE_DIR=joinpath(@__DIR__,"..","..","src","clicker_calibration")
include(joinpath(SOURCE_DIR,"support.jl"))
include(joinpath(SOURCE_DIR,"prepare.jl"))
include(joinpath(SOURCE_DIR,"physical_fit.jl"))
include(joinpath(SOURCE_DIR,"srp_review.jl"))
include(joinpath(SOURCE_DIR,"method_refinement.jl"))
include(joinpath(SOURCE_DIR,"workflow.jl"))
include(joinpath(SOURCE_DIR,"event_projection.jl"))

function fixture(parent;refined=true,method="gcc")
    original=joinpath(parent,"original");fit=joinpath(parent,"fit")
    mkpath(original);mkpath(fit)
    events=DataFrame(event_id=["e1","e2","e3","e4"],frame_id=["f1","f2","f3",""],
        device=fill("B",4),fs=fill(192000.0,4),recording_time_s=[3.,30.5,63.,93.],
        sample_index=[576000,5856000,12096000,17856000],audio_time_s=[3.,30.5,63.,93.],
        audio_file=fill("synthetic.wav",4),window_status=fill("complete",4),window_source_segments=fill("[]",4),
        sync_status=["mapped","mapped","mapped","outside_video"],
        width=fill(1280,4),height=fill(960,4),snippet_sha256=["snippet_$i" for i in 1:4])
    cache=nothing
    if method=="srp"
        events.snippet_file=fill("event_snippets.f64",4)
        events.snippet_byte_offset=[0,1,2,4].*(1001*3*8) # Last row intentionally has no available payload.
        events.snippet_samples=fill(1001,4);events.snippet_channels=fill(3,4);events.snippet_peak_index=fill(300,4)
        cache_path=joinpath(original,"event_snippets.f64")
        open(cache_path,"w") do io
            for i in 1:4
                snippet=zeros(1001,3)
                i==3 || (snippet[301,:].=1.0)
                bytes=_prep_snippet_bytes(snippet)
                events.snippet_sha256[i]=bytes2hex(SHA.sha256(bytes));write(io,bytes)
            end
        end
        cache=Dict("path"=>cache_path,"sha256"=>sha256_file(cache_path),"bytes"=>filesize(cache_path),
            "dtype"=>"Float64","endianness"=>"little","layout"=>"column-major")
    end
    for estimator in ("minmax2","gcc")
        events[!,estimator*"_d1"]=[0.,3.,900.,1.]
        events[!,estimator*"_d2"]=[0.,1.,0.,2.]
        events[!,estimator*"_d3"]=zeros(4)
    end
    frames=DataFrame(frame_id=["f1","f2","f3"],width=fill(1280,3),height=fill(960,3))
    annotations=DataFrame(frame_id=frames.frame_id,x=[640.,NaN,900.],y=[480.,NaN,500.],
        visibility=["visible","absent","uncertain"],medium=["underwater","air","underwater"],
        uncertainty_px=[2.,NaN,30.],annotator=fill("synthetic",3))
    CSV.write(joinpath(original,"events.csv"),events);CSV.write(joinpath(original,"frames.csv"),frames)
    CSV.write(joinpath(original,"annotations.csv"),annotations)
    write_json(joinpath(original,"annotations.complete.json"),Dict("frames_sha256"=>sha256_file(joinpath(original,"frames.csv")),
        "annotations_sha256"=>sha256_file(joinpath(original,"annotations.csv"))))
    prepared=Dict{String,Any}("status"=>"complete",
        "request"=>Dict("source_sha256"=>sha256_file(joinpath(SOURCE_DIR,"prepare.jl"))),
        "outputs"=>[Dict("path"=>joinpath(original,n),"sha256"=>sha256_file(joinpath(original,n))) for n in ("events.csv","frames.csv")])
    if cache!==nothing
        prepared["snippet_cache"]=cache
        push!(prepared["outputs"],Dict("path"=>cache["path"],"sha256"=>cache["sha256"]))
    end
    write_json(joinpath(original,"prepare_manifest.json"),prepared)
    input_hashes=Dict(key=>sha256_file(joinpath(original,name)) for (name,key) in
        (("events.csv","events_sha256"),("frames.csv","frames_sha256"),("annotations.csv","annotations_sha256"),("annotations.complete.json","annotations_complete_sha256")))
    q=_pf_initial(1280,960,"centered_pinhole")
    baseline_cal=_pf_calibration(q,1280,960,"centered_pinhole";device="B",method="minmax2")
    baseline_cal["hashes"]=input_hashes
    baseline=DataFrame(event_id=events.event_id,frame_id=events.frame_id,
        split=["train","purged","ineligible","ineligible"],eligibility_status=["eligible","boundary_purge","annotation_uncertain","no_annotation"])
    CSV.write(joinpath(fit,"predictions.csv"),baseline)
    write_json(joinpath(fit,"calibration.json"),baseline_cal)
    write_json(joinpath(fit,"heldout_evaluation.json"),Dict("hashes"=>input_hashes,"evaluation_calibration"=>baseline_cal))
    write_json(joinpath(fit,"fit_report.json"),Dict("calibration_sha256"=>sha256_file(joinpath(fit,"calibration.json")),
        "predictions_sha256"=>sha256_file(joinpath(fit,"predictions.csv")),
        "heldout_evaluation_sha256"=>sha256_file(joinpath(fit,"heldout_evaluation.json"))))
    model=fit
    if refined
        model=joinpath(fit,"refinement_v1");mkpath(model)
        bindings=Dict(name=>sha256_file(joinpath(original,name)) for name in ("events.csv","frames.csv","annotations.csv","annotations.complete.json"))
        for name in ("predictions.csv","heldout_evaluation.json")
            bindings[name]=sha256_file(joinpath(fit,name))
        end
        if method=="srp"
            review=DataFrame(event_id=["e1"],frame_id=["f1"],fs=[192000.],method=["srp"],
                snippet_sha256=[events.snippet_sha256[1]],d1=[0.],d2=[0.],d3=[0.],estimator_valid=[true])
            CSV.write(joinpath(fit,"estimator_review.csv"),review)
            write_json(joinpath(fit,"estimator_review.json"),Dict("estimator_review_csv_sha256"=>sha256_file(joinpath(fit,"estimator_review.csv")),
                "hashes"=>Dict(name=>bindings[name] for name in ("events.csv","frames.csv"))))
            for name in ("estimator_review.csv","estimator_review.json");bindings[name]=sha256_file(joinpath(fit,name));end
            bindings["event_snippets.f64"]=cache["sha256"]
        end
        cal=_pf_calibration(q,1280,960,"centered_pinhole";device="B",method=method)
        cal["estimator_spec"]=_mr_estimator_spec(method;fs=192000.)
        cal["hashes"]=bindings
        write_json(joinpath(model,"evaluation.json"),Dict("evaluation_calibration"=>cal))
        cal["evaluation_sha256"]=sha256_file(joinpath(model,"evaluation.json"))
        write_json(joinpath(model,"calibration.json"),cal)
        write_json(joinpath(model,"selection.json"),Dict("frozen"=>true))
        CSV.write(joinpath(model,"predictions.csv"),baseline[1:1,:])
        write_json(joinpath(model,"report.json"),Dict("hashes"=>bindings,
            "outputs_sha256"=>Dict(name=>sha256_file(joinpath(model,name)) for name in ("calibration.json","predictions.csv","evaluation.json","selection.json"))))
    end
    return (;original,fit,model)
end

function frozen_files(root)
    Dict(joinpath(directory,name)=>sha256_file(joinpath(directory,name)) for (directory,_,names) in walkdir(root)
        for name in names if !(name in ("inference_events.csv","inference_events.json","expected_errors.jsonl")))
end

@testset "Every candidate projects before original truth is joined" begin
    mktempdir() do root
        CALIBRATION_ERROR_LOG[]=joinpath(root,"expected_errors.jsonl")
        f=fixture(root);before=frozen_files(root)
        result=project_calibration_events(f.original,f.model)
        table=CSV.read(joinpath(f.model,"inference_events.csv"),DataFrame;stringtype=String)
        @test result["schema"]=="clicker-all-event-projection-v1"
        @test result["counts"]["n_original_events"]==4
        @test table.event_id==["e1","e2","e3","e4"]
        @test table.final_valid==[true,true,false,true]
        @test table.split[2]=="purged"
        @test table.visibility[2]=="absent"
        @test table.final_valid[2] # Neither absent nor boundary purging filters inference.
        @test table.sync_status[4]=="outside_video" && table.final_valid[4]
        @test table.final_status[3]=="impossible_tdoa" && isnan(table.final_x[3])
        @test table.original_visibility[3]=="uncertain"
        @test table.original_uncertainty_px[3]==30.0
        @test isequal(table.final_x,table.predicted_x)
        @test result["semantics"]["refit"]===false
        @test result["semantics"]["annotation_selects_or_changes_ray"]===false
        @test frozen_files(root)==before
        # Exercise the explicit phase boundary: arbitrary GT edits cannot enter
        # the ray producer, change validity, or remove/reorder an event.
        input=_ep_inputs(f.original,f.model)
        original=_ep_project(input.events,input.calibration,input.evaluation_calibration,input.saved_srp)
        changed=deepcopy(input.annotations)
        changed.x.=[100.,200.,300.];changed.y.=[50.,60.,70.]
        changed.visibility.= "absent";changed.medium.="air";changed.uncertainty_px.=[99.,88.,77.]
        joined=_ep_join_truth!(deepcopy(original),merge(input,(annotations=changed,)))
        @test joined.annotation_x[1]==100.0
        @test joined.final_error_px[1]!=hypot(original.final_x[1]-640,original.final_y[1]-480)
        for field in (:event_id,:final_x,:final_y,:final_valid,:final_status,:final_in_frame,:evaluation_x,:baseline_x)
            @test isequal(joined[!,field],original[!,field])
        end
        @test frozen_files(root)==before
    end
end

@testset "Resume, stale inputs, tampering, and immutable frozen products" begin
    mktempdir() do root
        CALIBRATION_ERROR_LOG[]=joinpath(root,"expected_errors.jsonl")
        f=fixture(root)
        first=project_calibration_events(f.original,f.model)
        manifest=joinpath(f.model,"inference_events.json");digest=sha256_file(manifest)
        @test project_calibration_events(f.original,f.model)==first
        @test sha256_file(manifest)==digest
        open(joinpath(f.model,"inference_events.csv"),"a") do io;write(io,"\n");end
        @test_throws ErrorException project_calibration_events(f.original,f.model)
        @test project_calibration_events(f.original,f.model;resume=false)["status"]=="complete"
        path=joinpath(f.original,"annotations.csv");bytes=read(path)
        open(path,"a") do io;write(io,"\n");end
        @test_throws ErrorException project_calibration_events(f.original,f.model;resume=false)
        write(path,bytes)
        path=joinpath(f.model,"predictions.csv");bytes=read(path)
        open(path,"a") do io;write(io,"\n");end
        @test_throws ErrorException project_calibration_events(f.original,f.model;resume=false)
        write(path,bytes)
        previous=read_json(manifest);previous["source_hashes"]["event_projection.jl"]="tampered"
        write_json(manifest,previous)
        @test_throws ErrorException project_calibration_events(f.original,f.model)
        @test project_calibration_events(f.original,f.model;resume=false)["status"]=="complete"
        # Model/estimator source tampering cannot be blessed by sidecar regeneration.
        calpath=joinpath(f.model,"calibration.json");cal=read_json(calpath)
        cal["estimator_spec"]["prepare_source_sha256"]="wrong"
        write_json(calpath,cal)
        reportpath=joinpath(f.model,"report.json");report=read_json(reportpath)
        report["outputs_sha256"]["calibration.json"]=sha256_file(calpath);write_json(reportpath,report)
        @test_throws ErrorException project_calibration_events(f.original,f.model;resume=false)
    end
end

@testset "Baseline mode and fixed SRP inference beyond frozen cohort" begin
    mktempdir() do root
        CALIBRATION_ERROR_LOG[]=joinpath(root,"expected_errors.jsonl")
        f=fixture(joinpath(root,"baseline");refined=false)
        result=project_calibration_events(f.original,f.model)
        @test result["method"]=="minmax2"
        @test result["counts"]["n_original_events"]==4
        f=fixture(joinpath(root,"srp");method="srp")
        before=frozen_files(root)
        result=project_calibration_events(f.original,f.model)
        table=CSV.read(joinpath(f.model,"inference_events.csv"),DataFrame)
        @test nrow(table)==4
        @test table.final_valid==[true,true,false,false]
        @test table.td_source[1:3]==["cached_review","fixed_estimator_cached_snippet","fixed_estimator_cached_snippet"]
        @test table.split[2]=="purged" && table.visibility[2]=="absent"
        @test table.final_status[3]=="fixed_srp_silent_channel"
        @test table.final_status[4]=="invalid_cached_snippet"
        @test all(isnan,table.final_x[3:4])
        # Exact fixed estimator, no fitting or coordinate-dependent steering.
        input=_ep_inputs(f.original,f.model)
        snippet=read_event_snippet(f.original,input.events[2,:];verify=true)
        raw=review_srp_phat(snippet;fs=192000.)
        @test [table.d1_samples[2],table.d2_samples[2],table.d3_samples[2]]==raw.tdoas_samples
        @test frozen_files(root)==before
        @test result["counts"]["projection_status"]["invalid_cached_snippet"]==1
        # An absent cache preserves all events and cached-review predictions.
        cachepath=joinpath(f.original,"event_snippets.f64");bytes=read(cachepath);rm(cachepath)
        missing=project_calibration_events(f.original,f.model;resume=false)
        @test missing["counts"]["projection_status"]["missing_cached_snippet"]==3
        write(cachepath,bytes)
        # A present but modified cache is a provenance error, not an invalid ray.
        open(cachepath,"a") do io;write(io,UInt8(0));end
        @test_throws ErrorException project_calibration_events(f.original,f.model;resume=false)
    end
end
end # module
