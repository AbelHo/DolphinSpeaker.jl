# Native synthetic controls only. These do not execute the real A plan.
module MethodRefinementTests
using Test,CSV,DataFrames,JSON,SHA,Dates,Statistics,LinearAlgebra,Random,Optim,DSP,FFTW,WAV,Peaks
const SOURCE_DIR=joinpath(@__DIR__,"..","..","src","clicker_calibration")
include(joinpath(SOURCE_DIR,"support.jl"))
include(joinpath(SOURCE_DIR,"physical_fit.jl"))
include(joinpath(SOURCE_DIR,"prepare.jl"))
include(joinpath(SOURCE_DIR,"srp_review.jl"))
include(joinpath(SOURCE_DIR,"method_refinement.jl"))

function synthetic_wave(u;fs=192000.0,c=1500.0,offsets=[2e-6,-1e-6,0.0])
    direction=[u[1],u[2],sqrt(1-sum(abs2,u))]
    t=-(clicker_geometry()*direction)./c.*fs.+offsets.*fs;delays=t.-t[3]
    time=collect(0:1000).-300.0
    pulse=exp.(-(time./3.5).^2).*(cos.(2pi*0.16.*time).+0.45cos.(2pi*0.32.*time))
    spectrum=rfft(pulse);bins=collect(0:length(spectrum)-1)
    snippet=hcat([irfft(spectrum.*cis.(-2pi.*bins.*delay./1001),1001) for delay in delays]...)
    return (snippet=snippet,direction=direction,delays=delays)
end

@testset "Method-matched raw-delay inference and single later correction" begin
    wave=synthetic_wave([0.24,-0.18]);fs=192000.0
    q=_pf_initial(1280,960,"centered_pinhole");q[7:8]=[0.02,-0.01]
    for method in ("gcc","srp")
        calibration=_pf_calibration(q,1280,960,"centered_pinhole";method=method)
        calibration["estimator_spec"]=_mr_estimator_spec(method;fs=fs)
        delays=estimate_method_tdoas(calibration,wave.snippet,[301];fs=fs,window=-300:700)
        expected=method=="gcc" ? first(prep_gcc(wave.snippet;fs=fs,c_min=1400.0,reference=3)) : review_srp_phat(wave.snippet;fs=fs).tdoas_samples
        @test size(delays)==(1,3)
        @test vec(delays)==expected
        @test delays[1,3]==0.0
        changed=deepcopy(calibration);changed["sound_speed_m_s"]=1400.0
        changed["relative_delay_offsets_s"]=[8e-6,3e-6,0.0]
        changed["rotation_array_to_camera"]=_pf_rows(_pf_rotation([0.01,0.02,0.03]))
        @test estimate_method_tdoas(changed,wave.snippet,[301];fs=fs)==delays
        # Applying calibrated corrections is explicitly a separate operation.
        original_pixel=project_tdoas(calibration,vec(delays);fs=fs)
        changed_pixel=project_tdoas(changed,vec(delays);fs=fs)
        @test original_pixel.valid && changed_pixel.valid
        @test abs(original_pixel.x-changed_pixel.x)>1.0
        @test_throws ErrorException estimate_method_tdoas(calibration,wave.snippet,[301];fs=48000)
        @test_throws ErrorException estimate_method_tdoas(calibration,wave.snippet,[301];fs=fs,window=-50:50)
    end
    gcc=Dict("method"=>"gcc")
    clipped=estimate_method_tdoas(gcc,wave.snippet,[1];fs=fs,window=-300:700)
    @test vec(clipped)==first(prep_gcc(wave.snippet[1:701,:];fs=fs,c_min=1400.0,reference=3))
    @test_throws ArgumentError estimate_method_tdoas(gcc,wave.snippet,[301];fs=fs,ref_channel=1)
    @test_throws DimensionMismatch estimate_method_tdoas(gcc,ones(1001,5),[301];fs=fs)
    @test_throws ArgumentError estimate_method_tdoas(gcc,wave.snippet,[0];fs=fs)
    @test_throws ArgumentError estimate_method_tdoas(Dict("method"=>"minmax2"),wave.snippet,[301];fs=fs)
    @test size(estimate_method_tdoas(gcc,wave.snippet,Int[];fs=fs))==(0,3)
end

function write_fixture(directory;heldout_shift=0.0,with_plan=true)
    n=39;fs=192000.0;width,height=1280,960;rng=MersenneTwister(91304)
    q=_pf_initial(width,height,"centered_pinhole");q[1]=0.0;q[2:4]=[0.02,-0.03,0.01]./pi;q[7:8]=[0.02,-0.01]
    truth=_pf_calibration(q,width,height,"centered_pinhole";device="synthetic")
    blocks=vcat(fill(33,3),fill(35,14),fill(40,3),fill(41,3),fill(34,16))
    within=Dict(block=>0 for block in unique(blocks))
    times=Float64[]
    for block in blocks;push!(times,30block+3+1.2within[block]);within[block]+=1;end
    events=DataFrame(event_id=["synthetic_$i" for i in 1:n],device=fill("synthetic",n),
        audio_file=fill("synthetic_only.wav",n),sample_index=[1000i for i in 1:n],audio_time_s=[1000i/fs for i in 1:n],
        recording_time_s=times,fs=fill(fs,n),frame_id=["frame_$i" for i in 1:n],width=fill(width,n),height=fill(height,n),
        sync_status=fill("mapped",n),window_status=fill("complete",n),
        window_source_segments=fill(JSON.json([Dict("audio_file"=>"synthetic_only.wav","first_sample"=>700,"samples"=>1001)]),n),
        snippet_file=fill("event_snippets.f64",n),snippet_byte_offset=collect(0:n-1).*(1001*3*8),
        snippet_samples=fill(1001,n),snippet_channels=fill(3,n),snippet_peak_index=fill(300,n),snippet_sha256=fill("",n))
    for method in ("minmax2","minmax","gcc"),j in 1:3;events[!,method*"_d"*string(j)]=zeros(n);end
    annotations=DataFrame(frame_id=events.frame_id,x=zeros(n),y=zeros(n),visibility=fill("visible",n),medium=fill("underwater",n),
        uncertainty_px=fill(1.0,n),note=fill("independent synthetic ray; not an estimate or real annotation",n),annotator=fill("synthetic",n))
    cache=joinpath(directory,"event_snippets.f64")
    open(cache,"w") do io
        for i in 1:n
            wave=synthetic_wave([0.7(rand(rng)-0.5),0.5(rand(rng)-0.5)])
            pixel=_pf_project_direction(truth,wave.direction)
            annotations[i,:x]=pixel.x+(blocks[i]==34 ? heldout_shift : 0.0);annotations[i,:y]=pixel.y
            mm2=prep_minmax(wave.snippet);mm=prep_minmax(wave.snippet;corrected=false);gcc,_=prep_gcc(wave.snippet;fs=fs)
            for (method,delays) in (("minmax2",mm2),("minmax",mm),("gcc",gcc)),j in 1:3
                events[i,method*"_d"*string(j)]=delays[j]
            end
            bytes=_prep_snippet_bytes(wave.snippet);write(io,bytes);events[i,:snippet_sha256]=bytes2hex(SHA.sha256(bytes))
        end
    end
    CSV.write(joinpath(directory,"events.csv"),events)
    CSV.write(joinpath(directory,"frames.csv"),select(events,:frame_id,:width,:height))
    CSV.write(joinpath(directory,"annotations.csv"),annotations)
    write_json(joinpath(directory,"annotations.complete.json"),Dict("frames_sha256"=>sha256_file(joinpath(directory,"frames.csv")),
        "annotations_sha256"=>sha256_file(joinpath(directory,"annotations.csv")),"n_frames"=>n))
    hashes=Dict(name*"_sha256"=>sha256_file(joinpath(directory,file)) for (name,file) in (("events","events.csv"),("frames","frames.csv"),("annotations","annotations.csv"),("annotations_complete","annotations.complete.json")))
    frozen=Dict("schema_version"=>"clicker-heldout-evaluation/1.0","hashes"=>hashes,"evaluation_calibration"=>truth,
        "split"=>Dict("train_event_ids"=>collect(events.event_id[1:23]),"heldout_event_ids"=>collect(events.event_id[24:39])))
    write_json(joinpath(directory,"heldout_evaluation.json"),frozen)
    write_json(joinpath(directory,"calibration.json"),truth)
    write_json(joinpath(directory,"fit_report.json"),Dict("calibration_sha256"=>sha256_file(joinpath(directory,"calibration.json")),
        "heldout_evaluation_sha256"=>sha256_file(joinpath(directory,"heldout_evaluation.json"))))
    CSV.write(joinpath(directory,"predictions.csv"),DataFrame(event_id=events.event_id,baseline_sentinel=fill("must remain unchanged",n)))
    outputs=[Dict("path"=>joinpath(directory,name),"sha256"=>sha256_file(joinpath(directory,name))) for name in ("events.csv","event_snippets.f64")]
    write_json(joinpath(directory,"prepare_manifest.json"),Dict("status"=>"complete","outputs"=>outputs,
        "snippet_cache"=>Dict("path"=>cache,"sha256"=>sha256_file(cache),"bytes"=>filesize(cache),"dtype"=>"Float64","endianness"=>"little","layout"=>"column-major")))
    review_estimators(directory;resume=false)
    canonical=["calibration.json","fit_report.json","heldout_evaluation.json","predictions.csv","estimator_review.json","estimator_review.csv",
        "events.csv","frames.csv","annotations.csv","annotations.complete.json","event_snippets.f64"]
    plan=Dict("version"=>"clicker-method-refinement-v1","canonical_input_sha256"=>Dict(name=>sha256_file(joinpath(directory,name)) for name in canonical),
        "methods"=>["gcc","srp"],"model"=>"centered_pinhole","parameters"=>8,"sound_speed_bounds_m_s"=>[1400,1600],
        "nstarts"=>6,"maxiters"=>2500,"robust_delta_px"=>12.0,"train_blocks"=>[33,35,40,41],
        "selection"=>"leave-one-training-block-out; equal-block mean invalid-penalized RMSE; deterministic GCC tie priority",
        "evaluation"=>"reused descriptive evaluation after frozen training selection")
    if with_plan
        mkpath(joinpath(directory,"refinement_v1"));write_json(joinpath(directory,"refinement_v1","plan.json"),plan)
    end
    return (canonical=canonical,plan=plan,truth=truth)
end

@testset "Train-only selection, reused evaluation, immutable baseline and resume" begin
    mktempdir() do parent
        CALIBRATION_ERROR_LOG[]=joinpath(parent,"expected_test_errors.jsonl")
        first_dir=joinpath(parent,"first");second_dir=joinpath(parent,"heldout_shifted");mkpath(first_dir);mkpath(second_dir)
        first_fixture=write_fixture(first_dir)
        second_fixture=write_fixture(second_dir;heldout_shift=70.0,with_plan=false)
        @test !isfile(joinpath(second_dir,"refinement_v1","plan.json"))
        protected=Dict(name=>sha256_file(joinpath(first_dir,name)) for name in first_fixture.canonical)
        plan_sha=sha256_file(joinpath(first_dir,"refinement_v1","plan.json"))
        # Compare the actual arrays reaching selection before native fitting.
        a=_mr_load(first_dir,_mr_bindings(first_dir));b=_mr_load(second_dir,_mr_bindings(second_dir))
        generated=read_json(joinpath(second_dir,"refinement_v1","plan.json"))
        @test generated["train_blocks"]==[33,35,40,41]
        @test generated["canonical_input_sha256"]==second_fixture.plan["canonical_input_sha256"]
        @test length(a.train)==23 && length(a.heldout)==16
        @test sort(unique(a.input.data.block[a.train]))==[33,35,40,41]
        for method in ("gcc","srp")
            @test _pf_subset(a.methods[method],a.train).x==_pf_subset(b.methods[method],b.train).x
            @test _pf_subset(a.methods[method],a.train).tau==_pf_subset(b.methods[method],b.train).tau
        end
        report=refine_calibration(first_dir;resume=false)
        report2=refine_calibration(second_dir;resume=false)
        dest=joinpath(first_dir,"refinement_v1");dest2=joinpath(second_dir,"refinement_v1")
        selection=read_json(joinpath(dest,"selection.json"));selection2=read_json(joinpath(dest2,"selection.json"))
        evaluation=read_json(joinpath(dest,"evaluation.json"));calibration=read_json(joinpath(dest,"calibration.json"))
        @test report["status"]=="complete"
        @test all(isfile(joinpath(dest,name)) for name in ("selection.json","evaluation.json","calibration.json","predictions.csv","report.json"))
        @test !selection["uses_heldout_labels"]
        @test selection["selected_method"]==selection2["selected_method"]
        @test selection["full_training_q"]==selection2["full_training_q"]
        @test selection["candidates"]==selection2["candidates"]
        @test report["reused_evaluation_metrics"]["rmse_px"]!=report2["reused_evaluation_metrics"]["rmse_px"]
        @test report["n_training"]==23 && report["n_reused_evaluation"]==16 && report["n_deployment"]==39
        @test evaluation["selection_sha256"]==sha256_file(joinpath(dest,"selection.json"))
        @test calibration["evaluation_sha256"]==sha256_file(joinpath(dest,"evaluation.json"))
        @test calibration["method"]==selection["selected_method"]==report["selected_method"]
        @test calibration["estimator_spec"]["method"]==calibration["method"]
        @test !report["automatic_defaults_installed"]
        for method in ("gcc","srp")
            folds=selection["candidates"][method]["folds"]
            @test length(folds)==4
            @test only(fold["n_training"] for fold in folds if fold["validation_block"]==35)==9
            @test all(!(fold["validation_block"] in fold["training_blocks"]) for fold in folds)
            @test all(isempty(intersect(Set(fold["training_event_rows"]),Set(fold["validation_event_rows"]))) for fold in folds)
            @test selection["candidates"][method]["equal_block_mean_penalized_rmse_px"]≈mean(fold["penalized_validation_rmse_px"] for fold in folds)
        end
        scores=[selection["candidates"][method]["equal_block_mean_penalized_rmse_px"] for method in ("gcc","srp")]
        @test selection["selected_method"]==(scores[2]<scores[1] ? "srp" : "gcc")
        table=CSV.read(joinpath(dest,"predictions.csv"),DataFrame)
        @test nrow(table)==39
        @test count(==("heldout"),table.split)==16
        @test all(==(selection["selected_method"]),table.method)
        @test isequal(table.predicted_x,table.final_x) && isequal(table.predicted_y,table.final_y)
        @test all(in(names(table)),["baseline_x","baseline_y","baseline_status","baseline_error_px","evaluation_status","evaluation_error_px","final_status","final_error_px"])
        @test all(sha256_file(joinpath(first_dir,name))==digest for (name,digest) in protected)
        @test sha256_file(joinpath(dest,"plan.json"))==plan_sha
        before=sha256_file(joinpath(dest,"report.json"))
        @test refine_calibration(first_dir)["outputs_sha256"]==report["outputs_sha256"]
        @test sha256_file(joinpath(dest,"report.json"))==before
        # A modified canonical input is rejected, including on resume, while
        # baseline/selection/evaluation outputs are never silently replaced.
        baseline=read(joinpath(first_dir,"predictions.csv"))
        open(joinpath(first_dir,"predictions.csv"),"a") do io;write(io,"\n");end
        @test_throws ErrorException refine_calibration(first_dir)
        @test sha256_file(joinpath(dest,"report.json"))==before
        write(joinpath(first_dir,"predictions.csv"),baseline)
        @test all(sha256_file(joinpath(first_dir,name))==digest for (name,digest) in protected)
        @info "Synthetic method refinement" selected=selection["selected_method"] block_scores=scores reused_rmse=report["reused_evaluation_metrics"]["rmse_px"] all_eligible_rmse=report["deployment_metrics_descriptive"]["rmse_px"]
    end
end
end # module
