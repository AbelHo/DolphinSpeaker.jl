# Standalone native-Julia bounded tests. No real recording or annotation input.
# JULIA_DEPOT_PATH=/home/spin/Documents/GitHub/DolphinSpeaker.jl/calibration_runs/clicker_20260209/20261030_v1/julia_depot:/home/spin/.julia \
# /home/spin/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia \
# --startup-file=no --compiled-modules=no --project=/home/spin/Documents/GitHub/DolphinSpeaker.jl test_physical_fit.jl
module PhysicalFitTests
using Test, CSV, DataFrames, JSON, SHA, Dates, Statistics, LinearAlgebra, Random, Optim
include(joinpath(@__DIR__,"..","..","src","clicker_calibration","support.jl"))
include(joinpath(@__DIR__,"..","..","src","clicker_calibration","physical_fit.jl"))

function synthetic_data(;n=120,outliers=false,degenerate=false)
    rng=MersenneTwister(7043)
    width,height=1280,960
    q=_pf_initial(width,height,"centered_pinhole")
    q[1]=0.27;q[2:4]=[0.035,-0.045,0.022]./pi;q[7:8]=[0.025,-0.018]
    q[5]+=0.06;q[6]-=0.05
    truth=_pf_calibration(q,width,height,"centered_pinhole";device="synthetic")
    geometry=clicker_geometry();baseline=geometry[1:2,1:2].-permutedims(geometry[3,1:2])
    directions=zeros(n,3);tau=zeros(n,2);xs=zeros(n);ys=zeros(n)
    for i in 1:n
        ux,uy=degenerate ? (0.1,0.05) : (0.95*(rand(rng)-0.5),0.85*(rand(rng)-0.5))
        u=[ux,uy,sqrt(1-ux^2-uy^2)];directions[i,:]=u
        tau[i,:]=-(baseline*u[1:2])./truth["sound_speed_m_s"].+truth["relative_delay_offsets_s"][1:2]
        pixel=_pf_project_direction(truth,u)
        @assert pixel.valid
        xs[i]=pixel.x;ys[i]=pixel.y
        if outliers && i%17 == 0
            xs[i]+=150;ys[i]-=120
        end
    end
    return (data=(tau=tau,x=xs,y=ys,weight=ones(n),block=floor.(Int,(0:n-1)./8),
        event_row=collect(1:n),frame=["frame_$i" for i in 1:n],inverse_baseline=inv(baseline)),
        q=q,truth=truth,width=width,height=height,directions=directions)
end

function write_fixture(directory;heldout_shift=0.0)
    fixture=synthetic_data(n=120)
    # Eight distinct annotations/block; every fifth held out. Additional exact
    # boundary observations exercise purging independently of acoustic geometry.
    n=120;fs=192000.0;times=[30floor((i-1)/8)+3+3*((i-1)%8) for i in 1:n]
    times[1]=0.5;times[2]=28.0
    events=DataFrame(event_id=["event_$i" for i in 1:n],device=fill("synthetic",n),
        frame_id=fixture.data.frame,width=fill(1280,n),height=fill(960,n),fs=fill(fs,n),
        recording_time_s=times,sync_status=fill("mapped",n))
    for method in ("minmax2","minmax","gcc")
        events[!,method*"_d1"]=fixture.data.tau[:,1].*fs
        events[!,method*"_d2"]=fixture.data.tau[:,2].*fs
        events[!,method*"_d3"]=zeros(n)
    end
    events[3,:sync_status]="ambiguous"
    annotations=DataFrame(frame_id=fixture.data.frame,x=fixture.data.x,y=fixture.data.y,
        visibility=fill("visible",n),medium=fill("underwater",n),uncertainty_px=fill(2.0,n),note=fill("synthetic independent ground truth",n),annotator=fill("synthetic",n))
    annotations[4,:visibility]="uncertain"
    annotations[5,:medium]="air";annotations[6,:medium]="mixed";annotations[7,:medium]="unknown"
    for i in 1:n
        floor(Int,times[i]/30)%5 == 4 && (annotations[i,:x]+=heldout_shift)
    end
    CSV.write(joinpath(directory,"events.csv"),events)
    # Acoustic fields intentionally absent; annotation generation has no model.
    CSV.write(joinpath(directory,"frames.csv"),select(events,:frame_id,:width,:height))
    CSV.write(joinpath(directory,"annotations.csv"),annotations)
    write_json(joinpath(directory,"annotations.complete.json"),Dict("schema_version"=>1,
        "frames_sha256"=>sha256_file(joinpath(directory,"frames.csv")),
        "annotations_sha256"=>sha256_file(joinpath(directory,"annotations.csv")),"n_frames"=>n,"annotator"=>"synthetic"))
    return fixture
end

@testset "D3 physical projection conventions and invalid geometry" begin
    f=synthetic_data(n=12)
    @test eltype(clicker_geometry()) == Float64
    @test all(isapprox(norm(clicker_geometry()[i,:]-clicker_geometry()[j,:]),0.18;atol=1e-15) for (i,j) in ((1,2),(1,3),(2,3)))
    @test det(_pf_matrix(f.truth["rotation_array_to_camera"])) ≈ 1 atol=1e-14
    rays=project_tdoas(f.truth,f.data.tau;units="seconds")
    @test all(r.valid for r in rays)
    @test [r.x for r in rays] ≈ f.data.x atol=1e-10
    @test [r.y for r in rays] ≈ f.data.y atol=1e-10
    @test rays[1].direction_array ≈ f.directions[1,:] atol=1e-13
    samples=hcat(f.data.tau.*192000.0,zeros(12))
    @test project_tdoas(f.truth,samples;fs=fill(192000.0,12))[1].x ≈ rays[1].x
    u=f.directions[1,:];el=asin(u[2]);az=atan(u[1],u[3])
    @test project_angles(f.truth,[az,el]).x ≈ rays[1].x atol=1e-10
    @test project_angles(f.truth,rad2deg.([az,el]);units="degrees").y ≈ rays[1].y atol=1e-10
    @test _pf_d3_baseline([0.0,0.0],1280,960).x == 640.0
    @test _pf_d3_baseline([0.0,0.0],1280,960).y == 480.0
    @test !_pf_d3_baseline([1.0,1.0],1280,960).valid
    @test length(project_angles(f.truth,[az el;0.0 0.0])) == 2
    @test !project_tdoas(f.truth,[1.0,1.0];units="seconds").valid
    @test project_tdoas(f.truth,[1.0,1.0];units="seconds").status == "impossible_tdoa"
    @test isnan(project_tdoas(f.truth,[1.0,1.0];units="seconds").x)
    @test project_tdoas(f.truth,[0.0,0.0,1.0];fs=192000.).status == "invalid_reference"
    @test project_tdoas(f.truth,[NaN,0.0];units="seconds").status == "nonfinite_tdoa"
    @test project_angles(f.truth,[pi,0.0]).status == "back_array_branch"
    @test project_angles(f.truth,[0.0,pi]).status == "invalid_elevation"
    @test_throws ArgumentError project_tdoas(f.truth,[0.0,0.0])
    @test_throws DimensionMismatch project_tdoas(f.truth,samples;fs=[192000.0])
    reversed=deepcopy(f.truth);reversed["rotation_array_to_camera"]=_pf_rows(_pf_rotation([0.0,pi,0.0]))
    @test project_angles(reversed,[0.0,0.0]).status == "behind_camera"
    broken=deepcopy(f.truth);broken["geometry_m"]=_pf_rows(zeros(3,3))
    @test_throws ArgumentError project_tdoas(broken,[0.0,0.0];units="seconds")
end

@testset "Native Optim physical recoverability and robust loss" begin
    f=synthetic_data(n=120)
    fit=_pf_fit(f.data,f.width,f.height,"centered_pinhole";nstarts=4,maxiters=3500)
    cal=_pf_calibration(fit.q,f.width,f.height,"centered_pinhole")
    metrics=_pf_metrics(cal,f.data)
    @info "Synthetic recoverability" rmse_px=metrics["rmse_px"] median_px=metrics["median_px"] p90_px=metrics["p90_px"] c=cal["sound_speed_m_s"]
    @test metrics["coverage"] == 1.0
    @test metrics["rmse_px"] < 0.05
    @test abs(cal["sound_speed_m_s"]-1527.0) < 2.0
    @test maximum(abs.(cal["relative_delay_offsets_s"].-f.truth["relative_delay_offsets_s"])) < 2e-7
    @test maximum(abs.(_pf_matrix(cal["rotation_array_to_camera"]).-_pf_matrix(f.truth["rotation_array_to_camera"]))) < 0.002
    diagnostics=_pf_diagnostics(fit.q,f.data,f.width,f.height,"centered_pinhole";profiles=false)
    @test diagnostics["jacobian_rank"] == 8
    @test diagnostics["cluster_blocks"] == 15
    profiled=_pf_diagnostics(fit.q,f.data,f.width,f.height,"centered_pinhole";profiles=true,maxiters=600)
    @test length(profiled["training_sound_speed_profiles"]) == 3
    @test [row["fixed_sound_speed_m_s"] for row in profiled["training_sound_speed_profiles"]] == [1400.0,1500.0,1600.0]
    # Robust objective must prefer the true broad-direction physical model even
    # with large sparse pixel errors. Tests a real scientific failure condition.
    dirty=synthetic_data(n=120,outliers=true)
    robust=_pf_fit(dirty.data,dirty.width,dirty.height,"centered_pinhole";nstarts=3,maxiters=2800)
    robust_cal=_pf_calibration(robust.q,dirty.width,dirty.height,"centered_pinhole")
    clean_metrics=_pf_metrics(robust_cal,f.data)
    @info "Sparse outlier fit clean-ray metrics" clean_metrics
    @test clean_metrics["median_px"] < 5.0
    @test clean_metrics["coverage"] == 1.0
end

@testset "Explicitly degenerate identifiability and coverage" begin
    f=synthetic_data(n=40,degenerate=true)
    diagnostics=_pf_diagnostics(f.q,f.data,f.width,f.height,"centered_pinhole";profiles=false)
    @test diagnostics["weak_identifiability"]
    @test diagnostics["jacobian_rank"] <= 2
    @test diagnostics["jacobian_condition_number"] > 1e5
    @test all(isnan,diagnostics["cluster_sandwich_se_normalized"])
    impossible=merge(f.data,(tau=fill(1.0,40,2),))
    metrics=_pf_metrics(f.truth,impossible)
    @test metrics["coverage"] == 0.0
    @test metrics["rmse_px"] === nothing
    @test metrics["n_eligible"] == 40
    _,selection=_pf_select_model(f.data,f.width,f.height;allow_extended=true)
    @test selection["selected"] == "centered_pinhole"
    @test !selection["uses_heldout"]
end

@testset "Gross impossible TDOA cannot pull clean physical fit to c bound" begin
    f=synthetic_data(n=120);q=copy(f.q);q[1]=0.0;q[7:8].=0.0
    truth=_pf_calibration(q,f.width,f.height,"centered_pinhole")
    geometry=clicker_geometry();B=geometry[1:2,1:2].-permutedims(geometry[3,1:2])
    tau=-(f.directions[:,1:2]*B')./1500.0
    pixels=[_pf_project_direction(truth,collect(row)) for row in eachrow(f.directions)]
    clean=merge(f.data,(tau=tau,x=[p.x for p in pixels],y=[p.y for p in pixels]))
    contaminated=(tau=vcat(tau,[0.005 0.0]),x=vcat(clean.x,f.width/2),y=vcat(clean.y,f.height/2),
        weight=ones(121),block=vcat(clean.block,15),event_row=collect(1:121),
        frame=vcat(clean.frame,"gross_tdoa_960samples"),inverse_baseline=clean.inverse_baseline)
    @test _pf_definitely_impossible(0.005,0.0,clean.inverse_baseline)
    @test !any(i->_pf_definitely_impossible(tau[i,1],tau[i,2],clean.inverse_baseline),1:120)
    # Exactly the audit contamination: [960,0] samples at192kHz,120 clean rays.
    @test _pf_loss(q,contaminated,f.width,f.height,"centered_pinhole",12.0) == _pf_loss(q,clean,f.width,f.height,"centered_pinhole",12.0)
    bound=copy(q);bound[1]=-1.0
    @test _pf_loss(bound,contaminated,f.width,f.height,"centered_pinhole",12.0) > _pf_loss(q,contaminated,f.width,f.height,"centered_pinhole",12.0)
    fit=_pf_fit(contaminated,f.width,f.height,"centered_pinhole";nstarts=4,maxiters=3500)
    cal=_pf_calibration(fit.q,f.width,f.height,"centered_pinhole")
    clean_metrics=_pf_metrics(cal,clean);all_metrics=_pf_metrics(cal,contaminated)
    @info "Gross TDOA contamination recovery" sound_speed=cal["sound_speed_m_s"] clean_rmse_px=clean_metrics["rmse_px"] coverage=all_metrics["coverage"]
    @test abs(cal["sound_speed_m_s"]-1500.0) < 2.0
    @test abs(fit.q[1]) < 0.98
    @test clean_metrics["rmse_px"] < 0.05
    @test all_metrics["coverage"] == 120/121
    @test all_metrics["status_counts"]["impossible_tdoa"] == 1
    @test _pf_screen_summary(contaminated)["n_screened_from_objective"] == 1
    @test isnan(project_tdoas(cal,[960.0,0.0,0.0];fs=192000.).x)
end

@testset "CSV pipeline, hash invalidation, frozen holdout and no annotation leakage" begin
    mktempdir() do directory
        # Keep intentional test failures out of shared operational incidents.
        CALIBRATION_ERROR_LOG[]=joinpath(directory,"expected_test_errors.jsonl")
        write_fixture(directory)
        input=_pf_load(directory,"minmax2")
        @test input.split[1:4] == ["purged","purged","ineligible","ineligible"]
        @test input.reason[3] == "sync_ambiguous"
        @test input.reason[5:7] == ["medium_air","medium_mixed","medium_unknown"]
        @test all(input.split[input.data.event_row[i]] == (input.data.block[i]%5 == 4 ? "heldout" : "train") for i in eachindex(input.data.x))
        @test isempty(intersect(input.data.frame[findall(i->input.split[input.data.event_row[i]]=="train",eachindex(input.data.x))],input.data.frame[findall(i->input.split[input.data.event_row[i]]=="heldout",eachindex(input.data.x))]))
        annotation_sha=sha256_file(joinpath(directory,"annotations.csv"))
        frames_sha=sha256_file(joinpath(directory,"frames.csv"))
        report=fit_device(directory;nstarts=2,maxiters=1800,allow_extended=false,diagnostics_profiles=false)
        @test all(isfile(joinpath(directory,name)) for name in ("calibration.json","predictions.csv","fit_report.json","heldout_evaluation.json"))
        @test sha256_file(joinpath(directory,"annotations.csv")) == annotation_sha
        @test sha256_file(joinpath(directory,"frames.csv")) == frames_sha
        @test report["heldout_metrics"]["n_eligible"] == 24
        @test report["heldout_metrics"]["coverage"] == 1.0
        @test report["heldout_metrics"]["rmse_px"] < 1.0
        @test !report["model_selection"]["uses_heldout"]
        frozen=read_json(joinpath(directory,"heldout_evaluation.json"))
        @test frozen["evaluation_status"] == "frozen_before_final_all_eligible_refit"
        @test Set(keys(report["post_first_fit_method_comparisons_training_only"]["methods"])) == Set(["minmax2","minmax","gcc"])
        predictions=CSV.read(joinpath(directory,"predictions.csv"),DataFrame)
        @test nrow(predictions) == 120
        @test all(in(names(predictions)),["final_x","evaluation_x","final_error_px","evaluation_error_px","eligibility_status"])
        @test predictions.predicted_x == predictions.final_x
        @test predictions.baseline_y != predictions.evaluation_y
        @test report["uncalibrated_d3_baseline_spec"]["sound_speed_m_s"] == 1540.0
        @test report["paired_baseline_comparison_heldout"]["common_coverage"] == 1.0
        frozen_sha=sha256_file(joinpath(directory,"heldout_evaluation.json"))
        @test fit_device(directory;nstarts=2,maxiters=1800,allow_extended=false,diagnostics_profiles=false)["heldout_evaluation_sha256"] == frozen_sha
        @test sha256_file(joinpath(directory,"heldout_evaluation.json")) == frozen_sha
        # Changing only heldout annotations cannot affect any training data,
        # model-selection input or fitted evaluation parameters.
        write_fixture(directory;heldout_shift=10.0)
        changed=_pf_load(directory,"minmax2")
        train=findall(i->input.split[input.data.event_row[i]]=="train",eachindex(input.data.x))
        @test _pf_subset(input.data,train).x == _pf_subset(changed.data,train).x
        @test _pf_subset(input.data,train).tau == _pf_subset(changed.data,train).tau
        report2=fit_device(directory;nstarts=2,maxiters=1800,allow_extended=false,diagnostics_profiles=false)
        frozen2=read_json(joinpath(directory,"heldout_evaluation.json"))
        @test frozen["evaluation_calibration"] == frozen2["evaluation_calibration"]
        @test report2["heldout_metrics"]["rmse_px"] > 9.0
        @test report2["hashes"]["annotations_sha256"] != annotation_sha
        # An annotation edit without a matching completion receipt fails closed.
        open(joinpath(directory,"annotations.csv"),"a") do io;write(io,"\n");end
        @test_throws ArgumentError fit_device(directory;nstarts=1,maxiters=10)
        @test isfile(joinpath(directory,"expected_test_errors.jsonl"))
    end
end
end # module
