# Bounded native synthetic tests; no recordings, real frames, or dataset jobs.
module SRPReviewTests
using Test, CSV, DataFrames, JSON, SHA, Dates, Statistics, LinearAlgebra, Random, Optim, DSP, FFTW, WAV, Peaks
const SOURCE_DIR=joinpath(@__DIR__,"..","..","src","clicker_calibration")
include(joinpath(SOURCE_DIR,"support.jl"))
include(joinpath(SOURCE_DIR,"physical_fit.jl"))
include(joinpath(SOURCE_DIR,"prepare.jl"))
include(joinpath(SOURCE_DIR,"srp_review.jl"))

function synthetic_ray(ux=0.24,uy=-0.18;fs=192000.0,n=1001)
    u=[ux,uy,sqrt(1-ux^2-uy^2)]
    arrival=-(clicker_geometry()*u)./1540.0.*fs;delays=arrival.-arrival[3]
    time=collect(0:n-1).-300.0
    pulse=exp.(-(time./3.5).^2).*(cos.(2pi*0.16.*time).+0.45cos.(2pi*0.32.*time))
    spectrum=rfft(pulse);bins=collect(0:length(spectrum)-1)
    snippet=hcat([irfft(spectrum.*cis.(-2pi.*bins.*delay./n),n) for delay in delays]...)
    return (snippet=snippet,u=u,delays=delays,fs=fs)
end

module FrozenOriginalSRP
using DSP,FFTW,LinearAlgebra,Statistics
default_c=1540.0;default_window=-300:700;ref_channel=3;impulsive_band_pass=[1000.0,Inf]
distance_pairwise(rx)=[norm(rx[:,j]-rx[:,k]) for j in axes(rx,2),k in axes(rx,2)]
get_relevant_channels(rx)=collect(axes(rx,2))
end

@testset "Fractional physical SRP sign, original native formulation and invalidity" begin
    ray=synthetic_ray();plan=_srv_plan(ray.fs)
    estimate=review_srp_phat(ray.snippet;fs=ray.fs,plan=plan)
    @test estimate.valid
    @test estimate.status=="ok"
    @test ray.delays[1]<0 && ray.delays[2]>0
    @test estimate.tdoas_samples[1]<0 && estimate.tdoas_samples[2]>0
    @test maximum(abs.(estimate.tdoas_samples-ray.delays))<0.12
    @test norm(estimate.uv-ray.u[1:2])<0.006
    @test estimate.tdoas_samples[3]==0.0
    @test sum(abs2,estimate.uv)<=1.0
    @test estimate.score>0.75
    @test estimate.score≈mean(estimate.pair_scores) atol=1e-13
    @test any(abs.(ray.delays-round.(ray.delays)).>0.1)
    # Read only the explicit frozen localization fragment; never include the
    # original script coordinator, dsp.jl, or any import-time legacy workflow.
    root=normpath(joinpath(@__DIR__,"..",".."))
    original_path=joinpath(root,"calibration_runs","clicker_20260209","20261030_v1","source_snapshot","src","localization.jl")
    @test sha256_file(original_path)==_srv_spec()["original_source_sha256"]
    original=read(original_path,String)
    start=Base.first(findfirst("using FFTW, LinearAlgebra, Statistics",original))
    stop=Base.first(findfirst("# unclamped direction cosines",original))-1
    Base.include_string(FrozenOriginalSRP,original[start:stop],original_path*"#gcc_srp_only")
    reference,quality=Base.invokelatest(FrozenOriginalSRP.get_tdoa_srp,ray.snippet,[301];window=-300:700,
        ref_channel=3,fs=ray.fs,rx_vect=permutedims(clicker_geometry()),c=1540.0,
        band=[1000.0,0.98ray.fs/2],interp=8,grid_deg=1.0,refine=10,taper=0.0,return_quality=true)
    @test estimate.tdoas_samples≈vec(reference) atol=1e-10
    @test estimate.uv≈vec(quality.uv) atol=1e-12
    @test estimate.score≈only(quality.srp) atol=1e-12
    @info "Native fractional physical SRP" true_delays=ray.delays estimated_delays=estimate.tdoas_samples max_delay_error_samples=maximum(abs.(estimate.tdoas_samples-ray.delays)) score=estimate.score
    @test !review_srp_phat(zeros(1001,3);fs=ray.fs,plan=plan).valid
    @test review_srp_phat(zeros(1001,3);fs=ray.fs,plan=plan).status=="silent_channel"
    missing=copy(ray.snippet);missing[:,2].=0
    @test !review_srp_phat(missing;fs=ray.fs,plan=plan).valid
    nonfinite=copy(ray.snippet);nonfinite[5,1]=NaN
    @test review_srp_phat(nonfinite;fs=ray.fs).status=="nonfinite_snippet"
    @test review_srp_phat(ones(3,3);fs=ray.fs).status=="too_short"
    @test_throws DimensionMismatch review_srp_phat(ones(1001,2);fs=ray.fs)
    @test_throws ArgumentError review_srp_phat(ray.snippet;fs=48000,plan=plan)
    @test_throws ArgumentError _srv_plan(1000)
end

function write_review_fixture(directory)
    n=16;fs=192000.0;width,height=1280,960
    q=_pf_initial(width,height,"centered_pinhole");q[1]=0.4
    evaluation=_pf_calibration(q,width,height,"centered_pinhole";device="synthetic")
    final=deepcopy(evaluation)
    final["rotation_array_to_camera"]=_pf_rows(_pf_rotation([0.0,0.25,0.0]))
    events=DataFrame(event_id=["event_$i" for i in 1:n],device=fill("synthetic",n),
        audio_file=fill("synthetic_source.wav",n),sample_index=[1000i for i in 1:n],
        audio_time_s=[1000i/fs for i in 1:n],recording_time_s=[(i<=8 ? 0 : 120)+3+2mod(i-1,8) for i in 1:n],
        fs=fill(fs,n),frame_id=["frame_$i" for i in 1:n],width=fill(width,n),height=fill(height,n),
        sync_status=fill("mapped",n),window_status=fill("complete",n),
        window_source_segments=fill(JSON.json([Dict("audio_file"=>"synthetic_source.wav","first_sample"=>700,"samples"=>1001)]),n),
        snippet_file=fill("event_snippets.f64",n),snippet_byte_offset=collect(0:n-1).*(1001*3*8),
        snippet_samples=fill(1001,n),snippet_channels=fill(3,n),snippet_peak_index=fill(300,n),snippet_sha256=fill("",n))
    events[1,:window_status]="cross_file_context"
    events[1,:window_source_segments]=JSON.json([Dict("audio_file"=>"part1.wav","first_sample"=>9900,"samples"=>100),Dict("audio_file"=>"part2.wav","first_sample"=>0,"samples"=>901)])
    for method in ("minmax2","minmax","gcc"),j in 1:3;events[!,method*"_d"*string(j)]=zeros(n);end
    annotations=DataFrame(frame_id=events.frame_id,x=zeros(n),y=zeros(n),visibility=fill("visible",n),medium=fill("underwater",n),
        uncertainty_px=fill(1.0,n),note=fill("independent synthetic physical ray",n),annotator=fill("synthetic",n))
    cache_path=joinpath(directory,"event_snippets.f64")
    open(cache_path,"w") do io
        for i in 1:n
            ray=synthetic_ray(0.1+0.012i,-0.2+0.009i)
            pixel=_pf_project_direction(evaluation,ray.u)
            annotations[i,:x]=pixel.x;annotations[i,:y]=pixel.y
            for method in ("minmax2","minmax","gcc"),j in 1:3
                events[i,method*"_d"*string(j)]=ray.delays[j]+(j==3 || method=="minmax2" ? 0.0 : method=="minmax" ? 0.5 : 0.1)
            end
            snippet=i==16 ? zeros(1001,3) : ray.snippet
            bytes=_prep_snippet_bytes(snippet);write(io,bytes);events[i,:snippet_sha256]=bytes2hex(SHA.sha256(bytes))
        end
    end
    events[14,:gcc_d1]=NaN;events[15,:minmax_d1]=960.0
    CSV.write(joinpath(directory,"events.csv"),events)
    CSV.write(joinpath(directory,"frames.csv"),select(events,:frame_id,:width,:height))
    CSV.write(joinpath(directory,"annotations.csv"),annotations)
    write_json(joinpath(directory,"annotations.complete.json"),Dict("frames_sha256"=>sha256_file(joinpath(directory,"frames.csv")),"annotations_sha256"=>sha256_file(joinpath(directory,"annotations.csv")),"n_frames"=>n))
    hashes=Dict(name*"_sha256"=>sha256_file(joinpath(directory,file)) for (name,file) in (("events","events.csv"),("frames","frames.csv"),("annotations","annotations.csv"),("annotations_complete","annotations.complete.json")))
    frozen=Dict("schema_version"=>"clicker-heldout-evaluation/1.0","hashes"=>hashes,"evaluation_calibration"=>evaluation,
        "split"=>Dict("train_event_ids"=>collect(events.event_id[1:8]),"heldout_event_ids"=>collect(events.event_id[9:16])))
    write_json(joinpath(directory,"heldout_evaluation.json"),frozen)
    write_json(joinpath(directory,"calibration.json"),final)
    write_json(joinpath(directory,"fit_report.json"),Dict("calibration_sha256"=>sha256_file(joinpath(directory,"calibration.json")),"heldout_evaluation_sha256"=>sha256_file(joinpath(directory,"heldout_evaluation.json"))))
    outputs=[Dict("path"=>joinpath(directory,name),"sha256"=>sha256_file(joinpath(directory,name))) for name in ("events.csv","event_snippets.f64")]
    write_json(joinpath(directory,"prepare_manifest.json"),Dict("status"=>"complete","outputs"=>outputs,
        "snippet_cache"=>Dict("path"=>cache_path,"sha256"=>sha256_file(cache_path),"bytes"=>filesize(cache_path),"dtype"=>"Float64","endianness"=>"little","layout"=>"column-major")))
    return (events=events,evaluation=evaluation,final=final)
end

@testset "Frozen evaluation, shared cache/population, paired coverage and no fit writes" begin
    mktempdir() do directory
        CALIBRATION_ERROR_LOG[]=joinpath(directory,"expected_test_errors.jsonl")
        fixture=write_review_fixture(directory)
        protected=["calibration.json","heldout_evaluation.json","fit_report.json","events.csv","annotations.csv","frames.csv","annotations.complete.json","prepare_manifest.json","event_snippets.f64"]
        before=Dict(name=>sha256_file(joinpath(directory,name)) for name in protected)
        report=review_estimators(directory;resume=false)
        @test report["status"]=="complete"
        @test report["n_population"]==16
        @test report["hashes"]["calibration.json"]==before["calibration.json"]
        @test report["hashes"]["heldout_evaluation.json"]==before["heldout_evaluation.json"]
        @test all(sha256_file(joinpath(directory,name))==before[name] for name in protected)
        @test report["specification"]["sound_speed_m_s"]==1540.0
        @test report["summaries"]["heldout"]["methods"]["minmax2"]["coverage"]==1.0
        @test report["summaries"]["heldout"]["methods"]["minmax2"]["rmse_px"]<1e-8
        @test report["summaries"]["heldout"]["methods"]["gcc"]["coverage"]==7/8
        @test report["summaries"]["heldout"]["methods"]["minmax"]["coverage"]==7/8
        @test report["summaries"]["heldout"]["methods"]["srp"]["coverage"]==7/8
        @test report["summaries"]["heldout"]["all_four_common_coverage"]==5/8
        @test length(report["summaries"]["heldout"]["all_four_common_event_ids"])==5
        @test length(report["summaries"]["heldout"]["paired_common_events"])==6
        @test report["summaries"]["heldout"]["methods"]["srp"]["srp_score_summary"]["n"]==7
        table=CSV.read(joinpath(directory,"estimator_review.csv"),DataFrame)
        @test nrow(table)==64
        @test all(count(==(id),table.event_id)==4 for id in unique(table.event_id))
        @test table[(table.event_id.=="event_16").&(table.method.=="srp"),:estimator_status]==["silent_channel"]
        @test table[(table.event_id.=="event_15").&(table.method.=="minmax"),:projection_status]==["impossible_tdoa"]
        @test all(==("cross_file_context"),table[table.event_id.=="event_1",:window_status])
        @test all(==(fixture.events[1,:window_source_segments]),table[table.event_id.=="event_1",:window_source_segments])
        @test Set(table[table.split.=="heldout",:event_id])==Set(fixture.events.event_id[9:16])
        review_sha=sha256_file(joinpath(directory,"estimator_review.json"))
        @test review_estimators(directory)["hashes"]==report["hashes"]
        @test sha256_file(joinpath(directory,"estimator_review.json"))==review_sha
        # Frozen fit tampering must fail before any new review output is written.
        calibration=read(joinpath(directory,"calibration.json"))
        open(joinpath(directory,"calibration.json"),"a") do io;write(io,"\n");end
        @test_throws ErrorException review_estimators(directory)
        @test sha256_file(joinpath(directory,"estimator_review.json"))==review_sha
        write(joinpath(directory,"calibration.json"),calibration)
        # Whole-cache integrity remains mandatory even on a resume path.
        cache=read(joinpath(directory,"event_snippets.f64"));bad=copy(cache);bad[10]⊻=0x01
        write(joinpath(directory,"event_snippets.f64"),bad)
        @test_throws ErrorException review_estimators(directory)
        write(joinpath(directory,"event_snippets.f64"),cache)
        annotation=read(joinpath(directory,"annotations.csv"))
        open(joinpath(directory,"annotations.csv"),"a") do io;write(io,"\n");end
        @test_throws ErrorException review_estimators(directory)
        write(joinpath(directory,"annotations.csv"),annotation)
        @test all(sha256_file(joinpath(directory,name))==before[name] for name in protected)
        @test isfile(joinpath(directory,"expected_test_errors.jsonl"))
    end
end
end # module
