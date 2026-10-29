using Test, Random
using DataFrames: nrow
# Run with Julia --project=/path/to/DolphinSpeaker.jl test/test_prepare.jl.
module PreparationTest
using CSV,DataFrames,JSON,SHA,Dates,Statistics,LinearAlgebra,Random,DSP,FFTW,WAV,Peaks
include(joinpath(@__DIR__,"../../src/clicker_calibration/support.jl"))
include(joinpath(@__DIR__,"../../src/clicker_calibration/prepare.jl"))
end
using .PreparationTest
const P=PreparationTest

try

if get(ENV,"CLICKER_PREPARE_TEST_SCOPE","all")=="all"
@testset "corrected polarity and reference-relative delays" begin
    x=zeros(128,3)
    for (j,d) in enumerate((4,-3,0))
        x[64+d,j]=-2; x[90,j]=0.5
    end
    @test P.prep_minmax(x)==[4,-3,0]
    @test P.prep_minmax(x;corrected=false)==[0,0,0]
end
@testset "all-pairs fractional GCC, signs, physical limits" begin
    n=1001; t=collect(1:n); shifts=(5.4,-3.2,0.0)
    x=hcat([exp.(-((t.-500 .-d)./8).^2).*cos.(2pi*0.13.*(t.-500 .-d)) for d in shifts]...)
    delay,closure=P.prep_gcc(x)
    @test maximum(abs.(delay.-collect(shifts)))<0.2
    @test delay[3]==0.0
    @test closure<0.15
    @test abs(P.prep_gcc_delay(x[:,1],x[:,3];fs=192000,maxlag=2.0)[1])<=2.0
    @test all(isfinite,P.prep_gcc(zeros(1001,3))[1])
    constrained=P.prep_reconcile_pairs([20.0,20.0,20.0],20.0)
    @test maximum(constrained)-minimum(constrained)<=20.0
end
@testset "isolated click detection and silence" begin
    Random.seed!(123);n=192000; x=1e-5randn(n,3)
    x[50000,:].+=1.0
    d=P.prep_detect(x)
    @test length(d.peaks)==1
    @test abs(only(d.peaks)-50000)<=2
    @test d.threshold>0
    @test isempty(P.prep_detect(zeros(20000,3)).peaks)
end
@testset "waveform offset, silence, repeated ambiguity" begin
    Random.seed!(71); y=randn(4096);x=0.01randn(16000); x[6001:10096].+=y
    hit=P.prep_wave_match(x,y;fs=2000)
    @test hit.time_s≈3.0
    @test hit.score>0.99
    @test hit.ambiguity_ratio>2
    @test P.prep_wave_match(zeros(16000),zeros(4096);fs=2000).score==0
    x[10001:14096].+=y
    @test P.prep_wave_match(x,y;fs=2000).ambiguity_ratio<1.15
end
@testset "robust drift fit, outlier rejection and source anchors" begin
    anchors=[Dict{String,Any}("recording_time_s"=>Float64(t),"video_global_time_s"=>-4.2+1.00012t+(t==300 ? 0.3 : 0),"accepted"=>true) for t in 0:100:1000]
    s=P.prep_fit_sync(anchors)
    @test s["status"]=="valid"
    @test s["offset_s"]≈-4.2 atol=1e-9
    @test s["drift_ppm"]≈120 atol=1e-6
    @test count(a->a["fit_inlier"],s["anchors"])==10
    @test P.prep_fit_sync(anchors[1:2])["status"]=="insufficient_waveform_anchors"
end
@testset "exact frame mapping and chapter edges" begin
    v=[Dict{String,Any}("path"=>"a.mp4","pts"=>[0.0,0.04,0.081],"global_start_s"=>0.0,"duration_s"=>0.12),
       Dict{String,Any}("path"=>"b.mp4","pts"=>[0.0,0.04],"global_start_s"=>0.12,"duration_s"=>0.08)]
    @test P.prep_match_frame(v,0.075).index==2
    @test P.prep_match_frame(v,0.12).index==0
    @test P.prep_match_frame(v,0.12).video["path"]=="b.mp4"
    @test isnothing(P.prep_match_frame(v,-0.01))
    @test isnothing(P.prep_match_frame(v,0.21))
    @test !P.prep_frame_ambiguous(v,0.01,0.001)
    @test P.prep_frame_ambiguous(v,0.02,0.001)
    @test P.prep_frame_ambiguous(v,0.12,0.001)
    @test P.prep_frame_ambiguous(v,0.199,0.003)
    sync=Dict("max_residual_s"=>0.002,"anchor_min_time_s"=>10.0,"anchor_max_time_s"=>110.0)
    @test P.prep_sync_uncertainty(sync,50.0)>=0.002+1/48000
    @test P.prep_sync_uncertainty(sync,200.0)>P.prep_sync_uncertainty(sync,50.0)
end
@testset "PCM boundary identity, no inferred duplicate removal" begin
    mktempdir() do dir
        a=reshape(collect(1.0:5000.0)./10000,1000,5);b=reverse(a;dims=1)
        pa=joinpath(dir,"a.WAV");pb=joinpath(dir,"b.WAV")
        P.WAV.wavwrite(a,pa;Fs=192000,nbits=24,compression=1);P.WAV.wavwrite(b,pb;Fs=192000,nbits=24,compression=1)
        infos=P.prep_wav_info.([pa,pb]);infos[1]["bext_time_reference"]=UInt64(100);infos[2]["bext_time_reference"]=UInt64(1099)
        boundary=only(P.prep_boundaries(infos))
        @test boundary["bext_minus_pcm_samples"]==-1
        @test boundary["samples_removed"]==0
        @test infos[2]["start_sample"]==1000
        segments=P._prep_window_sources(infos,1000,-300:700)
        @test segments[1]["first_sample"]==700
        @test segments[1]["samples"]==300
        @test segments[2]["first_sample"]==0
        @test segments[2]["samples"]==701
        @test size(P.prep_read_selected(pa,infos[1];channels=[3,5],samples=100))==(100,2)
        out=joinpath(dir,"result")
        m=P.prepare_device(dir,out;device="synthetic",detection_only=true,audio_limit_s=0.005)
        @test m["status"]=="bounded_or_unmapped"
        @test isfile(joinpath(out,"events.csv"))
        @test isfile(joinpath(out,"frames.csv"))
    end
end
@testset "authenticated resume and interrupted output" begin
    mktempdir() do dir
        path=joinpath(dir,"out.csv");write(path,"abc")
        manifest=joinpath(dir,"receipt.json");request=Dict("test"=>1)
        P.write_json(manifest,Dict("status"=>"complete","request"=>request,"outputs"=>P._prep_inventory([path])))
        @test !isnothing(P._prep_resume_valid(manifest,request))
        @test isnothing(P._prep_resume_valid(manifest,Dict("test"=>2)))
        write(path,"xyz")
        @test isnothing(P._prep_resume_valid(manifest,request))
        P.write_json(manifest,Dict("status"=>"incomplete","request"=>request,"outputs"=>P._prep_inventory([path])))
        @test isnothing(P._prep_resume_valid(manifest,request))
    end
end

@testset "exact source-frame extraction fixture" begin
    mktempdir() do dir
        video=joinpath(dir,"encoded.mp4")
        P.run_logged(`/usr/bin/ffmpeg -v error -nostdin -y -f lavfi -i testsrc2=size=320x180:rate=30000/1001 -t 3 -c:v libx264 -pix_fmt yuv420p $video`;operation="synthetic video fixture")
        frames=P._prep_table(P.PREPARE_FRAME_COLUMNS)
        for n in (0,7,60,89)
            push!(frames,("f$n",video,n,round(n*1001/30000;digits=6),round(n*1001/30000;digits=6),joinpath(dir,"f$n.png"),320,180))
        end
        P.prep_extract_frames(frames;batch_size=2,max_batch_span_s=0.5)
        @test all(isfile,frames.image_path)
        for row in eachrow(frames)
            expected=joinpath(dir,"expected$(row.frame_index).png")
            filter="select=eq(n\\,$(row.frame_index))"
            P.run_logged(`/usr/bin/ffmpeg -v error -nostdin -y -i $video -vf $filter -frames:v 1 $expected`;operation="independent expected frame")
            @test P.sha256_file(expected)==P.sha256_file(row.image_path)
        end
    end
end


@testset "coherent waveform segment mismatch remains diagnostic" begin
    anchors=[Dict{String,Any}("recording_time_s"=>Float64(t),"video_global_time_s"=>2.0+t+(t>=800 ? 0.04 : 0.0),
        "accepted"=>true,"audio_file"=>"rec$(t÷800)","video_file"=>"chapter$(t÷500)") for t in 0:100:1500]
    P.prep_fit_sync(anchors)
    diagnostic=P.prep_sync_diagnostics(anchors)
    @test diagnostic["coherent_segment_mismatch"]
    @test abs(diagnostic["best_chronological_step"]["offset_step_s"])≈0.04 atol=1e-9
    @test !isempty(diagnostic["source_groups"])
end


end # broad controls; set CLICKER_PREPARE_TEST_SCOPE=snippet for affected controls only

if get(ENV,"CLICKER_PREPARE_TEST_SCOPE","all") in ("all","warp")
@testset "known offset and drift removes half-window bias" begin
    fs=48000;duration=6.0;offset=0.3;scale=1.0002
    f(t)=sinpi(t/duration)^2*(sin(2pi*(180t+37t^2))+0.8sin(2pi*(490t+57t^2))+0.5sin(2pi*(800t+68t^2)))
    template=[f((i-1)/fs) for i in 1:round(Int,duration*fs)]
    camera=[begin t=((i-1)/fs-offset)/scale;0<=t<duration ? f(t) : 0.0 end for i in 1:round(Int,6.8fs)]
    old=P.prep_wave_match(camera,template;fs=fs)
    warped=P.prep_warp_template(template,scale)
    fixed=P.prep_wave_match(camera,warped;fs=fs)
    @test abs(old.time_s-offset)>0.0002
    @test fixed.time_s≈offset atol=1/fs
    @test fixed.score>0.999
    @test length(warped)==ceil(Int,length(template)*scale)
    @test_throws ErrorException P.prep_warp_template(template,1.1)
end

@testset "multi-WAV and multi-chapter waveform synchronization" begin
    mktempdir() do dir
        Random.seed!(414)
        lowfs=2000; seconds=50; raw=0.1randn(seconds*lowfs)
        # Recorder starts at camera2.5s; camera clock is 200ppm faster.
        offset=2.5; scale=1.0002
        camera=zeros(54lowfs)
        for k in eachindex(camera)
            pos=((k-1)/lowfs-offset)/scale*lowfs+1
            if 1<=pos<length(raw)
                j=floor(Int,pos);w=pos-j
                camera[k]=(1-w)*raw[j]+w*raw[j+1]
            end
        end
        r192=P.prep_resample_exact(raw,192000,2000)
        paths=String[]
        for j in 1:2
            path=joinpath(dir,"rec$j.WAV");a=r192[(j-1)*25*192000+1:j*25*192000]
            P.WAV.wavwrite(repeat(a,1,5),path;Fs=192000,nbits=24,compression=1)
            push!(paths,path)
        end
        infos=P.prep_wav_info.(paths);P.prep_boundaries(infos)
        camwav=joinpath(dir,"camera.wav")
        P.WAV.wavwrite(P.prep_resample_exact(camera,48000,2000),camwav;Fs=48000,nbits=16,compression=1)
        videos=Dict{String,Any}[]
        for j in 1:3
            path=joinpath(dir,"chapter$j.mp4")
            P.run_logged(`/usr/bin/ffmpeg -v error -nostdin -y -f lavfi -i color=size=320x180:rate=1 -ss $((j-1)*18) -i $camwav -t 18 -map 0:v:0 -map 1:a:0 -c:v libx264 -c:a aac -b:a 192k $path`;operation="synthetic camera chapter")
            push!(videos,Dict("path"=>path,"duration_s"=>18.0))
        end
        sync=P.prep_synchronize(infos,videos)
        @test sync["status"]=="valid"
        @test sync["offset_s"]≈offset atol=2/48000
        @test sync["drift_ppm"]≈200 atol=15
        @test length(sync["covered_audio_files"])==2
        @test length(sync["covered_video_files"])==3
        @test sync["max_residual_s"]<0.003
        @test length(sync["refinement_history"])==2
        @test sync["clock_scale_uncertainty"]>0
        @test P.prep_sync_uncertainty(sync,25.0)>=sync["max_residual_s"]+1/48000+sync["template_scale_bias_bound_s"]
    end
end


end

if get(ENV,"CLICKER_PREPARE_TEST_SCOPE","all") in ("all","snippet")
@testset "authenticated exact prepared snippet cache" begin
    mktempdir() do dir
        Random.seed!(819);x=1e-5randn(40000,5);x[18000,:].+=0.9
        source=joinpath(dir,"input.WAV")
        P.WAV.wavwrite(x,source;Fs=192000,nbits=24,compression=1)
        out=joinpath(dir,"prepared")
        result=P.prepare_device(dir,out;device="cache-control",detection_only=true)
        events=P.CSV.read(joinpath(out,"events.csv"),P.DataFrame)
        @test nrow(events)>=1
        row=events[1,:];snip=P.read_event_snippet(out,row)
        @test size(snip)==(row.snippet_samples,3)
        @test eltype(snip)==Float64
        @test P.prep_minmax(snip)==[row.minmax2_d1,row.minmax2_d2,row.minmax2_d3]
        @test row.snippet_peak_index==300
        @test P.prep_gcc(snip)[1]==[row.gcc_d1,row.gcc_d2,row.gcc_d3]
        @test result["snippet_cache"]["sha256"]==P.sha256_file(joinpath(out,"event_snippets.f64"))
        @test result["snippet_cache"]["endianness"]=="little"
        bad=Dict{Symbol,Any}(pairs(NamedTuple(row)));bad[:snippet_byte_offset]=-8
        @test_throws ErrorException P.read_event_snippet(out,bad)
        bad[:snippet_byte_offset]=1
        @test_throws ErrorException P.read_event_snippet(out,bad)
        bad[:snippet_byte_offset]=filesize(joinpath(out,"event_snippets.f64"))+8
        @test_throws ErrorException P.read_event_snippet(out,bad)
        bad=Dict{Symbol,Any}(pairs(NamedTuple(row)));bad[:snippet_sha256]=repeat("0",64)
        @test_throws ErrorException P.read_event_snippet(out,bad)
        bad=Dict{Symbol,Any}(pairs(NamedTuple(row)));bad[:event_id]="unbound-event"
        @test_throws ErrorException P.read_event_snippet(out,bad)
        cache=joinpath(out,"event_snippets.f64");original=read(cache);changed=copy(original);changed[1]⊻=0x01;write(cache,changed)
        @test_throws ErrorException P.read_event_snippet(out,row)
        # Even an updated whole-cache digest cannot bypass the per-payload digest.
        manifest=P.read_json(joinpath(out,"prepare_manifest.json"));manifest["snippet_cache"]["sha256"]=P.sha256_file(cache)
        for entry in manifest["outputs"]
            entry["path"]==cache && (entry["sha256"]=P.sha256_file(cache))
        end
        P.write_json(joinpath(out,"prepare_manifest.json"),manifest)
        @test_throws ErrorException P.read_event_snippet(out,row)
        @test P.read_event_snippet(out,row;verify=false)!=snip
    end
end

end # focused snippet scope

catch exception
    if exception isa Test.TestSetException
        for result in exception.errors_and_fails
            P.log_failure("synthetic prepare test detail",ErrorException(sprint(show,result)),catch_backtrace())
        end
    end
    P.log_failure("synthetic prepare test suite",exception,catch_backtrace())
    rethrow()
end
