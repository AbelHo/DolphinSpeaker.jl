using Test
module SyncEligibilityTest
using CSV,DataFrames,JSON,SHA,Dates,Statistics,LinearAlgebra,Random,DSP,FFTW,WAV,Peaks
include(joinpath(@__DIR__,"../../src/clicker_calibration/support.jl"))
include(joinpath(@__DIR__,"../../src/clicker_calibration/prepare.jl"))
include(joinpath(@__DIR__,"../../src/clicker_calibration/sync_eligibility.jl"))
end
const SE=SyncEligibilityTest

function fixture(dir; rounded=false)
    mkpath(dir)
    device=rounded ? "B-synthetic" : "A-synthetic"
    videos=Any[];hashes=Dict{String,String}()
    exact=[0.0,40.0,80.0];stored=copy(exact)
    if rounded
        exact=[0.0,40.000333,80.000666];stored=[0.0,40.0005,80.001]
    end
    for i in 1:3
        path=joinpath(dir,"chapter$i.mp4");hashes[path]=repeat(string(i),64)
        push!(videos,Dict{String,Any}("path"=>path,"pts"=>collect(0:0.01:39.99),"global_start_s"=>stored[i],
            "duration_s"=>40.0,"audio_duration_s"=>rounded ? 40.000333 : 40.0,"audio_start_s"=>0.0,"video_start_s"=>0.0,
            "width"=>2704,"height"=>1520))
    end
    true_scale=1.00006;true_offset=-1.25;anchors=Any[]
    for i in 1:3,localtime in (6.0,18.0,30.0)
        t=(exact[i]+localtime-true_offset)/true_scale
        push!(anchors,Dict{String,Any}("recording_time_s"=>t,"video_global_time_s"=>stored[i]+localtime,
            "video_file"=>videos[i]["path"],"audio_file"=>joinpath(dir,t<60 ? "rec1.WAV" : "rec2.WAV"),
            "accepted"=>true,"fit_inlier"=>true,"template_clock_scale"=>true_scale,"window_s"=>6.0))
    end
    x=Float64[a["recording_time_s"] for a in anchors];y=Float64[a["video_global_time_s"] for a in anchors]
    slope=sum((x.-SE.mean(x)).*(y.-SE.mean(y)))/sum(abs2,x.-SE.mean(x));offset=SE.mean(y)-slope*SE.mean(x)
    residual=y.-offset.-slope.*x;maxres=max(maximum(abs,residual),1e-8)
    scale_error=max(abs(slope-true_scale),2*(maxres+1/48000)/(maximum(x)-minimum(x)))
    sync=Dict{String,Any}("status"=>"valid","anchors"=>anchors,"offset_s"=>offset,"scale"=>slope,
        "max_residual_s"=>maxres,"anchor_min_time_s"=>minimum(x),"anchor_max_time_s"=>maximum(x),
        "clock_scale_uncertainty"=>scale_error,"template_scale_bias_bound_s"=>3scale_error)
    metadata=Dict("videos"=>videos)
    request=Dict("config"=>Dict("sync_residual_limit_s"=>0.015),"source_sha256"=>repeat("a",64),
        "inputs"=>[Dict("path"=>v["path"],"sha256"=>hashes[v["path"]]) for v in videos])
    manifest=Dict{String,Any}("status"=>"complete","event_count"=>0,"frame_count"=>0,"mapped_count"=>0,
        "request"=>request,"sync_status"=>"valid","outputs"=>Any[],"completed_utc"=>"2026-01-01T00:00:00")
    events=SE._prep_table(SE.PREPARE_COLUMNS)
    function add_event(t,status,id)
        t=round(Int,t*192000)/192000
        oldtime=offset+slope*t;h=SE._syncelig_frame(videos,stored,oldtime,device,hashes)
        row=Dict{Symbol,Any}(k=>missing for k in SE.PREPARE_COLUMNS)
        payload=SE._prep_snippet_bytes(reshape(Float64.(1:21),7,3)./100)
        eventid=device*"-"*id
        merge!(row,Dict(:event_id=>eventid,:device=>device,:audio_file=>joinpath(dir,"recø.WAV"),
            :sample_index=>round(Int,t*192000),:audio_time_s=>t,:recording_time_s=>t,:fs=>192000,
            :sync_error_s=>SE._syncelig_uncertainty(sync,t),:sync_status=>status,
            :window_status=>"complete; retained\nmultiline,\"quoted\"",:window_source_segments=>SE.JSON.json([Dict("audio_file"=>"recø.WAV","first_sample"=>0,"samples"=>7)]),
            :snippet_file=>"event_snippets.f64",:snippet_byte_offset=>nrow(events)*length(payload),:snippet_samples=>7,:snippet_channels=>3,
            :snippet_peak_index=>3,:snippet_sha256=>SE._syncelig_digest(payload),:minmax2_d1=>1.2345678901234567,:minmax2_d2=>-0.0,:minmax2_d3=>0.0,
            :minmax_d1=>3.0,:minmax_d2=>4.0,:minmax_d3=>0.0,:gcc_d1=>1.125,:gcc_d2=>-0.375,:gcc_d3=>0.0,:gcc_closure_samples=>0.25))
        if !isnothing(h)
            merge!(row,Dict(:video_file=>h["video_file"],:video_time_s=>oldtime-stored[findfirst(v->v["path"]==h["video_file"],videos)],
                :video_global_time_s=>oldtime,:frame_index=>h["frame_index"],:frame_pts_s=>h["frame_pts_s"],
                :frame_id=>h["frame_id"],:image_path=>joinpath(dir,h["frame_id"]*".png"),:width=>2704,:height=>1520))
        else
            row[:sync_error_s]=missing
        end
        push!(events,row)
        eventid
    end
    stable=add_event(20.127,"mapped","stable")
    preserved=add_event(24.127,"ambiguous","original-ambiguous")
    outside=add_event(130.0,"outside_video","original-outside")
    target=nothing
    if rounded
        # Independent analytic counterfactual true clock gives the relative local
        # shift. Choose a formerly-safe side of a half-frame boundary so that
        # the union interval crosses after declared-base correction.
        for i in 1:3,mid in (8.005,16.005,24.005,32.005),sign in (-1.0,1.0)
            t0=(stored[i]+mid-offset)/slope
            err=SE._syncelig_uncertainty(sync,t0)
            local_shift=(true_offset+true_scale*t0-exact[i])-mid
            local_shift*sign<0 || continue
            pad=min(abs(local_shift)/4,2e-5)
            t=round(Int,((stored[i]+mid+sign*(err+pad)-offset)/slope)*192000)/192000
            oldtime=offset+slope*t;oldhit=SE._syncelig_frame(videos,stored,oldtime,device,hashes)
            SE._syncelig_ambiguous(videos,stored,oldtime,SE._syncelig_uncertainty(sync,t),device,hashes) && continue
            newtime=true_offset+true_scale*t
            if SE._syncelig_ambiguous(videos,exact,newtime,SE._syncelig_uncertainty(sync,t),device,hashes)
                target=add_event(t,"mapped","new-boundary");break
            end
        end
        isnothing(target) && error("Synthetic fixture failed to construct timing-boundary event")
    end
    # A rejected saved window is never silently admitted into the refit.
    push!(anchors,Dict("accepted"=>false,"fit_inlier"=>false,"recording_time_s"=>99999.0))
    SE.CSV.write(joinpath(dir,"events.csv"),events)
    # Lexical patching must preserve CRLF plus quoted multiline/unicode fields.
    bytes=read(joinpath(dir,"events.csv"));bytes=Vector{UInt8}(codeunits(replace(String(bytes),"\n"=>"\r\n")))
    write(joinpath(dir,"events.csv"),bytes)
    frames=SE._prep_table(SE.PREPARE_FRAME_COLUMNS)
    for e in eachrow(events)
        ismissing(e.frame_id) && continue
        push!(frames,(e.frame_id,e.video_file,e.frame_index,e.frame_pts_s,e.video_global_time_s,e.image_path,e.width,e.height))
        write(e.image_path,"synthetic immutable image")
    end
    SE.CSV.write(joinpath(dir,"frames.csv"),frames)
    SE.write_json(joinpath(dir,"synchronization.json"),sync);SE.write_json(joinpath(dir,"prepare_metadata.json"),metadata)
    payload=SE._prep_snippet_bytes(reshape(Float64.(1:21),7,3)./100)
    write(joinpath(dir,"event_snippets.f64"),repeat(payload,nrow(events)))
    paths=[joinpath(dir,n) for n in ["events.csv","frames.csv","synchronization.json","prepare_metadata.json","event_snippets.f64"]]
    append!(paths,String.(frames.image_path))
    manifest["outputs"]=SE._prep_inventory(paths);manifest["event_count"]=nrow(events);manifest["frame_count"]=nrow(frames)
    manifest["mapped_count"]=count(==("mapped"),events.sync_status)
    manifest["snippet_cache"]=Dict("path"=>joinpath(dir,"event_snippets.f64"),"sha256"=>SE.sha256_file(joinpath(dir,"event_snippets.f64")),
        "bytes"=>filesize(joinpath(dir,"event_snippets.f64")),"dtype"=>"Float64","endianness"=>"little","layout"=>"column-major")
    SE.write_json(joinpath(dir,"prepare_manifest.json"),manifest)
    (;stable,preserved,outside,target,manifest,events,paths)
end

using .SyncEligibilityTest: nrow,eachrow
try
@testset "zero-addition A-style timing review preserves all canonical bytes" begin
    mktempdir() do dir
        f=fixture(dir)
        paths=vcat(f.paths,[joinpath(dir,"prepare_manifest.json")])
        before=Dict(p=>(SE.sha256_file(p),stat(p).mtime,stat(p).inode) for p in paths)
        q=SE.review_sync_eligibility(dir)
        @test q["status"]=="complete"
        @test q["additional_exclusion_count"]==0
        @test !q["canonical_changes_applied"]
        @test q["mapped_count_after"]==f.manifest["mapped_count"]
        @test !ispath(joinpath(dir,"preparation_raw"))
        @test all(p->(SE.sha256_file(p),stat(p).mtime,stat(p).inode)==before[p],paths)
        table=SE.CSV.read(joinpath(dir,"events.csv"),SE.DataFrame)
        @test only([r for r in eachrow(table) if r.event_id==f.preserved]).sync_status=="ambiguous"
        @test only([r for r in eachrow(table) if r.event_id==f.outside]).sync_status=="outside_video"
        qcsha=SE.sha256_file(joinpath(dir,"sync_eligibility.json"))
        @test SE.review_sync_eligibility(dir)==q
        @test SE.sha256_file(joinpath(dir,"sync_eligibility.json"))==qcsha
        @test all(p->(SE.sha256_file(p),stat(p).mtime,stat(p).inode)==before[p],paths)
    end
end
@testset "B-style uncertainty union, byte-exact originals, cache and idempotence" begin
    mktempdir() do dir
        f=fixture(dir;rounded=true)
        originalevents=read(joinpath(dir,"events.csv"));originalmanifest=read(joinpath(dir,"prepare_manifest.json"))
        others=Dict(p=>SE.sha256_file(p) for p in f.paths if basename(p)!="events.csv")
        table=SE.CSV.read(joinpath(dir,"events.csv"),SE.DataFrame)
        payloads=Dict(String(r.event_id)=>SE.read_event_snippet(dir,r) for r in eachrow(table))
        q=SE.review_sync_eligibility(dir)
        @test q["status"]=="complete"
        @test q["additional_excluded_event_ids"]==[f.target]
        @test q["additional_exclusion_count"]==1
        @test q["mapped_count_after"]==f.manifest["mapped_count"]-1
        @test read(joinpath(dir,"preparation_raw/events.csv"))==originalevents
        @test read(joinpath(dir,"preparation_raw/prepare_manifest.json"))==originalmanifest
        updated=SE.CSV.read(joinpath(dir,"events.csv"),SE.DataFrame)
        @test only([r for r in eachrow(updated) if r.event_id==f.target]).sync_status=="ambiguous"
        @test only([r for r in eachrow(updated) if r.event_id==f.preserved]).sync_status=="ambiguous"
        @test only([r for r in eachrow(updated) if r.event_id==f.outside]).sync_status=="outside_video"
        @test only([r for r in eachrow(updated) if r.event_id==f.stable]).sync_status=="mapped"
        @test all(p->SE.sha256_file(p)==others[p],keys(others))
        @test isequal(SE._syncelig_check_table_change(originalevents,read(joinpath(dir,"events.csv")),Set([f.target])),updated)
        m=SE.read_json(joinpath(dir,"prepare_manifest.json"));lineage=only(m["postprocessing"])
        @test m["request"]==f.manifest["request"]
        @test m["completed_utc"]==f.manifest["completed_utc"]
        @test m["mapped_count"]==q["mapped_count_after"]
        @test lineage["raw_manifest_sha256"]==SE._syncelig_digest(originalmanifest)
        @test lineage["raw_events_sha256"]==SE._syncelig_digest(originalevents)
        @test lineage["review_source_sha256"]==SE.sha256_file(joinpath(@__DIR__,"../../src/clicker_calibration/sync_eligibility.jl"))
        @test lineage["events_sha256_after"]==SE.sha256_file(joinpath(dir,"events.csv"))
        @test SE._syncelig_non_event_outputs(m,joinpath(dir,"events.csv"))==SE._syncelig_non_event_outputs(f.manifest,joinpath(dir,"events.csv"))
        @test all(r->SE.read_event_snippet(dir,r)==payloads[String(r.event_id)],eachrow(updated))
        canonical=Dict(n=>SE.sha256_file(joinpath(dir,n)) for n in ["events.csv","prepare_manifest.json","sync_eligibility.json","preparation_raw/events.csv","preparation_raw/prepare_manifest.json"])
        @test SE.review_sync_eligibility(dir)==q
        @test all(n->SE.sha256_file(joinpath(dir,n))==canonical[n],keys(canonical))
        # Recovery from the only non-atomic interval: table committed, producer
        # manifest still raw, originals already safe. No second raw copy allowed.
        write(joinpath(dir,"prepare_manifest.json"),originalmanifest)
        rm(joinpath(dir,"sync_eligibility.json"))
        recovered=SE.review_sync_eligibility(dir)
        @test recovered["status"]=="complete"
        @test recovered["additional_excluded_event_ids"]==[f.target]
        @test read(joinpath(dir,"preparation_raw/events.csv"))==originalevents
        @test read(joinpath(dir,"preparation_raw/prepare_manifest.json"))==originalmanifest
        @test all(r->SE.read_event_snippet(dir,r)==payloads[String(r.event_id)],eachrow(SE.CSV.read(joinpath(dir,"events.csv"),SE.DataFrame)))
        # Complete-cache verification precedes resume acceptance.
        cache=joinpath(dir,"event_snippets.f64");bytes=read(cache);changed=copy(bytes);changed[1]⊻=0x01;write(cache,changed)
        @test_throws ErrorException SE.review_sync_eligibility(dir)
    end
end
catch exception
    if exception isa Test.TestSetException
        for result in exception.errors_and_fails
            SE.log_failure("sync eligibility test detail",ErrorException(sprint(show,result)),catch_backtrace())
        end
    end
    SE.log_failure("sync eligibility synthetic tests",exception,catch_backtrace())
    rethrow()
end
