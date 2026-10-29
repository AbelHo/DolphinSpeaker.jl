using Test
include(joinpath(@__DIR__,"../../src/clicker_calibration/ClickerCalibration.jl"))
using .ClickerCalibration
const C=ClickerCalibration

@testset "shared command and durable error contract" begin
    mktempdir() do dir
        C.CALIBRATION_ERROR_LOG[]=joinpath(dir,"errors.jsonl")
        @test C.run_logged(`printf abc`;operation="bytes",capture=true)==UInt8[0x61,0x62,0x63]
        @test_throws ErrorException C.run_logged(`sh -c 'echo expected-fixture-error >&2; exit 3'`;operation="expected fixture subprocess failure",capture=true)
        lines=readlines(C.CALIBRATION_ERROR_LOG[])
        @test length(lines)==1
        @test occursin("expected-fixture-error",lines[1])
        @test occursin("expected fixture subprocess failure",lines[1])
        path=joinpath(dir,"safe.json")
        C.write_json(path,Dict("a"=>NaN,"b"=>Inf,"c"=>[1.0,2.0]))
        @test C.read_json(path)["a"]===nothing
        @test C.read_json(path)["b"]===nothing
        @test C.read_json(path)["c"]==[1.0,2.0]
    end
end

@testset "external annotation handoff and resumable empty stage" begin
    mktempdir() do dir
        source=joinpath(dir,"input");out=joinpath(dir,"output");mkpath(joinpath(source,"A"));dest=joinpath(out,"A");mkpath(dest)
        write(joinpath(dest,"frames.csv"),"frame_id,video_file,frame_index,frame_pts_s,video_global_time_s,image_path,width,height\n")
        state=C.run_clicker_calibration(source;devices=(:A,),output_dir=out,stage=:annotate,
            annotator=:external,notify=:none,wait_annotations=false,error_log=joinpath(dir,"errors.jsonl"))
        @test state["status"]=="needs_review"
        @test state["devices"]["A"]["status"]=="no_frames_to_annotate"
        write(joinpath(dest,"frames.csv"),"frame_id,video_file,frame_index,frame_pts_s,video_global_time_s,image_path,width,height\nA-f1,v.mp4,1,0.04,0.04,x.png,320,180\n")
        request=C.request_annotations(dest;annotator=:external)
        @test request["status"]=="awaiting_external_annotation"
        @test request["frame_count"]==1
        prompt=read(request["prompt_path"],String)
        @test occursin("gpt-6.1-terra",prompt)
        @test occursin("Do not open events.csv",prompt)
        @test occursin("medium: underwater/air/mixed/unknown",prompt)
        previous=Dict("status"=>"launched","frames_sha256"=>C.sha256_file(joinpath(dest,"frames.csv")),"agent_name"=>"old-annotation-agent")
        C.write_json(joinpath(dest,"annotation_request.json"),previous)
        spawned=Any[]
        live=C.request_annotations(dest;agent_lookup=name->Dict("name"=>name,"agent_status"=>"working"),
            agent_spawn=(name,prompt)->push!(spawned,(name,prompt)))
        @test live["agent_name"]=="old-annotation-agent"
        @test isempty(spawned)
        replacement=C.request_annotations(dest;agent_lookup=name->nothing,
            agent_spawn=(name,prompt)->push!(spawned,(name,prompt)))
        @test length(spawned)==1
        @test replacement["agent_name"]!="old-annotation-agent"
        @test C.read_json(joinpath(dest,"annotation_request.previous.json"))==previous
        state=C.run_clicker_calibration(source;devices=(:A,),output_dir=out,stage=:fit,
            annotator=:external,notify=:none,wait_annotations=false,error_log=joinpath(dir,"errors.jsonl"))
        @test state["status"]=="awaiting_annotations"
        @test_throws ErrorException C.run_clicker_calibration(source;output_dir=joinpath(source,"outputs"),notify=:none)
    end
end
