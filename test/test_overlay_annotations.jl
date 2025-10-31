using Test
using Dates
# no external FilePaths dependency
using DataFrames
using VideoIO
using Images
# Provide a noop @animate macro if Makie isn't available so the file can be parsed
macro animate(expr...)
    return esc(:(begin end))
end
using Printf

include(joinpath(@__DIR__, "..", "src", "video.jl"))

function make_synthetic_video(path::AbstractString; nframes=8, w=160, h=120, fps=4)
    # create temporary frames and use pic2vid to encode (avoids relying on VideoIO writer API)
    frames_dir = joinpath(dirname(path), "frames_tmp")
    mkpath(frames_dir)
    for i in 1:nframes
        img = fill(RGB{N0f8}(i / nframes, 0.2, 0.2), h, w)
        save(joinpath(frames_dir, @sprintf("%d.png", i)), img)
    end
    pic2vid(frames_dir, path; counter="%d.png", fps=fps, auto_mode=false)
    rm(frames_dir; force=true, recursive=true)
    return path
end

@testset "overlay_annotations_on_video" begin
    tmp = mktempdir()
    invid = joinpath(tmp, "in_test.mkv")
    outvid = joinpath(tmp, "out_test_stream.mkv")
    outvid2 = joinpath(tmp, "out_test_frames.mkv")

    # create synthetic video
    make_synthetic_video(invid; nframes=6, w=80, h=60, fps=4)

    # create synthetic annotations DataFrame: points moving diagonally
    df = DataFrame(frame = 1:6, px = collect(10:10:60), py = collect(5:10:55))

    # stream mode (default)
    res = overlay_annotations_on_video([df], invid, outvid; mode=:stream)
    @test isfile(outvid) && filesize(outvid) > 0

    # frames mode (write frames then encode)
    res2 = overlay_annotations_on_video([df], invid, outvid2; mode=:frames)
    @test isfile(outvid2) && filesize(outvid2) > 0

    # cleanup
    rm(tmp; force=true, recursive=true)
end

println("overlay_annotations_on_video tests completed")
