#!/usr/bin/env julia

# pic2vid.jl
# Convert a folder of images to a video using ffmpeg.
# Usage:
#   pic2vid.jl INFILEPATH OUTFILEPATH [COUNTER_PATTERN] [-r FPS] [--auto]
#
# Behavior mirrors the provided pic2vid.sh:
# - If --auto is given, builds a concat-compatible list.txt with all common image
#   files (png,jpg,jpeg,bmp,gif) sorted by name and runs:
#       ffmpeg -r FPS -f concat -safe 0 -i list.txt -vcodec libx264 -b 5m -vf format=yuv420p OUT
# - Otherwise uses a printf-style COUNTER_PATTERN (default "%d") and assumes
#   .png extension, running:
#       ffmpeg -framerate FPS -i INDIR/COUNTER_PATTERN.png -vcodec libx264 -b 5m -vf format=yuv420p OUT
#
# This script avoids external Julia package deps and uses Base only.

using Dates
using FFMPEG

function usage()
    println("usage: pic2vid.jl INFILEPATH OUTFILEPATH [COUNTER_PATTERN] [-r FPS] [--auto]")
    println("  INFILEPATH         Path to input image folder")
    println("  OUTFILEPATH        Output video file path")
    println("  COUNTER_PATTERN    (optional) printf-style pattern for image numbering, e.g. %06d (default: %d)")
    println("  -r FPS             (optional) set output frame rate (default: 4)")
    println("  --auto             (optional) automatically process all image files in folder, sorted by name")
end

function is_flag(s)
    startswith(s, "-")
end

function build_filelist(in_folder::String)
    exts = Set([".png", ".jpg", ".jpeg", ".bmp", ".gif"])
    all = readdir(in_folder)
    files = [f for f in all if lowercase(splitext(f)[2]) in exts]
    sort!(files)
    return files
end

function run_ffmpeg_concat(listfile::String, fps::Int, out_video::String)
    cmd = `ffmpeg -r $fps -f concat -safe 0 -i $listfile -vcodec libx264 -b 5m -vf format=yuv420p $out_video`
    println("Running: ", cmd)
    @ffmpeg_env run(cmd)
end

function run_ffmpeg_pattern(pattern::String, fps::Int, out_video::String)
    cmd = `ffmpeg -framerate $fps -i $pattern -vcodec libx264 -b 5m -vf format=yuv420p $out_video`
    println("Running: ", cmd)
    @ffmpeg_env run(cmd)
end

function pic2vid(input_folder::AbstractString, output_vidpath::AbstractString; counter::AbstractString="%d", fps::Integer=4, auto_mode::Bool=false)
    # Validate input folder
    if !isdir(input_folder)
        throw(ArgumentError("input folder does not exist: $input_folder"))
    end

    if auto_mode
        files = build_filelist(input_folder)
        if isempty(files)
            throw(ArgumentError("No image files found in $input_folder"))
        end
        listfile = joinpath(input_folder, "list.txt")
        open(listfile, "w") do io
            for f in files
                p = joinpath(input_folder, f)
                # escape single quote
                p2 = replace(p, "'" => "\\'")
                println(io, "file '", p2, "'")
            end
        end
        try
            run_ffmpeg_concat(listfile, fps, output_vidpath)
        finally
            isfile(listfile) && rm(listfile; force=true)
        end
    else
        # assume png extension unless pattern already contains an extension
        pattern = counter
        if !occursin(r"\.[a-zA-Z0-9]+$", pattern)
            pattern *= ".png"
        end
        pattern = joinpath(input_folder, pattern)
        run_ffmpeg_pattern(pattern, fps, output_vidpath)
    end
end

function main()
    args = copy(ARGS)
    if length(args) < 2
        usage()
        return
    end

    in_folder = args[1]
    out_video = args[2]
    counter = "%d"
    fps = 4
    auto_mode = false

    # Remove the first two positional args
    deleteat!(args, 1:2)

    # If next arg exists and doesn't start with '-', treat as counter pattern
    if !isempty(args) && !is_flag(args[1])
        counter = args[1]
        deleteat!(args, 1)
    end

    # parse remaining flags
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--auto"
            auto_mode = true
            i += 1
        elseif a == "-r" || a == "--fps"
            if i == length(args)
                println("Error: $a requires an argument")
                return
            end
            fps = tryparse(Int, args[i+1])
            if fps === nothing
                println("Error: invalid fps value: ", args[i+1])
                return
            end
            i += 2
        else
            println("Unknown option: ", a)
            usage()
            return
        end
    end

    try
        pic2vid(in_folder, out_video; counter=counter, fps=fps, auto_mode=auto_mode)
    catch e
        println("Error: ", e)
    end
end

# Execute main only when this file is the entry-point script (called from CLI)
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
