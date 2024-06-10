using FilePaths, Glob
# @info (ARGS[1], ARGS[2])
# "julia /home/sousa/GitHub/DolphinSpeaker.jl/src/mkv2mp4_frameno.jl"
in_dir = "/run/user/1000/gvfs/smb-share:server=10.246.128.21,share=data/Concretecho/data/temp/try8"
out_dir = "/run/user/1000/gvfs/smb-share:server=10.246.128.21,share=results/Concretecho/vid3"

function convert_vid(src::String, dst::String)
    if endswith(src, ".mkv")
        dst = joinpath(dst, replace(basename(src), ".mkv" => ".mp4"))
        @info (src, dst)
        run(`ffmpeg -i $src -filter_complex "drawtext=text='%{n}': x=10: y=35: fontsize=48: fontcolor=red" $dst -loglevel error`)
    end
end

function convert_videos_old(src_dir::String, dst_dir::String)
    mkv_files = glob("*.mkv", src_dir)

    for file in mkv_files
        rel_path = fileparts(file)[1][length(src_dir)+1:end]
        dst_file = joinpath(dst_dir, replace(rel_path, ".mkv" => ".mp4"))
        mkpath(dirname(dst_file))
        run(`ffmpeg -i $file -filter_complex "drawtext=text='%{n}': x=10: y=35: fontsize=48: fontcolor=red" $dst_file -loglevel error`)
    end
end

function convert_videos(src_dir::String, dst_dir::String)
    mkv_files = walkdir(src_dir) |> collect |> x -> x[2] |> x -> filter(i -> endswith(i, ".mkv"), x)

    for file in mkv_files
        rel_path = file[length(src_dir)+1:end]
        dst_file = joinpath(dst_dir, replace(rel_path, ".mkv" => ".mp4"))
        mkpath(dirname(dst_file))
        @info (file, dst_file)
        run(`ffmpeg -i $file -filter_complex "drawtext=text='%{n}': x=10: y=35: fontsize=48: fontcolor=red" $dst_file -loglevel error`)
    end
end

# convert_videos(ARGS[1], ARGS[2])