# Run examples for overlay_annotations_on_video
# Edit paths below to point to actual video files on your system before running.

using CSV, DataFrames

# Adjust these to real files on your machine
video_input = "/media/spin/anas2/data/calf/upload/results/good_stuff/results5_20251031_6e__Single_dolphin_target/17_10_2025_Tutti/combined__17_10_2025_Tutti.mp4.mp4"
output_dir = dirname(video_input)

mkpath(output_dir)

csv_path = joinpath(@__DIR__, "sample_annotations_all_features.csv")

println("Sample CSV: $csv_path")

# Example A: Per-row properties read from CSV columns
# This uses r,g,b,a,radius,shape columns in the CSV
outA = joinpath(output_dir, "output_in_annotations.mkv")
println("Running per-row (:in_annotations) example -> $outA")
# Make sure src/video.jl is included (it defines overlay_annotations_on_video)
include(joinpath(@__DIR__, "..", "src", "video.jl"))

# Call using :in_annotations for color/alpha/radius/shape
try
    overlay_annotations_on_video(csv_path, video_input, outA;
        default_color = :in_annotations,
        default_alpha = :in_annotations,
        radius = :in_annotations,
        default_shape = :in_annotations,
        flag_dryrun = false # set to false to actually write video
    )
    println("Per-row example completed (dry run). Frames written to tmpdir when not dry-run.")
catch err
    println("Per-row example failed: $err")
end

# Example B: Array-of-annotation-sets mode
# Read the CSV twice and shift the second set horizontally to simulate a second object
df1 = CSV.read(csv_path, DataFrame)
# Create df2 as a shifted copy
df2 = deepcopy(df1)
df2.px .= df2.px .+ 40  # shift x positions to the right

outB = joinpath(output_dir, "output_array_mode.mp4")
println("Running array-of-sets example -> $outB")
try
    overlay_annotations_on_video([df1, df2], video_input, outB;
        default_color = ["red@0.6", "blue@0.6"],
        default_alpha = [0.6, 0.8],
        radius = [25, 40],
        default_shape = [:circle, :rect],
        flag_dryrun = true # set to false to actually write video
    )
    println("Array-of-sets example completed (dry run).")
catch err
    println("Array-of-sets example failed: $err")
end

println("\nDone. To actually create videos, set flag_dryrun=false and ensure ffmpeg/VideoIO are configured and video_input points to a real file.")
