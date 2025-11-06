using VideoIO
using FFMPEG
using Images
using Printf
using Plots
using ProgressMeter
include("media_info.jl")
include("synchronization.jl")
include("pic2vid.jl")
#  @time open_video_out(newvidname, img2, framerate=get_fps(vidfname), encoder_options=encoder_options) do writer
    
# resolution=(1080,720)
# fig = Figure(;resolution=result_resolution)

# v1="/Users/abel/Documents/data/concretecho/2023-12-04/cam_uw1/2023-12-04_12.33.03_uw1.mkv"
# v2="/Users/abel/Documents/data/concretecho/2023-12-04/cam_topview/2023-12-04_12.33.03_topview.mkv"
# au="/Users/abel/Documents/data/concretecho/2023-12-04/acoustic/2023-12-04_12.33.03.ogg"
function combine_2v1a(v1,v2,au,output_file; sync_type=:new, interval = [-.2 .5])
    @info(v1,v2,au,output_file)
    isdir(output_file) && (output_file = joinpath(output_file, splitext(basename(au))[1]*"_norm.mp4" ))
    # output_file_norm = splitext(output_file)[1]*"_norm.mp4"

    if sync_type == :old
        v1_trigger = findVidAudioBlip(v1; argmax_len=0, plot_window_inS=nothing)
        v2_trigger = findVidAudioBlip(v2; argmax_len=0, plot_window_inS=nothing)
        au_trigger = findAudioBlip(au; argmax_len=0, plot_window_inS=nothing)
    else
        v1_trigger = findVidAudioBlip(v1; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5, flag_savefig=dirname(output_file))
        v2_trigger = findVidAudioBlip(v2; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5, flag_savefig=dirname(output_file))
        au_trigger = findAudioBlip(au; argmax_len=0, plot_window_inS=interval, threshold_percentMAX=0.5, flag_savefig=dirname(output_file))
    end
    # mini = min(v1_trigger, v2_trigger, au_trigger)
    # @info mini
    v1_trigger_n = v1_trigger - au_trigger
    v2_trigger_n = v2_trigger - au_trigger
    # au_trigger -= mini
    @info "delays: " * string(v1_trigger) * " " * string(v2_trigger) * " " * string(au_trigger)

    # width = get_whatever(v1, 0, "v"; entries_custom="width") |> Int
    try
        println(`ffmpeg -itsoffset $v1_trigger_n -i "$v1" -itsoffset $v2_trigger_n -i "$v2" -i "$au" -filter_complex "[2:a]loudnorm[a];[a]asplit[a0][a1]; [a1]showwaves=s=1920x400:mode=cline:colors=red:rate=25,format=yuv420p[vau]; [0:v][1:v][vau]vstack=inputs=3[v]" -map "[v]" -map "[a0]" $output_file -hide_banner`)
        @ffmpeg_env run(`ffmpeg -itsoffset $v1_trigger_n -i "$v1" -itsoffset $v2_trigger_n -i "$v2" -i "$au" -filter_complex "[2:a]loudnorm[a];[a]asplit[a0][a1]; [a1]showwaves=s=1920x400:mode=cline:colors=red:rate=25,format=yuv420p[vau]; [0:v][1:v][vau]vstack=inputs=3[v]" -map "[v]" -map "[a0]" $output_file -hide_banner`)
    catch err
        @error "FAILED!:    $au"
    end
    # @ffmpeg_env run(`ffmpeg -i "$v1" -i "$v2" -i "$au" -filter_complex "[0:v][1:v]vstack=inputs=2[v];[2:a]loudnorm[a]" -map "[v]" -map "[a]" $output_file`)
    # @ffmpeg_env run(`ffmpeg -i $output_file `) #-af loudnorm=I=-16:LRA=11:TP=-1.5
    # ffmpeg -i "$v1" -i "$v2" -i "$au" -filter_complex "[0:v][1:v]vstack=inputs=2[v];[2:a]anull[a]" -map "[v]" -map "[a]" output.mp4
    return v1_trigger, v2_trigger, au_trigger
end

function combine_2v1a_auto(aufolder, outfolder; filetype=".ogg")
    mkpath(outfolder)
    dname = dirname(aufolder)
    for fname in readdir(aufolder)|>skiphiddenfiles
        fname_split = splitext(fname)
        thisfiletype = fname_split[2]
        fname_split = fname_split[1]
        if thisfiletype != filetype
            continue
        end
        @info fname

        try
            combine_2v1a(joinpath(dname,"cam_topview",fname_split*"_topview.mkv"), joinpath(dname,"cam_uw1",fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
        catch err
            combine_2v1a(joinpath(aufolder,fname_split*"_topview.mkv"), joinpath(aufolder,fname_split*"_uw1.mkv"), joinpath(aufolder,fname), joinpath(outfolder,fname_split*"_norm.mp4"))
        end
    end
end

function plot_signal2vid(data,fs, outvidname; fps=25, kwargs...)
    sig = signal(data,fs)
    t_width = 1/fps
    min_max = extrema(data)
    anim = @animate for t in 0:1/fps:(size(data,1) / fs)
        plot(sig; kwargs...)
        plot!([t, t+t_width, t+t_width, t], [min_max[1], min_max[1], min_max[2], min_max[2]], fill=false, c=:red, alpha=0.8)
    end

    outgifname = splitext(outvidname)[1]*".gif"
    gif(anim,  outgifname; fps = fps)
    # outvidname = splitext(outgifname)[1]*".mp4"
    @ffmpeg_env run(`$ffmpeg -i $outgifname -pix_fmt yuv420p $outvidname -hide_banner -y`)
end

"""
# Example usage:
```
vidfname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/20250227_10.16.59_log.mp4"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/new_2025_freeplay_report/video/tarsier/6s"
make_video_clips_ffmpeg(df, vidfname; outdir=res_dir, min_duration=6.0)
```
"""
function make_video_clips_ffmpeg(df::DataFrame, input_video::AbstractString; outdir="clips", min_duration=-1,
    subject_col=:Subject, behavior_col=:Behavior, start_col=Symbol("Start (s)"), stop_col=Symbol("Stop (s)"), duration_col=Symbol("Duration (s)"))
    
    mkpath(outdir)
    for row in eachrow(df)
        subject   = row[subject_col]
        behavior  = row[behavior_col]
        start     = row[start_col]
        stop      = row[stop_col]
        duration  = row[duration_col]
        # If duration is less than min_duration, extend stop time
        if duration < min_duration
            stop = start + min_duration
        end
        # Sanitize filename
        fname = "$(behavior)_$(subject)_$(start)_$(stop).mp4"
        fname = basename(input_video) *"_"* replace(fname, r"[^\w\.\-]" => "_")
        outfile = joinpath(outdir, fname)
        # ffmpeg command: -ss (start), -to (stop), -i (input), -c copy (no re-encoding)
        cmd = `ffmpeg -y -ss $start -to $stop -i $input_video -c copy $outfile`
        @debug cmd
        @ffmpeg_env run(cmd)
    end
end

"""
    split_video_by_duration(vidfilename::AbstractString, dur::Real; outdir="clips")

Split the video into multiple clips, each of length `dur` seconds.
Output files are named as: basename_start_end.mp4

## Example usage:
    split_video_by_duration("input.mp4", 5.0; outdir="clips")
"""
function split_video_by_duration(vidfilename::AbstractString, dur::Real; outdir="clips")
    mkpath(outdir)
    total_duration = get_duration(vidfilename)
    base = splitext(basename(vidfilename))[1]
    start = 0.0
    clip_idx = 1
    while start < total_duration
        stop = min(start + dur, total_duration)
        outname = joinpath(outdir, "$(base)_$(round(start,sigdigits=4))_$(round(stop,sigdigits=4)).mp4")
        cmd = `ffmpeg -y -ss $start -to $stop -i $vidfilename -c copy $outname`
        @debug cmd
        @ffmpeg_env run(cmd)
        start += dur
        clip_idx += 1
    end
end


# Example usage:
# overlay_boxes_on_video(
#     "detection_pixels.csv",
#     "/media/spin/anas2/data_res/dolphin/calf/temp/delete/1/combined__1.GoPro_Clicker.MP4.mp4",
#     "output_with_circles_fill.mkv";
#     radius=100
# )
function overlay_boxes_on_video(csv_path::String, video_path::String, output_path::String; radius::Int=100,
    flag_dryrun=false)
    df = CSV.read(csv_path, DataFrame)
    filters = String[]
    for row in eachrow(df)
        frame = round(Int, row.frame)
        x = round(Int, row.px)
        y = round(Int, row.py)
        push!(filters,
            "drawbox=x=$(x-radius÷2):y=$(y-radius÷2):w=$radius:h=$radius:color=red@0.5:t=fill:enable='eq(n,$frame)'"
        )
    end

    # If the combined filter string is very long, the shell/OS can hit ARG_MAX
    # and raise E2BIG. To avoid that, write the filtergraph to a temporary file
    # and pass it to ffmpeg using -filter_complex_script (or -vf script for a
    # single-input video). This keeps the command-line short.
    filter_str = join(filters, ",")

    # Create a temporary file for the filtergraph
    tmp = tempname()
    # ffmpeg expects a plain text file; use UTF-8
    open(tmp, "w") do io
        write(io, filter_str)
    end

    # Build ffmpeg command using -filter_complex_script when possible. For a
    # single video input, -vf script can be used but -filter_complex_script is
    # acceptable and general.
    cmd = `ffmpeg -y -i $video_path -filter_complex_script $tmp -codec:a copy $output_path`
    # cmd = `ffmpeg -y -i $video_path -/filter_complex $tmp -codec:a copy $output_path`
    println(cmd)
    flag_dryrun && (rm(tmp); return)

    try
        @ffmpeg_env run(cmd)
    finally
        # Ensure temp file is removed
        isfile(tmp) && rm(tmp)
    end
end

function overlay_boxes_on_video_imageonly(csv_path::String, video_path::String, output_path::String; radius::Int=100,
    flag_dryrun=false)
    fps = get_fps(video_path)
    duration = get_duration(video_path)
    vid = VideoIO.openvideo(video_path)

    df = CSV.read(csv_path, DataFrame)
    mkpath(output_path)

    # For each CSV row: seek to the corresponding frame, read the image, overlay points
    for row in eachrow(df)
        frame = round(Int, row.frame)
        time = (frame - 1) / fps
        if time > duration
            @warn "Frame $frame at time $time exceeds video duration $duration. Ending!......."
            break
        end
        img = nothing
        # read the frame at the requested time
        try
            seek(vid, time)
            img = read(vid)
        catch err
            @error "Failed to read frame $frame at time $time" err
            continue
        end

        x = round(Int, row.px)
        y = round(Int, row.py)

        # Prepare extra_arg in the format expected by overlay_points!
        # overlay_points!(img, counter, extra_arg) expects extra_arg to be an
        # iterable of tuples (pind_vidframes, p_pixels, colour, ptsize).
        pind_vidframes = [frame]
        p_pixels = reshape([x, y], 1, 2)   # 1 x 2 matrix: rows are points, cols are x,y
        colour = [1.0, 0.0, 0.0]
        ptsize = radius
        extra_arg = [(pind_vidframes, p_pixels, colour, ptsize)]

        try
            overlay_points!(img, frame, extra_arg)
        catch err
            @error "overlay_points! failed on frame $frame" err
        end

        # Save the overlaid image named "<frame>_<time-in-seconds>.png"
        time_str = string(round(time, digits=3))
        fname = joinpath(output_path, "$(frame)_$(time_str)s.png")
        @debug "Saving overlaid image to $fname"
        flag_dryrun && continue
        save(fname, img)
    end
end

"""
    overlay_annotations_on_video(annotations, video_path, output_path; kwargs...)

Overlay annotations onto a video using pure-Julia drawing (no ffmpeg CLI for overlay).

Arguments
- `annotations`: A vector where each element is either
    - a String path to a CSV (must contain columns `frame`, `px`, `py`), or
    - a `DataFrame` with columns `frame`, `px`, `py`, or
    - a NamedTuple / Dict with keys `:df` or `:csv` (DataFrame or CSV path) and optional `:color`, `:radius`, `:alpha`, `:shape`.
- `video_path`: path to input video
- `output_path`: path to output video

Keyword arguments
- `tmpdir`: temporary folder to write frames (default: created with `mktempdir()`)
- `fps`: output FPS for encoder (defaults to video's fps)
- `default_radius`: Can be a single Int, an array of Int (one per annotation set), or `:in_annotations` to read from "radius" column
- `default_color`: Can be a single String/Tuple, an array (one per annotation set), or `:in_annotations` to read from "r", "g", "b" columns
- `default_alpha`: Can be a single Real, an array (one per annotation set), or `:in_annotations` to read from "a" column
- `default_shape`: Can be a single Symbol, an array (one per annotation set), or `:in_annotations` to read from "shape" column
- `frame_col`, `x_col`, `y_col`: column names in CSV/DataFrame
- `clean_tmp`: remove temporary frames after encoding

Notes
- Annotation drawing is done in Julia by directly modifying pixel values.
- After frame images are written, the function calls `pic2vid` (which uses ffmpeg) only to encode the frames into a video. The overlay itself is performed in Julia.
- When using arrays for default_color/alpha/radius/shape, the i-th element applies to the i-th annotation set.
- When using `:in_annotations`, each row in the DataFrame/CSV can have its own properties.
"""
function overlay_annotations_on_video(annotations, video_path::AbstractString, output_path::AbstractString;
    tmpdir::AbstractString = mktempdir(), fps=nothing, radius=OVERLAY_RADIUS, default_color="red@0.5",
    default_alpha=OVERLAY_DEFAULT_ALPHA, default_shape::Union{Symbol,AbstractVector}=:circle, 
    frame_col::Symbol=:frame, x_col::Symbol=:px, y_col::Symbol=:py,
    clean_tmp::Bool=true, flag_dryrun::Bool=false, mode::Symbol = :stream, encoder_options=(crf=23, preset="ultrafast"), max_frames::Union{Nothing,Int}=nothing,
    kwargs...)

    # Accept a single CSV path or DataFrame directly for convenience
    if annotations isa AbstractString || annotations isa DataFrame || annotations isa Dict || annotations isa NamedTuple
        annotations = [annotations]
    elseif !(annotations isa AbstractVector)
        throw(ArgumentError("annotations must be a Vector of items, a single CSV path, or a DataFrame"))
    end

    # Determine if we're reading properties from annotations
    use_in_annotations_color = (default_color === :in_annotations)
    use_in_annotations_alpha = (default_alpha === :in_annotations)
    use_in_annotations_radius = (radius === :in_annotations)
    use_in_annotations_shape = (default_shape === :in_annotations)

    # Convert scalar defaults to arrays for uniform handling
    num_annotations = length(annotations)
    
    # Handle color array
    if use_in_annotations_color
        default_colors = fill("red@0.5", num_annotations)  # placeholder, will be overridden per-row
    elseif default_color isa AbstractVector
        default_colors = default_color
        if length(default_colors) != num_annotations
            throw(ArgumentError("Length of default_color array ($(length(default_colors))) must match number of annotations ($num_annotations)"))
        end
    else
        default_colors = fill(default_color, num_annotations)
    end

    # Handle alpha array
    if use_in_annotations_alpha
        default_alphas = fill(OVERLAY_DEFAULT_ALPHA, num_annotations)  # placeholder
    elseif default_alpha isa AbstractVector
        default_alphas = default_alpha
        if length(default_alphas) != num_annotations
            throw(ArgumentError("Length of default_alpha array ($(length(default_alphas))) must match number of annotations ($num_annotations)"))
        end
    else
        default_alphas = fill(default_alpha, num_annotations)
    end

    # Handle radius array
    if use_in_annotations_radius
        default_radii = fill(25, num_annotations)  # placeholder
    elseif radius isa AbstractVector
        default_radii = radius
        if length(default_radii) != num_annotations
            throw(ArgumentError("Length of radius array ($(length(default_radii))) must match number of annotations ($num_annotations)"))
        end
    else
        default_radii = fill(radius, num_annotations)
    end

    # Handle shape array
    if use_in_annotations_shape
        default_shapes = fill(:circle, num_annotations)  # placeholder
    elseif default_shape isa AbstractVector
        default_shapes = default_shape
        if length(default_shapes) != num_annotations
            throw(ArgumentError("Length of default_shape array ($(length(default_shapes))) must match number of annotations ($num_annotations)"))
        end
    else
        default_shapes = fill(default_shape, num_annotations)
    end

    # helper: parse a simple color spec like "yellow@0.8" or "#rrggbb@a" or a 3-tuple
    color_map = Dict(
        "red"=> (1.0,0.0,0.0), "green"=> (0.0,1.0,0.0), "blue"=> (0.0,0.0,1.0),
        "yellow"=> (1.0,1.0,0.0), "white"=> (1.0,1.0,1.0), "black"=> (0.0,0.0,0.0),
        "cyan"=> (0.0,1.0,1.0), "magenta"=> (1.0,0.0,1.0), "orange"=> (1.0,0.5,0.0)
    )

    parse_color(s) = begin
        if s isa Tuple || s isa Vector && length(s) >= 3
            r,g,b = float(s[1]), float(s[2]), float(s[3])
            a = length(s) >= 4 ? float(s[4]) : default_alpha
            return (r,g,b,a)
        elseif s isa AbstractString
            parts = split(s, '@')
            col = parts[1]
            a = length(parts) == 2 ? parse(Float64, parts[2]) : default_alpha
            if startswith(col, '#') && length(col) in (4,7)
                # parse #rgb or #rrggbb
                hex = col
                if length(hex) == 4
                    r = parse(Int, repeat(string(hex[2]),2); base=16) / 255
                    g = parse(Int, repeat(string(hex[3]),2); base=16) / 255
                    b = parse(Int, repeat(string(hex[4]),2); base=16) / 255
                else
                    r = parse(Int, hex[2:3]; base=16) / 255
                    g = parse(Int, hex[4:5]; base=16) / 255
                    b = parse(Int, hex[6:7]; base=16) / 255
                end
                return (r,g,b,a)
            elseif haskey(color_map, lowercase(col))
                r,g,b = color_map[lowercase(col)]
                return (r,g,b,a)
            else
                # fallback: try parse numbers separated by commas
                parts2 = split(col, ',')
                if length(parts2) >= 3
                    r = parse(Float64, parts2[1]); g = parse(Float64, parts2[2]); b = parse(Float64, parts2[3])
                    return (r,g,b,a)
                else
                    # default red
                    return (1.0,0.0,0.0,a)
                end
            end
        else
            return (1.0,0.0,0.0,default_alpha)
        end
    end

    # normalize an annotation item into per-row properties or uniform properties
    # Returns: (frames_vec, Nx2 Int matrix of pixels, per_row_properties::Bool, color_or_colors, alpha_or_alphas, radius_or_radii, shape_or_shapes)
    function normalize_item(item, item_idx)
        df = nothing
        color = default_colors[item_idx]
        alpha = default_alphas[item_idx]
        rad = default_radii[item_idx]
        shape = default_shapes[item_idx]
        
        if item isa AbstractString
            # CSV path
            df = CSV.read(item, DataFrame)
        elseif item isa DataFrame
            df = item
        elseif item isa Dict || item isa NamedTuple
            if haskey(item, :csv) || haskey(item, :"csv")
                p = get(item, :csv, get(item, "csv", nothing))
                df = CSV.read(p, DataFrame)
            elseif haskey(item, :df) || haskey(item, :"df")
                df = get(item, :df, get(item, "df", nothing))
            end
            # Item-specific overrides (backward compatibility)
            color = get(item, :color, get(item, "color", color))
            rad = get(item, :radius, get(item, "radius", rad))
            alpha = get(item, :alpha, get(item, "alpha", alpha))
            shape = get(item, :shape, get(item, "shape", shape))
        else
            throw(ArgumentError("Unsupported annotation item type: $(typeof(item))"))
        end

        if df === nothing
            throw(ArgumentError("Annotation contains no dataframe or csv path"))
        end

        # Extract columns
        if !(frame_col in propertynames(df) || haskey(df, frame_col))
            # try symbol/string variations
        end
        frames = round.(Int, df[!, frame_col])
        xs = round.(Int, df[!, x_col])
        ys = round.(Int, df[!, y_col])
        pts = hcat(xs, ys)
        
        # Check if we should read per-row properties from DataFrame
        per_row_properties = false
        colors_per_row = nothing
        alphas_per_row = nothing
        radii_per_row = nothing
        shapes_per_row = nothing
        
        if use_in_annotations_color && (:r in propertynames(df) || :g in propertynames(df) || :b in propertynames(df))
            per_row_properties = true
            # Read r, g, b columns (and optionally 'a' if available)
            r_vals = hasproperty(df, :r) ? df[!, :r] : fill(1.0, nrow(df))
            g_vals = hasproperty(df, :g) ? df[!, :g] : fill(0.0, nrow(df))
            b_vals = hasproperty(df, :b) ? df[!, :b] : fill(0.0, nrow(df))
            a_vals = hasproperty(df, :a) ? df[!, :a] : fill(alpha, nrow(df))
            colors_per_row = [(Float64(r_vals[i]), Float64(g_vals[i]), Float64(b_vals[i]), Float64(a_vals[i])) for i in 1:nrow(df)]
        end
        
        if use_in_annotations_alpha && hasproperty(df, :a)
            per_row_properties = true
            alphas_per_row = Float64.(df[!, :a])
        end
        
        if use_in_annotations_radius && hasproperty(df, :radius)
            per_row_properties = true
            radii_per_row = Int.(round.(df[!, :radius]))
        end
        
        if use_in_annotations_shape && hasproperty(df, :shape)
            per_row_properties = true
            shapes_per_row = Symbol.(df[!, :shape])
        end
        
        # If using per-row properties, return them; otherwise return uniform properties
        if per_row_properties
            # Fill in any missing per-row arrays with uniform values
            if colors_per_row === nothing
                colt = parse_color(color)
                colors_per_row = fill(colt, nrow(df))
            end
            if alphas_per_row === nothing
                alphas_per_row = fill(alpha, nrow(df))
            end
            if radii_per_row === nothing
                radii_per_row = fill(rad, nrow(df))
            end
            if shapes_per_row === nothing
                shapes_per_row = fill(shape, nrow(df))
            end
            return (frames, pts, true, colors_per_row, alphas_per_row, radii_per_row, shapes_per_row)
        else
            # Uniform properties for all rows
            colt = parse_color(color)
            return (frames, pts, false, colt, alpha, rad, shape)
        end
    end

    # Build normalized list
    normalized = Vector{Any}(undef, length(annotations))
    for (i,item) in enumerate(annotations)
        normalized[i] = normalize_item(item, i)
    end

    
    flag_dryrun && (return)

    # Open video and read frames
    vid = VideoIO.openvideo(video_path)
    fps_vid = isnothing(fps) || fps === missing ? get_fps(video_path) : fps
    # Ensure an integer fps for encoders that expect it
    fps_vid_int = try
        Int(round(Float64(fps_vid)))
    catch
        25
    end
    num_frames = get_number_frames(video_path)
    if num_frames === missing || num_frames === nothing
        # fallback: estimate from duration
        dur = get_duration(video_path)
        num_frames = Int(ceil(dur * Float64(fps_vid)))
    end

    # Respect optional testing limit
    if !(max_frames === nothing)
        num_frames = min(num_frames, max_frames)
    end

    # read frames in sequence and draw
    @info "Processing $num_frames frames (fps=$fps_vid)"
    img = nothing
    try
        img = read(vid)
    catch err
        @warn "Failed to read first frame, retrying: $err"
        seekstart(vid)
        img = read(vid)
    end

    # Helper: clamp
    clamp1(a, lo, hi) = max(lo, min(hi, a))

    h_img, w_img = size(img,1), size(img,2)

    

    function blend_pixel!(img, x, y, overlay_rgba)
        # x,y are integer pixel coordinates (1-based) where indexing is img[y,x]
        if x < 1 || x > w_img || y < 1 || y > h_img
            return
        end
        base = img[y,x]
        # get base r,g,b in 0..1 float
        base_r = float(red(base)); base_g = float(green(base)); base_b = float(blue(base))
        or_, og, ob, oa = overlay_rgba
        a = oa
        nr = a*or_ + (1-a)*base_r
        ng = a*og + (1-a)*base_g
        nb = a*ob + (1-a)*base_b
        img[y,x] = typeof(img[1])(nr, ng, nb)
    end

    function draw_circle!(img, cx, cy, r, overlay_rgba)
        x0 = Int(clamp1(round(cx - r), 1, w_img))
        x1 = Int(clamp1(round(cx + r), 1, w_img))
        y0 = Int(clamp1(round(cy - r), 1, h_img))
        y1 = Int(clamp1(round(cy + r), 1, h_img))
        rr = r*r
        for yy in y0:y1, xx in x0:x1
            dx = xx - cx; dy = yy - cy
            if dx*dx + dy*dy <= rr
                blend_pixel!(img, xx, yy, overlay_rgba)
            end
        end
    end

    function draw_rect!(img, cx, cy, wrect, overlay_rgba)
        half = wrect/2
        x0 = Int(clamp1(round(cx-half), 1, w_img))
        x1 = Int(clamp1(round(cx+half), 1, w_img))
        y0 = Int(clamp1(round(cy-half), 1, h_img))
        y1 = Int(clamp1(round(cy+half), 1, h_img))
        for yy in y0:y1, xx in x0:x1
            blend_pixel!(img, xx, yy, overlay_rgba)
        end
    end

    # If the user requests explicit VideoIO mode, use the `open_video_out` helper
    # which mirrors the usage in `src/test_make_video.jl`. This attempts to open
    # a VideoIO writer and write frames via that writer. If it fails we fall
    # through to the existing stream/frame logic below.
    # encoder_options = (preset="medium", crf=23)
    if mode == :VideoIO
        @info "Attempting VideoIO mode via open_video_out for $output_path (fps=$fps_vid)"
        try
            open_video_out(output_path, img, framerate=fps_vid, encoder_options=encoder_options) do writer #codec_name = "h264_nvenc",
                seekstart(vid)
                @showprogress "Encoding video frames..." for frame_idx in 0:num_frames-1
                    try
                        read!(vid, img)
                    catch err
                        @warn "Failed to read frame $frame_idx: $err -- filling blank frame"
                        img = zeros(eltype(img), size(img))
                    end

                    # iterate all annotation sets
                    for norm_item in normalized
                        frames, pts, per_row, color_data, alpha_data, radius_data, shape_data = norm_item
                        inds = findall(==(frame_idx), frames)
                        if isempty(inds)
                            continue
                        end
                        @debug frame_idx
                        @debug inds
                        @debug color_data[inds]
                        @debug per_row
                        for i in inds
                            x = pts[i,1]; y = pts[i,2]
                            if per_row
                                overlay_rgba = color_data[i]
                                rad = radius_data[i]
                                shp = shape_data[i]
                            else
                                overlay_rgba = color_data
                                rad = radius_data
                                shp = shape_data
                            end
                            @debug rad
                            @debug overlay_rgba
                            if shp == :circle || shp == :dot
                                draw_circle!(img, x, y, rad, overlay_rgba)
                            else
                                draw_rect!(img, x, y, rad, overlay_rgba)
                            end
                        end
                    end

                    # write via VideoIO writer
                    write(writer, img)
                end
            end
            close(vid)
            @info "Wrote VideoIO-mode video to $output_path"
            return output_path
        catch err
            @warn "VideoIO mode failed: $err -- falling back to other modes"
            # fall through to existing logic (stream / frames)
        end
    end

    # main loop
    seekstart(vid)
    # If streaming mode, open a video writer. If VideoIO writer fails, fall back
    # to streaming via an ffmpeg stdin pipe (no intermediate frames on disk).
    writer = nothing
    writer_type = nothing  # :videoio or :ffmpegpipe

    # helper: convert image to rgb24 bytes in row-major order (ffmpeg expects rows)
    function frame_to_rgb24_bytes(img)
        h = size(img,1); w = size(img,2)
        ch = channelview(convert.(RGB{N0f8}, img))  # 3 x h x w array of N0f8
        buf = Vector{UInt8}(undef, w*h*3)
        idx = 1
        for y in 1:h
            for x in 1:w
                r = ch[1,y,x]; g = ch[2,y,x]; b = ch[3,y,x]
                buf[idx] = UInt8(clamp(round(Int, 255*float(r)), 0, 255)); idx += 1
                buf[idx] = UInt8(clamp(round(Int, 255*float(g)), 0, 255)); idx += 1
                buf[idx] = UInt8(clamp(round(Int, 255*float(b)), 0, 255)); idx += 1
            end
        end
        return buf
    end

    if mode == :stream
        @info "Attempting to open VideoIO writer for $output_path (fps=$fps_vid)"
        try
            writer = VideoIO.openvideo(output_path; framerate=fps_vid_int)
            writer_type = :videoio
        catch err
            @warn "VideoIO.openvideo failed: $err -- attempting ffmpeg stdin pipe writer"
            # Build ffmpeg command to accept raw rgb24 stdin
            sz = string(w_img)*"x"*string(h_img)
            cmd = `$(FFMPEG.ffmpeg) -y -f rawvideo -pix_fmt rgb24 -s $sz -r $fps_vid_int -i - -c:v libx264 -pix_fmt yuv420p $output_path`
            try
                writer = open(cmd, "w")
                writer_type = :ffmpegpipe
            catch err2
                @warn "ffmpeg stdin pipe failed: $err2 -- falling back to frame files mode"
                writer = nothing
                writer_type = nothing
                mode = :frames
            end
        end
    end
    
    # prepare tmpdir
    mkpath(tmpdir)
    @info "Created temporary directory: $tmpdir"
    @info "Writing temporary frames to $tmpdir"
    @showprogress "Encoding video frames..." for frame_idx in 0:num_frames-1
        try
            read!(vid, img)
        catch err
            @warn "Failed to read frame $frame_idx: $err -- filling blank frame"
            img = zeros(eltype(img), size(img))
        end

        # iterate all annotation sets
        for norm_item in normalized
            frames, pts, per_row, color_data, alpha_data, radius_data, shape_data = norm_item
            inds = findall(==(frame_idx), frames)
            if isempty(inds)
                continue
            end
            for i in inds
                x = pts[i,1]; y = pts[i,2]
                if per_row
                    overlay_rgba = color_data[i]
                    rad = radius_data[i]
                    shp = shape_data[i]
                else
                    overlay_rgba = color_data
                    rad = radius_data
                    shp = shape_data
                end
                if shp == :circle || shp == :dot
                    draw_circle!(img, x, y, rad, overlay_rgba)
                else
                    draw_rect!(img, x, y, rad, overlay_rgba)
                end
            end
        end

        if mode == :stream && writer !== nothing
            try
                if writer_type == :videoio
                    write(writer, img)
                elseif writer_type == :ffmpegpipe
                    buf = frame_to_rgb24_bytes(img)
                    write(writer, buf)
                else
                    error("Unknown writer_type: $writer_type")
                end
            catch err
                @error "Failed to write frame $frame_idx to video writer: $err"
                rethrow(err)
            end
        else
            # Write png frame named as %d.png (matches pic2vid default counter)
                fname = joinpath(tmpdir, @sprintf("%d.png", frame_idx))
            try
                save(fname, img)
            catch err
                @error "Failed to save frame $frame_idx -> $fname" err
                rethrow(err)
            end
        end
    end

    close(vid)

    # If streaming writer was used, close it and return
    if mode == :stream && writer !== nothing
        try
            close(writer)
        catch err
            @warn "Failed to close video writer cleanly: $err"
        end
        # If we used an ffmpeg stdin pipe, wait for the process to finish so the
        # output file is fully flushed and available to callers. `open(cmd, "w")`
        # returns a `Base.Process` which we should wait on.
        if writer_type == :ffmpegpipe
            try
                wait(writer)
            catch err
                @warn "Waiting on ffmpeg process failed: $err"
            end
        end
        @info "Wrote streamed video to $output_path"
        return output_path
    end

    # encode with helper pic2vid (uses ffmpeg under the hood to encode, but overlay done)
    counter = "%d.png"
    pic2vid(tmpdir, output_path; counter=counter, fps=fps_vid_int, auto_mode=false)

    if flag_dryrun
        @info "Dryrun: frames written to $tmpdir (not encoding)"
        return tmpdir
    end

    if clean_tmp
        try
            rm(tmpdir; force=true, recursive=true)
        catch err
            @warn "Failed to remove temporary frames: $err"
        end
    end

    return output_path
end

# """
# # Example: generate one frame per CSV row with overlay drawbox
# using CSV, DataFrames, FileIO, ImageIO, Images, ColorTypes, Colors, Printf
# """
# function generate_frames_from_csv(csv_path::AbstractString, base_image_path::AbstractString, out_dir::AbstractString)
#     df = CSV.File(csv_path) |> DataFrame
#     img = load(base_image_path)                           # load base image (height x width x channels)
#     mkpath(out_dir)

#     # helper to draw rectangle border (works on Images arrays)
#     function draw_rect!(img, x::Int, y::Int, w::Int, h::Int, color, thickness::Int=3)
#         h_img, w_img = size(img, 1), size(img, 2)
#         x1 = clamp(x, 1, w_img); y1 = clamp(y, 1, h_img)
#         x2 = clamp(x + w - 1, 1, w_img); y2 = clamp(y + h - 1, 1, h_img)
#         for t in 0:thickness-1
#             top_y = clamp(y1 + t, 1, h_img);       bottom_y = clamp(y2 - t, 1, h_img)
#             left_x = clamp(x1 + t, 1, w_img);      right_x = clamp(x2 - t, 1, w_img)
#             # top and bottom horizontal lines
#             img[top_y, left_x:right_x] .= color
#             img[bottom_y, left_x:right_x] .= color
#             # left and right vertical lines
#             img[y1:y2, left_x] .= color
#             img[y1:y2, right_x] .= color
#         end
#     end

#     for (i, row) in enumerate(eachrow(df))
#         frame = copy(img)
#         # adjust these to your CSV column names:
#         x = Int(round(row.x))           # left coordinate (1-based)
#         y = Int(round(row.y))           # top coordinate (1-based)
#         w = Int(round(row.w))           # width in pixels
#         h = Int(round(row.h))           # height in pixels
#         time_seconds = hasproperty(row, :time_seconds) ? row.time_seconds : (hasproperty(row, :time) ? row.time : i)

#         color = RGBA{N0f8}(1, 0, 0, 1)  # red box
#         draw_rect!(frame, x, y, w, h, color, 3)

#         fname = joinpath(out_dir, @sprintf("frame%04d_time-%.3f.png", i, Float64(time_seconds)))
#         save(fname, frame)
#     end
# end


# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=ow=1920:oh=ih*ow/iw:mode=line:rate=25,format=yuv420p[bottom]" -map "[top]" -map "[bottom]" -y output.mp4
# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=1920x480:mode=line:rate=25,format=yuv420p[bottom]" -map "[top]" -map "[bottom]" -y output.mp4

# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=ow=1920:oh=ih*ow/iw:mode=line:rate=25[bottom_audio]" -map "[top]" -map "[bottom_audio]" -c:a aac -strict experimental -y output.mp4

# input_resolution=$(ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of csv=s=x:p=0 input1.mp4)
# ffmpeg -i input1.mp4 -i input2.mp4 -i audio.ogg -filter_complex "[0:v][1:v]vstack=inputs=2[top];[2:a]showwaves=s=ow=${input_resolution%:*}:oh=${input_resolution#*:}/iw:mode=line:rate=25[bottom_audio]" -map "[top]" -map "[bottom_audio]" -c:a aac -strict experimental -y output.mp4


# "[2:a]showwaves=s=1920x400:mode=cline:colors=blue:rate=25[waves];[0:v][1:v][waves]vstack=inputs=3[v];[2:a]loudnorm[a];"