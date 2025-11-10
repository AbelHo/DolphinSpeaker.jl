"""
cam_calib.jl

Provide a small utility to calibrate pixel -> annotated positions using
paired CSVs produced by the pipeline. The primary function is

    calibrate_from_annotations(res_dir; vidfname=nothing, detection_pixels="detection_pixels.csv",
                                annotated_csv=nothing, annotated_frame_col=:frame_index,
                                out_prefix="combined__1.GoPro_Clicker.MP4_overlaidIMG",
                                write_csv=true, show_plots=true)

It returns a NamedTuple with fitted models, fitted closures, corrected DataFrame,
and basic statistics.

The script is lightweight and intentionally avoids heavy dependencies beyond
what the project already uses (CSV, DataFrames, GLM, RobustModels, Plots).

Usage as a script (example):
    julia src/cam_calib.jl /path/to/res_dir --annotated=annotated.csv

"""

# module CamCalib

using CSV
using DataFrames
using GLM
using RobustModels
using Statistics
using Printf
using Dates
using Plots

export calibrate_from_annotations

"""
    calibrate_from_annotations(res_dir; vidfname=nothing, detection_pixels="detection_pixels.csv",
                               annotated_csv=nothing, annotated_frame_col=:frame_index,
                               out_prefix="combined__1.GoPro_Clicker.MP4_overlaidIMG",
                               write_csv=true, show_plots=true)

Read detection pixels and annotated CSV, join on frame index, fit linear and robust
regressions mapping pixel x->world x and pixel y->world y, compute corrected
coordinates, optionally save corrected CSV and show simple diagnostic plots.

Returns a NamedTuple:
  (model_x, model_x_r, model_y, model_y_r, fitted_func_x, fitted_func_y, corrected_df, stats)

Notes:
- Expects the detection CSV to contain columns `frame` (or `frame_index`) and `px`, `py`.
- Expects the annotated CSV to contain columns `frame` (or `frame_index`) and `x`, `y`.

Example:
    res = calibrate_from_annotations("/path/to/res_dir"; annotated_csv="annotated.csv", vidfname="video.mp4")
    println("Calibration stats: ", res.stats)


annotated_csv = "/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1/combined__1.GoPro_Clicker.MP4_overlaidIMG_20251018_122735.csv"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1"
vidfname = joinpath(res_dir, "combined__1.GoPro_Clicker.MP4.mp4")
a = calibrate_from_annotations(res_dir, vidfname, annotated_csv)

annotated_csv = "/media/spin/anas2/data_res/dolphin/calf/calibration/D3/clicker/Calibration_05_11_2025/combined__GH011401.MP4_overlaidIMG_20251109_235724.csv"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/calibration/D3/clicker/Calibration_05_11_2025"
vidfname = "/media/spin/anas2/data_res/dolphin/calf/calibration/D3/clicker/Calibration_05_11_2025/combined__GH011401.MP4_overlaid.mkv_DYnormalized-audio.mp4" #"/media/spin/anas2/data_res/dolphin/calf/calibration/D3/clicker/Calibration_05_11_2025/combined__GH011401.MP4_overlaidIMG.mp4"
b = calibrate_from_annotations(res_dir, vidfname, annotated_csv)
"""
function calibrate_from_annotations(
    res_dir::AbstractString, vidfname::Union{Nothing,String}=nothing, annotated_csv::Union{Nothing,String}=nothing;
    detection_pixels::AbstractString="detection_pixels.csv",
    annotated_frame_col::Union{Symbol,String}=:frame_index,
    write_csv::Bool=true,
    show_plots::Bool=true,
    flag_save_csv::Bool=false)

    # build paths
    dp_path = joinpath(res_dir, detection_pixels)
    if !isfile(dp_path)
        error("detection pixels file not found: $dp_path")
    end

    # read detection pixels
    dp = CSV.read(dp_path, DataFrame)

    dur = get_duration(vidfname)
    fps = get_fps(vidfname)
    dp = dp[ 0 .<= dp.frame .<= dur*fps, :]
    dp.frame_index = 0:(length(dp.frame)-1)
  

    # load annotated CSV
    # annotated_csv = annotated_csv |> x -> isabspath(x) ? x : joinpath(res_dir, annotated_csv)
    d_annotated = CSV.read(annotated_csv, DataFrame)

    rename!(d_annotated, Dict("frame" => "frame_index"))
    # remove redundant columns if exist
    col_remove = intersect(names(d_annotated) .|> Symbol, [:r, :g, :b, :a, :radius])
    select!(d_annotated, Not(col_remove))  # remove extra columns if exist
    # dropmissing!(d_annotated, [:frame_index, :x, :y])  # drop rows with missing required columns
    if "tag" ∈ names(d_annotated)
        @info "Removing annotated outliers tagged as 'o' in the tag column"
        filter!(row -> ismissing(row.tag) || row.tag != "o", d_annotated)
    end
    # # normalize annotated frame column name to :frame_index if it uses `frame`
    # if annotated_frame_col in names(d_annotated) && annotated_frame_col != :frame_index
    #     rename!(d_annotated, Dict(annotated_frame_col => :frame_index))
    # elseif :frame in names(d_annotated) && :frame_index ∉ names(d_annotated)
    #     rename!(d_annotated, Dict(:frame => :frame_index))
    # end

    # join and clean
    dd = leftjoin(dp, d_annotated; on=:frame_index)
    sort!(dd, :frame_index)
    # remove rows with missing required columns
    required = [:px, :py, :x, :y]
    # miss = setdiff(required, names(dd))
    # if !isempty(miss)
    #     error("Missing required columns after join: $(miss)")
    # end
    dropmissing!(dd, required)

    # Fit linear and robust models for x<-px and y<-py
    form_x = @formula(x ~ px)
    form_y = @formula(y ~ py)
    model_x = fit(LinearModel, form_x, dd)
    model_x_r = rlm(form_x, dd, MEstimator{TukeyLoss}(); initial_scale=:L1, ridgeλ=1.0)
    model_y = fit(LinearModel, form_y, dd)
    model_y_r = model_y
    try
        model_y_r = rlm(form_y, dd, MEstimator{TukeyLoss}(); initial_scale=:L1, ridgeλ=1.0)
    catch err
        @warn "Robust model fitting for y failed: $err. Falling back to linear model."
    end
    # build fitted closures (robust fits preferred)
    cx = coef(model_x_r)
    cy = coef(model_y_r)
    fitted_func_x(px) = cx[1] .+ cx[2] .* px
    fitted_func_y(py) = cy[1] .+ cy[2] .* py

    # compute corrected coordinates
    dd[!, :x_new] = fitted_func_x.(dd.px)
    dd[!, :y_new] = fitted_func_y.(dd.py)

    # diagnostics
    resid = sqrt.((dd.x .- dd.x_new).^2 .+ (dd.y .- dd.y_new).^2)
    stats = Dict(
        :n_points => nrow(dd),
        :rms => sqrt(mean(resid .^2)),
        :median_err => median(resid),
        :mean_err => mean(resid)
    )

    if show_plots
        try
            vid_info = get_media_info(vidfname)
            vid_height = get(vid_info["streams"][1], "height", "")
            vid_width = get(vid_info["streams"][1], "width", "")

            p1 = Plots.scatter(dd.px, dd.x; alpha=0.4, label="x vs px", xlabel="Pixel px", ylabel="Annotated x")
            plot!(p1, dd.px, predict(model_x); label="linear fit")
            if model_x_r !== model_x
                plot!(p1, dd.px, fitted_func_x.(dd.px); label="robust fit")
            end
            p2 = Plots.scatter(dd.py, dd.y; alpha=0.4, label="y vs py", xlabel="Pixel py", ylabel="Annotated y")
            plot!(p2, dd.py, predict(model_y); label="linear fit")
            if model_y_r !== model_y
                plot!(p2, dd.py, fitted_func_y.(dd.py); label="robust fit")
            end
            display(plot(p1, p2; layout=(1,2), size=(1000,400)))

            # plot required correction vectors
            ## naive all
            Plots.quiver(dd.px, dd.py, quiver=(dd.x .- dd.px, dd.y .- dd.py); aspect_ratio=:equal, color=:black, label="error vectors",
                xlabel="Pixel px", ylabel="Pixel py", title="Annotated vs Estimated Positions",
                xlims=(0, vid_width), ylims=(0, vid_height)) |> display
            # ## remove high errrors as outlier
            dist = [dd.x'; dd.y'] - [dd.px'; dd.py'] |> eachcol .|> norm
            threshold = quantile(dist, 0.9)
            @info size(dist), threshold
            inds = findall(dist .<= threshold)
            Plots.quiver(dd.px[inds], dd.py[inds], quiver=(dd.x[inds] .- dd.px[inds], dd.y[inds] .- dd.py[inds]);
                aspect_ratio=:equal, color=:blue, label="error vectors (90% quantile)",
                xlabel="Pixel px", ylabel="Pixel py", title="Annotated vs Estimated Positions (90% quantile)",
                xlims=(0, vid_width), ylims=(0, vid_height)) |> display

            
            p3 = Plots.scatter(dd[:, [:px,:x,:x_new]] |> Matrix, dd[:,[:py,:y,:y_new]]|>Matrix; aspect_ratio=:equal,
                xlabel="X (pixels)", ylabel="Y (pixels)", title="Corrected vs Annotated Positions", label=["Estimated" "Annotated" "Corrected"],
                markerstrokewidth=0)#, markerstrokecolor=:transparent)
            Plots.quiver!(p3, dd.px, dd.py, quiver=(dd.x_new .- dd.px, dd.y_new .- dd.py); aspect_ratio=:equal, color=:black, label="correction vectors")
            # Plots.plot!(p3, (dd[:, [:px,:x]] |> Matrix)', (dd[:,[:py,:y]]|>Matrix)'; aspect_ratio=:equal, color=:black)
            # restrict to video height and width only
            xlims!(p3, 0, vid_width)
            ylims!(p3, 0, vid_height)
            display(p3)
        catch err
            @warn "Plotting failed: $err"
        end
    end

    # write corrected CSV if requested
    if flag_save_csv
        out_csv = joinpath(res_dir, out_prefix * "__annotated_corrected.csv")
        if write_csv
            try
                CSV.write(out_csv, dd)
                @info "Wrote corrected annotations to $out_csv"
            catch err
                @warn "Failed writing corrected CSV: $err"
            end
        end
    end

    # Prepare a minimal dd_new expected by overlay helper used in the repository
    dd_new = select(dd, :frame_index, :x_new, :y_new)
    rename!(dd_new, Dict(:frame_index => :frame, :x_new => :x, :y_new => :y))

    # Try to call overlay_annotations_on_video if available (optional)
    # overlay_out = nothing
    # try
    #     if isdefined(Main, :overlay_annotations_on_video)
    #         video_file = vidfname === nothing ? joinpath(res_dir, out_prefix * ".mp4") : vidfname
    #         overlay_out = try
    #             overlay_annotations_on_video(dd_new, video_file, replace(video_file, r"\.mp4$" => "_annotated_corrected.mkv"); x_col=:x, y_col=:y, default_color="green@0.4", default_radius=10, mode=:stream)
    #         catch e
    #             @warn "overlay_annotations_on_video failed: $e"
    #             nothing
    #         end
    #     end
    # catch _
    #     # ignore overlay errors
    # end

    return (;cx, cy, x_corrected=dd[!, :x_new], y_corrected=dd[!, :y_new], 
        model_x=model_x, model_x_r=model_x_r, model_y=model_y, model_y_r=model_y_r,
        fitted_func_x=fitted_func_x, fitted_func_y=fitted_func_y,
        corrected_df=dd, stats=stats)
end


# function calibrate_from_annotations_complex(res_dir::AbstractString; vidfname::Union{Nothing,String}=nothing,
#         detection_pixels::AbstractString="detection_pixels.csv",
#         annotated_csv::Union{Nothing,String}=nothing,
#         annotated_frame_col::Union{Symbol,String}=:frame_index,
#         out_prefix::AbstractString="combined__1.GoPro_Clicker.MP4_overlaidIMG",
#         write_csv::Bool=true,
#         show_plots::Bool=true)

#     # build paths
#     dp_path = joinpath(res_dir, detection_pixels)
#     if !isfile(dp_path)
#         error("detection pixels file not found: $dp_path")
#     end

#     # read detection pixels
#     dp = CSV.read(dp_path, DataFrame)

#     # If df contains a `frame` column (sample code used `frame`), create frame_index
#     if :frame in names(dp) && :frame_index ∉ names(dp)
#         # keep frame within video duration if vidfname provided
#         if vidfname !== nothing && isfile(vidfname)
#             try
#                 dur = get_duration(vidfname)
#                 fps = get_fps(vidfname)
#                 # filter same as original code: 0 <= frame <= dur*fps
#                 dp = dp[ 0 .<= dp.frame .<= dur*fps, :]
#             catch err
#                 @warn "Could not get video duration/fps for $vidfname: $err"
#             end
#         end
#         dp.frame_index = 0:(length(dp.frame)-1)
#     end

#     # load annotated CSV
#     if annotated_csv === nothing
#         # try common name used in pipeline
#         ann_path = joinpath(res_dir, out_prefix * "_20251018_122735.csv")
#         if isfile(ann_path)
#             annotated_csv = ann_path
#         else
#             # try a generic overlaid csv
#             candidate = joinpath(res_dir, out_prefix * ".csv")
#             annotated_csv = isfile(candidate) ? candidate : error("annotated CSV not provided and not found in $res_dir")
#         end
#     else
#         annotated_csv = annotated_csv |> x -> isabspath(x) ? x : joinpath(res_dir, annotated_csv)
#     end

#     d_annotated = CSV.read(annotated_csv, DataFrame)

#     # normalize annotated frame column name to :frame_index if it uses `frame`
#     if annotated_frame_col in names(d_annotated) && annotated_frame_col != :frame_index
#         rename!(d_annotated, Dict(annotated_frame_col => :frame_index))
#     elseif :frame in names(d_annotated) && :frame_index ∉ names(d_annotated)
#         rename!(d_annotated, Dict(:frame => :frame_index))
#     end

#     # join and clean
#     dd = leftjoin(dp, d_annotated, on=:frame_index)
#     sort!(dd, :frame_index)
#     # remove rows with missing required columns
#     required = [:px, :py, :x, :y]
#     miss = setdiff(required, names(dd))
#     if !isempty(miss)
#         error("Missing required columns after join: $(miss)")
#     end
#     dropmissing!(dd, required)

#     # Fit linear and robust models for x<-px and y<-py
#     form_x = @formula(x ~ px)
#     form_y = @formula(y ~ py)
#     model_x = fit(LinearModel, form_x, dd)
#     model_x_r = rlm(form_x, dd, MEstimator{TukeyLoss}(); initial_scale=:L1, ridgeλ=1.0)
#     model_y = fit(LinearModel, form_y, dd)
#     model_y_r = rlm(form_y, dd, MEstimator{TukeyLoss}(); initial_scale=:L1, ridgeλ=1.0)

#     # build fitted closures (robust fits preferred)
#     cx = coef(model_x_r)
#     cy = coef(model_y_r)
#     fitted_func_x(px) = cx[1] .+ cx[2] .* px
#     fitted_func_y(py) = cy[1] .+ cy[2] .* py

#     # compute corrected coordinates
#     dd[!, :x_new] = fitted_func_x.(dd.px)
#     dd[!, :y_new] = fitted_func_y.(dd.py)

#     # diagnostics
#     resid = sqrt.((dd.x .- dd.x_new).^2 .+ (dd.y .- dd.y_new).^2)
#     stats = Dict(
#         :n_points => nrow(dd),
#         :rms => sqrt(mean(resid .^2)),
#         :median_err => median(resid),
#         :mean_err => mean(resid)
#     )

#     if show_plots
#         try
#             p1 = Plots.scatter(dd.px, dd.x; alpha=0.4, label="x vs px", xlabel="Pixel px", ylabel="Annotated x")
#             plot!(p1, dd.px, fitted_func_x.(dd.px); label="robust fit")
#             p2 = Plots.scatter(dd.py, dd.y; alpha=0.4, label="y vs py", xlabel="Pixel py", ylabel="Annotated y")
#             plot!(p2, dd.py, fitted_func_y.(dd.py); label="robust fit")
#             display(plot(p1, p2; layout=(1,2), size=(1000,400)))
#         catch err
#             @warn "Plotting failed: $err"
#         end
#     end

#     # write corrected CSV if requested
#     out_csv = joinpath(res_dir, out_prefix * "__annotated_corrected.csv")
#     if write_csv
#         try
#             CSV.write(out_csv, dd)
#             @info "Wrote corrected annotations to $out_csv"
#         catch err
#             @warn "Failed writing corrected CSV: $err"
#         end
#     end

#     # Prepare a minimal dd_new expected by overlay helper used in the repository
#     dd_new = select(dd, :frame_index, :x_new, :y_new)
#     rename!(dd_new, Dict(:frame_index => :frame, :x_new => :x, :y_new => :y))

#     # Try to call overlay_annotations_on_video if available (optional)
#     overlay_out = nothing
#     try
#         if isdefined(Main, :overlay_annotations_on_video)
#             video_file = vidfname === nothing ? joinpath(res_dir, out_prefix * ".mp4") : vidfname
#             overlay_out = try
#                 overlay_annotations_on_video(dd_new, video_file, replace(video_file, r"\.mp4$" => "_annotated_corrected.mkv"); x_col=:x, y_col=:y, default_color="green@0.4", default_radius=10, mode=:stream)
#             catch e
#                 @warn "overlay_annotations_on_video failed: $e"
#                 nothing
#             end
#         end
#     catch _
#         # ignore overlay errors
#     end

#     return (model_x=model_x, model_x_r=model_x_r, model_y=model_y, model_y_r=model_y_r,
#             fitted_func_x=fitted_func_x, fitted_func_y=fitted_func_y,
#             corrected_df=dd, dd_new=dd_new, stats=stats, overlay_result=overlay_out, out_csv=out_csv)
# end

########## small CLI when executed directly ##########
function _print_usage()
    println("Usage: julia src/cam_calib.jl /path/to/res_dir [--annotated=annot.csv] [--detection=detection_pixels.csv] [--vid=video.mp4]")
end

function _main(argv)
    if length(argv) < 1
        _print_usage(); return
    end
    res_dir = argv[1]
    kwargs = Dict{String,String}()
    for a in argv[2:end]
        if occursin("=", a)
            k,v = split(a, "=", limit=2)
            kwargs[k] = v
        end
    end
    annotated = get(kwargs, "--annotated", nothing)
    detection = get(kwargs, "--detection", "detection_pixels.csv")
    vid = get(kwargs, "--vid", nothing)

    res = calibrate_from_annotations(res_dir; detection_pixels=detection, annotated_csv=annotated, vidfname=vid)
    @info "Calibration finished: $(res.stats)"
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end

# end # module
