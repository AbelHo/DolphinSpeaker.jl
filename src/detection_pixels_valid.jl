# NaN-safe detection-pixel CSV for the calibrated overlay (contract T6).
#
# Included from src/run_example_2.jl (T4), inside module DolphinSpeaker, so CSV and DataFrame
# are already in scope there. This file defines only write_valid_detection_pixels.
#
# Why: overlay_annotations_on_video (src/video.jl:436-438) runs round.(Int, ...) on the frame,
# px and py columns. round(Int, NaN) throws InexactError, so one NaN row aborts the whole
# overlay. This function removes non-finite rows before the overlay reads the CSV.
#
# Coordinates: px and py come from angle2pxd (src/localization.jl:246), which maps the field of
# view to the continuous range [0, width] x [0, height]. imsize is (height, width), as elsewhere
# in the repo. A row outside that range is kept in the valid CSV, not clamped. The overlay
# clamps its drawing box to the image itself (src/video.jl:561-566).

"""
    write_valid_detection_pixels(csv_path; out_path, status_path, imsize) -> NamedTuple

Read a detection-pixel CSV (header `frame,px,py,radius,r,g,b,a,shape`) and write:

- `out_path` (default `*_valid.csv`): every row whose `frame`, `px` and `py` are finite.
  Rows with `outside_image` status are kept, not clamped.
- `status_path` (default `*_status.csv`): one row per input row with columns `row,frame,status`.
  `row` is the 1-based data-row index in `csv_path`. `status` is `ok`, `nonfinite` (any of
  frame, px, py is NaN, Inf or missing) or `outside_image` (finite, but px<0, px>width,
  py<0 or py>height, with `imsize = (height, width)`).

Returns `(; out_path, status_path, n_input, n_written, n_ok, n_outside_image, n_nonfinite)`.
"""
function write_valid_detection_pixels(csv_path::AbstractString;
        out_path::AbstractString = replace(csv_path, ".csv" => "_valid.csv"),
        status_path::AbstractString = replace(csv_path, ".csv" => "_status.csv"),
        imsize = imsize)
    df = CSV.read(csv_path, DataFrame)
    for col in ("frame", "px", "py")
        Symbol(col) in propertynames(df) || throw(ArgumentError("write_valid_detection_pixels: column '$col' missing in $csv_path"))
    end
    height, width = imsize[1], imsize[2]
    isfin(x) = !ismissing(x) && isfinite(x)

    n = nrow(df)
    status = Vector{String}(undef, n)
    for i in 1:n
        f = df.frame[i]; x = df.px[i]; y = df.py[i]
        if !(isfin(f) && isfin(x) && isfin(y))
            status[i] = "nonfinite"
        elseif x < 0 || x > width || y < 0 || y > height
            status[i] = "outside_image"
        else
            status[i] = "ok"
        end
    end

    keep = status .!= "nonfinite"
    valid_df = df[keep, :]
    CSV.write(out_path, valid_df)

    status_df = DataFrame(row = collect(1:n), frame = df.frame, status = status)
    CSV.write(status_path, status_df)

    n_ok = count(==("ok"), status)
    n_outside = count(==("outside_image"), status)
    n_nonfinite = count(==("nonfinite"), status)
    return (; out_path, status_path, n_input = n, n_written = count(keep),
              n_ok, n_outside_image = n_outside, n_nonfinite)
end
