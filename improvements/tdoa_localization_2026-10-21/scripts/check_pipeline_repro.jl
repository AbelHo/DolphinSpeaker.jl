# Does the harness reproduce the pipeline's own output? Compare its minmax + NelderMead pixels with the
# Impulse CSV written by run_syncVid.jl (process_detections) for the same file.
using CSV, DataFrames, DelimitedFiles, Printf
res_dir = ARGS[1]
tag = "251006_001_0001"
pipe_dir = "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_nooverlay2/New_device"
fs = 192_000
lab = readdlm(joinpath(pipe_dir, "$(tag)___t490.376461598132_d77__cps60.0.txt"), '\t')
pind_pipe = round.(Int, lab[:, 1] .* fs) .+ 1
pix = CSV.read(joinpath(pipe_dir, "$(tag)_Impulse_t490.376461598132_d77.csv"), DataFrame)
@assert nrow(pix) == length(pind_pipe)
pc = CSV.read(joinpath(res_dir, "per_click_$tag.csv"), DataFrame)
mm = pc[pc.method.=="minmax [Impulse, overlay]", :]
row_of = Dict(p => i for (i, p) in enumerate(pind_pipe))
found = [haskey(row_of, p) for p in mm.pind]
idx = [row_of[p] for p in mm.pind[found]]
fov = [80.0, 54.46086642417594]
imsize = (1520, 2740)
px = mm.az_nm[found] ./ (fov[1] / 2) .* (imsize[2] / 2) .+ imsize[2] / 2
py = mm.el_nm[found] ./ (-fov[2] / 2) .* (imsize[1] / 2) .+ imsize[1] / 2
dpx = abs.(px .- pix.px[idx])
dpy = abs.(py .- pix.py[idx])
@printf("harness clicks found in the pipeline's click list: %d/%d (pipeline has %d)\n", count(found), nrow(mm), length(pind_pipe))
@printf("pixel difference vs pipeline Impulse CSV: max |Δpx| %.3g, max |Δpy| %.3g; within 0.01 px: %.2f%%\n",
    maximum(dpx), maximum(dpy), 100 * count((dpx .< 0.01) .& (dpy .< 0.01)) / length(dpx))
offscreen = .!((0 .<= pix.px .<= imsize[2]) .& (0 .<= pix.py .<= imsize[1]))
@printf("pipeline Impulse detections drawn outside the %dx%d frame: %.1f%%\n", imsize[2], imsize[1], 100 * count(offscreen) / nrow(pix))
