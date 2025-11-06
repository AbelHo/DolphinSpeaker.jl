# using SignalAlignment
include("audio.jl")
include("dsp.jl")
using PlotlyJS
include("plotting.jl")
using LinearAlgebra

aufname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2/250519_002_0001.WAV"
vidfname = "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2/GH019986.MP4"

aufname = "/media/spin/anas2/data/calf/upload/Double_ball_S1/250403_002_0001.WAV"
vidfname = "/media/spin/anas2/data/calf/upload/Double_ball_S1/GH019888.MP4"


#~ get sync time between video and audio file #######################################
data, fs = readAudio(aufname)
data_v, fs_v = get_videos_audiodata(vidfname)

# @info "duration of audio file: $(size(data,1)/fs)s"
# @info "duration of video file: $(size(data_v,1)/fs_v)s"

data_down = resample(@view(data[1:fs*120,1]), fs_v/fs)
template_start = fs_v*60
template_end = fs_v*2*60
template_sig = data_down[template_start:template_end]
# max_time = fs2*60
delays, delay_conf = finddelay2(@view(data_v[:,1]), template_sig)
# @info "Delay between video and audio signal: $((template_start - delays)/fs2)s, confidence: $delay_conf"
@info "Delay of audio against video signal: $((delays - template_start)/fs_v)s, confidence: $delay_conf"
###################################################

#~ Organized split video and audio files
data_dir = "/media/spin/anas2/data/calf/upload/Double_ball_S1"
data_dir = "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/temp/delete"
occursin.( Ref(Regex(join(vidtypes, '|'))), readdir(data_dir))

vidlist = filter( x -> occursin(Regex(join(vidtypes, '|')), x|>lowercase), readdir(data_dir; join=true))
audlist = filter( x -> occursin(Regex(join(autypes, '|')), x|>lowercase), readdir(data_dir; join=true))

# combine all video files into one file
mkpath(res_dir)
temp_filelist = joinpath(res_dir, "temp_filelist.txt")
write(temp_filelist, join(["file '$v'" for v in vidlist], "\n"))

cmd = `ffmpeg -f concat -safe 0 -i $temp_filelist -c copy $res_dir/combined__$(join(basename.(vidlist), '_')).mp4`
print(cmd)
try
    @ffmpeg_env run(cmd)
    rm(temp_filelist)
catch e
    @error "FFmpeg command failed: $e"
end

find_vid_vs_audio_syncdiff_timesegment(vidlist[1], audlist[1]; flag_verbose=true)

#~ only need to run this
include("run_example.jl")

res_dir = "/media/spin/anas2/data_res/dolphin/calf/temp/delete"
res_dir = ""
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Double_ball_S1"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Double_ball_S1"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_4males_19_09_2025"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_4males_19_09_2025"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Calibration_19_09_2025"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Calibration_19_09_2025"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_3males_19_09_2025"; res_dir=res_dir, auto_segment_len=120, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Two-way_3males_19_09_2025"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

run_analysis_split_vidau("/media/spin/anas2/data/calf/upload/Double_ball_S1"; res_dir=res_dir)

delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/1"; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/2"; res_dir="", auto_segment_len=250, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/3"; res_dir="", auto_segment_len=250, flag_norm_rms=true)
delays, conf = run_analysis_split_vidau("/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/4"; res_dir="", auto_segment_len=250, flag_norm_rms=true)

aufname = "/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/1/1.F6_Clicker.WAV"
vidfname = "/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/1/1.GoPro_Clicker.MP4"

####################################################################################################################################################################################
#~ run analysis on folder with a set of video and audio files
# already streamed line the the block below, ignore this
include("run_example.jl")
include("video.jl")
set_device__ophk_acoustic_D3()

folname = "/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/1"
res_dir = "/media/spin/anas2/data_res/dolphin/calf/calibration/D3/localization0"
delays, conf, output_vidname, vidlist, audlist, res_dir2 = run_analysis_split_vidau(folname; res_dir=res_dir, auto_segment_len=250, flag_norm_rms=true)

results = process_detections.(audlist, Ref(vidlist[1]); res_dir=res_dir2)
# res[1][1] = res[1][1] .+ (delays*get_fps(vidfname))

detection_pixels = joinpath(res_dir2, "detection_pixels.csv")
detection_type = 1  # 1: impulsive, 2: tonal, 3: boat
cum_duration = 0.0
for (ind, res) in enumerate(results)
    write_mode = ind==1 ? "w" : "a"
    open( detection_pixels, write_mode) do io
        writedlm(io, ["frame" "px" "py"], ',')
        writedlm(io, [res[detection_type][1] .+ ( (cum_duration + delays)*get_fps(vidlist[1])) res[detection_type][2]], ',') # add sync delay and cumulative duration
    end
    cum_duration += get_duration(audlist[ind])
end

# vidpath = "/media/spin/anas2/data_res/dolphin/calf/temp/delete/1/combined__1.GoPro_Clicker.MP4.mp4"
overlay_boxes_on_video(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaid.mp4"; radius=100)
overlay_boxes_on_video_imageonly(detection_pixels, output_vidname, splitext(output_vidname)[1]*"_overlaidIMG"; radius=100)


##########################################################
## run the function
include("run_example.jl")
include("video.jl")
set_device__ophk_acoustic_D3()

set_device__ophk_acoustic_D3_clicker()
folders = "/media/spin/anas2/data/calf/Calibration/Housingv3_test_28_03_2025/" .* string.(1:4)
run_contiguous_folders.(folders; res_dir = "/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2",
    overlay_radius=25)#, flag_overlayvideo=false)

folder = "/media/spin/anas2/data/calf/upload/Two-way_2males_06_10_2025"
folder = "/media/spin/anas2/data/calf/Calibration/device3_calibration/Calibration_19_09_2025"
run_contiguous_folders(folder; res_dir = "/media/spin/anas2/data/calf/upload/results", flag_overlayvideo=false)

folders = ["/media/spin/anas2/data/calf/upload/Two-way_2males_06_10_2025",
"/media/spin/anas2/data/calf/upload/Two-way_3males_19_09_2025", 
    "/media/spin/anas2/data/calf/upload/Two-way_4males_15_10_2025",
    "/media/spin/anas2/data/calf/upload/Two-way_4males_19_09_2025",
    "/media/spin/anas2/data/calf/Calibration/device3_calibration/Calibration_19_09_2025"]
run_contiguous_folders.(folders; res_dir = "/media/spin/anas2/data/calf/upload/results3", flag_overlayimages=false)

#######
include("pic2vid.jl")
pic2vid("/media/spin/anas2/data_res/dolphin/calf/calibration/D3/localization5/4/combined__4.GoPro_Swimming_Clicker_SDrivers.MP4_overlaidIMG",
    "/media/spin/anas2/data_res/dolphin/calf/calibration/D3/localization5/4/combined__4.GoPro_Swimming_Clicker_SDrivers.MP4_overlaidIMG_2.mp4";
    auto_mode=true)

##########################################
#~ calibrate with manually annotated pixel points
res_dir = "/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1"
vidfname = joinpath(res_dir, "combined__1.GoPro_Clicker.MP4.mp4")
dur = get_duration(vidfname)
fps = get_fps(vidfname)
d = load("/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1/1.F6_Clicker_t466.4552869078375_d15000__cps0.375_angles.jld2")
d["ang_impulsive"]
dp = CSV.read("/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1/detection_pixels.csv", DataFrame)
dp = dp[ 0 .<= dp.frame .<= dur*fps, :]
dp.frame_index = 0:length(dp.frame)-1

d_annotated = CSV.read("/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1/combined__1.GoPro_Clicker.MP4_overlaidIMG_20251018_122735.csv", DataFrame)
rename!(d_annotated, Dict("frame" => "frame_index"))

dd = leftjoin(dp, d_annotated, on=:frame_index)
sort!(dd, :frame_index)
# remove rows with missing values
dropmissing!(dd, [:px, :py, :x, :y])

Plots.scatter(dd.px, dd.x; alpha=0.3, markerstrokewidth=0, label="x vs px", xlabel="Pixel position", ylabel="Annotated position")
Plots.scatter!(dd.py, dd.y; alpha=0.3, markerstrokewidth=0, title="Pixel vs annotated positions", label=["py vs y"], xlabel="Pixel position", ylabel="Annotated position")
# regression fit for px to x0
using GLM#, StatsModels
using RobustModels
form = @formula(x ~ px)
model_x = fit(LinearModel, form, dd)
model_x_r = rlm(form, dd, MEstimator{TukeyLoss}(); initial_scale=:L1, ridgeλ=1.0)
Plots.scatter(dd.px, dd.x; alpha=0.3, markerstrokewidth=0, label="x vs px", xlabel="Pixel position", ylabel="Annotated position")
Plots.plot!(dd.px, predict(model_x); label="linear fit line")
Plots.plot!(dd.px, predict(model_x_r); label="robust fit line")

form = @formula(y ~ py)
model_y = fit(LinearModel, form, dd)
model_y_r = rlm(form, dd, MEstimator{TukeyLoss}(); initial_scale=:L1, ridgeλ=1.0)
Plots.scatter(dd.py, dd.y; alpha=0.3, markerstrokewidth=0, label="y vs py", xlabel="Pixel position", ylabel="Annotated position")
Plots.plot!(dd.py, predict(model_y); label="linear fit line")
Plots.plot!(dd.py, predict(model_y_r); label="robust fit line")

# use camera calibration matrix parameters to get the correction function instead to predict the new data points
using Optim
include("camera2.jl")
# distort_pixels(vcat([dd.px dd.py]', ones(length(dd.px),)'))
pixels = vcat(dd.px', dd.py', ones(length(dd.px),)')
annotated_pts = vcat(dd.x', dd.y', ones(length(dd.x),)')
# optimize the camera parameters to minimize the residuals between annotated and predicted points
r2 = optimize(input -> sum(abs2.(annotated_pts .- distort_pixels(pixels, input))), [0.0,0.0,1.0,1.0,0.0,0.0,0.0],
    LBFGS())#; autodiff = :forward) #BFGS() #SimulatedAnnealing())#LBFGS()
correction_vals = Optim.minimizer(r2)
pixels_corrected = distort_pixels(pixels, correction_vals)
annotated_pts - pixels_corrected

# px_range = range(minimum(dd.px), stop=maximum(dd.px), length=100)
# py_range = range(minimum(dd.py), stop=maximum(dd.py), length=100)

# # predict with new fitted models
# pred_x = predict(model_x_r, DataFrame(px=px_range))
# pred_y = predict(model_y_r, DataFrame(py=py_range))

# create function with fitted models
coef(model_x_r)
coef(model_y_r)

fitted_func_x(px) = coef(model_x_r)[1] .+ coef(model_x_r)[2] .* px
fitted_func_y(py) = coef(model_y_r)[1] .+ coef(model_y_r)[2] .* py

Plots.scatter(dd.px, dd.x; alpha=0.3, markerstrokewidth=0, label="x vs px", xlabel="Pixel position", ylabel="Annotated position")
Plots.plot!(dd.px, fitted_func_x.(dd.px); label="fitted line x vs px")
Plots.scatter(dd.py, dd.y; alpha=0.3, markerstrokewidth=0, label="y vs py", xlabel="Pixel position", ylabel="Annotated position")
Plots.plot!(dd.py, fitted_func_y.(dd.py); label="fitted line y vs py")

dd.x_new = fitted_func_x.(dd.px)
dd.y_new = fitted_func_y.(dd.py)
Plots.scatter(dd.px, dd.py; alpha=0.3, markerstrokewidth=0, label="pixel positions", xlabel="px", ylabel="py")
Plots.scatter(dd.x_new, dd.y_new; alpha=0.3, markerstrokewidth=0, label="corrected positions", xlabel="x_new", ylabel="y_new")
Plots.scatter!(dd.x, dd.y; alpha=0.3, markerstrokewidth=0, label="annotated positions", xlabel="x", ylabel="y")

CSV.write(joinpath(res_dir, "combined__1.GoPro_Clicker.MP4_overlaidIMG__annotated_corrected.csv"), dd)

# keep only frame_indx, x_new, y_new, change the names to frame, x, y
dd_new = copy(dd)
dd_new = select(dd_new, :frame_index, :x_new, :y_new)
rename!(dd_new, Dict("frame_index"=>"frame", "x_new" => "x", "y_new" => "y"))
# dd_new = select(dd_new, :frame_index, :x, :y)
overlay_annotations_on_video(
       dd_new,
       "/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1/combined__1.GoPro_Clicker.MP4_overlaidIMG__annotated.mp4",
       "/media/spin/anas2/data_res/dolphin/calf/calibration/D2/localization13_fovH-80_small2/1/combined__1.GoPro_Clicker.MP4_overlaidIMG__annotated_corrected.mkv";
       x_col=:x, y_col=:y,
       default_color="green@0.4",
       default_radius=10,
       mode=:stream
       )



"""
    pixel_to_world_correction(K, dist, R, t; z_plane=0.0)

Create a correction closure that maps image pixel coordinates (px, py) to world coordinates (x, y)
assuming the world points lie on a plane with z = z_plane.

Inputs:
- K: 3x3 camera intrinsic matrix (Float64).
- dist: distortion coefficients vector (k1,k2,p1,p2,k3, ...). Missing values are assumed zero.
- R: 3x3 rotation matrix that transforms world -> camera coordinates.
- t: 3-element translation vector that transforms world -> camera coordinates.

Returns: a function f(px, py) -> (x, y) for scalar inputs or (xs, ys) for array inputs.

Notes:
- This uses a simple iterative undistortion (Brown–Conrady) and a pinhole projection model.
- Assumes R and t follow X_cam = R*X_world + t. The function computes X_world from the
  camera ray intersecting the plane z = z_plane.
"""
function pixel_to_world_correction(K::AbstractMatrix, dist::AbstractVector, R::AbstractMatrix, t::AbstractVector; z_plane::Real=0.0)
    # unpack intrinsics
    fx = K[1,1]
    fy = K[2,2]
    cx = K[1,3]
    cy = K[2,3]

    # distortion coefficients with safe defaults
    k1 = length(dist) >= 1 ? dist[1] : 0.0
    k2 = length(dist) >= 2 ? dist[2] : 0.0
    p1 = length(dist) >= 3 ? dist[3] : 0.0
    p2 = length(dist) >= 4 ? dist[4] : 0.0
    k3 = length(dist) >= 5 ? dist[5] : 0.0

    # small helper to map a single point
    function _map_single(px::Real, py::Real)
        # convert to normalized image coordinates
        x = (px - cx) / fx
        y = (py - cy) / fy

        # iterative undistortion (approximate inverse of Brown-Conrady)
        xu = x
        yu = y
        for _ in 1:6
            r2 = xu^2 + yu^2
            radial = 1.0 + k1*r2 + k2*r2^2 + k3*r2^3
            x_dist = xu*radial + 2.0*p1*xu*yu + p2*(r2 + 2.0*xu^2)
            y_dist = yu*radial + p1*(r2 + 2.0*yu^2) + 2.0*p2*xu*yu
            # update estimate by small correction (fixed-point like)
            xu += x - x_dist
            yu += y - y_dist
        end

        # direction vector in camera coordinates
        cam_dir = SVector(xu, yu, 1.0)

        # Solve for scale s such that the back-projected world point has z = z_plane.
        # Assuming X_cam = s * cam_dir, and X_cam = R * X_world + t => X_world = R'*(X_cam - t)
        # So X_world = s*(R'*cam_dir) - R'*t. Let a = R'*cam_dir, b = R'*t, then s = (z_plane + b[3]) / a[3]
        a = R' * collect(cam_dir)
        b = R' * collect(t)
        s = (z_plane + b[3]) / a[3]
        Xw = s .* a .- b
        return Xw[1], Xw[2]
    end

    # returned closure accepts scalars or arrays
    function mapper(px, py)
        if isa(px, AbstractArray) && isa(py, AbstractArray)
            n = length(px)
            xs = Vector{Float64}(undef, n)
            ys = Vector{Float64}(undef, n)
            for i in 1:n
                xs[i], ys[i] = _map_single(px[i], py[i])
            end
            return xs, ys
        elseif isa(px, AbstractArray) && !isa(py, AbstractArray)
            n = length(px)
            ys = Vector{Float64}(undef, n)
            xs = Vector{Float64}(undef, n)
            for i in 1:n
                xs[i], ys[i] = _map_single(px[i], py)
            end
            return xs, ys
        elseif !isa(px, AbstractArray) && isa(py, AbstractArray)
            n = length(py)
            xs = Vector{Float64}(undef, n)
            ys = Vector{Float64}(undef, n)
            for i in 1:n
                xs[i], ys[i] = _map_single(px, py[i])
            end
            return xs, ys
        else
            return _map_single(px, py)
        end
    end

    return mapper
end


# Rodrigues rotation vector -> rotation matrix
function rodrigues(r::AbstractVector{<:Real})
    θ = norm(r)
    if θ == 0.0
        return I(3)
    end
    k = r / θ
    K = [  0.0  -k[3]  k[2];
          k[3]   0.0  -k[1];
         -k[2]  k[1]   0.0 ]
    return Matrix(I(3) + sin(θ)*K + (1-cos(θ))*(K*K))
end


# internal helper: build K, dist, R, t from parameter vector
function _unpack_camera_params(p::AbstractVector)
    # p layout: [fx, fy, cx, cy, k1, k2, p1, p2, k3, r1, r2, r3, tx, ty, tz]
    fx, fy, cx, cy = p[1], p[2], p[3], p[4]
    k1, k2, p1, p2, k3 = p[5], p[6], p[7], p[8], p[9]
    rvec = p[10:12]
    t = p[13:15]
    K = [fx 0.0 cx; 0.0 fy cy; 0.0 0.0 1.0]
    dist = [k1, k2, p1, p2, k3]
    R = rodrigues(rvec)
    return K, dist, R, t
end


# model that maps parameter vector p and input (px,py) -> vector of predicted [x; y]
function _camera_model(p, xdata)
    # xdata is expected as a 2×N or vector concatenated [px; py]
    if isa(xdata, AbstractMatrix) && size(xdata,1) == 2
        px = vec(xdata[1,:])
        py = vec(xdata[2,:])
    else
        # assume concatenated vector [px; py]
        n = length(xdata) ÷ 2
        px = xdata[1:n]
        py = xdata[n+1:2n]
    end
    K, dist, R, t = _unpack_camera_params(p)
    mapper = pixel_to_world_correction(K, dist, R, collect(t); z_plane=0.0)
    xs, ys = mapper(px, py)
    return vcat(xs, ys)
end

using Random, LsqFit, StaticArrays
"""
    fit_camera_from_point_pairs(px, py, x, y; ransac_iters=150, sample_size=8, inlier_thr=0.05)

Robustly fit camera parameters so that pixel coords (px,py) map to world (x,y) using the
`pixel_to_world_correction` model. Returns a NamedTuple with fields:
 - params: fitted parameter vector [fx,fy,cx,cy,k1,k2,p1,p2,k3,r1,r2,r3,tx,ty,tz]
 - inlier_mask: boolean vector marking inliers
 - rms: root-mean-square error on inliers
 - residuals: per-point residual distances

This uses RANSAC to reject outliers, and LsqFit (Levenberg-Marquardt) for local refinement.
"""
function fit_camera_from_point_pairs(px::AbstractVector, py::AbstractVector, x::AbstractVector, y::AbstractVector; 
        ransac_iters::Int=150, sample_size::Int=8, inlier_thr::Real=0.05, verbose::Bool=false)
    # using Random, LsqFit

    N = length(px)
    @assert length(py) == N && length(x) == N && length(y) == N

    # initial parameter guess
    fx0 = 1000.0
    fy0 = fx0
    cx0 = mean(px)
    cy0 = mean(py)
    k10 = 0.0; k20 = 0.0; p10 = 0.0; p20 = 0.0; k30 = 0.0
    r0 = zeros(3)
    t0 = [0.0, 0.0, 1.0]
    p0 = vcat([fx0, fy0, cx0, cy0, k10, k20, p10, p20, k30], r0, t0)

    best = nothing
    best_inliers = falses(N)
    best_count = 0
    best_rms = Inf

    xydata = vcat(px, py)
    yobs = vcat(x, y)

    for iter in 1:ransac_iters
        idx = rand(1:N, sample_size)
        pxs = px[idx]; pys = py[idx]; xs = x[idx]; ys = y[idx]
        xdat_sample = vcat(pxs, pys)
        yobs_sample = vcat(xs, ys)
        try
            fit = curve_fit(_camera_model, xdat_sample, yobs_sample, p0; autodiff = :forward)
        catch err
            continue
        end
        p_est = coef(fit)
        ypred = _camera_model(p_est, xydata)
        xs_pred = ypred[1:N]; ys_pred = ypred[N+1:2N]
        resid = sqrt.( (x .- xs_pred).^2 .+ (y .- ys_pred).^2 )
        inliers = resid .<= inlier_thr
        cnt = count(inliers)
        if cnt > best_count
            # refine using all inliers
            xdat_in = vcat(px[inliers], py[inliers])
            yobs_in = vcat(x[inliers], y[inliers])
            try
                fit2 = curve_fit(_camera_model, xdat_in, yobs_in, p_est; autodiff=:forward)
            catch err2
                continue
            end
            p_ref = coef(fit2)
            ypred_ref = _camera_model(p_ref, xydata)
            xs_pred_ref = ypred_ref[1:N]; ys_pred_ref = ypred_ref[N+1:2N]
            resid_ref = sqrt.( (x .- xs_pred_ref).^2 .+ (y .- ys_pred_ref).^2 )
            inliers_ref = resid_ref .<= inlier_thr
            rms = sqrt(mean(resid_ref[inliers_ref].^2))
            best = (params = p_ref, residuals = resid_ref)
            best_inliers = inliers_ref
            best_count = count(inliers_ref)
            best_rms = rms
            if verbose
                @info "RANSAC iter=$iter best_count=$best_count rms=$best_rms"
            end
        end
    end

    if best === nothing
        # fallback: fit all points
        try
            fit_all = curve_fit(_camera_model, xydata, yobs, p0; autodiff=:forward)
            p_all = coef(fit_all)
            ypred_all = _camera_model(p_all, xydata)
            xs_pred = ypred_all[1:N]; ys_pred = ypred_all[N+1:2N]
            resid_all = sqrt.( (x .- xs_pred).^2 .+ (y .- ys_pred).^2 )
            inlier_mask = resid_all .<= inlier_thr
            rms = sqrt(mean(resid_all[inlier_mask].^2))
            return (params=p_all, inlier_mask=inlier_mask, rms=rms, residuals=resid_all)
        catch err
            error("Camera fit failed: $err")
        end
    end

    return (params = best.params, inlier_mask = best_inliers, rms = best_rms, residuals = best.residuals)
end

#############
include("run_example.jl")
# include("video.jl")
set_device__ophk_acoustic_D3()
#~ run on dolphin recordings
folder = "/media/spin/anas2/data/calf/upload/Single_dolphin_target"
run_contiguous_folders.(readdir(folder, join=true); res_dir = "/media/spin/anas2/data/calf/upload/results/results8_20251031_033500", flag_overlayimages=false, flag_rm_oldfile=true)

folders=["/media/spin/anas2/data/calf/upload/Double_4males_24_10_2025/New_acoustics_Double_4males_24_10_2025",
    "/media/spin/anas2/data/calf/upload/Double_ball_S1",
    "/media/spin/anas2/data/calf/new_2025_freeplay_report/Ball_3way_S2"
]
run_contiguous_folders.(folders;
    res_dir = "/media/spin/anas2/data/calf/upload/results/good_stuff/results5_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_rm_oldfile=true)


## run on one set
### small
run_contiguous_folders("/media/spin/anas2/data/calf/upload/Single_dolphin_target/17_10_2025_Tutti";
    res_dir = "/media/spin/anas2/data/calf/upload/results/res_old/results5_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_rm_oldfile=true)
run_contiguous_folders("/media/spin/anas2/data/calf/upload/Single_dolphin_target/24_10_2025_Tutti"; res_dir = "/media/spin/anas2/data/calf/upload/results5_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_overlayimages=false)

### normal
run_contiguous_folders("/media/spin/anas2/data/calf/upload/Two-way/Two-way_3males_19_09_2025";
    res_dir = "/media/spin/anas2/data/calf/upload/results/results6_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_rm_oldfile=false)
res = run_contiguous_folders("/media/spin/anas2/data/calf/upload/Two-way/Two-way_3males_19_09_2025";
    res_dir = "/media/spin/anas2/data/calf/upload/results/results6_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_rm_oldfile=true,
    overlay_radius=:in_annotations, default_color=:in_annotations, default_alpha=:in_annotations, default_shape=:in_annotations,
    flag_return=true)

folder = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_3males_19_09_2025"
folder = "/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025"
res = run_contiguous_folders(folder;
    res_dir = "/media/spin/anas2/data/calf/upload/results/results7_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_rm_oldfile=false,
    overlay_radius=:in_annotations, default_color=:in_annotations, default_alpha=:in_annotations, default_shape=:in_annotations,
    flag_return=true, flag_overlayvideo=false, force_extension_type=".flac"
    )
res = run_contiguous_folders(folder;
    res_dir = "/media/spin/anas2/data/calf/upload/results/results7_$(Dates.format(now(), "yyyymmdd_HHMMSS"))", flag_rm_oldfile=false,
    overlay_radius=:in_annotations, default_color=:in_annotations, default_alpha=:in_annotations, default_shape=:in_annotations, detection_types = [1,4],
    flag_return=true, flag_overlayvideo=true, force_extension_type=".flac"
    )



# clicker calibration
set_device__ophk_acoustic_D3_clicker()
run_contiguous_folders("/media/spin/anas2/data/calf/Calibration/device3_calibration/Calibration_19_09_2025";
    res_dir = "/media/spin/anas2/data_res/dolphin/calf/calibration/D3/clicker", flag_rm_oldfile=true, flag_overlayimages=true)


## get clips
res_dir = "/media/spin/anas2/data/calf/upload/results/results7_20251105_132919/Two-way_2males_06_10_2025"
auda_file = "/media/spin/anas2/data/calf/upload/results/results7_20251105_132919/Two-way_2males_06_10_2025/combined_impulse__audacity.txt"
aufname = "/media/spin/anas2/data/calf/upload/results/results7_20251105_132919/Two-way_2males_06_10_2025/combined__251006_001_0001.WAV_251006_001_0002.WAV_251006_001_0003.WAV.flac"
detection_px_file = "/media/spin/anas2/data/calf/upload/results/results7_20251105_132919/Two-way_2males_06_10_2025/detection_pixels.csv"
delay_file = "/media/spin/anas2/data/calf/upload/results/results7_20251105_132919/Two-way_2males_06_10_2025/sync_delay.csv"
delays = CSV.read(delay_file, DataFrame)[1, :delay_s]

data, fs = readAudio(aufname)
data = data[:, get_relevant_channels(rx_vect)]

df = CSV.read(auda_file, DataFrame; header=false)
t = round.(Int, df[!, 1] .* fs)
windows = t .+ [window_impulsive[1] window_impulsive[end]] .+ 1

data_filt = filter_simple(data, impulsive_band_pass; fs=fs); data=nothing;
clips = extract_clips(data_filt, windows)
clipsm = extract_clips_matrix(data_filt, windows)

angs = detection2angle(data_filt, t, rx_vect[:,get_relevant_channels(rx_vect)]; fs=fs, window=window_impulsive, return_residual=true,
    getTDOA_func=default_getTDOA_func)
angs_xcorr = detection2angle(data_filt, t, rx_vect[:,get_relevant_channels(rx_vect)]; fs=fs, window=window_impulsive, return_residual=true,
    getTDOA_func=get_tdoa_raw_MaxEnergyRefChannel)
p_pixels = angle2px(angs[1][1], fov_angle)
p_pixels_xcorr = angle2px(angs_xcorr[1][1], fov_angle)

Plots.scatter(angs[1][2]; alpha=0.01, markerstrokewidth=0)
ylims!(0, 2e-8)

df_px = CSV.read(detection_px_file, DataFrame)
df_px[:, :frame] = round.(Int, df_px[:, :frame])
df_px_im = filter(row -> row.shape =="circle", df_px)
# df_px_to = filter(row -> row.shape !="circle", df_px)

# try different tdoa methods
tdoa_funcs = [get_tdoa_envelope, get_tdoa_envelope_filtered, get_tdoa_findTrigger, get_tdoa_max, get_tdoa_min, get_tdoa_minmax, get_tdoa_raw, get_tdoa_raw_MaxEnergyRefChannel, get_tdoa_raw_MaxEnergyRefChannel_resample, get_tdoa_raw_MaxPeakRefChannel]#, get_tdoa_raw_flexi]
irange = 1:size(df_px_im, 1)
# for (i, func) in enumerate(tdoa_funcs)
collections = Array{String}(undef, length(tdoa_funcs))
df_new = Array{DataFrame}(undef, length(tdoa_funcs))
Threads.@threads for i in 1:length(tdoa_funcs)
    func = tdoa_funcs[i]
    angs_tmp = detection2angle(data_filt, t, rx_vect[:,get_relevant_channels(rx_vect)]; fs=fs, window=window_impulsive, return_residual=true,
        getTDOA_func=func)
    p_pixels_tmp = angle2px(angs_tmp[1][1], fov_angle)
    # export csv for visualization
    # irange = 1:size(angs_tmp[1][1], 1)
    csv_fname = joinpath(res_dir, "impulse_$(string(func)).csv")
    df_new[i] = DataFrame([df_px_im[irange, 1] p_pixels_tmp[irange, :] repeat([15+i*2 'a'+i string(func) "square"], length(irange)) angs_tmp[1][2][irange]],
        ["frame", "x", "y", "radius", "tag", "tag_name", "type", "residuals"])
    CSV.write(csv_fname, Tables.table([df_px_im[irange, 1] p_pixels_tmp[irange, :] repeat([15+i*2 'a'+i string(func) "square"], length(irange)) angs_tmp[1][2][irange]]);
        header=["frame", "x", "y", "radius", "tag", "tag_name", "type", "residuals"]
    )
    collections[i] = csv_fname
end
# df_new = [CSV.read(fname, DataFrame) for fname in collections]
# df_new = vcat(df_new...)
CSV.write(joinpath(res_dir, "impulse_all_tdoa_methods.csv"), vcat(df_new...))
# get best method based on residuals
df_best_residuals = similar(df_new[end-1])
df_best = copy(df_new[end-1])
Threads.@threads for i in 1:size(df_best, 1)
    # least residual
    residuals = [df_new[j][i, :residuals] for j in 1:length(df_new)]
    best_idx = argmin(residuals)
    df_best_residuals[i, :] = df_new[best_idx][i, :]

    # median pixel position
    xs = [df_new[j][i, :x] for j in 1:length(df_new)]
    ys = [df_new[j][i, :y] for j in 1:length(df_new)]
    # Plots.scatter(xs; labels="x"); hline!([median(xs) mode(xs)]; labels=["median" "mode"])
    # Plots.scatter!(ys; labels="y"); hline!([median(ys) mode(ys)]; labels=["median" "mode"])
    df_best[i, :x] = mode(xs)
    df_best[i, :y] = mode(ys)
end
CSV.write(joinpath(res_dir, "impulse_best-residuals_tdoa_method.csv"), df_best_residuals)

df_best.radius .= 50; df_best.tag .= "m"; df_best.tag_name .= "mode"; df_best.type .= "square"
CSV.write(joinpath(res_dir, "impulse_best-mode_pixel_position.csv"), df_best)

i=63382 #40341;#39947; #39895; #106355 #106344 #40800
# angs[1][1] # angles, angs[2][1] # residuals, angs[2] #tdoa

plot(clips[i])
tdoa = angs[2][i, :]
angs[1][1][i, :]
angs[1][2][i] # residual
p_pixels[i, :]

# plot(data_filt[121789310 .+ (-1000:1000),:]) # buzz
# plot(data_filt[256399551:256600000,:])

tdoa_diff = maximum(tdoa) .- tdoa
p=plot();
for j in 1:length(tdoa)
    plot!(p,[zeros(tdoa_diff[j],); clips[i][:,j]])
end
display(p)

savepath = "temp/minusing_31"; mkpath(savepath)
err_sum = 0
for i = 1:100:length(clips)
    er = false
    @info i
    tdoa = angs[2][i, :]
    tdoa_min = minimum(tdoa)
    tdoa = tdoa .- tdoa_min .+ 31
    p=plot();
    for j in 1:length(tdoa)
        # @info tdoa[j]-30:tdoa[j]+30
        try 
            plot!(p,clips[i][tdoa[j]-30:tdoa[j]+30,j])
        catch e
            @warn "Error plotting clip $i channel $j"#: $e"
            @info "tdoa: $tdoa"
            er = true; 
        end
    end
    if er
        err_sum += 1
    end
    plot(p; title=string(i)) |> display
    savefig(p, joinpath(savepath, @sprintf("%010i.png", i)))
    # readline()
end
err_sum

irange = 106344:106783; length(irange)
irange = 109606:109609
i = irange[2]
plot(clipsm[:,1,irange])
plot(clipsm[:,:,irange[4]])
[irange df_px_im[irange,:frame] p_pixels[irange, :] angs[1][2][irange] p_pixels_xcorr[irange, :] angs_xcorr[1][2][irange]] |> showall

[df_px_im[irange,:frame] ]

#~ export csv for visualization
irange = 1:size(df_px_im, 1)
CSV.write("temp/impulse_mimax.csv", Tables.table([df_px_im[!, :frame] p_pixels repeat([15 "a" "min_max" "square"], length(irange))]);
    header=["frame", "x", "y", "radius", "tag", "tag_name", "type"]
)
CSV.write("temp/impulse_rawMaxEnergy.csv", Tables.table([df_px_im[!, :frame] p_pixels_xcorr repeat([15 "b" "raw_MaxEnergyRefChannel" "square"], length(irange))]);
    header=["frame", "x", "y", "radius", "tag", "tag_name", "type"]
)


###########################################################################
#~ run analysis on split video and audio files

data_dir = "/media/spin/anas2/data/calf/Single_ball"
result_directory = "/media/spin/anas2/data_res/dolphin/calf/Single_ball/res_20250910"

data_folders = sort(readdir(data_dir; join=true); lt=natural) #sort(readdir(data_dir; join=true) |> filter(isdir) |> filter(x -> occursin("Single", x)); lt=natural)
ftype = ".flac"
mkpath(result_directory)

# tabulate data
summary_filepath = joinpath(result_directory, "summary.csv")
# tabulate_data("/media/spin/anas2/data/marecet/datai/deployment_18012024_18042024"; output=summary_filepath, filetype="wav", fname2timestamp_func=fname2dt_soundtrap, flag_filepath=true)
tabulate_data(data_dir; output=summary_filepath, filetype="flac", fname2timestamp_func=DEFAULT_fname2timestamp_func, flag_filepath=true)
# df = CSV.read(summary_filepath, DataFrame)
# df = sort(df, :datetime)
# CSV.write(summary_filepath, df)

for aufol in data_folders
    @info "--- Processing folder: $aufol"
    files = readdir(aufol; join=true) |> filter(endswith(ftype))
    res_dir = joinpath(result_directory, Dates.format(DEFAULT_fname2timestamp_func(files[1]), "yyyymmdd_HHMMSS"))
    isnothing(res_dir) || mkpath(res_dir); cp("web/index.html", joinpath(res_dir, "index.html"))
    for file in files
        fname, ext = splitext(basename(file))
        @info "Processing file: $(dirname(file))/\033[32m$fname\033[0m$ext"
        try
            detect_impulseNtonal(file, res_dir)
        catch e
            @error "Error processing file $file: $e"
        end
        GC.gc() # run garbage collector to free memory
    end
end














asdf
#~ old code below ##########################################################
template_window = (-fs2÷2+1:fs2÷2) #(-fs2÷50+1:fs2÷50)
max_time = argmax(d_z)


@info "Audio sync time: $(max_time/fs2)s"
win2 = template_window .+ max_time
template_sig = d_z[win2]
template_sig |> PlotlyJS.plot
# m = mfilter(template_sig, data2[:,1])


m = mfilter(template_sig, @view(data2[:,1]))
video_sync_time = argmax(abs.(m))
correl_sign = m[video_sync_time] |> sign

[data2[video_sync_time .+ (1:length(template_window))] template_sig*correl_sign] |> plot_norm
@info "Video sync time: $(video_sync_time/fs2)s"
@info "Audio sync time: $(max_time/fs2)s"
vid_aud_sync_diff = video_sync_time - max_time
@info "Video-Audio sync time difference: $(vid_aud_sync_diff/fs2)s"


aud_vid_sync_diff = audio_sync_time - max_time
@info "Audio-Video sync time difference: $(aud_vid_sync_diff/fs2)s"
@info "Audio-Video sync time difference(samples): $(aud_vid_sync_diff)"
win_test = (1:2000) .+ max_time
plot_norm([d_z[win_test .+ aud_vid_sync_diff] data2[win_test,1]])#; title="Test plot around video max correlation point")

########################################
# Further analysis
########################################
data2[]
d_z[argmax(m) .+ (-fs2÷25+1:fs2÷25)] |> PlotlyJS.plot
#  |> PlotlyJS.plot

delays = align_signals([d_z[1:fs2*5], data2[1:fs2*5,1]], Delay(delay_method=DTWDelay()))

plot(d_z[win2[1]+delays[1][1]:win2[1]+delays[1][1]+1000,1])
plot(data2[win2[1]:win2[1]+1000,1])

win = (23*fs + 87843):(23*fs + 91223)

# sum(data[win,3:5]; dims=2) |> plot
a=resample(data[win,1], fs2/fs)
win2 = 19*fs2+5000+58:19*fs2+5845+58
[a -data2[19*fs2+5000+58:19*fs2+5845+58,1]] |> plot_norm
d_z=resample(data[:,1], fs2/fs)

# sanity check
[sum(data[win,3:5]; dims=2) data[win,1]] |> plot_norm

(data2[:,1] - data2[:,2]) .|> abs |> mean # both channels of video audio are identical
( norm_max(sum(data[:,3:5]; dims=2)) - norm_max(data[:,1])) .|> abs |> mean # sum of last 3 channels of audio are identical to the first two channels