# Bounded source-level reproductions for audio_to_video_projection_review_2026-10-30.md.
# Run from the repository root with: julia --startup-file=no --project=. reviews/audio_to_video_projection_checks_2026-10-30.jl
# Selected definitions are read from the reviewed Git commit, not the moving
# working tree. Do not include the workflow scripts, which have top-level
# processing and media-writing side effects.
using LinearAlgebra, Statistics, Random, Printf
using Optim, LsqFit, StaticArrays, ForwardDiff

const REPO = normpath(joinpath(@__DIR__, ".."))
const REVIEW_COMMIT = "d22056c39b949ee5bead6ddd3b9644e4261ab267"
function load_definition(relative_path, marker)
    source = read(Cmd(["git", "-C", REPO, "show", "$REVIEW_COMMIT:$relative_path"]), String)
    location = findfirst(marker, source)
    isnothing(location) && error("Definition not found: $marker")
    expression, _ = Meta.parse(source, first(location))
    Core.eval(Main, expression)
end

default_c = 1540.0
ref_channel = 3
default_window = -5000:10000
default_tdoa2dir_solver = NelderMead
imsize = (1520, 2740)
fov_angle = [80.0, 54.46086642417594]
for marker in ["function cost_tdoa2ang(", "mod_sign(a,b) =", "function tdoa2dir(",
               "function angle2pxd(", "function get_tdoa_minmax(data, peaks;",
               "_clip_window(win, n) =", "function get_tdoa_minmax2(",
               "function tdoa2uv(", "function tdoa2dir_ls("]
    load_definition("src/localization.jl", marker)
end
load_definition("src/dsp.jl", "function extrema_and_indices(arr::Array{T, 2}")
for marker in ["function pixel_to_world_correction(", "function rodrigues(",
               "function _unpack_camera_params(", "function _camera_model(",
               "function fit_camera_from_point_pairs("]
    load_definition("src/dev2_sync.jl", marker)
end

function expected_error(label, exception_type, f)
    try
        f()
    catch exception
        @assert exception isa exception_type
        # Keep output short, without huge concrete ForwardDiff type expansions.
        println(label, ": ", nameof(typeof(exception)))
        return
    end
    error("Expected an error for $label")
end

println("Source commit: ", REVIEW_COMMIT)
println("Julia ", VERSION, "; Optim ", pkgversion(Optim), "; LsqFit ", pkgversion(LsqFit),
        "; ForwardDiff ", pkgversion(ForwardDiff))

# A negative direct pulse and a smaller positive lobe with different peak times.
waveforms = zeros(30, 3)
for channel in 1:3
    waveforms[9 + channel, channel] = -10
    waveforms[17 + 2channel, channel] = 1
end
old_peaks = get_tdoa_minmax(waveforms, [15]; window=-14:15)
fixed_peaks = get_tdoa_minmax2(waveforms, [15]; window=-14:15)
@assert vec(old_peaks) == [19, 21, 23]
@assert vec(fixed_peaks) == [10, 11, 12]
println("Polarity example: old peak indices = ", vec(old_peaks), "; fixed = ", vec(fixed_peaks))
expected_error("Old minmax at file boundary", BoundsError,
    () -> get_tdoa_minmax(waveforms, [2]; window=-5:5))

# Compare the actual linear-angle function to a specified, hypothetical pinhole
# camera with the same FOV and optical axes. This is not a measured camera error.
fx = imsize[2] / (2tand(fov_angle[1] / 2))
fy = imsize[1] / (2tand(fov_angle[2] / 2))
for (az, el) in [(20.0, 0.0), (30.0, 0.0), (30.0, 20.0)]
    linear = vec(angle2pxd([az el], fov_angle, imsize))
    pinhole = [imsize[2]/2 + fx*tand(az), imsize[1]/2 - fy*tand(el)/cosd(az)]
    @printf("Projection az=%.0f el=%.0f: linear=(%.3f, %.3f), pinhole=(%.3f, %.3f), delta=(%.3f, %.3f) px\n",
            az, el, linear..., pinhole..., (linear - pinhole)...)
end

rx_complex = 0.18/sqrt(3) .* exp.(im .* deg2rad.([-30 90 210]))
rx = [real(rx_complex); imag(rx_complex); zeros(1, 3)]
sample_rate = 192000
angles = deg2rad.([30.0, 10.0])
direction = [sin(angles[1])*cos(angles[2]), sin(angles[2]), cos(angles[1])*cos(angles[2])]
arrival_samples = vec(-direction' * rx / default_c * sample_rate)
reflected_angles = [pi - angles[1], angles[2]]
front_cost = cost_tdoa2ang(arrival_samples, angles, rx, sample_rate)
back_cost = cost_tdoa2ang(arrival_samples, reflected_angles, rx, sample_rate)
@assert front_cost < 1e-20 && back_cost < 1e-20
println("Planar front/back cost: ", front_cost, ", ", back_cost)

# The c variable in config.jl is not the default_c used by the cost function.
c = 1480.0
implicit_cost = cost_tdoa2ang(arrival_samples, angles, rx, sample_rate)
explicit_cost = cost_tdoa2ang(arrival_samples, angles, rx, sample_rate; c=c)
@assert implicit_cost < 1e-20 && explicit_cost > 0.0
@printf("After c=1480: default cost=%.6f, explicit c=1480 cost=%.6f samples^2\n", implicit_cost, explicit_cost)

impossible = reshape([100.0, 0.0, 0.0], 1, :)
uv = tdoa2uv(impossible, rx, sample_rate)
bad_angle, bad_residual = tdoa2dir(impossible, rx, sample_rate; return_residual=true)
@assert norm(uv) > 1 && all(isfinite, bad_angle) && bad_residual[1] > 1.0
@printf("Impossible delays: uv norm=%.3f; optimizer still returns finite angles, residual=%.3f samples^2\n",
        norm(uv), bad_residual[1])
@printf("D3 0.18 m baseline at 192 kHz: max lag=%.3f samples; broadside one-sample angle=%.3f degrees\n",
        0.18*sample_rate/default_c, asind(default_c/(0.18*sample_rate)))

# Demonstrate the two time/index transformations without opening a video.
fps = 30.0
export_frame_at_zero = round(Int, 0.0*fps) + 1
@assert export_frame_at_zero == 1
@printf("t=0: exporter frame=%d; zero-based renderer uses frame 1 at %.6f s; image exporter uses %.6f s\n",
        export_frame_at_zero, export_frame_at_zero/fps, (export_frame_at_zero-1)/fps)
local_time, sync_offset = 0.49/fps, 0.49/fps
two_rounds = round(Int, (round(Int, local_time*fps)+1) + sync_offset*fps)
one_round = round(Int, (local_time+sync_offset)*fps)+1
@assert two_rounds == 1 && one_round == 2
println("Early frame rounding example: two rounds = ", two_rounds, "; one round = ", one_round)

source_frames = [100, 100, 200]
montage_frames = sort(unique(source_frames)) # same names overwrite, then filenames are sorted
wrong_second_match = source_frames[2]
@assert montage_frames[2] == 200 && wrong_second_match == 100
println("Duplicate-frame montage: source rows = ", source_frames, "; montage = ", montage_frames,
        "; ordinal join maps montage index 1 to source frame ", wrong_second_match)

# Actual camera functions, avoiding all top-level calls in dev2_sync.jl.
p0 = [1000.0, 1000.0, 1370.0, 760.0, zeros(8)..., 0.0, 0.0, 1.0]
@assert length(p0) == 15
input_x = [800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 1100.0, 1500.0]
input_y = [400.0, 800.0, 500.0, 900.0, 700.0, 600.0, 1000.0, 300.0]
xdata = vcat(input_x, input_y)
observations = _camera_model(p0, xdata)
reversed_result = _camera_model(xdata, p0) # how LsqFit actually invokes the callback
@assert length(observations) == 16 && length(reversed_result) == 14
println("Camera callback: correct output length = ", length(observations),
        "; LsqFit callback order yields length = ", length(reversed_result))
expected_error("Direct camera curve_fit", DimensionMismatch,
    () -> curve_fit(_camera_model, xdata, observations, p0; autodiff=:forward))
expected_error("ForwardDiff after correcting callback order", MethodError,
    () -> ForwardDiff.jacobian(p -> _camera_model(p, xdata), p0))
rotation_ad = ForwardDiff.jacobian(r -> vec(Matrix(rodrigues(r))), zeros(3))
rotation_fd = hcat([(vec(Matrix(rodrigues(1e-6 .* Matrix{Float64}(I, 3, 3)[:, j]))) -
                     vec(Matrix(rodrigues(-1e-6 .* Matrix{Float64}(I, 3, 3)[:, j])))) ./ 2e-6 for j in 1:3]...)
@assert norm(rotation_fd) > 2.4
@assert !all(isfinite, rotation_ad) || norm(rotation_ad - rotation_fd) > 1.0
@printf("Rodrigues Jacobian at zero: ForwardDiff norm=%.6f; central-difference norm=%.6f\n",
        norm(rotation_ad), norm(rotation_fd))
Random.seed!(20261030)
expected_error("Full camera fitter on exact synthetic pairs", ErrorException,
    () -> fit_camera_from_point_pairs(input_x, input_y, observations[1:8], observations[9:16]; ransac_iters=1))
println("All bounded review checks reproduced the reported behavior. No field-accuracy test was run.")
