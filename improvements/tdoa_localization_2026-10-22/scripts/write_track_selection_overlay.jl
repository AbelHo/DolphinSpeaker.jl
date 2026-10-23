# Viewer CSV of the track-rule selection (select_candidates: the candidate on the confirmed track with the smallest
# track-level arrival lag) from a round-2 feature table, without rerunning. Click colours and frames come from the
# Impulse6 file; circle = track rule says direct, triangle = only lagging (echo) tracks or no confirmed track.
# usage: julia --project=. improvements/tdoa_localization_2026-10-22/scripts/write_track_selection_overlay.jl sic|sicw [res_folder]
using CSV, DataFrames, Statistics, Printf
tag = length(ARGS) >= 1 ? ARGS[1] : "sicw"
res_folder = length(ARGS) >= 2 ? ARGS[2] :
    "/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device"
res = joinpath(@__DIR__, "..", "results")
folder = basename(rstrip(res_folder, '/')); sfx = folder == "New_device" ? "" : "_$folder"   # output names for other recordings get the folder suffix
f = CSV.read(joinpath(res, "data", "click_features_$(tag)_$folder.csv"), DataFrame)
w6 = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__Impulse6_srpPhat.csv"), DataFrame)
N = nrow(f)
W, H = 2740, 1520; fovh, fovv = 80.0, 54.46086642417594      # set_device__ophk_acoustic_D3: imsize (1520, 2740), fov_angle
azel2px(az, el) = (az / (fovh / 2) * (W / 2) + W / 2, el / (-fovv / 2) * (H / 2) + H / 2)   # angle2pxd
@assert all(f.frame .== w6.frame[1:N]) "click order differs from the Impulse6 file"
k = f.selected_track
az = [f[i, "az_$(k[i])"] for i in 1:N]; el = [f[i, "el_$(k[i])"] for i in 1:N]
px = azel2px.(az, el)
imp = DataFrame(frame=w6.frame[1:N], px=first.(px), py=last.(px), radius=w6.radius[1:N], r=w6.r[1:N], g=w6.g[1:N], b=w6.b[1:N],
    a=w6.a[1:N], shape=[s == "direct" ? "circle" : "triangle" for s in f.track_status])
imp.type = imp.shape
out = vcat(imp, w6[N+1:end, names(imp)])
fname = joinpath(res_folder, "detection_pixels__color255_webui__Impulse9$(tag)t_trackSelected.csv")
CSV.write(fname, out)
@printf("%s: %d impulse rows (%.1f%% circle = direct), %d tonal rows; %.1f%% inside the frame\n", basename(fname), N,
    100mean(imp.shape .== "circle"), nrow(out) - N, 100mean((0 .<= imp.px .<= W) .& (0 .<= imp.py .<= H)))
