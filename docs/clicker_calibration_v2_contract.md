# Clicker calibration v2: shared implementation contract

Working contract for the agents that implement
`reviews/plans/read-reviews-camera-calib-review-2026-11-jaunty-nova.md`.
Interfaces change only through the lead. Rules that apply to everyone:

- Never edit existing functions. Add new functions in new files. The only existing-file edits allowed are the ones listed under "Owned files".
- Agents never commit, stage, stash, reset or checkout. Report the exact paths you changed in the final `DONE:` message.
- Never modify these files, which hold others' pending work: `.CondaPkg/meta`, `src/config.jl`, `src/run_hydromoth.jl`, `src/run_syncVid.jl`, `src/click_features.jl`, `tools/clicktrain_agent_template/`. You may append entries to `docs/error.log` (long-task skill), but never rewrite it.
- Source data is read-only. Calibration and recording data come only from `/media/spin/anas/data/calf`, plus the repo's ignored `calibration_runs/`.
- Ask the lead before any run that takes longer than about 10 minutes, uses more than about 8 GB of RAM, or decodes whole videos.

## Schema `clicker-physical-calibration/2.0` (JSON)

v1 field names are kept wherever the existing projection code reads them (`project_tdoas`, `project_angles`, `_pf_project_direction`, `calibrated_tdoa2dir`, `calibrated_angle2px`).

| key | type | rule |
|---|---|---|
| `schema_version` | string | `"clicker-physical-calibration/2.0"` |
| `device` | string | `"A"` or `"B"` |
| `method` | string | `"agnostic"`; v2 is not tied to a TDOA estimator |
| `estimator_spec` | object | `{"method":"agnostic","sample_rate_hz":192000.0}`; `_clicker_sample_rate` reads `sample_rate_hz` |
| `physical_channels` | [3,4,5] | as v1 |
| `reference_channel` / `reference_physical_channel` | 3 / 5 | as v1 |
| `geometry_m` | 3×3 rows | must equal `clicker_geometry()` (D3, 0.18 m) |
| `relative_delay_offsets_s` | [0,0,0] | exactly zero in v2 |
| `sound_speed_m_s` | number | the value used during fitting; the runtime overrides it (see below) |
| `sound_speed_role` | string | `"runtime_environment"` |
| `rotation_array_to_camera` | 3×3 | proper SO(3) |
| `rotation_euler_xyz_rad` | [3] | same convention as v1 (`_pf_rotation`, Rz·Ry·Rx) |
| `width`, `height` | int | 2704, 1520 |
| `camera` | object | `fx, fy, cx, cy, k1, k2`: always present, the pinhole or radial view; exact for the `pinhole` family |
| `lens` | object | `{"family": "pinhole" \| "fisheye" \| "superview", "params": {...}, "source": "<intrinsics.json path + sha256>" \| "clicker_fit_only"}` |
| `lens_mode` | string or null | GPMF lens mode (e.g. `"W"`, `"S"`, `"L"`) when known |
| `camera_chain` | string | e.g. `"gamma_1265-2051"` (A) or `"beta_0059-0133"` (B) |
| `valid_from`, `valid_until` | ISO date or null | the mounting validity window; the runner warns outside it |
| `model_ladder_choice` | string | `"M0"` to `"M3"` (Step 2) |
| `metrics` | object | angular metrics (median/p80/p90 degrees, % within 5°/10°, coverage) and fold stability |
| `provenance` | object | input hashes, software, git commit |

`lens.params` per family:
- `pinhole`: `{}`; uses `camera`.
- `fisheye`: OpenCV equidistant model, `{"fx","fy","cx","cy","k1","k2","k3","k4"}`.
- `superview`: `{"inner": {<pinhole or fisheye params>}, "inner_family": "pinhole"|"fisheye", "sensor_width_px", "remap_coeffs": [a1,a3,a5,...]}`. The inner model projects onto a 4:3 sensor frame, then the odd polynomial `x' = Σ a_k x^k` (normalized horizontal coordinate) maps it to output pixels. The exact definition belongs to T2/T3 and must be written into `tools/checkerboard_intrinsics.py` and `constrained_fit.jl`.

## Julia types and dispatch (inside `module ClickerCalibration`)

- `struct CalibrationV2 <: AbstractDict{String,Any}` wraps `Dict{String,Any}`. It forwards `getindex`, `get`, `haskey`, `keys`, `values`, `length`, `iterate` and `setindex!`. All existing code that takes `AbstractDict` keeps working unchanged.
- `load_calibration_2(path; device=nothing) -> CalibrationV2` validates every rule in the table and throws `ErrorException` with a clear message on failure.
- Lens-aware projection is a new **method**, `_pf_project_direction(calibration::CalibrationV2, direction)`. The untyped original stays untouched. For `pinhole` it applies exactly the original formula. For `fisheye` and `superview` it calls `_v2_lens_project(Val(family), lens_params, Xc/Zc, Yc/Zc)`, which T3 implements. T4 provides a pinhole-only version and an `ArgumentError` stub for the other families. Because `project_tdoas` and `project_angles` call `_pf_project_direction(calibration, …)`, they become lens-aware with no edits.
- Inverse (for the fit and for undistorting annotations): `_v2_lens_unproject(Val(family), params, x_px, y_px) -> (xn, yn)`, owned by T3.

## Runtime sound speed

The default site value is `CLICKER_SITE_SOUND_SPEED_M_S = 1530.0` (warm seawater; overridable). Activation copies the profile and replaces `sound_speed_m_s` with the caller's `c_m_s` (1400 to 1600). It also sets the module globals `c` and `default_c` to `c_m_s` inside the scoped block, so `tdoa2dir_ls` agrees.

## Owned files and interfaces

| Task | Owns (new files) | Allowed edit to existing file |
|---|---|---|
| T3 | `src/clicker_calibration/constrained_fit.jl`; the fisheye/superview methods of `_v2_lens_project`/`_v2_lens_unproject` | `ClickerCalibration.jl`: one `include("constrained_fit.jl")` line plus exports |
| T4 | `src/clicker_calibration/calibration_v2.jl` (type, `load_calibration_2`, pinhole `_v2_lens_*`); `src/clicker_calibration/preset_bridge_2.jl` (activation, scoped block, runner); `src/run_example_2.jl` (`run_contiguous_folders_2`) | `ClickerCalibration.jl`: one `include("calibration_v2.jl")` line plus exports; `preset_bridge.jl`: one line appended at the very end, `include(joinpath(@__DIR__,"preset_bridge_2.jl"))`; `src/DolphinSpeaker.jl`: one `include("run_example_2.jl")` plus exports |
| T5 | `src/clicker_calibration/folder_sync.jl` | `ClickerCalibration.jl`: one `include("folder_sync.jl")` line plus exports |
| T6 | `src/detection_pixels_valid.jl` | `src/run_example_2.jl` is T4's, so T6 does NOT edit it. T4 includes `detection_pixels_valid.jl` from `run_example_2.jl`. |

`ClickerCalibration.jl` include order: the existing includes, then `calibration_v2.jl`, `folder_sync.jl`, `constrained_fit.jl`. Each agent appends only its own line.

### T5: `fit_folder_sync_affine`
```julia
fit_folder_sync_affine(vidlist::AbstractVector{<:AbstractString}, audlist::AbstractVector{<:AbstractString};
                       out_path=nothing, kwargs...) -> Dict{String,Any}
```
- Reuses `prep_wav_info`, `prep_video_info` and `prep_synchronize` (`prepare.jl`).
- Recorder timeline: the concatenated `audlist` (the sum of each file's samples/fs).
- Video timeline: the concatenated chapters, as `prep_synchronize` defines it (`global_video_timeline`).
- Returns `"status"` (`"valid"`, `"insufficient_anchors"` or `"failed"`), `"scale"`, `"offset_s"`, `"drift_ppm"`, `"rms_residual_s"`, `"max_residual_s"`, `"n_anchors"`, `"anchor_min_time_s"`, `"anchor_max_time_s"`, `"video_timeline"`, `"recorder_timeline"`, `"inputs"`.
- Writes JSON when `out_path` is given.

Helper: `video_time_from_recorder(t_rec_s, sync) = sync["scale"]*t_rec_s + sync["offset_s"]`. The convention is `t_video = scale·t_rec + offset`. In the existing pipeline, the constant-delay equivalent is `scale = 1` and `offset = delays` (`frame = (t + cum + delays)·fps`).

### T6: `write_valid_detection_pixels`
```julia
write_valid_detection_pixels(csv_path; out_path=replace(csv_path, ".csv"=>"_valid.csv"),
                             status_path=replace(csv_path, ".csv"=>"_status.csv"), imsize=imsize) -> NamedTuple
```
- Input header: `frame,px,py,radius,r,g,b,a,shape`.
- The valid CSV keeps only rows with finite `frame`, `px` and `py`.
- The status CSV has `row,frame,status`, where status is `ok`, `nonfinite` or `outside_image`. `outside_image` rows stay in the valid CSV, so nothing is clamped.
- T6 also reports how the existing overlay renderer treats NaN rows, and whether frame origins are consistent (I7).

### T4: `run_contiguous_folders_2` and the calibrated runner
- `run_contiguous_folders_2(folname; sync=:affine | :constant, pixel_filter=true, kwargs...)` is a NEW COPY of `run_contiguous_folders`:
  - frames come from exact detection times (`pind_good_inS` in the returned `res`) through the sync model;
  - with `sync=:affine` it calls `_clicker_call(:fit_folder_sync_affine, vidlist, audlist; out_path=joinpath(res_dir2, "sync_affine.json"))`, and falls back to the constant delay with `@warn` when the status is not `valid`;
  - the overlay reads `detection_pixels_valid.csv` when `pixel_filter=true`.
- `run_contiguous_folders_calibrated(folder, device::Symbol; res_dir, c_m_s=CLICKER_SITE_SOUND_SPEED_M_S, calibration_path=default_v2_path(device), detection=:dolphin, tdoa_estimator=nothing, sync=:affine, kwargs...)`:
  - runs `run_contiguous_folders_2` inside `with_clicker_calibration_2`;
  - asserts the calibration is active before and after;
  - checks the video resolution against the profile;
  - checks the GPMF lens mode against the profile when both are known;
  - warns outside the validity window;
  - writes `calibration_used.json` into the results folder.
- `default_v2_path(device)` is `src/clicker_calibration/calibrations/clicker_20260209_$(device)_v2.json`. That file is created by T3; until it exists, tests use a fixture under `test/clicker_calibration/fixtures/`.

### T4: `with_clicker_calibration_2`
```julia
with_clicker_calibration_2(f, device::Union{Symbol,AbstractString};
    c_m_s=CLICKER_SITE_SOUND_SPEED_M_S, detection=:dolphin, calibration_path=default_v2_path(device), tdoa_estimator=nothing)
```
- Validates before changing anything.
- Snapshots and restores globals with `_clicker_global_snapshot`/`_clicker_changed_names`.
- `detection=:dolphin` calls `set_device__ophk_acoustic_D3()`; `detection=:clicker` calls `set_device__ophk_acoustic_D3_clicker()`.
- Installs `ACTIVE_CLICKER_CALIBRATION` (a `CalibrationV2` copy carrying `c_m_s`) and `ACTIVE_CLICKER_TOKEN`, plus `imsize`, `fov_angle` (from fx/fy), `ref_channel = 3`, `FILE_device_ID = "ophk_acoustic_D3_$(device)2"`, and `c`/`default_c = c_m_s`.
- TDOA estimator: if `tdoa_estimator !== nothing`, use it. Otherwise, if `default_getTDOA_func === _clicker_tdoas` (left by a v1 activation), use `CLICKER_UNCALIBRATED_TDOA`. Otherwise keep the current one. v2 never installs `_clicker_tdoas`.
- `active_clicker_calibration()`, `tdoa2dir` and `angle2px` then apply v2 through the existing dispatch.
