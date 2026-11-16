# Clicker camera–audio calibration: fix the A/B fits and integrate them for real recordings

## Context

`reviews/camera_calib_review_2026-11-15.md` says the clicker calibration (`src/clicker_calibration/`) is wired into the pipeline through `preset_bridge.jl`. The A (GCC) and B (SRP) profiles are committed but provisional. Goal: find what is wrong with the calibration and make A and B usable on real dolphin recordings.

Constraints:
- existing data only, all under `/media/spin/anas/data/calf`; search nowhere else;
- explicit device argument;
- typical site sound speed;
- lens mode and medium unknown, so infer them from the data;
- existing functions untouched; fixes go into new functions and files;
- implementation by several Herdr agents (`/herdr-agent`) with the model and effort given per task below.

## Findings: calibration errors (most severe first)

1. **B's profile is degenerate; do not use it.** fx = 1647 px but **fy = 204 px**, which is a 150° vertical FOV (at its bound).
   - c = 1401 m/s, also at its bound.
   - The offsets of +26.3/−30.6 µs (±5–6 samples, about 25% of the 22.4-sample maximum TDOA) move the acoustic ray by about 26°. The −25.9°/−12.7° rotation turns it back by the same amount. They cancel; neither is physical.
2. **The 8-parameter model has built-in parameter trade-offs (structural non-identifiability)** (`physical_fit.jl:_pf_decode`/`_pf_pixels`). Near the image centre the pixel is ≈ f·c·B⁻¹(τ−b):
   - the offsets b act like a pitch/yaw rotation;
   - c and f only appear as a product;
   - fx and fy are free.
   - Evidence: in the same water on the same day, A and B fit c = 1562 and 1401 m/s; the c objective is flat (2.75%); A's pose varies by 22° across folds; A's offsets changed between refits (13.6/2.6 µs → 1.1/0.05 µs).
3. **The loss lets the fit give up on bad events.** `_pf_loss` caps residuals at about 155 px, and an invalid ray costs the same cap. That is how A has a 44 px median but a 1755 px RMSE. Pixel RMSE is dominated by near-90° rays; use angular error.
4. **Shrinkage (errors-in-variables) bias.** Pixel regression on noisy TDOAs shrinks c, fx and fy toward the image centre, consistent with B's c and fy.
5. **Lens model is probably wrong.** `Housingv3_test_28_03_2025/Notes.rtf` says "GoPro9 – 2.7K, SV, 60FPS". SuperView is a non-linear horizontal stretch that a centred pinhole with k1 = k2 = 0 cannot represent. It may also be why the checkerboard fit is imperfect. A's fx/fy of 1.22 could be partly real under SuperView.
6. **Weak truth and coverage.** Terra estimates only ("uncertain", 35–50 px); one short, spatially narrow sequence per device; an already-inspected evaluation block.
7. **Effects the model does not include:** camera–array parallax (depends on range), the front/back mirror ambiguity of the planar array, flat-port refraction, and near-field curvature (minor).

## Integration gaps for real recordings

- **I1. The calibration forces the clicker detector settings.** Activation calls `set_device__ophk_acoustic_D3_clicker()` (1 kHz high-pass, `dist_impulsive = 15000`, an enforced −300:700 window). Dolphin settings cannot run calibrated, and a later `set_device__ophk_acoustic_D3()` call switches the calibration off without warning.
- **I2. Offsets are tied to one estimator.** Tonal, Impulse2–4 and `run_multitdoa` (`tdoa2dir_ls`, `default_c`) get the wrong correction.
- **I3. Sound speed is stored in the rig profile.** It should be an input for each deployment.
- **I4. No drift correction** (see "Clock drift: where the numbers come from" below).
- **I5. NaN pixels.** NaN pixels go to `detection_pixels.csv`; `video.jl` has no NaN handling (to verify).
- **I6. No guard against mismatched videos, silent deactivation, no record of what ran.** No resolution or lens-mode check against the profile, silent deactivation, no provenance file. Report only, don't edit: the legacy D3 setters use `imsize = (1520, 2740)`, a typo for 2704.
- **I7. Possible frame off-by-one.** Frames are 1-origin plus a fractional delay·fps (`aspod.jl`, `run_example.jl:403`); the calibration uses 0-origin.
- **I8. Profiles may expire.** `upload/README.md` says B's housing went "back for reprint due to the leakage issue". If B was remounted after 2026-02-09, its rotation does not carry over, so profiles need validity windows.

### Clock drift: where the numbers come from
- The 61.4 ppm (A) and 65.3 ppm (B) values were not computed in this session. They are outputs of the calibration preparation (`prep_synchronize`, `prepare.jl:356`), stored in `calibration_runs/clicker_20260209/20261030_v1/{A,B}/synchronization.json` and quoted in `docs/clicker_calibration.md` §Device A/B preparation.
- **Method:** multi-window normalized waveform correlation of the GoPro's 48 kHz audio track against Zoom physical channel 1 (the mixdown). That gives 17 anchors from 74 s to 1349 s across all chapters. A robust affine fit `t_video = scale·t_recorder + offset` then gives:

  | | A | B |
  |---|---|---|
  | scale | 1.0000614 | 1.0000653 |
  | offset | −10.036 s | −9.641 s |
  | RMS / max residual | 56 / 121 µs | 75 / 146 µs |
  | scale uncertainty | about 0.2 ppm | about 0.2 ppm |
- **Caveat:** this is one session of about 21 minutes per device. Whether the drift is stable across sessions or temperatures is unknown, so the plan **measures drift per recording** instead of assuming 61–65 ppm.
- **The pipeline today:** `find_vid_vs_audio_syncdiff_timesegment` (`synchronization.jl:194`) correlates one 250 s segment near the loudest video-audio peak and returns a single constant delay. The error grows by about 0.06 ms per second away from that segment: roughly 110 ms, or 6–7 frames at 59.94 fps, after 30 minutes.

## Data attribution (corrected with your device dates; all under `/media/spin/anas/data/calf`)

Device A has been in use since about Sept 2025 (GoPro counter chain γ, 1265→2051). Device B has been in use since about Feb 2026 (chain β, 0059→0133). The older recorder is `.mat` + `.mkv`.

| Device | Calibration data | Real recordings for rollout |
|---|---|---|
| **A** | **Checkerboard:** `Calibration/device3_calibration/checkerboard_2025-09-19` (15 MP4s), `checkerboard_2025-09-22` (9 MP4s). **Clicker:** `device3_calibration/Calibration_19_09_2025` (GoPro + F6), `upload/Calibration_05_11_2025` (human annotations in the D3/clicker results), `Clicktest/clicker_20260209/A`. **Other:** `device3_calibration/3d_modeling/camera_params.json` (per-view intrinsics, a sanity check only) | `upload/TT-shape/TDA_06022026/A`, `upload/Single_dolphin_target/13_11_2026_Tutti/A` |
| **B** | **Clicker:** `Clicktest/clicker_20260209/B` only. **No B-era checkerboard found.** | `upload/TT-shape/TDA_06022026/B`, `upload/Single_dolphin_target/13_11_2026_Tutti/B` (subject to the remount question) |
| neither (Mar–May 2025, before A and B) | `temp/delete/from_mac/cam_calib/Calibration18032025` (38 MP4s, plus `SV.MP4`/`W.MP4`), `Calibration/device_2_calibration`, `Housingv3_test_*` | Chain α (98xx–9986) wraps plausibly into B's 0059. The physical camera may have moved into B; it could serve only as a camera-intrinsics prior for B, and only if you confirm that. It sits under `temp/delete/`, so preserve it before any cleanup. |

**OLD device (inventoried; out of scope this round).** IP camera, HEVC 3840×2160 at 20 fps with an 8 kHz µ-law audio track, plus NI acoustics in `.mat`/`.flac`/`.bin`.
- **Checkerboard videos:** `Calibration/Checkerboard_Calibration/20191128` (with `Procedure.txt`; file name says "air"), `2023-03-31` (352 PNG frames), `2023-04-21` (copy in `old/cam_calib`), `Calibration/checkerboard_2023`.
- **Clicker:** `old/clicker_calib` (2019–2023, `.mkv` + `.flac`).
- **FOV ruler measurements:** `Calibration/cam_fov` (2018).
- **Other:** `pinger/` (2018), `sync_test/` (2018).
- **Old-device sync clip:** `Calibration/device_new/Housingv3_test28_02_2025/Old_device` (2 kHz tone on both video and audio).
- **What is missing:** no intrinsic result for the old camera, and no checkerboard board or square size recorded anywhere.
- **Code:** neither calf nor this repo contains a checkerboard solver (`cam_calib.jl` is a pixel-to-annotation fit only).
- **Security:** `old/old/Calf_Project/runVid67.py` holds a plaintext RTSP camera credential. Move it out of the file or rotate it; the plan does not touch it.

## Plan

### Step 0: Housekeeping (lead)
- Copy the approved plan into `reviews/plans/`, with the first draft as `_v0`.
- Save memories:
  - the device timeline (A since ~Sept 2025, B since ~Feb 2026, old = `.mat`/`.mkv`);
  - calibration data lives only under `/media/spin/anas/data/calf`;
  - inventories and broad file searches go to a Haiku subagent; the lead does only very small searches itself;
  - implement with Herdr agents and give each task a model and effort, preferring Claude opus/sonnet over Codex Sol for coding;
  - commit progressively, staging only this session's own paths.

### Step 1: Checkerboard intrinsics for each device (improving the checkerboard code)
- **Code to start from.** No checkerboard solver exists in this repo or in calf. If your current checkerboard script lives elsewhere, give T2 its path to improve. Otherwise T2 writes `tools/checkerboard_intrinsics.py` (OpenCV 4.10 is installed).
- **Board size.** Intrinsics need only the inner-corner count, which the tool detects automatically. The unknown square size affects only metric extrinsics, which this plan does not use.
- **Lens-mode evidence (T1).** For every checkerboard, clicker and rollout clip, extract the GPMF `VFOV`/`ZFOV` tags where present and record them in `calibration_runs/intrinsics/lens_manifest.csv`. For B, also search `calf` for Feb 2026+ checkerboard clips under any filename (by `creation_time`).
- **Improvements (T2):**
  - blur rejection (Laplacian variance);
  - `findChessboardCornersSB` with automatic board-size detection and sub-pixel refinement;
  - per-view outlier rejection;
  - a frame-coverage map, requiring corners near the edges and corners of the frame;
  - three lens families: (a) pinhole with k1, k2, p1, p2, k3; (b) `cv2.fisheye`; (c) SuperView-style, a 4:3 sensor model followed by an odd-polynomial horizontal remap;
  - selection by reprojection error on held-out **whole videos**.
- **Per-device runs:**
  - **A:** fit on 2025-09-19 and validate on 2025-09-22, then the reverse. That is a cross-session stability check; the final fit uses both sessions.
  - **B:** use a B-era checkerboard if T1 finds one. Otherwise B gets A's lens family and shape as a prior (same HERO9, and the same mode if GPMF agrees), and records a B checkerboard as future data. The α set is used only if you confirm it is the same camera.
- **Report:** RMS error, per-view errors, coverage, fx/fy ratio, the implied FOV, and an air-vs-water verdict (a flat port reduces the in-air FOV by roughly 1.33×); for A, also a comparison with `camera_params.json`.
- **Output:** `calibration_runs/intrinsics/{A,B}/intrinsics.json` in the `cam_projection.jl:get_cam_fov_fromparam` format, plus a `model` block.
- **Visual QA (T2v):** contact sheets showing detected corners and undistorted frames.

### Step 2: Constrained, identifiable refit on the existing clicker data
New `src/clicker_calibration/constrained_fit.jl`, with one `include` and an export in `ClickerCalibration.jl`.
- **Reuse:**
  - `_pf_load` (`physical_fit.jl:613`), with the GCC/SRP delays from `estimator_review.csv` as `_mr_load` does (`method_refinement.jl:156`);
  - the frozen blocks and purge;
  - `clicker_geometry`, `_pf_rotation`, `project_tdoas`;
  - inputs `calibration_runs/clicker_20260209/20261030_v2_fit/A` and `.../20261030_v2/B`.
- **Physics fixed:**
  - b = 0 (simultaneous F6 sampling; hydrophone phase 0.07–0.09 rad);
  - one c shared by A and B, site default 1530 m/s, overridable;
  - lens shape from Step 1, falling back to a square-pixel pinhole.
- **Rays:**
  - acoustic: u = [uv, √(1−|uv|²)] with uv = −c·B⁻¹τ;
  - camera: the annotation pixel undistorted through the Step 1 model;
  - impossible TDOAs are invalid and counted, never clamped.
- **Fit:** rotation R by robust weighted Wahba/SVD with IRLS on the angular residual; scale f by a bounded 1-D search. No free cap.
- **Model ladder:**
  - M0: R only;
  - M1: R + f;
  - M2: M1 + one distortion scale;
  - M3: M1 + offsets with a σ = 2 µs prior (diagnostic only).
  - Choose model and estimator by the existing leave-one-training-block-out rule: median angular error, ties to the simpler model. Replay the evaluation block only after selection.
- **Optional, pre-declared:** a temporal-consensus flag (the `select_candidates` idea from `multiarrival.jl`). Report with and without it; truth is never relabelled.
- **Independent check for A:** prepare `Calibration_05_11_2025` with `prepare_device`. Convert the human CSV into `annotations.csv`, with `o`-tagged rows excluded. Use it as untouched validation for v1 and v2, never for fitting (provided A's mounting was unchanged between Nov 2025 and Feb 2026).
- **Outputs:**
  - `calibrations/clicker_20260209_{A,B}_v2.json`, schema `clicker-physical-calibration/2.0`, with the same projection fields and offsets `[0,0,0]`;
  - new fields: `estimator: "agnostic"`, `sound_speed_role: "runtime_environment"`, `lens_model`, `lens_mode`, `camera_chain`, `valid_from`/`valid_until`, fold stability, angular metrics;
  - new projection functions for lens families (b)/(c) next to `_pf_project_direction`, which stays untouched.
- **Report**, v1 vs v2 on identical events: median/p80/p90 angle, % within 5°/10°, coverage, pixel error for in-frame events only, and fold spread. Saved to `calibration_runs/clicker_20260209/constrained_v2/`.

### Step 3: v2 bridge and explicit-device runner
New `src/clicker_calibration/preset_bridge_2.jl`, included by one line appended to the end of `preset_bridge.jl`. That keeps `src/config.jl`, which has others' pending edits, untouched. Definitions are guarded like `preset_bridge.jl`'s, so the repeated includes of `config.jl` stay safe.
- `load_calibration_2(path; device)`: the `workflow.jl:load_calibration` checks for schema 2.0, plus b = 0 and lens-block validation.
- `with_clicker_calibration_2(f, device::Symbol; c_m_s=1530, detection=:dolphin, calibration_path, tdoa_estimator=nothing)`:
  - snapshot and restore globals (`_clicker_global_snapshot`/`_clicker_changed_names`);
  - apply `set_device__ophk_acoustic_D3()` (or the clicker preset);
  - install the profile with c = `c_m_s`, plus `imsize`/`fov_angle`/`FILE_device_ID` (`..._A2`/`_B2`);
  - inside the block, set `c`/`default_c`;
  - keep the preset's TDOA estimator, or use the optional one (for example a short-window wrapper around `prep_gcc`, `prepare.jl:140`);
  - the existing `tdoa2dir`/`angle2px` dispatch then applies v2 unchanged.
- `run_contiguous_folders_calibrated(folder, device; res_dir, c_m_s, kwargs...)`:
  - asserts the calibration is active before and after;
  - refuses a mismatched resolution or GPMF lens mode;
  - warns outside the validity window;
  - writes `calibration_used.json` (device, SHA-256, c, preset, lens mode, git commit).

### Step 4: NaN safety and frame origin
- Check how the renderer handles NaN. If needed, add a post-step that writes `detection_pixels_valid.csv` (finite rows only) plus a status sidecar for the calibrated overlay.
- Verify I7 with a synthetic impulse.

### Step 5: Drift-aware sync for each recording
- New function: fit offset + drift per folder with `prep_synchronize`/`prep_wave_match`/`prep_fit_sync` (`prepare.jl:356/239/265`) and write `sync_affine.json`.
- Frames are mapped as `(t·scale + offset)·fps`.
- Fall back to the constant delay with a warning; full-folder runs log to `docs/error.log`.

### Step 6: Tests (new files only)
- `test/clicker_calibration/test_constrained_fit.jl`
- `test/clicker_calibration/test_presets_2.jl`
- synthetic-board tests for each lens family

The cases are as listed in the Step 2–5 contracts.

### Step 7: Real-recording rollout
Run `run_contiguous_folders_calibrated` on short excerpts of `TT-shape/TDA_06022026/{A,B}` (four days before the clicker take) and then `13_11_2026_Tutti/{A,B}`. Compare v2 with legacy: % valid, % in-frame, direction stability along trains, and a visual overlay review.

## Agent assignment (Herdr)

The lead is this session (Claude Opus 5.5, xhigh effort). It runs `hagent init cal-lead` once, spawns each agent with `hagent spawn <claude|codex> <name> --model M --effort E -- "TASK"`, and reviews every `DONE:` message. Heavy jobs are announced to you before launch and run through `tools/run_clicker_job.py`. Names must be lowercase and unique.

| ID | Task | Agent / model / effort | Depends on | Deliverable |
|---|---|---|---|---|
| T1 | Lens-mode GPMF manifest; checkerboard clip list per device; B-era checkerboard search (calf only) | claude / **haiku** / medium | none | `lens_manifest.csv`, clip lists |
| T2 | Checkerboard tool improvements + A cross-session fits + B fit or prior | claude / **opus** / high (sonnet / high for the follow-up fixes) | T1 | `tools/checkerboard_intrinsics.py`, `intrinsics/{A,B}` |
| T2v | Visual QA of corners and undistortion contact sheets | codex / **gpt-6.1-terra** / medium | T2 | QA notes, flagged views |
| T3 | `constrained_fit.jl`, model ladder, v2 JSONs, comparison report | claude / **opus** / xhigh | T2 (lens model); code can start in parallel with a pinhole default | v2 profiles + report |
| T3r | Independent review of the T3 statistics: bias, CV leakage, fixed parameters held | codex / **gpt-6-astra** / high | T3 | review notes |
| T4 | `preset_bridge_2.jl`: loader, scoped block, runner, provenance, resolution and lens checks | claude / **sonnet** / high | none (uses a mock v2 JSON) | bridge + runner |
| T5 | Drift-aware sync per recording | claude / **sonnet** / high | none | sync function + `sync_affine.json` |
| T6 | NaN-safe overlay CSV, frame-origin check | claude / **haiku** / high | none | post-step + finding |
| T7 | Tests for T3–T6 and the checkerboard tool, written by a different agent than the implementers | claude / **sonnet** / medium | T2–T6 | new test files passing |
| T8 | Run the jobs: refits, `Calibration_05_11_2025` preparation, rollout excerpts; log failures to `docs/error.log` | claude / **haiku** / medium | T3–T7 | job receipts, results |
| T9 | Visual review of rollout overlays (legacy vs v1 vs v2) for A and B | codex / **gpt-6.1-terra** / high | T8 | review gallery notes |
| T10 | Docs: `docs/clicker_calibration.md` v2 section, a review note in `reviews/` | claude / **sonnet** / medium | T3, T8 | docs |
| T11 | Final code review of the whole change (`/code-review high`) | claude / **fable** / high | all | findings, fixed by owners |

Waves: (1) T1, T4, T5, T6 → (2) T2, T2v, and T3 coding → (3) T3 fits, T3r, T7 → (4) T8, T9, T10 → (5) T11, then the lead integrates.

## Git: commit progressively, this session's work only
Branch `dev2`. Before this session started, the working tree already had uncommitted changes by others: `.CondaPkg/meta`, `docs/error.log`, `src/config.jl` (HydroMoth dugong preset), `src/run_hydromoth.jl` and `src/run_syncVid.jl`, plus untracked files `reviews/camera_calib_review_2026-11-15.md`, `src/click_features.jl` and `tools/clicktrain_agent_template/`. Never stage these.
- **Only the lead commits.** Herdr agents never commit, and they report the exact paths they changed in their `DONE:` message.
- **One commit per reviewed task or wave**, made after its tests pass. The first commit (Step 0) holds the plan copies.
- **Stage by explicit path only** (`git add <paths>`, never `-A`, `.` or `-a`), then check `git diff --cached --stat` against the task's file list before committing.
- **Avoid files others have touched.** `preset_bridge_2.jl` is included by one line appended to `preset_bridge.jl` (clean in git), not from `src/config.jl`. Long-run failure entries go to `docs/error.log` as the long-task skill requires. That file is not committed, because it already holds others' pending lines; if our entries ever need committing, stage only our hunk via `git diff docs/error.log > patch`, cut the patch down to our lines, and `git apply --cached`.
- **Commit messages** use the repo style (`feat(clicker): …`, `doc(plans): …`) and end with the session's attribution lines. No pushing unless you ask.

## Future data that would help most
- Per rig: clicker on a 3×3 grid plus the corners, at 1, 2 and 4 m, holding still for several clicks, with silent gaps between positions.
- A separate validation take, never inspected.
- Human or LED-confirmed positions.
- Water temperature and salinity.
- A measured camera-to-array offset.
- An underwater checkerboard per camera in the deployment lens mode, on the same day as the clicker take; B needs one most.

## Verification
- Run the new tests and the existing `test/clicker_calibration/test_presets.jl` (to show v1 is unchanged), using `julia --project`. Run the Python tool tests.
- Step 1: A's cross-session reprojection error and parameter agreement; lens family consistent with GPMF.
- Step 2: B's v2 has no degenerate aspect ratio and a rotation within a few degrees of A's unless the mount differs. Check the v1/v2 angular table and the A check on `Calibration_05_11_2025`.
- Step 7 on one excerpt per rig: the overlay video, `calibration_used.json` and `sync_affine.json`, and NaN-free `detection_pixels_valid.csv`.
