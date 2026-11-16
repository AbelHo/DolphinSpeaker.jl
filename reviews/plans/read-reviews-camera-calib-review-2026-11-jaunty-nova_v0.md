# Clicker camera–audio calibration — fix the A/B fits and integrate them for real recordings

## Context

`reviews/camera_calib_review_2026-11-15.md` reports the clicker calibration (`src/clicker_calibration/`) is wired into the pipeline through `preset_bridge.jl`, with committed A (GCC) and B (SRP) profiles, but provisional. The user wants to know what is wrong with the calibration and how to make A and B usable on real dolphin recordings. Constraints for this round: **existing 2026-02-09 clicker data only** (new-data ideas listed at the end), **explicit device argument**, **typical site sound speed**. Checkerboard intrinsics are not reliable yet (inventory pending, Step 2). Existing functions stay untouched; fixes go in new functions/files (memory rule).

## Findings — calibration errors (ranked)

1. **B profile is degenerate — do not use.** `clicker_20260209_B.json`: fx = 1647 px but **fy = 204 px** (vertical FOV 150°, at bound) — an 8:1 pixel aspect no GoPro has; c = 1401 m/s (at bound); delay offsets +26.3/−30.6 µs (±5–6 samples, ~25% of the 22.4-sample max TDOA); rotation −25.9°/−12.7°. The offsets shift the acoustic ray by ≈26° (≈−12° az, +23° el) and the rotation turns it back by the same amount — the two cancel, so neither is a physical measurement.
2. **Structural non-identifiability of the 8-parameter model** (`physical_fit.jl:_pf_decode`, `_pf_pixels`). Near the image centre pixel ≈ f·c·B⁻¹(τ−b): (a) the 2 delay offsets b are equivalent to a small pitch/yaw rotation; (b) c and focal lengths appear only as a product; (c) fx and fy are fitted independently (no square-pixel constraint). Evidence: same-day, same-water A/B fit c = 1562 vs 1401 m/s; c objective changes only ~2.75% over 1400–1600; A fold pose varies 22°; A offsets differ between refits (13.6/2.6 µs earlier, 1.1/0.05 µs now).
3. **The robust loss lets the fit "give up"**: `_pf_loss` caps any residual at ~155 px (`max(8δ, 0.05·diag)`) and charges invalid rays the same cap, so the optimizer gains nothing by fixing large failures and can push events into "invalid". Result: A median 44 px but RMSE 1755 px. Pixel RMSE is also dominated by near-90° rays (X/Z blow-up); angular error is the meaningful metric.
4. **Errors-in-variables shrinkage**: fitting scale parameters (c, fx, fy) by pixel regression on noisy TDOAs biases predictions toward the image centre — consistent with B's c → lower bound and fy collapse, and A's fy/fx = 0.82.
5. **Weak truth and narrow coverage**: all fitted A points are Terra "uncertain" estimates (35–50 px radius), B uses unconfirmed Terra estimates; one short sequence per device (A frames 459–701, B 140–318) in a narrow image region; the evaluation block was already inspected.
6. **Unmodelled**: camera–array parallax (range-dependent, unknown range); planar-array front/back ambiguity (sources behind the rig mirror into the image); flat-port refraction / GoPro lens distortion (k1 = k2 = 0 in both); near-field curvature (minor, < 1 sample beyond 0.5 m).

## Integration gaps for real recordings

- **I1 Calibration forces clicker detector settings**: `_activate_clicker_calibration!` calls `set_device__ophk_acoustic_D3_clicker()` (HP 1 kHz, `dist_impulsive = 15000` = 78 ms min spacing, window −300:700), and `estimate_method_tdoas` throws if the window changes. Dolphin settings (`set_device__ophk_acoustic_D3`: 5 kHz, dist 77, window −50:50) cannot run calibrated; calling that setter afterwards silently deactivates the calibration.
- **I2 Method-matched offsets**: offsets belong to GCC/SRP; tonal (`get_tdoa_raw_flexi`), Impulse2–4 and `run_multitdoa` (`tdoa2dir_ls` with `default_c`) apply the wrong correction (serious with B's 5–6-sample offsets).
- **I3 Sound speed baked into the rig profile** — should be a per-deployment input.
- **I4 No drift correction**: `run_contiguous_folders` applies one constant delay (`find_vid_vs_audio_syncdiff_timesegment`), but preparation measured 61–65 ppm Zoom/GoPro drift ≈ 110 ms (≈3 frames) per 30 min — ~100+ px for a dolphin moving 2 m/s at 3 m.
- **I5 Invalid pixels**: calibrated paths return NaN pixels, written to `detection_pixels.csv`; `video.jl` has no NaN handling (to verify).
- **I6 Fragile state / provenance**: silent deactivation; no record in `res_dir` of the calibration used; `imsize` not checked against the actual video resolution/lens mode; legacy D3 setters use `imsize = (1520, 2740)` (typo for 2704 — report, don't edit).
- **I7 Frame indexing**: pipeline uses 1-origin frames plus fractional delay·fps (`aspod.jl`, `run_example.jl:403`); calibration contract is 0-origin — verify for an off-by-one.

## Plan

### Step 0 — Record the plan
Copy this approved plan unchanged to `reviews/plans/`, first draft as `_v0`.

### Step 1 — Constrained, identifiable refit on existing data
New file `src/clicker_calibration/constrained_fit.jl` (included from `ClickerCalibration.jl`, one added `include` + export).
- **Reuse** `_pf_load` (`physical_fit.jl:613`) and the per-method delays from `estimator_review.csv` as in `_mr_load` (`method_refinement.jl:156`), the frozen train/held-out blocks and purge, `clicker_geometry`, `_pf_rotation`, `project_tdoas`/`_pf_project_direction` for final projection. Inputs: the device output dirs the v1 JSONs were fitted from (A: `calibration_runs/clicker_20260209/20261030_v2_fit/A`, B: `…/20261030_v2/B`) — no re-preparation, no new annotations.
- **Model (physics-fixed)**: b = 0 (Zoom F6 simultaneous sampling; round-5 hydrophone phase 0.07–0.09 rad); c fixed and shared by A and B (site value, default 1530 m/s, overridable); square pixels fx = fy = f, centred principal point. Acoustic ray u = [uv, √(1−|uv|²)] with uv = −c·B⁻¹τ; camera ray v = normalize([(x−cx)/f, −(y−cy)/f, 1]). Impossible TDOAs → invalid, counted, never clamped.
- **Fit**: rotation R by weighted robust Wahba/Kabsch (SVD) with IRLS (Tukey/Huber on angular residual ∠(v, Ru)); weights from annotation uncertainty converted to angle (σ_px/f). No saturating cap, no "free" invalid rays. f by 1-D bounded search (HFOV 60–90°) wrapped around the inner Wahba.
- **Model ladder** M0 (R only, f from a lens prior), M1 (R + f), M2 (M1 + k1), M3 (M1 + offsets with tight prior σ = 2 µs, diagnostic only), each with GCC and SRP delays; select model × estimator with the existing leave-one-training-block-out rule (equal block weights), median angular error, pre-declared tie-break to the simpler model. Replay the reused evaluation block only after selection.
- **Optional pre-declared robustification**: temporal consensus — flag a click whose direction departs strongly from its ±1 s neighbours (multipath), using the `multiarrival.jl` `select_candidates` track-level idea; report with and without, never relabel truth.
- **Outputs**: `calibrations/clicker_20260209_A_v2.json`, `_B_v2.json`, schema `clicker-physical-calibration/2.0`, keeping the existing field names read by `project_tdoas`/`project_angles` (so projection code is reused), offsets `[0,0,0]`, new fields `estimator: "agnostic"`, `sound_speed_role: "runtime_environment"`, `lens_mode`/`resolution`/`camera_files` (rig identity), fold stability and angular metrics.
- **Comparison report** (v1 vs v2 on identical events): median/p80/p90 angular error, % within 5°/10°, coverage (valid / in-frame), pixel error for in-frame events only, parameter spread across folds. Written to `calibration_runs/clicker_20260209/constrained_v2/` and summarised in `docs/clicker_calibration.md`.

### Step 2 — Camera intrinsics
_Pending the data inventory (Haiku agent) of `/media/spin/anas/data/calf`; updated below once it reports._

### Step 3 — v2 bridge, explicit-device runner, site sound speed
New file `src/clicker_calibration/preset_bridge_2.jl` (included right after `preset_bridge.jl` in `config.jl`); existing setters and `with_clicker_calibration` unchanged.
- `load_calibration_2(path; device)` — validates schema 2.0 (same checks as `workflow.jl:load_calibration` plus fx = fy, b = 0).
- `with_clicker_calibration_2(f, device::Symbol; c_m_s=SITE_C, detection=:dolphin, calibration_path)` — snapshot/restore every changed global (reuse `_clicker_global_snapshot`/`_clicker_changed_names`); applies `set_device__ophk_acoustic_D3()` (dolphin detection) or the clicker preset when `detection=:clicker`, then installs the profile with `sound_speed_m_s = c_m_s`, `imsize`, `fov_angle`, `FILE_device_ID` (`…_A2`/`_B2`); sets `c`/`default_c = c_m_s` inside the block only (so `tdoa2dir_ls` paths agree); keeps the preset's `default_getTDOA_func` (estimator-agnostic) with an optional `tdoa_estimator=` (e.g. a dolphin-window GCC built on `prep_gcc`, `prepare.jl:140`). The existing `active_clicker_calibration()`/`tdoa2dir`/`angle2px` dispatch then applies v2 unchanged, because `calibrated_tdoa2dir`/`calibrated_angle2px` only read JSON fields.
- `run_contiguous_folders_calibrated(folder, device; res_dir, c_m_s, kwargs...)` — wraps `run_contiguous_folders` in that block; asserts the calibration is active before and after; refuses a video whose resolution differs from the profile (via `media_info.jl`); writes `calibration_used.json` (device, JSON sha256, c, preset, git commit) into `res_dir`.

### Step 4 — Pipeline consistency and invalid-pixel safety
- With b = 0 the tonal, Impulse2–4 and multitdoa paths become consistent automatically; document which outputs are "calibrated" (Impulse1, tonal/Impulse2–4 under v2).
- Check the overlay renderer's handling of NaN rows; if it draws or crashes, add a new post-step writing `detection_pixels_valid.csv` (drop non-finite rows, add a `status` column in a sidecar) and point the calibrated runner's overlay at it.
- Check frame-origin/rounding consistency (I7) with a known synthetic impulse.

### Step 5 — Drift-aware audio–video sync for the calibrated runner
New function reusing `prep_synchronize`/`prep_wave_match`/`prep_fit_sync` (`prepare.jl:356/239/265`) to fit offset + drift per folder; writes `sync_affine.json`; the calibrated runner maps frames with `(t·(1+drift)+offset)·fps`. Falls back to the old constant delay with a logged warning if too few anchors (`docs/error.log` per the long-task skill for full-folder runs).

### Step 6 — Tests (new files, existing suites unchanged)
- `test/clicker_calibration/test_constrained_fit.jl`: synthetic rays with known R, f, c + TDOA noise → R recovered without bias, f bias bounded, invalid rows counted, fixed parameters really fixed.
- `test/clicker_calibration/test_presets_2.jl`: v2 loader rejects fx≠fy / b≠0; scoped activation with dolphin detection; c override; estimator-agnostic dispatch; restoration on exception; video-resolution refusal; provenance file; NaN rows excluded from the overlay CSV.

### Step 7 — Real-recording rollout per device
For one short real folder per rig: `run_contiguous_folders_calibrated(folder, :A)` / `(folder, :B)`; compare against the legacy run: % valid TDOAs, % in-frame, direction stability along click trains, visual review of the overlay. If the old D3 human-labelled `Calibration_05_11_2025` (165 labels) is confirmed to be a different rig, keep it out of A/B validation.

## Future data that would help most (not this round)
Per rig: clicker moved over a 3×3 image grid + corners at 1, 2 and 4 m, holding still for several clicks, with silent pauses; a separate never-inspected validation take; human-confirmed positions or an LED on the clicker; water temperature/salinity; tape-measured camera-to-array offset; underwater checkerboard in the deployment lens mode.

## Verification
- `julia --project -e 'include("test/clicker_calibration/test_constrained_fit.jl"); include("test/clicker_calibration/test_presets_2.jl")'` plus the existing `test_presets.jl` (v1 unchanged).
- Run the Step 1 refit for A and B (≈1–2 min each on cached delays); expect B v2 with fy = fx and |R| within a few degrees of A unless the mounting really differs; compare the v1/v2 angular-error table.
- Step 7 runs on one real folder per rig; inspect the overlay video, `calibration_used.json`, `sync_affine.json`, and `detection_pixels_valid.csv` (no NaN rows).
