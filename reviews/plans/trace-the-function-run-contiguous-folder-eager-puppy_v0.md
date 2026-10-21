# Better TDOA estimation for dolphin-click localization (ophk_acoustic_D3)

## Context
`run_syncVid.jl` → `run_contiguous_folders` turns each detected click into an image pixel. We want TDOA estimates that are more accurate and more robust than the current ones, built as new functions in `src/localization.jl` and compared on real dolphin recordings. There is no ground truth yet. The annotated clicker calibration data comes later (see the end of this plan).

## Trace (as called from run_syncVid.jl)
1. `set_device__ophk_acoustic_D3()` ([config.jl:166](src/config.jl#L166)) calls `set_device__ophk_acoustic_D3_clicker` first, then sets:
   - Array (`rx_vect`): 3 hydrophones in a **planar equilateral triangle**, side 0.18 m, z = 0. Columns 1–2 are NaN, so `get_relevant_channels` keeps 3 channels.
   - `ref_channel = 3`, `window_impulsive = -50:50`, `impulsive_band_pass = [5000, Inf]`, `dist_impulsive = 77` (comment implies fs = 192 kHz).
2. `run_contiguous_folders` ([run_example.jl:374](src/run_example.jl#L374)) runs `run_analysis_split_vidau` (A/V sync), then `process_detections` per audio file (threaded), then writes `detection_pixels.csv`.
3. `process_detections` ([aspod.jl:165](src/aspod.jl#L165)):
   - `readAudio` → `detect_impulseNtonal(...; return_datafilt=true)`.
   - That gives `res.res_impulse.data_filt` (band-passed, relevant channels) and `res.res_impulsetrain.pind_good` (peak index of the Hilbert envelope on `ref_channel`).
   - Then `detection2angle` → `getTDOA_func` → `tdoa2dir` (NelderMead on `cost_tdoa2ang`) → `angle2px`.
   - It runs 4 impulse variants:

| tag | getTDOA_func | window |
|---|---|---|
| Impulse | `get_tdoa_minmax` (default) | -50:50 |
| Impulse2 (used for colouring) | `get_tdoa_raw_MaxEnergyRefChannel` | -130:130 |
| Impulse3 | `get_tdoa_envelope_filtered` | -50:50 |
| Impulse4 | `get_tdoa_raw_MaxPeakRefChannel` | -50:50 |

## Why the current estimates are limited
- **Integer-sample resolution.** The largest possible TDOA is 0.18 m / 1540 m/s × 192 kHz ≈ **22 samples**. Rounding to ±0.5 sample is about 2–3° on axis and worse off axis. All 4 methods return `Int`.
- **Cycle skipping.** Band-passed dolphin clicks (roughly 60–130 kHz) have a period of about 1.5–3 samples. Picking peaks (`minmax`) or raw xcorr can jump a whole cycle between channels. `minmax` also picks the max or min polarity by majority vote.
- **No physical lag limit.** `finddelay` searches every lag in the window, up to ±130 samples. Multipath or noise can produce |τ| ≫ 22, which no direction can explain.
- **Only 2 of the 3 pairs are used** (each channel against the reference). The third pair would give a free consistency check (τ12+τ23+τ31 = 0).
- **Bugs to fix along the way:**
  - `envelope_lpfiltered` hard-codes `fs=500_000`, so the low-pass cutoff is wrong at 192 kHz.
  - The 3-arg aliases are miswired: `get_tdoa_minmax(…,ref)` calls `get_tdoa_max`, and `get_tdoa_envelope_filtered(…,ref)` calls `get_tdoa_envelope`.
  - The envelope methods do not clip the window at file edges.
- **The solver is overkill.** For a planar array (z = 0), `cost_tdoa2ang` is *linear* in u = sin(az)·cos(el) and v = sin(el): Δτ_j = −(Δx_j·u + Δy_j·v)/c·fs. A closed-form least-squares solve is exact, fast and deterministic. It also flags infeasible TDOAs when u²+v² > 1.

## Proposed new functions (in src/localization.jl)
All TDOA functions keep the existing signature `f(data, peaks; window, ref_channel)` and return an N×m `Float64` matrix. Convention: **positive = channel arrives later than the reference**, the same as `finddelay(data[:,ch], data[:,ref])`. The existing `tdoa2dir` / `cost_tdoa2ang` accept them unchanged. A `return_quality=false` kwarg can also return a per-click quality matrix (normalized peak height and closure error).

Shared helpers:
- `_clip_window(win, n)`: factor out the repeated edge-clipping code.
- `tau_max_samples(rx_vect, fs; c, margin=2)`: per-pair physical lag limit.
- `subsample_peak(s, k)`: parabolic interpolation of a correlation peak.
- `gcc(x, y; weighting=:phat|:scot|:none, band=nothing, fs, interp=8, maxlag)`: FFT cross-spectrum (FFTW is already a dependency).
  - Optional band mask, so PHAT does not whiten out-of-band noise.
  - Zero-padded IFFT gives upsampled correlation.
  - Searches only |lag| ≤ maxlag, then refines with `subsample_peak`.
  - Returns (τ, peak value).
- Optional Tukey taper on the snippet, to cut edge effects and later reflections.

Methods to compare:
1. **`get_tdoa_xcorr_subsample`**: today's raw xcorr with the MaxEnergy reference, plus the lag limit and parabolic sub-sample refinement. This is the "minimal fix" baseline and shows how much comes from sub-sample resolution alone.
2. **`get_tdoa_gcc_phat`**: band-limited GCC-PHAT (default band from `impulsive_band_pass`, top capped at 0.98·fs/2), with upsampled correlation and the lag limit. This is the standard robust TDOA estimator for broadband transients.
3. **`get_tdoa_envelope_subsample`**: Hilbert-envelope xcorr with the fs bug fixed, the lag limit and sub-sample refinement. It cannot cycle-skip but is less precise.
4. **`get_tdoa_hybrid`**: a coarse τ from the envelope (method 3), then the GCC-PHAT peak (method 2) nearest to it within ±1 carrier period. This aims for envelope robustness with phase precision.
5. **`get_tdoa_srp`** (joint, all 3 pairs): compute the 3 pairwise GCC-PHAT curves, then grid-search (u,v) over the unit disk. The coarse grid is about 1°, refined locally. It maximizes the sum of the GCC curves interpolated at the predicted lags. It returns the TDOAs predicted at the optimum (consistent by construction) and the SRP peak as quality. This is the most robust option against one bad channel or a multipath peak.

Solver additions:
- **`tdoa2dir_ls(tdoas, rx_vect, fs; c, ref_channel, use_all_pairs=true)`**: closed-form least squares for planar arrays.
  - Solves (u,v), then el = asin(v) and az = asin(u/cos el), using the forward hemisphere as in the current solver.
  - Returns the same `mod_sign` angle format, plus the residual and an infeasibility flag.
  - Non-planar arrays fall back to `tdoa2dir`.
- Add a backward-compatible `tdoa2dir_func=tdoa2dir` kwarg to `detection2angle` ([localization.jl:224](src/localization.jl#L224)).

## Comparison harness: new `src/compare_tdoa.jl` (not wired into the pipeline yet)
- `include("run_example.jl"); set_device__ophk_acoustic_D3()`.
- Take 1–3 audio files from the run_syncVid folder (`extract_vid_au`). Run `readAudio` → `detect_impulseNtonal(...; return_datafilt=true)`, exactly as `process_detections` does.
- Run every method (old 4 + new 5) × solver (`tdoa2dir` NelderMead, `tdoa2dir_ls`) on `res.res_impulse.data_filt`, `res.res_impulsetrain.pind_good`.
- Metrics, all without ground truth. Each is reported overall and by SNR quartile (`ppeak`):
  - **Closure error** |τ12+τ23+τ31| (samples), for the pairwise methods.
  - **Feasibility rate**: share of clicks with |τ_ij| ≤ τmax_ij and u²+v² ≤ 1.
  - **Solver residual.**
  - **Within-train angular jitter**: median |Δaz| and |Δel| between consecutive clicks with ICI < `click_train_check_interval`. Lower is better, since an animal barely moves within a few ms.
  - **Integer fraction**: share of TDOAs on integer values (sanity check).
  - **Pairwise method agreement**: median angular difference between methods.
  - **Runtime per click.**
- **Semi-synthetic test** (known truth from real data):
  1. Take high-SNR clicks from the reference channel.
  2. Make 3 channels with FFT fractional delays for a grid of (az,el).
  3. Add real noise snippets from click-free parts of the same file, at several SNRs.
  4. Report RMSE of τ and of angle per method.
- Outputs to a `res_dir`:
  - CSV of per-click results.
  - A summary table printed and saved.
  - Plots: az/el vs time per method; closure and jitter histograms; error vs SNR for the synthetic test.

## Later: using the clicker calibration data
- Invert `angle2pxd` to turn annotated clicker pixels into reference (az,el). Then compute per-method angular error directly.
- Fit the extrinsic and global parameters to that data: array roll/yaw relative to the camera, `fov_angle`, sound speed c, and a constant per-channel delay offset. Per-channel offsets from ADC or cable skew show up as a constant TDOA bias, and calibration can remove them.

## Files
- Modify: `src/localization.jl`.
  - New functions and helpers.
  - `tdoa2dir_func` kwarg.
  - Bug fixes: `envelope_lpfiltered` takes `fs`, and the miswired aliases are fixed.
- New: `src/compare_tdoa.jl` (harness).
- `aspod.jl`: no change for now. After the comparison, the best method can become `default_getTDOA_func` or an "Impulse5" entry.

## Verification
1. Unit sanity in the REPL:
   - Delay a real click by known fractional delays. Each new estimator should recover them to ≲0.1 sample at high SNR.
   - Its sign must match `get_tdoa_raw`.
2. `tdoa2dir_ls` should match NelderMead `tdoa2dir` on feasible TDOAs (|Δangle| < 0.1°).
3. Run `julia --project src/compare_tdoa.jl` on one file from `Two-way_2males_06_10_2025/New_device`. Check the summary table and plots.
4. Confirm that `run_syncVid.jl` still runs unchanged (the defaults are untouched).
