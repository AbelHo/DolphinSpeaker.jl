# Better TDOA estimation for dolphin-click localization (ophk_acoustic_D3)

## Context
`run_syncVid.jl` → `run_contiguous_folders` turns each detected click into an image pixel. We want TDOA estimates that are more accurate and more robust than the current ones. They will be new functions in `src/localization.jl`, compared on real dolphin recordings (`Two-way_2males_06_10_2025/New_device`). There is no ground truth yet. The annotated clicker calibration data comes later (see the end of this plan).

## Trace (as called from run_syncVid.jl)
1. `set_device__ophk_acoustic_D3()` ([config.jl:166](src/config.jl#L166)) first calls `set_device__ophk_acoustic_D3_clicker`, then sets:
   - Array (`rx_vect`): 3 hydrophones in a **planar equilateral triangle**, side 0.18 m, z = 0. Columns 1–2 are NaN, so `get_relevant_channels` keeps 3 of the 5 channels.
   - `ref_channel = 3`, `window_impulsive = -50:50`, `impulsive_band_pass = [5000, Inf]`.
   - The WAVs are **192 kHz, 5 ch, 24-bit** (checked with ffprobe on `251006_001_0001.WAV`).
2. `run_contiguous_folders` ([run_example.jl:374](src/run_example.jl#L374)):
   - `run_analysis_split_vidau` (A/V sync).
   - `process_detections` per audio file (threaded; a throw loses that whole file's results).
   - Writes `detection_pixels.csv` from `results[i][1]` (impulse) and `[2]` (tonal), because `detection_types = 1:2`.
3. `process_detections` ([aspod.jl:165](src/aspod.jl#L165)):
   - `readAudio` → `detect_impulseNtonal(...; return_datafilt=true)`. This gives `res.res_impulse.data_filt` (band-passed, 3 channels) and `res.res_impulsetrain.pind_good`, the Hilbert-envelope peak on `ref_channel`. Clicks outside trains are already dropped.
   - Then `detection2angle` → `getTDOA_func` → `tdoa2dir` (NelderMead on `cost_tdoa2ang`) → `angle2px`.

| tag | getTDOA_func | window | where it ends up |
|---|---|---|---|
| Impulse | `get_tdoa_minmax` (default) | -50:50 | **detection_pixels.csv → video overlay** |
| Impulse2 | `get_tdoa_raw_MaxEnergyRefChannel` | -130:130 | own CSV only |
| Impulse3 | `get_tdoa_envelope_filtered` | -50:50 | own CSV only |
| Impulse4 | `get_tdoa_raw_MaxPeakRefChannel` | -50:50 | own CSV only |

## Why the current estimates are limited
- **Integer-sample resolution.** The largest possible TDOA per pair is 0.18 m / 1540 m/s × 192 kHz ≈ **22.4 samples**.
  - One sample step ≈ **2.5° on axis ≈ 87 px horizontally** (2740 px over 80° FOV), while the overlay radius is 20 px. Detections snap to a visible grid.
  - All 4 pipeline methods return `Int`. The only sub-sample attempt, `get_tdoa_raw_MaxEnergyRefChannel_resample` (×10 resampling), is unused.
- **The overlay method is plain argmax.** In `get_tdoa_minmax` the polarity vote tests `max_val > min_val`, which is always true, so the min branch never runs. The result is the argmax of the signed waveform per channel.
- **Cycle skipping.** Clicks band-limited to under 96 kHz (Nyquist) have carrier periods of about 2–6 samples. Peak picking or raw xcorr can jump a whole cycle between channels.
- **No physical lag limit.** `finddelay` searches every lag in the window, up to ±130 samples, so multipath or noise can produce |τ| ≫ 22.
- **Only 2 of the 3 pairs are used** (each against the reference). The third pair gives a free consistency check (closure τ12+τ23+τ31 = 0).
- **Other bugs:**
  - `envelope_lpfiltered` hard-codes `fs=500_000`, so its 10 kHz low-pass actually cuts at about 3.8 kHz.
  - 3-arg aliases are miswired: `get_tdoa_minmax(…,ref)` calls `get_tdoa_max`, and `get_tdoa_envelope_filtered(…,ref)` calls `get_tdoa_envelope`.
  - `get_tdoa_minmax` and the envelope methods do not clip the window at file edges, so a click within 50 samples of either end throws.
- **The solver is overkill.** With z = 0, `cost_tdoa2ang` is *linear* in u = sin(az)·cos(el) and v = sin(el): Δτ_j = −(Δx_j·u + Δy_j·v)/c·fs.
  - With 3 sensors and one reference the system is exactly determined. So today's residual is non-zero only when the TDOAs are infeasible (u²+v² > 1).

## Proposed new functions (in src/localization.jl)
Every TDOA function keeps the signature `f(data, peaks; window, ref_channel)`. Each returns an N×m `Float64` matrix with the existing sign convention: **positive = channel arrives later than the reference**, as with `finddelay(data[:,ch], data[:,ref])`. `tdoa2dir` / `cost_tdoa2ang` therefore accept them unchanged. `return_quality=false` can also return per-click quality: normalized peak height and closure error.

Shared helpers:
- `_clip_window(win, n)`: replaces the repeated edge-clipping code.
- `tau_max_samples(rx_vect, fs; c, margin=2)`: per-pair physical lag limit.
- `subsample_peak(s, k)`: parabolic peak interpolation.
- `gcc(x, y; weighting=:phat|:scot|:none, band, fs, interp=8, maxlag)`: FFT cross-spectrum (FFTW is already a dependency).
  - Band mask, so PHAT does not whiten out-of-band noise.
  - Zero-padded IFFT for upsampled correlation.
  - Searches only |lag| ≤ maxlag, then refines with `subsample_peak`.
- `pairs_to_delays(τpair)`: least-squares projection of all m(m−1)/2 pairwise delays to per-channel delays, t_j = (1/m)·Σ_k τ_jk, then re-referenced. Its residual is the closure error.
- Optional Tukey taper on the snippet, to reduce edge effects and late reflections.

Already tried, so not proposed again:
- sub-sample xcorr (= `get_tdoa_raw_MaxEnergyRefChannel_resample`)
- envelope xcorr (= `get_tdoa_envelope`, `get_tdoa_envelope_filtered`)
- peak picking (`get_tdoa_max/min/minmax`), template xcorr (`_MaxPeakRefChannel`), threshold onset (`get_tdoa_findTrigger*`)

New methods (none of the existing functions does any of these):
1. **`get_tdoa_gcc`**: generalized cross-correlation, per pair.
   - What is new: frequency weighting, `weighting = :phat` (default), `:scot` or `:ml`. `:ml` weights each frequency by its SNR, using a noise spectrum taken just before the click.
   - Band mask (default `impulsive_band_pass`, top capped at 0.98·fs/2; the harness also tries a 20 kHz low edge to keep whistle energy out).
   - Lag limited to ±τmax, upsampled peak with parabolic refinement.
   - `all_pairs=true`: all 3 pairs, then `pairs_to_delays`, with closure as quality.
   - `guide=:envelope` option: when the GCC curve has several near-equal peaks, take the one nearest the envelope-correlation estimate. The envelope is used only as a coarse guide.
2. **`get_tdoa_srp`** (SRP-PHAT): estimates the direction jointly instead of TDOA pair by pair.
   - Grid search over (u,v) in the unit disk, about 1° coarse grid then local refinement.
   - Maximizes the sum of the 3 pairwise GCC-PHAT curves interpolated at the predicted lags.
   - Returns the predicted TDOAs (consistent by construction) and the SRP peak as quality.
   - Most robust against one bad channel or a multipath peak.

Solver:
- **`tdoa2dir_ls(tdoas, rx_vect, fs; c, ref_channel, return_residual, kwargs...)`**: closed-form (u,v), exact for m = 3 and least squares for m > 3.
  - el = asin(v), az = asin(u/cos el) in the forward hemisphere, which is where NelderMead from [0,0] lands.
  - Infeasible points are clamped to the unit circle and flagged.
  - Call signature and return format match `tdoa2dir`. The residual is `cost_tdoa2ang` at the solution, so it compares directly with NelderMead. Non-planar arrays fall back to `tdoa2dir`.
- Add a backward-compatible `tdoa2dir_func=tdoa2dir` kwarg to `detection2angle` ([localization.jl:224](src/localization.jl#L224)).

**Expected outcome:** I expect `get_tdoa_gcc` (PHAT, all pairs) or `get_tdoa_srp` to beat all existing methods. The existing `_resample` row shows how much of that gain is sub-sample resolution alone.

## Comparison harness: new `src/compare_tdoa.jl` (not wired into the pipeline)
- `include("run_example.jl"); set_device__ophk_acoustic_D3()`.
- Load 1–3 WAVs from the New_device folder, then `readAudio` → `detect_impulseNtonal(...; return_datafilt=true)`, exactly as in `process_detections`.
- Methods: existing functions run as-is as baselines (the 4 pipeline methods plus `_resample`; nothing re-implemented). New: `get_tdoa_gcc` (`:phat`, `:scot`, `:ml`, with and without `guide=:envelope`) and `get_tdoa_srp`. Each runs with both `tdoa2dir` and `tdoa2dir_ls`.
- Metrics (no ground truth), overall and by `ppeak` quartile (envelope peak, an SNR proxy):
  - **Closure error** |τ12+τ23+τ31|, computed from independent pairwise estimates.
  - **Feasibility rate**: share of clicks with |τ_ij| ≤ τmax and u²+v² ≤ 1. This also covers the solver residual.
  - **Within-train angular jitter**: median |Δaz| and |Δel| between consecutive clicks of the same train (membership from `train_start`/`train_end`). Medians keep it robust to the two males interleaving clicks.
  - **Method agreement**: median angular difference between each pair of methods.
  - **Runtime per click.**
- **Semi-synthetic benchmark** (known truth):
  1. Take high-SNR real clicks from the reference channel.
  2. Make 3 channels with FFT fractional delays over an (az,el) grid.
  3. Add real noise snippets from click-free stretches at several SNRs.
  4. Report τ and angle RMSE per method.
  - Caveat: every channel gets the same waveform, so this tests precision against noise, not differences in waveform between channels.
- Outputs to a `res_dir`:
  - Per-click CSV.
  - A summary table (printed and saved).
  - Plots: az/el vs time per method, closure and jitter histograms, synthetic error vs SNR.

## Later: clicker calibration data
- Invert `angle2pxd` to turn annotated clicker pixels into reference (az,el), then compute per-method angular error directly.
- Fit these parameters to the same data:
  - array roll/yaw relative to the camera;
  - `fov_angle`;
  - c;
  - constant per-channel delay offsets. If the ADC samples channels sequentially, each channel has a fixed fraction-of-a-sample skew, which matters once TDOAs are sub-sample.

## Files
- Modify `src/localization.jl`:
  - New helpers, methods and `tdoa2dir_ls`.
  - The `tdoa2dir_func` kwarg.
  - Fixes: `envelope_lpfiltered` takes `fs` (changes only the Impulse3 CSV), the 3-arg aliases are corrected, and window clipping is added.
  - `get_tdoa_minmax` is left as is, because it drives the overlay; whether to change it is decided after the comparison.
- New `src/compare_tdoa.jl`.
- `aspod.jl`: unchanged for now. Later the winner becomes the default or an "Impulse5" entry.

## Verification
1. REPL checks:
   - Apply known fractional delays to a real click. Every new estimator should recover them to ≲0.1 sample at high SNR, with the same sign as `get_tdoa_raw`.
   - `pairs_to_delays` gives zero closure on consistent input.
2. `tdoa2dir_ls` matches NelderMead `tdoa2dir` on feasible TDOAs (|Δangle| < 0.1°).
3. Run `compare_tdoa.jl` on `251006_001_0001.WAV` and inspect the summary table and plots.
4. `run_syncVid.jl` still runs, and the video overlay output is unchanged.
