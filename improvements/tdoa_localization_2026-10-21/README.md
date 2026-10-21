# Better TDOA estimation for click localization (ophk_acoustic_D3), 2026-10-21

This folder records one working session on how `run_contiguous_folders` turns dolphin clicks into video pixels: the trace, the new TDOA methods, the comparison against the existing ones, what we learned, and what to do next.

## Summary

- The new TDOA methods (GCC with frequency weighting, SRP-PHAT) are about **10× more precise** than the existing integer methods on clicks with a known answer. On the strongest real clicks, consecutive clicks agree within **about 0.6°**, against 3.7° for the current overlay method.
- The **current video overlay is mostly unusable**. Its method (`get_tdoa_minmax`) produces physically impossible TDOAs for about half the clicks, and **80% of the overlay points land outside the frame**.
- On real data, the **limit is not sub-sample precision** but weak detections: false detections, echoes, and overlapping clicks give the three hydrophones different signals. No method fixes those. The new methods' quality outputs can flag them.
- **Changing the solver** (NelderMead → closed form) gives the same directions. It is faster and simpler, not more accurate.

## What was done

1. **Traced the pipeline** as called from `src/run_syncVid.jl`:
   `set_device__ophk_acoustic_D3` → `run_contiguous_folders` (src/run_example.jl) → `run_analysis_split_vidau` (A/V sync) → `process_detections` (src/aspod.jl) → `detect_impulse` + `detect_impulsetrain` → `detection2angle` → `getTDOA_func` → `tdoa2dir` (NelderMead) → `angle2px`.
   `process_detections` runs four impulse TDOA methods. Only the first, `get_tdoa_minmax`, goes into `detection_pixels.csv` and the video overlay; the others are written to their own CSVs.

   | tag | getTDOA_func | window | used for |
   |---|---|---|---|
   | Impulse | `get_tdoa_minmax` (default) | -50:50 | overlay |
   | Impulse2 | `get_tdoa_raw_MaxEnergyRefChannel` | -130:130 | own CSV |
   | Impulse3 | `get_tdoa_envelope_filtered` | -50:50 | own CSV |
   | Impulse4 | `get_tdoa_raw_MaxPeakRefChannel` | -50:50 | own CSV |

2. **Added new functions to `src/localization.jl`.** These are additions only (258 lines inserted before `default_getTDOA_func`); no existing line was changed. All TDOA functions keep the `f(data, peaks; window, ref_channel, ...)` signature and the existing sign convention (positive = later than `ref_channel`). They also take `fs` as a required keyword, because the pipeline has no global `fs`.
   - `get_tdoa_gcc`: generalized cross-correlation.
     - Weighting: `:phat` (default), `:scot` (smoothed spectra), `:ml` (per-frequency SNR weights, noise taken before the click) or `:none`.
     - Band mask, a lag limit at what the array geometry allows, and an upsampled peak with parabolic refinement.
     - `all_pairs=true` estimates all three pairs and reconciles them by least squares (`pairs_to_delays`). `guide=:envelope` uses the envelope correlation to choose between near-equal peaks.
     - `return_quality=true` returns the peak heights and the closure.
   - `get_tdoa_srp`: SRP-PHAT. It searches directions (u, v) for the one where the three GCC curves agree best, then returns the TDOAs predicted there; these are consistent and feasible by construction.
   - `tdoa2dir_ls`, `tdoa2uv`: closed-form direction solver for the planar array. It has the same call and output as `tdoa2dir`; the front hemisphere is assumed.
   - Fixed copies of buggy existing functions, with the originals left untouched:
     - `get_tdoa_minmax2`: the polarity vote compares |max| with |min|.
     - `get_tdoa_envelope_filtered2`: passes the real `fs`.
   - Helpers: `gcc_curve`, `gcc_delay`, `gcc_weight`, `noise_psd`, `subsample_peak`, `curve_at`, `pairs_to_delays`, `tau_max_samples`, `_clip_window`, `_rx_for`, `default_gcc_band`.

3. **Comparison harness `src/compare_tdoa.jl`.**
   - It runs detection exactly as the pipeline does, then every method with both solvers, and computes metrics on real clicks plus a semi-synthetic benchmark. It writes a report, CSVs and plots.
   - `replot_comparison` rebuilds the plots from the CSVs.
   - `export_webui=true` writes one `detection_pixels__color255_webui__<file>__<method>.csv` per method, in the pipeline's webui format, for overlaying on the combined video.

4. **Unit tests `test/test_tdoa.jl`** (24 tests, synthetic data only).

5. **Plans:** `reviews/plans/trace-the-function-run-contiguous-folder-eager-puppy.md` (approved) and `..._v0.md` (first draft).

## How to reproduce

All commands run from the repository root. The data are on `/media/spin/anas2`.

```bash
# unit tests (about 10 s, no recordings needed)
julia --project=. -t auto test/test_tdoa.jl

# full comparison on 251006_001_0001.WAV (about 6 min, about 10 GB RAM); writes a new
# /media/spin/anas2/data/calf/upload/results/compare_tdoa_<timestamp>/ folder
julia --project=. -t 16 src/compare_tdoa.jl            # add --webui for the per-method overlay CSVs
julia --project=. -t 16 src/compare_tdoa.jl a.WAV b.WAV
```

From the REPL:

```julia
include("src/compare_tdoa.jl")
compare_tdoa([joinpath(COMPARE_TDOA_FOLDER, "251006_001_0001.WAV")];
    res_dir="/your/out", run_synthetic=true,
    export_webui=true, webui_solver=:nm,                       # :ls = closed-form solver
    webui_kwargs=(; sync_dir="/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_nooverlay2/New_device"))
replot_comparison("/your/out", "251006_001_0001")
```

For the webui CSVs, the A/V sync delay is taken from `sync_delay_s`, or else from `sync_dir/sync_delay.csv`, or else computed the way the pipeline does.

The diagnostic scripts in `scripts/` are one-off analyses that take a results folder:

| script | what it shows |
|---|---|
| `check_pipeline_repro.jl <res_dir>` | the harness's `minmax` pixels against the pipeline's own `Impulse` CSV |
| `diag_channels.jl <res_dir>` | SNR of all 5 channels at the clicks, waveform similarity between hydrophones, example clicks |
| `analyze_per_click.jl <per_click.csv>` | jitter, spread, repeated estimates and closure per method |
| `subset_metrics.jl <per_click.csv>` | metrics on the strongest detections and on closure-gated clicks |

`results/` holds the report, summary CSVs and figures of the final run (`compare_tdoa_20261021_170157`). The per-click CSV (about 250 MB) stays on the data drive.

## Results

Data: `Two-way_2males_06_10_2025/New_device/251006_001_0001.WAV`, a Zoom F6 recording at 192 kHz, 5 ch, 24-bit. Detection found 56,447 clicks in 153 trains, with auto threshold 490.376, the same as the pipeline. The real clicks' window SNR is about 15 / 19 / 29 dB (p10 / median / p90).

**Checks**
- The unit tests pass, 24 of 24.
- The harness reproduces the pipeline: all 56,447 overlay pixels match with 0 px difference.
- With `--webui`, the `minmax` file matches the pipeline's webui CSV: frame, px, py, radius and alpha are identical, and the colours differ by at most 2e-5 out of 255.

**Semi-synthetic benchmark.** Real clicks, delayed by the exact TDOAs of known directions inside the field of view and added to real noise. Values are the median angle error:

| method | 0 dB | 10 dB | 15 dB | 30 dB |
|---|---|---|---|---|
| minmax [Impulse, overlay] | 1.48° | 1.37° | 1.36° | 1.36° |
| xcorr MaxEnergyRef [Impulse2] | 1.77° | 1.21° | 1.20° | 1.19° |
| envelope_filtered [Impulse3] | 15.99° | 3.35° | 2.02° | 1.26° |
| xcorr MaxEnergyRef resample x10 | 0.29° | 0.12° | 0.11° | 0.11° |
| gcc phat (+guide is the same) | 0.35° | 0.16° | 0.12° | 0.08° |
| gcc ml | 0.23° | 0.06° | 0.03° | 0.01° |
| srp phat | 0.37° | 0.18° | 0.15° | 0.12° |

This is a best case: every channel gets the same waveform, so it measures precision against noise but not robustness to echoes.

**Real clicks** (no ground truth). Definitions:
- **Jitter:** the median angle between consecutive clicks of a train.
- **Ratio:** jitter divided by the same quantity for random click pairs. Lower is better, and 1 means consecutive clicks agree no better than random ones.
- **Feasible:** TDOAs that some plane wave could produce.

| method | feasible, all clicks | ratio, all | ratio, top 25% | strongest 10%: jitter / feasible |
|---|---|---|---|---|
| minmax [Impulse, overlay] | 52% | 0.84 | 0.42 | 3.7° / 77% |
| xcorr MaxEnergyRef [Impulse2] | 53% | 0.74 | 0.25 | 0.0°\* / 80% |
| envelope_filtered [Impulse3] | 93% | 0.58 | 0.30 | 4.0° / 85% |
| envelope_filtered2 (fs fixed) | 91% | 0.68 | 0.23 | 3.0° / 85% |
| gcc phat | 96% | 0.78 | 0.40 | 0.65° / 96% |
| gcc phat +guide | 95% | 0.74 | 0.29 | 0.63° / 96% |
| srp phat | 100% | 0.78 | 0.28 | 0.71° / 100% |

\* Integer TDOAs repeat between consecutive clicks; the 0.0° is rounding, not accuracy.

The full table (18 methods) is in `results/report_251006_001_0001.txt`. The figures are:
- `results/metrics_*.png`
- `results/synthetic_*.png`
- `results/closure_*.png`
- `results/az_vs_time_*.png`
- `results/example_clicks.png`

## Observations

- **`get_tdoa_minmax` bug.** The polarity vote tests `max_val > min_val`, which is always true, so the method is plain argmax of the signed waveform.
- **Infeasible TDOAs.** About 48% of clicks have impossible TDOAs under the current overlay method; they are drawn at az or el = ±90° (px ≈ −1712 or 4452), off-screen. That accounts for most of the 80% of off-frame overlay points.
- **`envelope_lpfiltered` assumes 500 kHz.** At 192 kHz its "10 kHz" low-pass cuts at about 3.8 kHz. It is correct for the 500 kHz D1 device.
- **Integer xcorr cycle-skips.** `finddelay` maximizes |xcorr|, so on band-limited clicks it can lock onto a negative lobe. In synthetic tests, 82 of 200 pair delays from `get_tdoa_raw` were off by half a 60 kHz period.
- **NelderMead mirrors some directions.** It sometimes reports the back-hemisphere angles (e.g. az 121° instead of 59°). A planar 3-hydrophone array cannot tell front from back, so both are equally valid answers; the closed-form solver assumes front. On the real data, switching solvers moves 0% of clicks into the frame, because all mirrored cases are outside the field of view.
- **Waveforms differ between hydrophones.** The median normalized cross-correlation between hydrophone pairs is only about 0.5. Below 0.5, GCC closure is 6–10 samples; above it, under 0.4. Many detection windows contain several arrivals per channel, 100–200 µs apart: echoes or overlapping clicks from the two males (`results/example_clicks.png`). A tighter -30:30 window did not help.
- **False detections** (pointed out during the session) also produce inconsistent arrivals. The closure < 1 sample test passes 64% of the strongest quarter of detections but only 23% of the weakest, so closure works as a true-click indicator.
- **Channels 1–2** are bit-identical and correlate about 1/√3 with each of channels 3–5. They are the Zoom F6 stereo mixdown, so marking them NaN in `rx_vect` is correct.
- **Weighting on real data.** `gcc ml` is best on synthetic clicks but not on real ones; its noise estimate before the click is probably contaminated by other arrivals. The envelope-guided variants and SRP are the most robust of the new methods on real clicks.

## Useful points from the discussion

- **Methods already tried.** The first plan included sub-sample xcorr and sub-sample envelope methods. They were dropped because they duplicate `get_tdoa_raw_MaxEnergyRefChannel_resample` and `get_tdoa_envelope(_filtered)`, which had already been tried. `get_tdoa_gcc(...; weighting=:none)` is kept only as an ablation row, to separate the effect of weighting from the effect of sub-sample refinement and the lag limit.
- **Bugs in existing functions** are handled by adding fixed copies (`*2`); the originals are left untouched because they drive the current pipeline.
- **Ground truth.** There was none for this session. The clicker calibration data (clicker pixel annotated) is the planned ground truth.
- **NelderMead vs closed form.** For three sensors, the direction is fully determined by the TDOAs, so any correct solver gives the same answer for feasible TDOAs. The solver was never the bottleneck.

## Recommendations

**Fixes**
1. **Stop drawing infeasible clicks.** Mark clicks with infeasible TDOAs (u² + v² > 1) as "no direction" instead of placing them at ±90°.
2. **Switch the overlay TDOA method** away from `get_tdoa_minmax` to `get_tdoa_srp` or `get_tdoa_gcc(...; guide=:envelope)`, with `tdoa2dir_ls`. For example, add an "Impulse5" entry in `process_detections` by passing a closure as `getTDOA_func`: `(d, p; kw...) -> get_tdoa_srp(d, p; fs=fs, kw...)`.
3. **Gate by quality.** Only draw clicks with GCC closure < about 1 sample (or a high SRP peak). This also removes many false detections. Tune the threshold on the calibration data.
4. **Decide whether to fix the originals** (`get_tdoa_minmax` vote, `envelope_lpfiltered` fs, the miswired 3-argument aliases, missing edge clipping) or retire them in favour of the `*2` copies.

**Improvements**
5. **Calibrate with the clicker data.**
   - Convert the annotated clicker pixels back to angles to get per-method angular error. This is the only way to measure accuracy rather than consistency.
   - Fit the array roll/yaw relative to the camera, `fov_angle`, the sound speed `c` and constant per-channel delay offsets. These bias every method equally.
6. **Handle multiple arrivals.** Possible approaches:
   - choose the earliest consistent arrival (direct path) among the GCC peaks;
   - accumulate SRP maps over consecutive clicks of a train, since an animal barely moves within milliseconds;
   - cluster the directions within a train, so the two males can be separated.
7. **Improve detection.** Use the closure / waveform similarity between hydrophones as a detection feature to reject false detections earlier.
8. **Front/back ambiguity** can only be resolved with a fourth, out-of-plane hydrophone.
9. **Run the comparison on files 0002 and 0003** (and other sessions); these conclusions are based on one file.
10. **Test infrastructure.** `test/runtests.jl` is still the PkgTemplates suite (see `reviews/repository_audit.md`). Replace it with a DolphinSpeaker runner that includes `test/test_*.jl`, so `test_tdoa.jl` runs under `Pkg.test()`.

## Follow-up 2026-10-22: one webui overlay CSV per TDOA method from the real pipeline

Plan: `reviews/plans/multitdoa-webui-csvs_2026-10-22.md` (first draft: `..._v0.md`). The aim is to compare the TDOA methods in the real-time overlay viewer, without rendering video.

**What was added** (no existing function was changed)
- `src/run_multitdoa.jl`:
  - `run_contiguous_folders_multitdoa(folder; res_dir, closure_max=1.0)` has the same flow as `run_contiguous_folders` (A/V sync, combined video, threaded detection, audacity labels) but writes one webui CSV per method and no overlay video.
  - `process_detections_multitdoa` calls the unchanged `process_detections` for Impulse1–4 and computes Impulse5/6 from the same filtered data.
  - `impulse_variant_rows` and `flag_clicks` are the pure helpers.
- `src/run_syncVid.jl` now includes it and calls the new runner on `New_device` (files 0001–0003 as one set). `res_dir` ends in `_multitdoa`, so the earlier results are not overwritten.
- `test/test_multitdoa.jl` (14 tests, synthetic only) and `scripts/check_multitdoa_outputs.jl`.

**Output files** in `<res_dir>/<folder name>/`. Each one is a drop-in replacement for `detection_pixels__color255_webui.csv` (that method's impulse rows plus the tonal rows):

| file suffix | method | notes |
|---|---|---|
| `Impulse1_minmax` | `get_tdoa_minmax` (current overlay) | NelderMead |
| `Impulse2_xcorrMaxEnergy` | `get_tdoa_raw_MaxEnergyRefChannel` | window -130:130 |
| `Impulse3_envelopeFiltered` | `get_tdoa_envelope_filtered` | |
| `Impulse4_xcorrMaxPeak` | `get_tdoa_raw_MaxPeakRefChannel` | |
| `Impulse5_gccPhatGuided` | `get_tdoa_gcc`, PHAT, `guide=:envelope` | `tdoa2dir_ls` |
| `Impulse6_srpPhat` | `get_tdoa_srp` | `tdoa2dir_ls` |
| `Impulse5g_…-gated`, `Impulse6g_…-gated` | same clicks as 5 and 6, flagged ones dropped | |

In the ungated Impulse5/6 files a click flagged as probably wrong has `shape = type = "triangle"`. A click is flagged when its GCC closure |τ12+τ23+τ31| exceeds `closure_max` (1 sample), or when its TDOAs are infeasible (u² + v² > 1), or when either is NaN. Impulse6 takes the closure from the Impulse5 GCC pass, because SRP's own TDOAs are consistent by construction. `multitdoa_summary.csv` lists clicks, flagged and inside-frame counts per file.

**Run it**
```bash
julia --project=. -t 12 src/run_syncVid.jl          # edit folder / res_dir / closure_max in the script
julia --project=. improvements/tdoa_localization_2026-10-21/scripts/check_multitdoa_outputs.jl <res_dir>/New_device [<earlier run>/detection_pixels__color255_webui.csv]
julia --project=. -t auto test/test_multitdoa.jl
```

**Results**, `New_device` 0001–0003 (135,200 clicks; share of points inside the 2740 × 1520 frame):

| file | flagged | inside frame |
|---|---|---|
| Impulse1 minmax | – | 21.9% |
| Impulse2 xcorrMaxEnergy | – | 22.1% |
| Impulse3 envelopeFiltered | – | 74.1% |
| Impulse4 xcorrMaxPeak | – | 24.9% |
| Impulse5 gccPhatGuided | 82,162 (61%) | 55.6% |
| Impulse6 srpPhat | 80,160 (59%) | 57.4% |
| Impulse5g gated | 0 (53,038 clicks left) | 51.8% |
| Impulse6g gated | 0 (55,040 clicks left) | 52.1% |

- **Checks passed.** The standard `detection_pixels__color255_webui.csv` is identical to the earlier pipeline run (0 difference in every column). The Impulse1 file equals the standard file, with colours within 7.6e-6 out of 255. All six ungated files have the same number of impulse rows.
- **"Inside frame" is not accuracy.** The 80% off-frame share of Impulse1 reproduces the first session's finding. Impulse3 is mostly inside the frame because its low-pass envelope TDOAs stay small (points cluster near the centre), not because they are right. Dolphins can also be outside the camera's field of view. Use the viewer, and later the clicker data, to judge accuracy.
- **The gate is strict.** About 60% of clicks have a closure above 1 sample, in line with the first session (64% of the strongest quarter pass, 23% of the weakest). The flagged clicks are not more often off-frame than the others (about 58% vs 52% inside), so the flag is not an off-frame detector. `closure_max` is a keyword; tune it on the clicker data.
- **Quick test on `Calibration_05_11_2025`** (see `results/multitdoa/`): the same checks pass. Note that `split_vid_au` reads sub-folders, so this run also took `not_so_good/251105_001.WAV` and `GH011400.MP4` (17,922 clicks from two recordings). For the ground-truth study (next item), use a folder with only `251105_002.WAV` and `GH011401.MP4`.

**Not verified here:** how the viewer draws `"triangle"`. If it ignores or renders it as a circle, change `flag_shape` (a keyword of `run_contiguous_folders_multitdoa`).

**Next** (from the recommendation list):
1. The D3 clicker calibration data already exists (`Calibration_05_11_2025`: 165 usable annotated pixels, 123 outliers tagged `o`). Use it for per-method angular error and the offset fit (recommendation 5); `src/cam_calib.jl` is the starting point.
2. Tune `closure_max` and choose the default overlay method from that error.
3. Multi-arrival handling and closure-based false-detection rejection (recommendations 6 and 7).
4. Housekeeping: other sessions, the `runtests.jl` runner (add `test/test_multitdoa.jl` to it as well), and the decision on the old functions (recommendations 9, 10 and 4).
