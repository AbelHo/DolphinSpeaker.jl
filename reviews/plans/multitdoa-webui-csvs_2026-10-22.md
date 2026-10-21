# Next steps after the TDOA study: per-TDOA webui overlays from the real pipeline (Impulse1–6)

## Context
The first session (`improvements/tdoa_localization_2026-10-21/README.md`) added better TDOA methods (`get_tdoa_gcc`, `get_tdoa_srp`, `tdoa2dir_ls`) and a comparison harness. Nothing from the README's recommendations has reached the real pipeline yet. This plan lists what is worth doing next, and details the item you chose first: run the pipeline the way `src/run_syncVid.jl` does, and get **one `detection_pixels__color255_webui*.csv` per TDOA method**, so the methods can be compared in your real-time viewer without rendering video.

Decisions from our Q&A:
- Start with the pipeline run, not with the calibration study.
- **Impulse5 = `get_tdoa_gcc` (PHAT, `guide=:envelope`) + `tdoa2dir_ls`; Impulse6 = `get_tdoa_srp` + `tdoa2dir_ls`.**
- Impulse5/6 each come as an ungated file and a gated twin. In the ungated file the clicks flagged as probably wrong are drawn with **shape = type = `"triangle"`**.
- The run covers the whole `New_device` folder (0001, 0002, 0003), as a set. `run_contiguous_folders` already treats the audio files of a folder as one contiguous recording.
- Existing functions stay untouched (memory rule). The new code goes in new functions.

## Recommended list (ordered)
| # | Item | Status |
|---|---|---|
| **1** | **Per-TDOA webui CSVs from the real pipeline, with Impulse5/6 (this plan)** | **now** |
| 2 | **Ground-truth error and offset fit on the D3 clicker data.** It already exists: `Calibration_05_11_2025/251105_002.WAV` (259 s) has 165 usable annotated pixels (123 more are tagged `o` = outlier), pipeline outputs and a sync delay of −5.07 s, in `data_res/dolphin/calf/calibration/D3/clicker/Calibration_05_11_2025/`. Build on `src/cam_calib.jl`: invert `angle2pxd` for per-method angular error, then fit roll/yaw, `fov_angle`, `c` and per-channel delay offsets. | next; the README thought this was "later" |
| 3 | Tune the quality gate (`closure_max`, now 1 sample) and pick the default overlay method from item 2's numbers. | after 2 |
| 4 | Multi-arrival handling: direct-path choice among GCC peaks, SRP accumulated over a click train, clustering the two males. Closure/waveform similarity as detection features to reject false detections. | research |
| 5 | Housekeeping: run the harness on other sessions (the README's conclusions rest on one file); replace `test/runtests.jl` (still PkgTemplates) with a runner that includes `test/test_*.jl`; decide whether to fix or retire the buggy originals. | any time |
| – | Front/back ambiguity needs a fourth, out-of-plane hydrophone. This is hardware, so it isn't planned here. | — |

## Item 1: design
Facts from reading the code:
- `process_detections` ([aspod.jl:165](src/aspod.jl#L165)) already computes Impulse1–4 and, with `flag_return=true`, returns them plus `res` (with `res.fs`, `res.res_impulse.data_filt` and `res.res_impulsetrain.pind_good`). Only Impulse1 and Impulse2 carry the radius/colour/alpha/shape columns. Impulse3 and Impulse4 are bare px/py.
- `run_contiguous_folders` ([run_example.jl:374](src/run_example.jl#L374)) writes only `res[1]` and `res[2]` (impulse, tonal) to `detection_pixels.csv`, and from it `detection_pixels__color255_webui.csv` (r, g, b × 255; `type` = copy of `shape`). Neither function has a hook for more methods, so a wrapper is needed.
- `get_tdoa_gcc(...; return_quality=true)` gives `pair_rms`, and for 3 hydrophones the closure |τ12+τ23+τ31| = 3·`pair_rms`. `tdoa2uv` gives u, v, and the TDOAs are infeasible when hypot(u, v) > 1.

### New file `src/run_multitdoa.jl` (nothing existing is edited)
1. `impulse_variant_rows(px, colour_cols, flagged; flag_shape="triangle")` is a pure helper. It returns the rows [px py radius r g b a shape] for one variant, with `shape` set to `flag_shape` where `flagged`. Unit-testable.
2. `process_detections_multitdoa(aufname, vidfname; res_dir, closure_max=1.0, kwargs...)`:
   - calls the existing `process_detections(...; flag_return=true)`;
   - takes the colour columns (radius, r, g, b, a) from Impulse1's rows. All variants have the same click set, so they can share them;
   - computes Impulse5 and Impulse6 from `res.res_impulse.data_filt`, `pind_good`, `window_impulsive`, `rx_vect[:, get_relevant_channels(rx_vect)]` and `res.fs`. Directions come from `tdoa2dir_ls`, pixels from `angle2px(ang, fov_angle)` as for Impulse1–4;
   - flags a click when closure > `closure_max` (from the GCC pass; Impulse6 reuses it) or its own TDOAs are infeasible;
   - returns the original fields plus `variants`: Impulse1 `minmax`, Impulse2 `xcorrMaxEnergy`, Impulse3 `envelopeFiltered`, Impulse4 `xcorrMaxPeak`, Impulse5 `gccPhatGuided`, Impulse6 `srpPhat`, and the gated twins `Impulse5g` and `Impulse6g` (flagged clicks dropped);
   - frees `data_filt` afterwards.
3. `run_contiguous_folders_multitdoa(folname; res_dir, closure_max=1.0, detection_types=1:2, kwargs...)` has the same flow as `run_contiguous_folders`: `run_analysis_split_vidau` (A/V sync and combined video), threaded processing per audio file, cumulative-duration plus sync-delay frame offset, and the audacity combination. It writes:
   - the standard `detection_pixels.csv` and `detection_pixels__color255_webui.csv`, unchanged;
   - **`detection_pixels__color255_webui__Impulse<N>_<name>.csv`** for each variant above. Each file holds that variant's impulse rows plus the tonal rows, so it is a drop-in replacement for the standard webui CSV;
   - `multitdoa_summary.csv` with, per variant, the number of clicks, the number flagged, and the share inside the video frame.
   - There is no video rendering; your viewer replaces it.

### Edit `src/run_syncVid.jl`
Include `run_multitdoa.jl` and call `run_contiguous_folders_multitdoa(folder; res_dir=..., closure_max=1.0)`. `folder` is `.../New_device`, so 0001–0003 run together. `res_dir` gets a new suffix (e.g. `..._multitdoa`), so the previous results aren't overwritten.

### Tests (`test/test_multitdoa.jl`, synthetic only)
- `impulse_variant_rows` marks exactly the flagged rows as `"triangle"`, and the gated twin drops them.
- The infeasible/closure flag logic on synthetic TDOAs: a consistent click is not flagged, a click with closure 3 samples is, and an infeasible one is.

## Verification
1. `julia --project=. -t auto test/test_tdoa.jl` and `test/test_multitdoa.jl` pass.
2. **Quick end-to-end run** on the single-WAV `Calibration_05_11_2025` folder (259 s) before the long run. This is the folder that holds the ground-truth annotation, so its outputs are reused in item 2.
3. **Full run** on `New_device` (3 WAVs, about 2 GB each, 94 GB RAM available, run with `-t 8` or more). Check:
   - all variants have the same number of impulse rows;
   - `detection_pixels__color255_webui.csv` equals the file from the earlier pipeline run (`.../tmp_20261007/Two-way_2males_06_10_2025_nooverlay2/New_device/`);
   - the Impulse1 file's impulse rows equal that standard file's impulse rows;
   - `multitdoa_summary.csv` agrees with the earlier findings: roughly 80% of Impulse1 points off-frame, and a much smaller share for Impulse5/6.
4. Open the Impulse1, 5 and 6 files in your viewer and look at the triangles. **Risk:** if the viewer does not draw `"triangle"`, the flagged clicks may vanish or fall back to circles. Tell me and I'll switch to another shape.

## Housekeeping when implementation starts (memory rule)
- Copy this plan into `reviews/plans/` **under a new name** (e.g. `multitdoa-webui-csvs_2026-10-22.md`). The name `trace-the-function-run-contiguous-folder-eager-puppy.md` is already taken there by the previous plan, which must not be overwritten. Save the first draft as `<that name>_v0.md`.
- Afterwards add a short "follow-up" section to `improvements/tdoa_localization_2026-10-21/README.md`: how to run, the output files, and the results of the three-file run.
