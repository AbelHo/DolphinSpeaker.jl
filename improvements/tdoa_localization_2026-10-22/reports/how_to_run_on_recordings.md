# How to run the multi-arrival methods on a real recording

This produces, for any ophk_acoustic_D3 recording folder, the same files as for New_device:
- viewer CSVs `detection_pixels__color255_webui__Impulse*.csv` next to the combined video;
- the per-click tables and summaries;
- optionally, Audacity label files.

All commands run from the repository root (`julia --project=.`, Julia 1.12). Each viewer CSV is a drop-in replacement for the pipeline's `detection_pixels__color255_webui.csv`: that method's impulse rows plus the tonal rows.

## Step 0 — the pipeline's result folder (Impulse1–6)
Every method starts from the folder written by `run_contiguous_folders_multitdoa` (`src/run_multitdoa.jl`). The steps:
1. Set `folder` (audio + video) and `res_dir` (a new folder each time, since ffmpeg does not overwrite) in `src/run_syncVid.jl`.
2. Run:
```bash
julia --project=. -t 12 src/run_syncVid.jl
```
3. The result folder is `<res_dir>/<folder name>/`. It holds:
   - `detection_pixels__color255_webui__Impulse1..6*.csv`, plus `Impulse5g` / `Impulse6g` (gated);
   - the standard `detection_pixels.csv` and `detection_pixels__color255_webui.csv`;
   - `combined_impulse__audacity.txt`, `combined_impulsetrain__audacity.txt`, `combined_tonal__audacity.txt`;
   - `multitdoa_summary.csv`;
   - the per-file `*__cps*.jld2` detections;
   - the combined video.

Details (memory-adaptive scheduling, options): the follow-up section of `improvements/tdoa_localization_2026-10-21/README.md`.

## Step 1 — all multi-arrival methods (rounds 1–4), optionally with Audacity labels
```bash
bash improvements/tdoa_localization_2026-10-22/scripts/run_on_recording.sh <res_folder> [--audacity] [--n 2] [--threads 12]
```
- `<res_folder>` is the step-0 result folder, e.g. `.../Two-way_2males_06_10_2025_multitdoa/New_device`.
- `--n` is the number of dolphins for the separation (default 2).
- `--audacity` also writes the label files (see Step 2).

On New_device (3 WAVs, 135,200 clicks, 24 threads) the steps take about 14, 8 + 14, 1, 1 and 1 minutes. Tested on New_device and on `Calibration_05_11_2025` (2 WAVs).

| round | step (script) | viewer CSVs written to `<res_folder>` (`detection_pixels__color255_webui__…csv`) | tables in `improvements/tdoa_localization_2026-10-22/results/` |
|---|---|---|---|
| 1 | `run_multiarrival_newdevice.jl` (SRP candidates + tracker) | `Impulse7_srpTracks` (tracks, coloured), `Impulse8_srpCandidates` (all candidates; rank > 1 triangles) | `data/click_features_<folder>.csv`, `data/tracks_<folder>.csv`, `multiarrival_summary_<folder>.md/.csv` |
| 2 | `run_multiarrival2_newdevice.jl sic` and `sicw` (SIC candidates; windows −50:50 / −110:50) | `Impulse7sic_tracks`, `Impulse8sic_candidates`, `Impulse9sic_selected` (per-click earliest arrival), and the same with `sicw` | `data/click_features_sic(w)_<folder>.csv`, `data/tracks_sic(w)_<folder>.csv`, `multiarrival_summary_sic(w)_<folder>.md/.csv` |
| 2 | `write_track_selection_overlay.jl sic` and `sicw` (track-lag rule) | `Impulse9sict_trackSelected`, `Impulse9sicwt_trackSelected` (triangle = echo-only / untracked) | – |
| 3 | `apply_rules_real.jl sicw "default" "longest track + followers inherit"` | `Impulse10sicw_default_longesttrack-inherit` (triangle = follower) | `real_rules_sicw_default_longesttrack-inherit[_<folder>].md/.csv` (all rules' proxies), `data/real_rules_…_<folder>.csv` (per-click directions) |
| 4 | `apply_sources_real.jl <n> tuned` (dolphin separation) | `Impulse11sicw_sources_n<n>_tuned` (dolphin 1 blue, 2 orange, …; unassigned grey; echoes triangles), `Impulse11sicw_sources_auto_tuned` (n estimated) | `real_sources_sicw_tuned[_<folder>].md` (per-dolphin statistics), `data/real_sources_sicw_tuned_<folder>.csv` (label, role, direction per click) |

Output names for recordings other than New_device carry a `_<folder>` suffix in `results/`, so nothing is overwritten.

**Which file to look at first:**
1. `Impulse10sicw_default_longesttrack-inherit` — the most consistent single direction per click.
2. `Impulse11sicw_sources_n2_tuned` — the dolphin separation.

The other files are for comparison.

## Step 2 — Audacity labels (also via `--audacity`)
```bash
julia --project=. improvements/tdoa_localization_2026-10-22/scripts/export_audacity.jl <res_folder> [n]
```
Labels go into `<res_folder>`, one point label per click, on the same combined timeline as the pipeline's `combined_impulse__audacity.txt`: files concatenated in name order, no A/V sync delay. The files:

| file | label text |
|---|---|
| `combined_Impulse9sicwt__audacity.txt` | `az el direct\|echo\|untracked` (track rule) |
| `combined_Impulse10sicw_longesttrack-inherit__audacity.txt` | `az el [follower]` |
| `combined_Impulse11sicw_sources_n<n>__audacity.txt` | `dolphin k\|unassigned role az el` (role: leader, own, repeat, echo) |

Import them in Audacity with *File → Import → Labels* next to the combined audio, or next to the pipeline's own label tracks.

## Optional diagnostics on the same recording
| script | output |
|---|---|
| `compare_selection_rules.jl sicw <folder>` | proxies of all selection rules (`selection_rules_sicw_<folder>.md`) |
| `diag_detection_timing.jl <res_folder>`, `diag_train_clusters.jl <res_folder>`, `diag_window_arrivals.jl <res_folder>` | round-1 diagnostics (follower share, lag histogram, clusters per train, arrivals inside the window); they write to the fixed `results/diag_*.md`, so New_device's copies are overwritten |
| hydrophone equalization (round 5) | not in the runner: it had no effect on New_device. To check another recording: `estimate_eq_real.jl` (New_device paths inside; change `wavdir` and the feature table) or, from the REPL, `estimate_hydrophone_eq_blind` + `apply_hydrophone_eq` (`src/hydrophone_eq.jl`) |

## Limits
- The scripts assume the ophk_acoustic_D3 geometry and camera: `set_device__ophk_acoustic_D3`, image 2740 × 1520, field of view 80° × 54.46° (hard-coded in the CSV writers).
- `SEPARATION_PARAMS_D3` (round 4) were tuned on the simulation; round 5 found its fixed-response mechanism unrealistic, so treat them as a starting point.
- Frames in the viewer CSVs are copied from the Impulse6 file of the same folder, so the A/V sync is the pipeline's.
