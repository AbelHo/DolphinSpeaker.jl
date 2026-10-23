# Item 4, multiple arrivals: per-click direction candidates + tracking, with diagnostics and a feature table

## Context
This is item 4 of `reviews/plans/multitdoa-webui-csvs_2026-10-22.md`. The multitdoa run on `Two-way_2males_06_10_2025/New_device` gives one direction per click (135,200 clicks; results in `/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device`). About 60% of clicks fail the closure gate. Within a train, consecutive clicks often point in unrelated directions.

Today's read-only analysis of that run (table below) shows why:
- multiple arrivals occur on two time scales: inside the ±50-sample TDOA window, and as extra detections 0.4–1.6 ms later;
- two animals click at the same time.

The aim of this round is to stop forcing one direction per click. Each click keeps several direction candidates, the candidates are linked over time into tracks, and the result goes to the viewer. Today's analyses are saved as scripts with their results, and a per-click feature table is written.

Decisions from the Q&A:
- First round: candidates + tracking, plus today's analyses as scripts/results and a per-click feature table.
- The array is portable. Its depth and its distances to walls drift slowly and are not recorded; the hydrophone–camera offset is fixed. So there is no fixed pool or image-source model, and echo handling is data-driven.
- `Calibration_05_11_2025` was recorded at a different place with a different mounting. Don't rely on it; more calibration sets will come later. Evaluation uses a synthetic multipath benchmark and internal-consistency metrics.
- All work of this round goes in the new folder `improvements/tdoa_localization_2026-10-22/` (scripts, results, scratch, README). Existing functions stay untouched; new code goes in new files.

## Current state (read-only analysis, 2026-10-22)
| finding | number | meaning |
|---|---|---|
| detections 0.4–1.6 ms after the previous one ("followers") | 51% of 135,200 | the detector (minimum spacing 77 samples = 0.4 ms) also picks up later arrivals |
| direction of a follower vs. the previous detection | median 43–55° apart (random pairs: 45.5°), even when both pass the gate | separate plane waves, not noisy copies |
| followers on the same bearing (\|Δaz\| < 3°) at a different elevation | 30% (random: 7%); Δel p25 −35°, p75 +41° | surface/bottom image pattern |
| consecutive detections 1.8–2.2 ms apart (burst pulses, ~500 clicks/s) | equal amplitude (median ratio 1.00); both unflagged: median 0.9°, 74% < 2° | real trains, very stable direction |
| clean clicks: isolated (gap ≥ 1.6 ms) and closure ≤ 1 | 23% | consecutive clean clicks agree within median 0.9–2.7° up to 30 ms apart |
| TDOA windows with ≥ 2 envelope peaks on some channel (sample of 3,764 clicks, file 0001) | 93%; flagged 61%, against 21% for single-arrival windows | overlap inside the window causes most closure failures |
| a weaker envelope peak 0.4–0.7 ms before the detection | 84% (76% of isolated clicks) | a window on that earliest peak is worse (closure pass 22% vs 32%; consecutive spread 37.5° vs 14.9°), so it is mostly reverberation, not a usable direct path |
| trains (≥ 20 clean clicks) with ≥ 2 direction clusters holding ≥ 10% | 159 of 215 | the second cluster is mostly not time-locked to the first (other male or movement); a small echo-locked share at 1.6–2.4 ms |
| flag rate by amplitude quartile | 74% (weakest) → 37% (strongest) | closure works as a true-click indicator |

## Options, by problem
Status: **now** = this round; next; later; ✗ = tested, no gain.

| problem | option | status |
|---|---|---|
| overlapping arrivals in the window | multi-peak SRP: top-K local maxima of the SRP-PHAT map per click, not only the global maximum | **now** |
| | beam arrival time per candidate (delay-and-sum to it, envelope peak), to order arrivals | **now** (feature) |
| | consistent-triplet search over several GCC peaks (zero-sum closure) | next, only if the SRP candidates miss arrivals |
| | iterative CLEAN / sparse multipath deconvolution (subtract the strongest beamformed arrival, rerun SRP) | later |
| | window on the earliest envelope peak (today); narrower window −30:30 (2026-10-21) | ✗ |
| echo detections (followers) | follower label by timing (gap < 1.6 ms) | **now** |
| | polarity-aware waveform similarity between follower and leader (surface echoes invert polarity) | next |
| | direct vs. image tracks: lead/lag between co-occurring tracks, same bearing with a vertical offset, amplitude | next |
| | self-calibrating image model: range, source depth and the slowly drifting array depth from direct + surface pairs (needs the array tilt relative to gravity) | later |
| integration over time | multi-target tracker over the candidates (Kalman random walk in (u, v), gating, nearest-neighbour assignment, birth/confirm/delete) | **now** |
| | SRP accumulated over a sliding block of clicks (≤ 30 ms) for weak clicks, feeding the tracker | next |
| | spatio-temporal DBSCAN on (t, u, v) (Clustering.jl) as a baseline | next (optional) |
| two males | group tracks into sources (tracks co-occurring at a fixed lag = one animal's direct + image paths) | next |
| | ICI-rhythm continuity per track | next |
| | spectral clustering with the existing HDBSCAN code (`src/hdbscan.jl`, `src/cluster_plot.jl`) as an extra feature | next (reuse) |
| | which track sits on which animal, judged from the video in your viewer | needs you |
| false detections | per-click feature table (closure, waveform similarity, GCC peaks, envelope arrivals, candidates, isolation, amplitude, track) | **now** |
| | rule-based score, then a classifier trained on calibration clicks plus weak labels from long tracks | later (calibration sets) |
| pipeline | move the validated parts into `run_contiguous_folders_multitdoa` as Impulse7/8 | later |

## This round

### Folder `improvements/tdoa_localization_2026-10-22/`
- `README.md`: findings, options, what was done, results, how to run, next steps.
- `scripts/`: diagnostics, benchmark and run drivers.
- `results/`: reports (.txt), small CSVs and figures, tracked in git.
- `results/data/` and `scratch/`: large tables and temporary files. A `.gitignore` in the folder keeps them out of git.

### New `src/multiarrival.jl`
Nothing existing is edited; it uses the helpers in `src/localization.jl`.
1. `get_srp_candidates(data, peaks; window, ref_channel, fs, rx_vect, c, weighting=:phat, band, interp=8, grid_deg=1.0, refine=10, K=3, min_sep_deg=8.0, min_rel=0.5, taper=0.0)`. Per click it returns up to K directions (u, v), each with its SRP score, TDOAs and beam arrival time. It also returns the GCC-PHAT closure, the waveform similarity `ncc_min` (smallest time-domain normalized cross-correlation over the pairs, taken at each pair's lag), the smallest pair peak height, and the GCC peak count.
   - It uses the same GCC curves, grid and refinement as `get_tdoa_srp`, so candidate 1 equals `get_tdoa_srp`.
   - Grid scores are kept as a 2D (u, v) lattice. Local maxima are found over 8 neighbours, then selected greedily by score: at least `min_sep_deg` apart and scoring at least `min_rel` × best. Each one is refined the way `get_tdoa_srp` refines its maximum.
   - Beam arrival time: shift the three channels by the candidate's TDOAs (FFT fractional delay), sum them, take the argmax of the Hilbert envelope. The result is in samples relative to the detection.
2. `envelope_arrivals(data, peaks; window, rel=0.3, min_sep=8)` returns the envelope peak count per channel (today's measure).
3. `follower_flags(samples, fs; echo_window_ms=1.6)` returns gap_prev and gap_next (ms) and the follower flag.
4. `track_candidates(t, uv, score; gate_deg=6, sigma_meas_deg=1.5, sigma_q_deg=20 (°/√s), n_confirm=3, t_confirm=0.05, t_delete=0.25, t_delete_tentative=0.03, birth_min_score=0.2)` returns:
   - a track id per candidate (0 = none) and the filtered direction of each association;
   - a track table: id, start/end, hits, mean az/el, spread, consecutive jitter, median score, confirmed.

   It runs online in time order:
   - predict: random walk, with variance growing with Δt;
   - gate by angle and by Mahalanobis distance;
   - assign greedily by distance: one candidate per track, but one click's candidates may go to different tracks (for example a direct path and its image);
   - Kalman update; start new tracks from unassigned candidates; confirm or delete tracks by hits and time.

   Defaults come from today's data (clean consecutive clicks agree within ~1–3°) and are tuned in the run.
5. `track_colours(tracks; palette)` gives each track the first colour of a fixed 8-colour palette that no concurrent track uses. Concurrent tracks never share a colour; colours are reused over time.
6. `run_multiarrival(res_folder; out_dir, viewer_dir=res_folder, K=3, kwargs...)` post-processes an existing multitdoa result folder. It does no re-detection and no video work.
   - Inputs: each `*__cps*.jld2` (clicks, amplitudes, trains, audio path, fs) and the Impulse6 webui CSV (the frame of each click, and the tonal rows).
   - Audio: for each train, reads blocks of at most 10 s (padded by 50 ms) with `wavread(...; subrange)` and filters them as `detect_impulse` does (`filter_simple`, `impulsive_band_pass`). Trains run on threads; memory stays at a few GB.
   - To `out_dir` (`improvements/tdoa_localization_2026-10-22/results/`):
     - `data/click_features_<folder>.csv`, one row per click. Columns: file, click, sample, t_s, frame, train, ppeak, gap_prev_ms, gap_next_ms, follower, closure_phat, ncc_min, gcc_peak_min, gcc_npeaks, env_peaks, Impulse5 flag, Impulse6 az/el, and for k = 1..3 az_k, el_k, score_k, tarr_k and track_k.
     - `data/tracks_<folder>.csv`.
     - `multiarrival_summary_<folder>.txt/.csv`.
   - To `viewer_dir`, next to the combined video and the Impulse1–6 files: drop-in webui CSVs that include the tonal rows.
     - `detection_pixels__color255_webui__Impulse7_srpTracks.csv`: one row per (click, confirmed track); px/py from the track-filtered direction; the track's colour; the click's alpha.
     - `detection_pixels__color255_webui__Impulse8_srpCandidates.csv`: every candidate; rank 1 as a circle, ranks 2–3 as triangles; the click's colour; alpha scaled by relative score.

### Scripts (`improvements/tdoa_localization_2026-10-22/scripts/`)
- `diag_detection_timing.jl <res_folder>`: interval and all-pair lag histograms, amplitude ratios, direction agreement by lag, follower geometry and flag rates → `results/diag_detection_timing.txt` plus a figure.
- `diag_train_clusters.jl <res_folder>`: direction map of clean clicks, clusters per train, timing and ICI of the second cluster → `results/diag_train_clusters.txt`.
- `diag_window_arrivals.jl <res_folder>`: envelope arrivals inside the TDOA window from WAV snippets, and the earliest-peak test → `results/diag_window_arrivals.txt`.
- `bench_multipath.jl`: semi-synthetic multipath benchmark.
  - Templates: clean real clicks (isolated, closure ≤ 0.3, high SRP).
  - Added arrivals: 1–2 per click, delayed 0.05–0.3 ms (inside the window), gain 0.5–1.25, polarity ±1. Directions are vertical images of the direct path or random azimuths. Real noise is added at the measured SNR.
  - Metrics: error of `get_tdoa_srp` and of guided `get_tdoa_gcc` against the direct path; recall@K (the direct path among the candidates, within 2°); its rank; how often the earliest beam arrival is the direct path.
  - Output: `results/bench_multipath.csv/.txt/.png`.
- `run_multiarrival_newdevice.jl` runs `run_multiarrival` on the New_device result folder.

Figures follow the dataviz guidance: a fixed categorical order and one axis per chart.

### Tests: `test/test_multiarrival.jl` (synthetic, no recordings)
- One plane wave: candidate 1 equals `get_tdoa_srp`, and no second candidate reaches `min_rel`.
- Two plane waves 30 samples apart in one window, from different directions: both directions appear among the candidates within 1°, and the beam arrival order is correct.
- `follower_flags` on hand-made times.
- Tracker, on two slowly moving sources plus images and clutter (fixed seed):
  - two confirmed source tracks plus their image tracks;
  - at least 95% of true associations correct;
  - clutter tracks are not confirmed;
  - a gap longer than `t_delete` ends a track.
- `track_colours`: concurrent tracks never share a colour.

### Reused, not changed
- `src/localization.jl`: `get_tdoa_srp`, `get_tdoa_gcc`, `gcc_curve`, `curve_at`, `subsample_peak`, `tau_max_samples`, `_rx_for`, `default_gcc_band`, `_clip_window`, `tdoa2uv`, `tdoa2dir_ls`, `angle2px`.
- `src/dsp.jl`: `filter_simple`.
- `src/run_multitdoa.jl`: `flag_clicks`, `MULTITDOA_COLUMNS`.
- `src/config.jl`: `set_device__ophk_acoustic_D3`.
- WAV.jl: `wavread(...; subrange)`.

## Verification
1. `julia --project=. -t auto test/test_multiarrival.jl` passes; `test/test_tdoa.jl` (24) and `test/test_multitdoa.jl` (38) still pass.
2. The three diagnostics scripts reproduce the numbers in the table above.
3. `run_multiarrival` on New_device:
   - The feature table has 135,200 rows, and its frames equal the Impulse6 CSV.
   - Candidate 1 matches Impulse6 for at least 99% of clicks within 0.5°; the rest are near-ties from block-wise filtering. This checks that the audio blocks and filtering reproduce the pipeline.
   - The summary reports the share of clicks with ≥ 2 candidates, coverage by confirmed tracks, track durations, per-track jitter, and the share inside the frame for Impulse7 vs. Impulse6/6g.
4. Benchmark: recall@3 against single-estimate accuracy, by echo delay and gain.
5. You open Impulse7/8 next to Impulse6 in the viewer: do the coloured tracks sit on the animals and on their reflections?

## Housekeeping
- Before implementing, copy this plan to `reviews/plans/multiarrival-candidates-tracking_2026-10-22.md` and the first draft to `..._v0.md`.
- After the run: the new folder's README, a one-line pointer in the "Next" list of the 2026-10-21 README, and updated project memory (followers, overlap, the portable array, calibration sets to come).
