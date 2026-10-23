# Multiple arrivals per click: direction candidates and tracking (ophk_acoustic_D3), 2026-10-22/23

Item 4 of `reviews/plans/multitdoa-webui-csvs_2026-10-22.md`. This round's plan is `reviews/plans/multiarrival-candidates-tracking_2026-10-22.md`; the first draft is `_v0`. The previous round is in `../tdoa_localization_2026-10-21/README.md`.

Data: `Two-way_2males_06_10_2025/New_device`, 3 WAVs, 135,200 in-train clicks. The input is the multitdoa run in `/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device`.

## Summary
- **Half the detections are later arrivals.** 51% of detections come 0.4–1.6 ms after another one. Their direction is unrelated to it: a median of 43–55° apart, the same as random pairs. 30% of them share its bearing but differ in elevation, which is the surface/bottom echo pattern.
- **Overlap inside the TDOA window causes most closure failures.** 93% of windows have two or more envelope peaks. Single-arrival windows are flagged 21% of the time, the others 61%.
- **The clean subset is very stable.** Clean clicks are isolated and pass closure; they are 23% of the total. Consecutive clean clicks agree within a median 0.9–2.7° up to 30 ms apart.
- **Several candidates per click do help.** In the synthetic benchmark the direct path is among the 3 SRP candidates in 69% of trials, against 56% for the single SRP estimate. When the echo is as strong as the direct path, the figures are 66% vs 43%.
- **Picking the right candidate is the open problem.** SRP sidelobes also produce a second candidate for 98% of clicks, synthetic ones included. The earliest beam arrival does not pick the direct path (50%, worse than taking candidate 1).
- **The first tracker works but fragments.** 84% of clicks sit on a confirmed track. But the tracks are short (median 0.19 s, 8 hits), and the per-track jitter is 3.1°.

## Findings from the existing run (read-only diagnostics)
Reports are in `results/`; each is produced by a script in `scripts/`.

| finding | number | report |
|---|---|---|
| followers (0.4–1.6 ms after the previous detection) | 51.4% | diag_detection_timing |
| follower vs. predecessor direction | median 43–55° (random 45.5°) | diag_detection_timing |
| follower on the same bearing (\|Δaz\| < 3°) | 30% (random 7%); Δel p25 −35°, p75 +41° | diag_detection_timing |
| all-pair lag histogram | peaks at 2, 4, 6, 8, 10 ms = burst pulses of ~500 clicks/s with equal amplitudes; excess at 0.4–0.6 ms | diag_detection_timing |
| clicks 1.8–2.2 ms apart, both unflagged | median 0.9°, 74% < 2° | diag_detection_timing |
| clean clicks | 23.3% | diag_detection_timing |
| flag rate by amplitude quartile | 74.5% (weakest) → 37.1% (strongest) | diag_detection_timing |
| windows with ≥ 2 envelope peaks (3,764 clicks) | 92.9%; flagged 61% vs 21% | diag_window_arrivals |
| window on the earliest envelope peak instead of the detection | worse: closure pass 22% vs 32%; consecutive spread 37.5° vs 14.9°. The earlier peaks are reverberation, not a direct path | diag_window_arrivals |
| trains with ≥ 2 direction clusters (≥ 10%) | 159 of 215 | diag_train_clusters |
| second cluster time-locked to the first? | mostly not: 54% come > 30 ms after the last cluster-1 click (random: 50%). Excess at 1.6–2.5 ms: 8.0% vs 3.8% | diag_train_clusters |

Interpretation: two kinds of multiple arrivals (inside the window, and later detections), plus a second animal or movement. Because the array is portable and its depth and distance to walls drift, the echo handling has to be data-driven.

## What was added
Nothing existing was changed.
- `src/multiarrival.jl`:
  - `get_srp_candidates`: the top-K local maxima of the SRP-PHAT map; candidate 1 equals `get_tdoa_srp`. It also returns the beam arrival time and, per click, the closure, `ncc_min`, the GCC peak height and the GCC peak count.
  - `envelope_arrivals`, `follower_flags`.
  - `track_candidates`: a Kalman random walk in (u, v), greedy gated assignment, birth/confirm/delete, and merging of converged tracks.
  - `track_colours`: concurrent tracks never share a colour (Okabe-Ito palette plus white).
  - `run_multiarrival`: post-processing of a multitdoa result folder. It reads padded per-train audio blocks of at most 10 s, filters them as `detect_impulse` does, and writes the feature table, the track table, the summary and the viewer CSVs.
- `test/test_multiarrival.jl`: 219 tests, synthetic.
- `scripts/`:
  - `common.jl`: the shared loader;
  - `diag_detection_timing.jl`, `diag_train_clusters.jl`, `diag_window_arrivals.jl`: the planning-phase analyses saved as scripts;
  - `run_multiarrival_newdevice.jl`;
  - `bench_multipath.jl`.

## Outputs
- `results/*.md`: diagnostics reports, `multiarrival_summary_New_device.md` (also as .csv), and `bench_multipath.md`.
- `results/bench_multipath.csv` and `.png`.
- `results/lag_histogram.csv`.
- `results/data/` (git-ignored, large):
  - `click_features_New_device.csv`, 135,200 rows. Columns: file, click, sample, t_s, frame, train, ppeak, gap_prev_ms, gap_next_ms, follower, closure_phat, ncc_min, gcc_peak_min, gcc_npeaks, env_peaks, impulse5_flagged, nfound; then az/el/score/tarr/track for each of the 3 candidates; then track_confirmed.
  - `tracks_New_device.csv`: one row per track. Columns: id, file, start/end, hits, confirmed, colour, mean az/el, jitter, median score.
- Viewer CSVs, written next to the combined video and the Impulse1–6 files on the data drive:
  - `detection_pixels__color255_webui__Impulse7_srpTracks.csv`: one row per (click, confirmed track), at the track-filtered direction, coloured by track. 241,688 rows, 63.8% inside the frame (Impulse6: 57.4%).
  - `detection_pixels__color255_webui__Impulse8_srpCandidates.csv`: every candidate. Rank 1 is a circle, ranks 2–3 are triangles; alpha is scaled by the score relative to rank 1. 398,859 rows.

## Results
**Real data** (`results/multiarrival_summary_New_device.md`):
- Candidate 1 matches Impulse6 for 100% of clicks (within 0.5°).
- Clicks with ≥ 2 candidates: 98.2%; ≥ 3: 96.8%. The median score ratio of candidate 2 to candidate 1 is 0.89.
- Tracks: 107,098 in total, of which 11,443 confirmed. Confirmed tracks last a median 0.19 s with 8 hits, at a jitter of 3.1°.
- Clicks on a confirmed track: 83.5% (isolated 85.4%, followers 81.7%).
- Run time: about 14 min with 12 threads; the threads are underused because of the many small per-train blocks.

**Benchmark** (`results/bench_multipath.md`; 3,000 trials, real clean clicks plus 1–2 echoes inside the window, SNR 19 dB). Share of trials with the direct path within 2°:

| group | SRP single | GCC guided | any of 3 candidates | earliest beam |
|---|---|---|---|---|
| all | 55.8 | 50.8 | 69.2 | 50.3 |
| echo gain 0.5 | 98.3 | 98.6 | 98.8 | 71.1 |
| echo gain 0.8 | 80.3 | 72.5 | 88.4 | 61.0 |
| echo gain 1.0 | 42.6 | 32.8 | 66.3 | 49.4 |
| echo gain 1.25 | 4.6 | 2.3 | 25.2 | 21.3 |

- Guided GCC (Impulse5) is slightly worse than SRP (Impulse6) under multipath.
- The candidates carry the direct path much more often than the single estimate, but nothing selects it yet.
- The earliest beam arrival is not a usable selector in its current form. A beam aimed at a sidelobe also peaks at the strongest arrival's time.

## Caveats
- Sidelobes: for the D3 geometry, SRP-PHAT sidelobes of a single broadband click reach 0.4–0.6 × the peak, so `min_rel=0.5` admits many of them. Most second and third candidates are sidelobes, not arrivals. Impulse8 therefore shows many triangles that are not real.
- The tracker links sidelobes too. The track count is high and the tracks are fragmented; the colours in Impulse7 change often.
- The benchmark puts echoes inside the window only, uses one click template channel for all hydrophones, and uses white-ish real noise. It measures selection under overlap, not detection-level echoes.
- The viewer rendering of Impulse7/8 is not checked here.

## Run it
```bash
julia --project=. -t auto test/test_multiarrival.jl
R=/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device
julia --project=. improvements/tdoa_localization_2026-10-22/scripts/diag_detection_timing.jl $R
julia --project=. improvements/tdoa_localization_2026-10-22/scripts/diag_train_clusters.jl $R
julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/diag_window_arrivals.jl $R
julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/run_multiarrival_newdevice.jl     # ~14 min
julia --project=. -t 12 improvements/tdoa_localization_2026-10-22/scripts/bench_multipath.jl 3000          # needs the feature table
```

## Round 2 (2026-10-23): sidelobe-free candidates and one direction per click
The full report, with figures and reproduction steps, is `reports/round2_sic_selection.md`; round 1's is `reports/round1_candidates_tracking.md`. The plan is `reviews/plans/multiarrival-round2-clean-selection_2026-10-23.md`.

- **Successive interference cancellation (SIC)** replaces local maxima: `get_srp_candidates_sic` subtracts each accepted arrival (beam waveform gated ±10 samples) from the channels and recomputes PHAT before searching again.
  - Benchmark v2 (3,000 trials, a third echo-free): 0.01–0.02 false candidates per echo-free click (round 1: 1.6–2.0); direct path among the candidates in 73.9% (window −50:50) and 85.7% (−110:50) of echo trials (round 1: 67.9%).
  - CLEAN on the GCC curves fails (PHAT is not additive).
- **In the benchmark** the earliest SIC arrival selects the direct path in 71.6% / 83.4% of echo trials, against 55.7% for the single SRP estimate.
- **On real data** 78–79% of clicks keep an accepted second arrival (reverberation, echoes, the other male), and per-click earliest selection is noisier than Impulse6.
- The **track lag rule** is the best real-data selector. With the −110:50 window:
  - consecutive isolated clicks 2.2–10 ms apart change by a median 19.2° (Impulse6: 24.3°), and 10–30 ms apart by 25.3° (33.6°);
  - burst pulses are unchanged (11.5° vs 10.4°);
  - 67.7% of points are inside the frame (Impulse6: 57.4%).
- Viewer file to check: `detection_pixels__color255_webui__Impulse9sicwt_trackSelected.csv` (circle = direct, triangle = echo-only or untracked). The others are listed in the round-2 report.
- New code: `get_srp_candidates_sic`, `get_srp_candidates_clean`, `srp_support`, `track_lags`, `select_candidates`, `run_multiarrival2` (`src/multiarrival.jl`); 238 tests.
- New scripts: `bench_common.jl`, `bench_multipath2.jl`, `run_multiarrival2_newdevice.jl`, `compare_selection_rules.jl`, `write_track_selection_overlay.jl`, `report_figures.jl`, `report_figures_round2.jl`.

## Round 3 (2026-10-24): reverberant pool simulation and an end-to-end benchmark
Full report: `reports/round3_pool_simulation.md`. Plan: `reviews/plans/multiarrival-round3-pool-simulation_2026-10-23.md`.

- `src/pool_sim.jl`: image-source pool simulator (382 tests). `scripts/bench_sequences.jl` runs the pipeline detector, SIC candidates, tracker and rules on 55 random pools and scores them against the truth. Amplitudes and follower share are calibrated to the real data.
- **Followers inherit the previous detection's direction:** about 68% → **about 80%** of detections at the direct path, for every rule. Tracker settings barely matter.
- **The simulation is cleaner than the real data:** 15% vs 78% of clicks with a second SIC candidate, so it cannot validate the real-data proxies.
- **Real data:** the **longest-track rule** is the most consistent (burst pulses 8.0° vs 10.4° for Impulse6; 2–10 ms 17.3° vs 24.3°; 10–30 ms 25.2° vs 33.6°). With followers inheriting, 69.3% of points are in frame (Impulse6: 57.4%).
- Recommended viewer file: `detection_pixels__color255_webui__Impulse10sicw_default_longesttrack-inherit.csv`.

## Round 4 (2026-10-24): realistic simulation, smarter inheritance, separating n dolphins
Full report: `reports/round4_realism_sources.md`. Plan: `reviews/plans/multiarrival-round4-realism-sources_2026-10-24.md`.

- **Realism.** `src/pool_sim.jl` gained hydrophone scattering, diffuse reverb and off-axis low-pass; the defaults reproduce round 3 exactly.
  - Per-hydrophone scattering, not reverb, reproduces the real windows: 78.5% vs 78% of clicks with a second SIC candidate, closure ≤ 1 in 33% vs 35%, followers 53% vs 51%.
  - On that simulation, one direction per click is within 3° of the direct path for only 35% (2 dolphins) to 54% (3 dolphins) of detections, with the best method (longest track + inheritance).
- **`src/source_separation.jl`** handles any number of dolphins (`n` given or estimated) and does smart inheritance (echo / repeat / the other dolphin's own click).
  - Tuned (`SEPARATION_PARAMS_D3`), pairwise F1 is 0.84 with 2 dolphins and 0.83 with 3 (held out; k-means 0.80).
  - A click-rhythm cost term did not help and was removed.
- **Real two-male data (`n` = 2, tuned).** 88% of clicks get a dolphin label, and direction continuity within a label is twice the unseparated stream's. But click intervals within a label are not more regular, so direction alone does not isolate the two males.
- Viewer file to check: `detection_pixels__color255_webui__Impulse11sicw_sources_n2_tuned.csv` (dolphin 1 blue, dolphin 2 orange, unassigned grey, echoes as triangles).

## Next
1. Viewer check of `Impulse11sicw_sources_n2_tuned`.
2. Equalize the per-hydrophone responses (the dominant error source per the simulation): estimate them from strong clicks, deconvolve before GCC/SRP, and test in the simulation, where the true responses are known.
3. Multi-hypothesis click-train separation (ICI continuity) for dolphins in similar directions; spectral identity from the existing HDBSCAN code.
4. Video labels of the males' positions as ground truth for the separation; validation on the calibration sets when they arrive.
