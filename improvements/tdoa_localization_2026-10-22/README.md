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

## Next
1. **Sidelobe suppression** before tracking: CLEAN, i.e. subtract the best candidate's predicted GCC pattern and rescore, or a sidelobe model. This should cut the false candidates and fragmentation, and make Impulse8 readable.
2. **Candidate selection** using the track context (prefer the candidate that continues a long track), then re-run the benchmark with selection.
3. Retune the tracker on the cleaned candidates (gate, `birth_min_score`, `merge_deg`), and parallelise over clicks instead of blocks.
4. Direct vs. image tracks; tracks → sources (the two males), using ICI and the existing HDBSCAN spectra; a polarity-aware follower ↔ leader similarity.
5. False-detection score from the feature table, once the calibration sets arrive.
