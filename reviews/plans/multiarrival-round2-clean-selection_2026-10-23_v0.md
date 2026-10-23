# Multi-arrival round 2: sidelobe suppression (CLEAN) and one selected direction per click

Follows `multiarrival-candidates-tracking_2026-10-22.md`. Round-1 report: `improvements/tdoa_localization_2026-10-22/reports/round1_candidates_tracking.md`. The user said "proceed to the next" (2026-10-23) and asked for a commit at every stage plus a mini report in `improvements/tdoa_localization_2026-10-22/reports`.

## Problem
Round 1's SRP candidates hold the direct path more often than the single estimate (69% vs 56% in the benchmark), but:
- sidelobes of the main peak (0.4–0.6 × the peak for this array) give a second candidate for 98% of clicks;
- the tracker links those sidelobes, so tracks are short (median 0.19 s);
- nothing selects the direct path; the per-click earliest beam scores only 50%.

## Approach (new functions only; round-1 functions stay unchanged)
1. `get_srp_candidates_clean` in `src/multiarrival.jl` (CLEAN on the GCC curves).
   - With PHAT, a single arrival's GCC curve is the band-mask kernel h(τ) = mean over band bins of cos(2πkτ/N), shifted to its lag.
   - Candidate 1 = the SRP maximum, identical to `get_tdoa_srp`.
   - Each accepted candidate's kernel is subtracted from every pair curve, with amplitude equal to the curve value at its lag. The next candidate is the residual SRP maximum, at least `min_sep_deg` from earlier ones.
   - Per candidate it returns the residual score, the per-pair support (`pair_min`) and the beam arrival time. Thresholds come from a benchmark sweep.
   - The grid lookup is precomputed, and a `threaded` keyword allows a parallel driver.
2. Benchmark `bench_multipath2.jl`.
   - Same generator as round 1, plus echo-free trials, which measure false candidates.
   - Per-candidate rows (error to the direct path and to each echo, score, residual score, pair support, arrival time) for plain and CLEAN candidates.
   - Sweep: recall@3 vs false-candidate rate per threshold.
   - Per-click selection rules: rank 1, earliest beam, most consistent.
3. Track-level selection.
   - `label_tracks` gives each track its median arrival delay relative to the earliest co-occurring confirmed-track candidate in the same click. Direct tracks are near 0; image tracks lag.
   - `select_candidates` picks per click the candidate on the confirmed track with the smallest delay; ties go to more hits, then score. A click whose only tracks lag is marked as an echo.
   - Benchmarked on synthetic click sequences: a slowly moving source plus a vertical image whose gain varies per click and is sometimes stronger than the direct path.
4. `run_multiarrival2`: the same outputs as `run_multiarrival`, parallel over blocks, using the CLEAN candidates, track labels and selection. New viewer files:
   - `Impulse7c_srpCleanTracks`;
   - `Impulse8c_srpCleanCandidates`;
   - `Impulse9_selected` (selected direction; clicks judged echo-only drawn as triangles).
5. Tests in `test/test_multiarrival.jl` (new testsets):
   - CLEAN finds no second candidate for a single wave;
   - it finds both waves of a pair;
   - its candidate 1 equals `get_tdoa_srp`;
   - `label_tracks` / `select_candidates` on hand-made tracks.
6. Report `reports/round2_clean_selection.md` with figures and reproduction commands. Commit after each stage.

## Verification
- Tests pass (old and new).
- Benchmark: false candidates on echo-free clicks drop well below round 1's 98%, while recall@3 stays at or above round 1's 69%. Sequence selection beats single SRP when echoes are as strong as the direct path.
- Real data: candidate 1 still equals Impulse6. Fewer, longer tracks than round 1. Compare consecutive-click agreement and the inside-frame share of Impulse9 with Impulse6.
