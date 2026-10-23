# Multi-arrival round 3: reverberant pool simulation and an end-to-end benchmark

Follows `multiarrival-round2-clean-selection_2026-10-23.md`; round-2 report: `improvements/tdoa_localization_2026-10-22/reports/round2_sic_selection.md`. The user said "proceed" (2026-10-23). As before: commit at every stage, and write a mini report in `improvements/tdoa_localization_2026-10-22/reports/` with figures and reproduction steps.

## Problem
Round 2's window benchmark had no reverberation, and its conclusion failed on real data:
- per-click "earliest arrival" picked the direct path in 72–83% of synthetic trials;
- but it made consecutive real clicks less consistent.

The track lag rule did best on real data, but only by proxies (consecutive consistency, share inside the frame), and those proxies are unvalidated. To tune anything further we need ground truth with realistic multipath: echoes of the same click, reverberation of earlier clicks, the other animal, and the pipeline's own detector producing follower detections.

## Approach (new code only)
1. **`src/pool_sim.jl`: image-source simulator** for a rectangular pool and the D3 array.
   - Geometry: surface (reflection −0.95), bottom (+0.6) and four walls (+0.5), images up to order 3. The array is level, faces along the pool length, and has a planar front/back ambiguity like the real one.
   - Sources: a directional click beam (Gaussian main lobe, −3 dB at about 12°, −30 dB floor; images use the mirrored heading), spherical spreading, real clean click templates, and real quiet noise (file 0001, 53–98 s).
   - Click schedules: burst pulses (ICI 1.8–2.5 ms) and echolocation trains (ICI 8–40 ms) from 1–2 moving sources.
   - Output: 5-channel data (channels 1–2 zero), as the pipeline expects, plus a truth table per arrival (source, click, image order, arrival sample, amplitude, az/el).
   - Pool size, array depth and array position are randomized per scenario, because the real array is portable.
2. **`scripts/bench_sequences.jl`: end-to-end benchmark on simulated scenarios.**
   - Chain: the pipeline detector (`detect_impulse` with the D3 settings and a threshold of 9 × the noise envelope's 75% quantile, as the auto-threshold) → `detect_impulsetrain` → SIC candidates (windows −50:50 and −110:50) → tracker → selection rules.
   - Each detection is matched to its arrival in the truth table (direct or image, and of which click).
   - Rules: rank 1, earliest, track lag, longest track, followers inheriting the leader's selection, and combinations.
   - Metrics: the selected direction within 3° of the emitting click's direct path, for all / direct-path / echo detections; the precision and recall of the echo label; per-source track purity.
   - **Proxy validation:** compute the real-data proxies on the simulated output and check that they rank the rules as the truth does.
   - Source level calibrated so that detection amplitudes / threshold match the real distribution (p10/p50/p90 = 1.09 / 1.79 / 7.0).
3. **Tuning** of the tracker (gate, birth score, delete time, process noise) and of the selection on the simulated ground truth.
4. **Real data:** run the best configuration on New_device; new viewer file(s); compare the proxies with Impulse6 and with round 2.
5. **Tests:** `test/test_pool_sim.jl`.
   - Image geometry: the surface image's delay and elevation match the closed form.
   - The direct arrival's TDOAs match the localization model.
   - Arrivals are placed at the right sample.
   - The truth table is consistent.
6. Report `reports/round3_pool_simulation.md`, README section, memory. Commit after each stage.

## Verification
- Tests pass (new, plus `test_multiarrival.jl` 238, `test_tdoa.jl` 24, `test_multitdoa.jl` 38).
- The simulated detections look like the real ones: amplitude/threshold distribution, follower share (real 51%), share of clicks with ≥ 2 SIC candidates (real 78%), lag-histogram structure. Mismatches are reported.
- The benchmark reports, per rule, accuracy against the truth and the proxies; the chosen configuration is the best on the truth.
- Real data: the chosen configuration's proxies, compared with Impulse6 and the round-2 track rule.

## Caveats
The pool, the reflection coefficients and the beam pattern are assumptions; the real recording's geometry is unknown. Conclusions are relative comparisons of methods under realistic multipath, not absolute accuracies.
