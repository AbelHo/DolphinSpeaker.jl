# Multi-arrival round 5: hydrophone equalization

Follows `multiarrival-round4-realism-sources_2026-10-24.md`; round-4 report: `improvements/tdoa_localization_2026-10-22/reports/round4_realism_sources.md`. The user said "proceed" (2026-10-24). As before: commit at every stage, write a mini report in `reports/` with figures and reproduction steps, and keep code general (any number of dolphins; here also any number of hydrophones, with a reference).

## Problem
Round 4 found that a fixed short response per hydrophone (mount and housing scattering, 20–150 µs) reproduces the real within-window statistics. With it, a per-click direction is within 3° of the truth for only 35–54% of detections. The response is fixed relative to the array, so it should be **estimable from the data and removable**.

## Approach
1. **New `src/hydrophone_eq.jl`.**
   - `estimate_hydrophone_eq(data, peaks; fs, rx_vect, ref, window, N, iters, …)`, from strong single-arrival clicks:
     - each click's direction comes from SRP on the currently equalized snippet;
     - the channels are aligned by the predicted plane-wave delays;
     - the relative transfer function `G_m(f) = Σ A_m conj(A_ref) / Σ |A_ref|²` and its coherence are accumulated;
     - the **linear-phase (pure delay) part is removed**: a constant delay per hydrophone cannot be told from a direction bias without ground truth, and is left for calibration data;
     - bins with low coherence are left uncorrected; magnitude correction is optional, with regularization;
     - iterate (directions improve as the equalization improves).
   - `apply_hydrophone_eq(data, eq)`: a zero-phase FIR (inverse of `G_m`) per hydrophone, applied by convolution; the reference channel is unchanged.
   - Works for any number of hydrophones.
2. **`src/pool_sim.jl`:** `simulate_scenario` also returns the hydrophone responses (`responses`; additive, the fingerprint test still holds), so an oracle equalizer can be built.
3. **Simulation benchmark (`scripts/bench_eq.jl`)** on the calibrated realistic simulation (2 dolphins).
   - Per scenario, equalizers are estimated from its own strong isolated clicks (4 s of data).
   - Compared: no equalization, the estimated equalization, and the oracle (true responses).
   - Metrics:
     - phase error of the estimated response vs the truth;
     - clicks with ≥ 2 SIC candidates, closure ≤ 1, waveform similarity;
     - accuracy within 3° (all / direct-path detections) for longest track + smart inheritance;
     - separation F1 with `SEPARATION_PARAMS_D3`.
4. **Real data (New_device).**
   - Estimate per file and check that the estimates agree across files (the evidence that the response is fixed).
   - Diagnostics: coherence before/after, phase/magnitude of `G_m`.
   - Rerun `run_multiarrival2` with an optional `eq` keyword (default unchanged), SIC, window −110:50, on equalized blocks.
   - Compare the within-window statistics, the proxies, the in-frame share and the dolphin separation (incl. ICI regularity within a dolphin) with round 4.
   - New viewer files.
5. **Tests:** `test/test_hydrophone_eq.jl`, synthetic:
   - known per-hydrophone responses are recovered (dispersive part) from clicks in random directions;
   - the equalizer applied to the data restores inter-channel coherence;
   - the reference channel is untouched;
   - identity when there is no distortion.
6. Report `reports/round5_hydrophone_eq.md`, README, memory.

## Verification
- Tests pass (new and existing).
- Simulation: the estimated equalizer reduces the phase error and recovers a large part of the oracle's gain in accuracy and separation.
- Real data: consistent estimates across files; better within-window coherence; proxies compared with round 4.

## Risks
- A real mount response may depend on direction; a single equalizer then removes only the average part.
- The estimate needs clicks whose direction is roughly right; iteration should converge, and the simulation shows whether it does.
