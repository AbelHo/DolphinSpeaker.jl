# Multi-arrival round 4: realistic simulation, smarter inheritance, separating n dolphins

Follows `multiarrival-round3-pool-simulation_2026-10-23.md`; round-3 report: `improvements/tdoa_localization_2026-10-22/reports/round3_pool_simulation.md`. The user agreed to round 3's three next steps (2026-10-24), with one added requirement: **the code must handle any number n of dolphins; testing with 2 is enough (3 optional).** As before: commit at every stage, and write a mini report in `reports/` with figures and reproduction steps.

## Targets from the real data (SIC −50:50 feature table, New_device)
| statistic | real | round-3 simulation |
|---|---|---|
| clicks with ≥ 2 SIC candidates | 78% | 15% |
| GCC-PHAT closure ≤ 1 sample (median) | 35% (4.6 samples) | – |
| waveform similarity between hydrophones `ncc_min`, median | 0.21 | – |
| energy explained by the first arrival `efrac_1`, median | 0.30 | – |
| windows with ≥ 2 envelope peaks | 95% | – |
| followers | 51% | 42% |
| detection amplitude / threshold, median | 1.79 | 2.11 |

## Approach
1. **Realism options in `src/pool_sim.jl`** (keywords; the defaults reproduce round 3 exactly):
   - `hyd_scatter`: a short random impulse response per hydrophone (mount and housing scattering; 3 taps, 4–30 samples). It makes the waveforms differ between hydrophones.
   - `reverb_level`, `reverb_t60`, `reverb_n`: a diffuse tail per click, i.e. random plane-wave arrivals with density growing as t² and exponential decay, from random directions.
   - `offaxis_lowpass`: the click is low-passed by its beam angle (off-axis clicks are duller), so direct paths and echoes differ in waveform.

   Calibration (`scripts/calibrate_realism.jl`): a sweep over scatter, reverb level and source level against the targets above; pick the closest setting.
2. **Separating n dolphins: new `src/source_separation.jl`.**
   - **Events:** each non-follower detection, with its followers, is one acoustic event.
   - **Source tracker over events,** for any n: one state per source (direction, Kalman random walk; last time; recent inter-click interval). Each event goes to the source with the lowest gated cost. A new source starts when no source matches, and with `n` given at most n sources are active; `n = nothing` estimates the number. Labels are reused like track colours.
   - **Smarter inheritance:** a follower takes its leader's source and direction, unless its own direction fits another active source better than its leader's, in which case it keeps its own direction and that source.
   - **Evaluation helpers:** pairwise precision/recall of "same source" for detections close in time (robust to label permutations and identity resets), adjusted Rand index (Clustering.jl), and the error in the number of sources.
3. **Benchmark (`scripts/bench_sources.jl`)** on the calibrated simulation:
   - 2-dolphin scenarios (main) and 3-dolphin scenarios (stress test);
   - selection rules with inherit-all vs smart inheritance;
   - source separation with n given vs estimated, against a k-means-on-direction baseline;
   - the round-3 proxies recomputed, to see whether they track the truth now.
4. **Real data (`scripts/apply_sources_real.jl`):** SIC −110:50 candidates from round 2, longest-track selection, smart inheritance, sources with n = 2 and n estimated.
   - Viewer CSVs coloured by dolphin; echo followers as triangles.
   - Per-source ICI and direction statistics as plausibility checks.
5. **Tests:**
   - `test/test_source_separation.jl`, with synthetic event streams of 2 and 3 sources, followers that are echoes or the other dolphin, and clutter;
   - `test/test_pool_sim.jl`, extended for the realism options.
6. Report `reports/round4_realism_sources.md`, README section, memory. Commit after each stage.

## Verification
- Tests pass (new and existing).
- Calibrated simulation statistics within a few points of the real targets; mismatches reported.
- Benchmark tables per rule / inheritance / separation method, with 2 and 3 dolphins.
- Real data: viewer files, per-source statistics, proxies.
