# Click localization with the ophk_acoustic_D3 array under multipath: methods, findings and lessons

Session report, 2026-10-21 to 2026-10-24 · Repository `DolphinSpeaker.jl`, branch `dev2`
Detailed reports:
- `improvements/tdoa_localization_2026-10-21/README.md` (TDOA study, pipeline integration);
- `improvements/tdoa_localization_2026-10-22/reports/round1…round5_*.md`;
- how to run on recordings: `how_to_run_on_recordings.md`.

## Abstract
The DolphinSpeaker pipeline overlays the direction of each detected dolphin click on synchronized video. On the ophk_acoustic_D3 device (three hydrophones in a 0.18 m triangle, 192 kHz), most overlay points of a two-male recording fell outside the frame. This session traced the pipeline, replaced its TDOA estimation, diagnosed why directions remain unreliable, and built and evaluated methods for multiple arrivals, temporal tracking, separating several dolphins, and equalizing the hydrophones. Methods were tested against:
- semi-synthetic benchmarks;
- an image-source pool simulator calibrated to the real data;
- internal-consistency proxies on 135,200 real clicks.

The main findings:
1. **Sub-sample TDOA precision is not the limit.** Robust GCC/SRP estimators are about 10× more precise than the original integer methods, yet real directions stay inconsistent.
2. **Real recordings are dominated by multiple arrivals that vary from click to click.**
   - Half of all detections are later arrivals of another click, often vertical-image echoes.
   - 93% of analysis windows hold overlapping arrivals.
   - The strongest arrival explains only about a third of a window's energy.
3. **Echo detections should take the direction of the click they follow,** and directions should be chosen with track context (the longest track). Together these give the most consistent real directions: median change between consecutive clicks 17° vs 24°, 69% of points in frame vs 57%.
4. **Separating dolphins by direction works in simulation (F1 about 0.83) but does not clearly isolate the two males in the real recording.**
5. **The real hydrophones have only a tiny fixed response, so the real incoherence is not a property of the array**, and equalization cannot fix it.

Absolute accuracy needs a calibration recording with the same mounting.

## 1. System and problem
- **Device.** ophk_acoustic_D3 = Zoom F6, 192 kHz, 5 channels.
  - Channels 1–2 are bit-identical: a stereo mixdown, correlation ≈ 1/√3 with each hydrophone. They are correctly excluded (`rx_vect` NaN).
  - Channels 3–5 are the hydrophones on a 0.18 m equilateral triangle; the largest TDOA is about 22.4 samples. The reference is column 3 of the hydrophone matrix (channel 5).
  - The array plane is perpendicular to the camera axis. A planar 3-element array measures only the direction cosines (u, v) = (sin az · cos el, sin el), so **front and back are indistinguishable**; this needs a fourth, out-of-plane hydrophone.
- **Pipeline (`src/run_syncVid.jl`).** `run_contiguous_folders` → A/V sync → `detect_impulse` (5 kHz high-pass, envelope peaks, minimum spacing 77 samples = 0.4 ms, auto threshold 9 × the 75% envelope quantile) → `detect_impulsetrain` (≥ 5 clicks in 0.1 s) → TDOA → direction → pixel (`angle2px`, field of view 80° × 54.46°). Only the first TDOA method fed the overlay.
- **Data.**
  - `Two-way_2males_06_10_2025/New_device`: 3 WAVs, about 28 min, 135,200 in-train clicks, two male bottlenose dolphins in a pool. No ground truth.
  - `Calibration_05_11_2025`: a clicker recording, but from a different place and mounting.

## 2. Methods developed (all new functions; existing pipeline functions untouched)
| area | method | where |
|---|---|---|
| TDOA | GCC with PHAT / SCOT / ML weighting, band mask, physical lag limit, upsampled sub-sample peak, all-pairs least squares, envelope guide, quality outputs (peak heights, closure); SRP-PHAT; closed-form direction `tdoa2dir_ls`; fixed copies of buggy originals | `src/localization.jl` |
| pipeline | one webui CSV per TDOA method (Impulse1–6, gated twins), flagging by closure > 1 sample or infeasibility; memory-adaptive parallel scheduling (portable probes) | `src/run_multitdoa.jl` |
| multiple arrivals | top-K SRP candidates; beam arrival time; GCC-kernel CLEAN (rejected); **time-domain successive interference cancellation (SIC)**; follower flags; per-click quality features | `src/multiarrival.jl` |
| time | Kalman random-walk tracker over candidates (gating, birth/confirm/delete, merging); track colours; selection rules (rank 1, earliest arrival, track lag, longest track; followers inherit) | `src/multiarrival.jl`, `scripts/apply_rules_real.jl` |
| dolphins | `separate_sources` for any number n (given or estimated): gated Kalman per source; echo / repeat / other-dolphin follower roles; pairwise and ARI scores; k-means baseline; tuned parameters `SEPARATION_PARAMS_D3` | `src/source_separation.jl` |
| simulation | image-source pool (surface, bottom, walls; order 3); directional click beam; moving sources with burst-pulse and echolocation schedules; real click templates and noise; options for hydrophone scattering, diffuse reverb and off-axis low-pass; truth tables | `src/pool_sim.jl` |
| array | hydrophone equalization: aligned self-calibration and a direction-free blind estimator; zero-phase FIR; `with_delays` for calibration delays | `src/hydrophone_eq.jl` |
| tests | 783 unit tests (TDOA 24, multitdoa 38, multiarrival 238, pool_sim 427, source separation 32, equalization 24) | `test/` |

## 3. Findings

### 3.1 TDOA estimation (2026-10-21)
- **The original overlay method was broken.**
  - `get_tdoa_minmax` has a polarity vote that is always true, so it is plain argmax.
  - It gave infeasible TDOAs for about 48% of clicks, which were drawn at ±90°; 80% of overlay points fell off-frame.
  - `envelope_lpfiltered` assumes 500 kHz (at 192 kHz its "10 kHz" low-pass cuts at about 3.8 kHz).
  - Integer xcorr cycle-skips on band-limited clicks.
- **The new estimators are about 10× more precise where the answer is known** (synthetic directions, real clicks + real noise): median error 0.12–0.16° at 10 dB, against 1.2–1.4°. On the strongest real clicks, consecutive clicks agree within about 0.6° vs 3.7°.
- **On real clicks, precision is not the limit.** Waveforms differ between hydrophones (pair correlation about 0.5, later measured 0.21–0.26), and GCC closure exceeds 1 sample for about 60% of clicks. The closure test passes 64% of the strongest quarter of clicks and 23% of the weakest, so it is an indicator of a true click.
- **The solver is irrelevant.** NelderMead and the closed form give the same directions; NelderMead sometimes returns the back-hemisphere mirror.
- **In the pipeline (Impulse1–6 on New_device):** inside-frame share 21.9% (original) → 55.6% (GCC) / 57.4% (SRP). Impulse3's 74% comes from TDOAs biased toward the centre, not accuracy.

### 3.2 What the real recording looks like (round 1 diagnostics)
| observation | number |
|---|---|
| detections 0.4–1.6 ms after another ("followers") | 51% |
| their direction vs. the previous detection | unrelated: median 43–55°, the same as random pairs |
| followers on the same bearing with an elevation offset (surface/bottom image) | 30% (random: 7%) |
| burst pulses (about 500 clicks/s, equal amplitude) | consecutive clicks 1.8–2.2 ms apart agree within 0.9° (both passing closure) |
| clean clicks (isolated, closure ≤ 1) | 23%; consecutive ones agree within 0.9–2.7° up to 30 ms apart |
| windows with ≥ 2 envelope peaks | 93% (closure failure 61% vs 21% for single-arrival windows) |
| a window placed on the earliest envelope peak | worse than on the detection: the earlier peaks are reverberation, not the direct path |
| trains with ≥ 2 direction clusters | 159 of 215 (the other male or movement; a small echo-locked share) |

### 3.3 Multiple candidates, cancellation and selection (rounds 1–2)
- **Top-3 SRP local maxima contain the direct path far more often than the single maximum** (69% vs 56% of synthetic echo trials). But SRP sidelobes of this array reach 0.4–0.6 of the peak, so 98% of clicks get a spurious second candidate.
- **CLEAN on the GCC curves fails:** PHAT normalizes each frequency, so arrivals do not add in the GCC domain. **Successive interference cancellation (SIC)** in the time domain works instead: beamform the accepted arrival, gate ±10 samples, subtract it per channel, recompute PHAT. It gives 0.01–0.02 false candidates per echo-free click, with the direct path among the candidates in 74% of echo trials (86% with a window reaching 110 samples back).
- **Per-click "earliest arrival" works in a quiet benchmark but fails on real data,** where reverberation of earlier clicks arrives first: burst-pulse consistency gets worse (21° vs 10°). **Selection with track context** (track lag, longest track) is the best real-data selector.

### 3.4 Simulation and the end-to-end benchmark (rounds 3–4)
- **The image-source pool simulator was calibrated to the real detection statistics** (source level A0 = 9000, then 6000 with realism; detection threshold of the real file):
  - amplitude/threshold p10/p50/p90 1.12/1.99/6.29 (real 1.09/1.79/7.0);
  - followers 52% (real 51%).
- **Followers inheriting the previous detection's direction** raises the share of detections drawn at the right direct path from about 68% to about 80% in the clean simulation, for every rule. Tracker settings barely matter (±0.2 points).
- **Real-data proxies** (median direction change of consecutive isolated clicks; inside-frame share), window −110:50:

| rule | burst pulses | 2–10 ms | 10–30 ms | inside frame (+ inherit) |
|---|---|---|---|---|
| Impulse6 (SRP, pipeline) | 10.4° | 24.3° | 33.6° | 57.4% |
| rank 1 | 11.2° | 23.6° | 31.3° | 61.3% (63.6%) |
| earliest arrival | 21.0° | 27.1° | 31.6° | 62.9% (65.6%) |
| track lag | 11.5° | 19.2° | 25.3° | 67.7% (69.2%) |
| **longest track** | **8.0°** | **17.3°** | **25.2°** | **68.0% (69.3%)** |

- **A simulation clean inside the windows cannot validate these proxies.** With realism added (round 4), the proxies roughly ranked the rules as the truth did.

### 3.5 Separating dolphins (round 4)
- `separate_sources` handles any number of dolphins. Tuned on the simulation: pairwise F1 0.84 with 2 dolphins, 0.83 with 3 (held out), vs 0.85 / 0.80 for k-means on directions.
- **Estimating n over-counts** (spurious short sources), so give n when it is known.
- **A rhythm (inter-click interval) cost term does not help:** a source absorbs both interleaved trains when they share a direction gate. Using rhythm needs multi-hypothesis train separation.
- **On the real two-male recording (n = 2):**
  - 88% of clicks are labelled, and direction continuity within a label is twice the unseparated stream's;
  - one label follows a coherent swim track;
  - but click intervals within a label are not more regular, and both labels share similar mean directions, so direction alone does not isolate the males.

### 3.6 Is the array the problem? (rounds 4–5)
- **Round 4 (hypothesis).** A strong fixed scattering response per hydrophone reproduced the real window statistics in the simulator:
  - second SIC candidate in 78.5% of clicks (real 78%);
  - closure ≤ 1 in 33% (real 35%);
  - diffuse reverberation did not reproduce them.
- **Round 5 (test).**
  - **In simulation, equalization would be decisive:** an oracle equalizer gives 54% → 91% within 3° and separation F1 0.81 → 0.97. The blind and aligned estimators restore coherence, but they shift every direction by a common unknown offset, because per-hydrophone pure delays are not identifiable without calibration.
  - **On real data, the fixed responses are tiny** (dispersive phase 0.07–0.09 rad, magnitude spread 1.3–1.4×, against 1.46 rad and 6.7× in the simulator) and consistent across the three files.
  - **Equalizing them changes nothing** (closure ≤ 1 42.8% → 42.9–43.3%).
  - **No in-window echo tied to the array exists:** the cepstrum and envelope autocorrelation look the same in every direction bin.
- **Conclusion.** The real incoherence varies from click to click: reverberation of earlier clicks 2 ms apart in burst pulses, source-dependent echoes, and the other animal. The round-4 calibration matched the statistics with the wrong mechanism.

## 4. Major learnings
1. **Measure the failure mode before improving precision.** The original problem looked like poor TDOA precision; it is multipath and overlap. Better estimators help only on clean clicks.
2. **For this planar array, treat (u, v) as the measurement.** TDOAs are linear in (u, v), the closed-form solver is exact, front/back is ambiguous, and constant per-hydrophone delays equal a constant (u, v) shift. **Without ground truth, absolute direction biases cannot be detected.**
3. **PHAT is nonlinear.** Subtract arrivals in the signal domain (SIC), not in the GCC domain.
4. **A detection is not a click.** About half the detections are echoes of an earlier click, or the other animal within 1.6 ms. Handle events (a leader with its followers), not detections. Inheriting the leader's direction is the single most effective correction found.
5. **Per-click decisions are noisy in reverberant pools; decisions over tracks are robust.** The longest-track rule was the best selector on real data.
6. **Benchmarks must contain the confound that matters.** Round 2's quiet-noise benchmark favoured "earliest arrival", which failed on real data. Round 3's clean simulation could not validate proxies. Round 4's simulation matched aggregate statistics with the wrong mechanism (round 5). Calibrate simulations against several statistics, including ones tied to the mechanism, such as the fixed hydrophone response and the cepstrum.
7. **Internal-consistency proxies are useful but not truth.** A method that consistently picks an echo looks consistent. Use several proxies (consecutive-click change at several gaps, inside-frame share, rhythm regularity) and validate them where truth exists.
8. **Engineering:**
   - Per-file processing peaks at about 11 GB of RAM per 745-s file. Memory-adaptive scheduling prevents the OOM kills.
   - Audio blocks per click train, filtered like the pipeline, reproduce the pipeline's SRP exactly (100% within 0.5°).
   - Spawned tasks must not capture loop variables.
   - Driver scripts must include files by absolute path.

## 5. Negative results (do not repeat without a new idea)
- An analysis window centred on the earliest envelope peak (reverberation precedes the direct path); a narrower −30:30 window.
- CLEAN on GCC-PHAT curves.
- Per-click earliest-arrival selection on real data.
- A click-rhythm cost term in direction-based separation.
- Fixed hydrophone equalization on New_device: the response is too mild to matter.
- Simulation realism via strong fixed hydrophone taps: it matches the window statistics but not the real mechanism.

## 6. Open problems and recommended next steps
1. **Ground truth.** This is the most valuable next input.
   - A calibration recording (known clicker positions) with the same mounting would give absolute accuracy, constant per-hydrophone delays (`with_delays`), field-of-view/orientation errors, and validation of the proxies.
   - Video labels of which male is where would validate the separation.
2. **A simulator with click-varying clutter**, calibrated to the real window statistics *and* the mild fixed response. Use it to redo the absolute accuracy and the separation tuning.
3. **Multi-hypothesis click-train separation** (ICI continuity plus direction), for animals in similar directions; spectral identity per click (the existing HDBSCAN code).
4. **A fourth, out-of-plane hydrophone** for front/back and better elevation (hardware).
5. **Integrate into the pipeline** once validated: SIC candidates, longest-track selection with inheritance and separation as Impulse variants inside `run_contiguous_folders_multitdoa` (now run as post-processing).

## 7. Reproducibility
- **Runbook for any recording:**
  - `reports/how_to_run_on_recordings.md`
  - `scripts/run_on_recording.sh <res_folder> [--audacity] [--n N]` writes the viewer CSVs `Impulse7`–`Impulse11` and optionally Audacity labels.
- **Per-round reproduction commands:** in each round report.
- **Plans:** `reviews/plans/` (`trace-the-function-…`, `multitdoa-webui-csvs_2026-10-22`, `multiarrival-candidates-tracking_2026-10-22`, `multiarrival-round2…round5_*`).
- **Commits:** the earlier commits (TDOA estimators, multitdoa pipeline, memory scheduling), then `4549348` … `7d08c44` (rounds 1–5) on `dev2`.
- **Large tables** (per-click features, tracks, benchmark detections) are in `results/data/` (git-ignored); viewer CSVs sit next to each recording's combined video.
