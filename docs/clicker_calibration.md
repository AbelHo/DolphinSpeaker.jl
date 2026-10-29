# A/B clicker calibration

Implementation is being restored and integrated directly in this Julia repository
following the 2026-10-30 machine restart. The user authorized progressive scoped
commits. Persistent source is under `src/clicker_calibration`; dedicated tests
are under `test/clicker_calibration`. Generated evidence, runtime depots, scratch
files, and outputs belong under ignored `calibration_runs/`, not `/tmp`.

The input drive currently mounts at `/media/spin/anas2`. The two independent
inputs are `data/calf/Clicktest/clicker_20260209/A` and `B` on that drive.
Source recordings remain read-only. Use physical channels 3–5, reference physical
channel 5, and the confirmed 0.18 m D3 equilateral array.

The planned native Julia entry point prepares exact timestamped click frames,
hands them to a gpt-6.1-terra/high Herdr annotator, produces review overlays and a
Telegram notification, then fits physical parameters independently for A and B.
Fitting continues during user review. Sound speed is fitted within 1400–1600 m/s;
the camera model includes alignment, focal lengths and supported optical terms,
plus relative channel delays. Only visible underwater annotations with valid
synchronization enter fitting. An empirical pixel correction is not requested.

## Recovery checkpoint

The lost temporary sources are recovered from persisted code-writing records,
without executing unrelated recorded commands. Recovery manifests in
`calibration_runs/clicker_20260209/recovery` identify source transcripts, exact
recovered SHA-256 values, and changed include paths. Preparation and physical-fit
sources match their last reported precrash hashes. Historical execution evidence
does not substitute for the remaining post-recovery integration checks.

Current work: finish integration, validate the invalid-TDOA influence correction,
complete the SRP comparator, then run A and B. Full real-data preparation, Terra
event annotation, and final A/B fitted presets have not yet been completed.
Native Julia had intermittent CSV/LLVM SIG11 failures before the machine crash;
record any recurrence durably and use bounded checks before heavy computation.

## Numerical and annotation contract

Public tables use zero-origin samples/frames and original-resolution pixels,
x right/y down; audio delays are Float64 samples relative to selected channel 3.
Array directions use Julia azimuth/elevation:
`[sin(az)*cos(el), sin(el), cos(az)*cos(el)]`, with +y up and +z forward.
Camera rotation is proper SO(3); pixel y reverses the vertical camera coordinate.

`frames.csv` lists unique source frame identities and image paths without
predictions. `annotations.csv` has
`frame_id,x,y,visibility,medium,uncertainty_px,note,annotator`.
Visibility is visible/uncertain/occluded/absent; medium is
underwater/air/mixed/unknown. The completion manifest binds both CSV hashes.

Every fifth 30-second block is held out, with a two-second boundary purge.
Optical model selection uses only training blocks. Evaluation is saved before
the final deployment fit. Reports distinguish held-out predictions, final-fit
predictions, the dimension-corrected original D3 baseline, and estimator coverage.
