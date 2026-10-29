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

Integration and the invalid-TDOA/SRP controls now pass. A real-data preparation
is complete and Terra annotation is active; B preparation follows independently.
Final A/B fitted parameters remain pending annotation and physical fitting.
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

## Run and resume

The standalone entry point avoids the legacy script include chain:

```julia
include("src/clicker_calibration.jl")
run_clicker_calibration(
    "/media/spin/anas2/data/calf/Clicktest/clicker_20260209";
    output_dir=joinpath(pwd(), "calibration_runs/clicker_20260209/20261030_v1"),
    devices=(:A, :B), annotator=:herdr, notify=:telegram,
)
```

Run from an initialized Herdr pane for automatic Terra handoff. The agent reads
every frame without acoustic predictions; completion is authenticated by image
identity and CSV hashes. A persisted launch record is checked against current
agent liveness on resume. `annotator=:external` instead prepares the prompt and
returns until a complete external annotation CSV is available. The helper
`annotation_tools.py seal DEVICE_OUTPUT` validates and seals such a CSV.

Supported stages are `:prepare`, `:annotate`, `:fit`, `:review`, and `:all`.
`wait_annotations=false` returns immediately at the annotation boundary.
Existing preparation is reused only when inputs/settings/software and outputs
match. Annotation edits invalidate fits. The Telegram notification contains the
review gallery and CSV paths; fitting continues without an approval pause.

`tools/run_clicker_calibration.jl INPUT OUTPUT --device=A --stage=prepare`
provides the same interface as a CLI. `tools/run_clicker_job.py --job-dir NEW_DIR
--wake -- COMMAND...` adds persistent stdout/stderr, sampled process-tree RSS,
timeout/28 GiB default RSS limits, native-crash logging, and completion wake-up.
It holds the execution session through watcher delivery and cleans up the owned
process group on failures or termination. Use a fresh job directory for each
attempt; keep the calibration output directory stable for stage resume.

This machine's successful post-recovery Julia checks used:
`--startup-file=no --compiled-modules=no --pkgimages=no --threads=1
--gcthreads=1 -O0 --inline=no`. The job helper disables core dumps, sets a64 MiB
stack and single-thread BLAS/OpenMP, and keeps Julia depot/scratch in this repo.
The low-optimization flags are a runtime workaround, not a repair of the earlier
native compiler failures. Generated annotation graphics need Python with Pillow;
`CLICKER_PYTHON` selects it. No dependency installation occurs on import.

## Apply a fitted device

`set_device__ophk_acoustic_D3_clicker_A(calibration_path=...)` and the equivalent
B setter are defined in `config.jl`. They configure the corrected peak method,
sound speed, relative channel delays, camera rotation/projection and measured
image dimensions. Localization and pixel conversion apply each correction once.
Calling another existing device setter restores the uncalibrated behavior.
The no-argument defaults point to the final A/B calibration JSONs beside the
implementation; those files will be installed only after the real fits exist.

## Checks after recovery

- Physical fitting:85 native assertions, including one impossible960-sample
  outlier among120 clean rays. Recovered c1500.000000621 m/s and clean RMSE
  3.536e-7 px; invalid-event coverage remains120/121. These are synthetic values.
- Exact filtered snippet cache:16 native assertions.
- SRP comparator:53 native assertions, including original frozen-source parity,
  fractional-delay sign, coverage, held-out model identity and resume/tamper.
- Shared workflow and calibrated presets:45 native assertions after repairing
  lazy-module world-age handling. The failed first integration remains recorded.
- Annotation renderer/validation:6 Python tests; launcher failure/descendant/
  termination lifecycle:3 Python tests.

Run receipts/logs are under `calibration_runs/clicker_20260209/20261030_v1`.
These checks establish the tested implementation controls, not field accuracy.

## Device A preparation (2026-10-30)

Job `prepare_A_01` completed with exit0 in6949.83 s; sampled process-tree peak
RSS was17,289,895,936 bytes and completion notification was delivered. It produced
791 candidate events and791 unique original-resolution frames:775 mapped and16
ambiguous at frame-selection boundaries. Every listed image exists.

Seventeen waveform anchors cover both WAVs and all three video chapters.
The fitted video-minus-recorder offset is-10.03635704395856 s, clock drift
61.4488198059 ppm, maximum anchor residual121.013 microseconds and RMS55.800
microseconds. Minimum normalized waveform score is0.8335 and minimum peak
uniqueness ratio6.584. Localization accuracy remains unmeasured until annotation
and fitting complete.

The bounded preparation review reproduced all event/frame/sample/cache identities.
Its exact-camera-base diagnostic found that coarse2 kHz camera-audio lengths
round chapter bases by167/334 microseconds. Replacing those bases and refitting
the existing anchors changes six central frames, all already ambiguous. None
of the775 mapped frames changes or becomes ambiguous, even preserving the old
conservative uncertainty bound. Keep the original16 exclusions and frozen A v1
artifacts. The diagnostic is recorded in
`calibration_runs/clicker_20260209/20261030_v1/diagnostics/A_exact_camera_bases.json`;
an exact-base code correction remains a documented follow-up.

Recorded artifact timing attributes90m25s to sequential audio detection/TDOA/
snippet generation and21m08s to video probing/synchronization/mapping before
the first PNG. Image writing took about78s. These are combined-stage intervals,
not isolated kernel benchmarks.

The first automated annotation handoff failed because Herdr requires lowercase
agent names. Name generation now conforms to that rule; the real second handoff
completed and launched `clickann-a-1029222739-tf01` with gpt-6.1-terra/high for
all791 frames. The annotation agent has no acoustic prediction inputs. The
gallery-ready Telegram message and fit will follow its complete validated CSV.
