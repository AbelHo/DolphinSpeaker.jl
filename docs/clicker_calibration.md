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
The landmark is consistently the center of the identifiable clicker housing.
A finger covering the button does not invalidate an otherwise clearly visible
housing; the uncertainty radius covers plausible emitter position within it.
A hand or cord alone cannot establish a hidden housing center. Mostly hidden or
indistinct clickers remain uncertain/occluded and excluded from physical fitting.

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

## A annotation and first-fit findings

Terra completed and sealed all791 frames:180 visible,207 uncertain,1 occluded,
403 absent. There are44 visible underwater clickers. After the unchanged timing
and block-boundary exclusions,39 events remain:23 training events in four blocks
and16 evaluation events in block34. The annotation SHA256 is
`ada438af58486074f8c46ed3a55c5a574a0b49e8ba69fb35272453cee3e0aa65`.
The review gallery is under A/annotation_review/ada438af5848; its location was
delivered through Telegram. Use its Visible underwater filter to inspect the
calibration landmarks. Annotation and image files remain unchanged during fitting.

The first physical fit is diagnostic and unsuitable as A's installed default.
Only5/16 evaluation events project validly, with1750.40 px RMSE and1991.96 px
median error. The other11 have impossible delays. Its sound speed reaches the
lower bound and vertical field of view the upper bound; a training-only local
sound-speed standard error of about480 m/s indicates weak physical identification.
The final all39-point refit also remains poor:18/39 valid projections,779.56 px
RMSE. These values are retained in the original fit/evaluation reports.

The subsequent estimator comparison uses the same frozen training-only physical
model for every method. It is not a method-specific recalibration:

| Training data | Valid projections | RMSE, px | Median, px |
|---|---:|---:|---:|
| Corrected extrema |13/23|1637.79|259.79|
| Guided GCC-PHAT |23/23|230.90|146.96|
| Fixed original SRP-PHAT |23/23|256.49|147.72|

Eight training extrema estimates are impossible even across the permitted
sound-speed/delay parameter box. All23 training GCC estimates are physically
feasible at the original1540 m/s. No sign, reference-channel, unit, or rotation
implementation fault was found in the bounded audit. This supports changing the
delay estimator before introducing extra optics; it does not validate a new
physical calibration yet.

A separate method-refinement plan is frozen at A/refinement_v1/plan.json,
SHA256`0931f5bfd68587af9efde414e6d1025f6f49ce1fc5c71212dd57646e137a45fa`.
It compares two predeclared physical fits (GCC and SRP), using leave-one-training-
block-out selection with equal block weights, the same eight-parameter centered
model and fixed optimization settings. The nine-point fold is reported explicitly.
Selection and its training fit must be frozen before replay of block34, which is
already inspected and therefore follow-up evidence rather than fresh confirmation.
The original first-fit artifacts remain untouched; refined results have separate
paths. B preparation continues independently.
