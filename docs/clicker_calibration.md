# A/B clicker calibration

The workflow is implemented directly in this Julia repository and was recovered
after the 2026-10-30 machine restart. Persistent source is under
`src/clicker_calibration`; dedicated tests
are under `test/clicker_calibration`. Generated evidence, runtime depots, scratch
files, and outputs belong under ignored `calibration_runs/`, not `/tmp`.

The input drive currently mounts at `/media/spin/anas2`. The two independent
inputs are `data/calf/Clicktest/clicker_20260209/A` and `B` on that drive.
Source recordings remain read-only. Use physical channels 3–5, reference physical
channel 5, and the confirmed 0.18 m D3 equilateral array.

The native Julia entry point prepares exact timestamped click frames,
hands them to a gpt-6.1-terra/high Herdr annotator, produces review overlays and a
Telegram notification, then fits physical parameters independently for A and B.
Fitting continues during user review. Sound speed is fitted within 1400–1600 m/s;
the camera model includes alignment, focal lengths and supported optical terms,
plus relative channel delays. The default policy admits only visible underwater
annotations with valid synchronization. Explicit provisional-coordinate policies
also admit localizable uncertain points with their original uncertainty weights;
the original Terra labels are never changed. An empirical pixel correction is
not applied.

## Recovery checkpoint

The lost temporary sources are recovered from persisted code-writing records,
without executing unrelated recorded commands. Recovery manifests in
`calibration_runs/clicker_20260209/recovery` identify source transcripts, exact
recovered SHA-256 values, and changed include paths. Preparation and physical-fit
sources match their last reported precrash hashes. Historical execution evidence
does not substitute for the remaining post-recovery integration checks.

Integration and the invalid-TDOA/SRP controls pass. Both real-data preparations
are complete. Both have sealed independent annotations, initial fits and
method-specific refinements. Their results remain diagnostically weak,
as detailed below; a successfully completed fit does not establish field accuracy.
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
Visibility is visible/uncertain/occluded/absent, or excluded for an explicit
human-selected scope; medium is
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

Supported stages are `:prepare`, `:sync`, `:annotate`, `:fit`, `:review`, `:refine`, and `:all`.
`wait_annotations=false` returns immediately at the annotation boundary.
Existing preparation is reused only when inputs/settings/software and outputs
match. Annotation edits invalidate fits. The Telegram notification contains the
review gallery and CSV paths; fitting continues without an approval pause.

`tools/run_clicker_calibration.jl INPUT OUTPUT --device=A --stage=prepare`
provides the same interface as a CLI. `tools/run_clicker_job.py --job-dir NEW_DIR
--wake -- COMMAND...` adds persistent stdout/stderr, sampled process-tree RSS,
timeout/memory monitoring, native-crash logging, and completion wake-up.
It holds the execution session through watcher delivery and cleans up the owned
process group on failures or termination. Use a fresh job directory for each
attempt; keep the calibration output directory stable for stage resume.

New runs use `--startup-file=no --compiled-modules=no --pkgimages=no
--threads=auto -O2`. Both the 139-assertion native policy suite and actual B fit
passed with these settings. The job helper disables core dumps, sets a 64 MiB
stack, uses all available CPU threads for Julia/BLAS/OpenMP, and retains runtime
files in this repo. Heavy jobs share a queue lock. Videos use qualified NVDEC
and NVENC on the RTX 4090, with all available CPU cores for supporting work.
Historical O0/single-thread crash-recovery flags are no longer the defaults.
Generated annotation graphics need Python with Pillow;
`CLICKER_PYTHON` selects it. No dependency installation occurs on import.
Saved-result plots additionally require NumPy and Matplotlib in
`CLICKER_PLOT_PYTHON` (default `python3`), independently of the Pillow interpreter.
Plots run automatically after each fit and refinement; `diagnostic_plots=false`
or CLI `--no-plots` disables only these plots. `render_diagnostics(DEVICE_OUTPUT;
model_dir=...)` regenerates them without fitting. Automatic outputs are under
`diagnostics/automatic`; explicitly selected waveform examples can use a separate
output directory. Plots show spatial/time errors, speed profiles, and authenticated
cached waveforms. Acoustic descriptors are descriptive and never exclude events.

## Apply a fitted device

`set_device__ophk_acoustic_D3_clicker_A(calibration_path=...)` and the equivalent
B setter are defined in `config.jl`. They configure the corrected peak method,
sound speed, relative channel delays, camera rotation/projection and measured
image dimensions. Localization and pixel conversion apply each correction once.
Calling another existing device setter restores the uncalibrated behavior.
Refined presets use their selected GCC or SRP estimator with the actual sample
rate and trained window, preserving the clicker1 kHz high-pass preprocessing.
The fitted delays and camera rotation are applied only in projection. Switching
devices also restores the previous default localization window.
The no-argument defaults load separate committed A/B GCC calibration JSONs from
this recording. These are provisional fits with substantial outliers, not a
validated accuracy guarantee or transferable parameters for another mounting.

## Current A/B results

Final delivery root: `/media/spin/anas2/data/calf/Clicktest/clicker_20260209/results`.
Each device index links the full Terra-only synchronized hydrophone review movie,
individual raw/Terra/comparison frames, the calibrated comparison movie when
complete, and preserved numerical artifacts. Green/amber markers are original
Terra points and their uncertainty circles; red markers are calibrated acoustic
predictions. Impossible or outside-image rays are retained as statuses and never
clamped to the image edge. Markers hold for 120 ms for review, without trajectory
interpolation. Movies retain every original video frame and exact chapter timing.

A uses the user's inclusive gallery selection 00459–00701: 233 uncertain point
estimates, 10 unlocalizable selected frames and 548 excluded frames. Its separate
fit-input adapter preserves original truth and explicitly admits those estimates.
B retains its original labels and uses the declared
`terra_estimated_underwater_coordinates` observation policy; these are provisional
visual estimates, not human-confirmed positions. The same timing and temporal
purge rules apply to both.

| Result | A | B |
|---|---:|---:|
| Candidate events retained | 791 | 469 |
| Fitting cohort (train / reused evaluation) | 202 (163 / 39) | 200 (160 / 40) |
| Selected estimator | GCC | GCC |
| Training-block GCC / SRP score (px; lower better) | 247.09 / 302.96 | 301.42 / 408.77 |
| Reused evaluation projected / eligible | 35 / 39 | 40 / 40 |
| Reused evaluation median / RMSE (px) | 176.74 / 2564.63 | 132.30 / 2007.23 |
| Deployment descriptive median / RMSE (px) | 44.31 / 1755.38 | 43.11 / 911.42 |
| Deployment effective sound speed (m/s) | 1562.18 | 1600.00 (bound) |
| Deployment horizontal / vertical FOV (degrees) | 75.89 / 56.17 | 86.82 / 63.78 |

Errors are original-resolution Euclidean pixels on projected rays; denominator
coverage is reported separately and includes invalid events. Reused evaluation
was already inspected, so it is descriptive evidence rather than an untouched
confirmation set. Both models report `fitted_weak_identifiability`; small median
errors do not cancel the severe outliers or establish a physical sound speed.

`project_calibration_events(original_dir, model_dir)` additionally writes
`inference_events.csv/json` for every original candidate, including purged and
unannotated events. It uses the frozen selected model and original delay data,
joins visual truth only after projection, and never changes fit metrics. Media
consumers authenticate these sidecars and retain the frozen cohort CSV separately.

The implemented improvement is method-matched GCC/SRP fitting selected only on
training blocks, with invalidity and coverage preserved. Further work should
collect independently confirmed clicks over more poses/heights, constrain camera
and array geometry or sound speed independently, and validate on a new recording.
Frequency-band GCC comparisons are a possible next experiment, requiring a new
training-only selection plan; residuals or spectral heuristics alone must not be
used to declare which events really came from the clicker.

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
all791 frames. The annotation agent had no acoustic prediction inputs. The
gallery-ready Telegram message was delivered and the fits below are complete.

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
paths. B annotation proceeds independently.

The refinement endpoint is `refine_calibration(DEVICE_OUTPUT)`, also available
as `stage=:refine`; `stage=:all` runs it after the first-fit estimator review.
For a new device it freezes the same fixed protocol and input identities before
fitting. Results stay under `refinement_v1`, leaving the first fit intact.
Refined galleries use the original annotation set with the separate refined
model/predictions directory. Invalid projections are stated in each caption.

Refinement synthetic controls passed67/67, including exact method-matched
inference, the nine-training-point fold, immutable baseline/resume checks, and
unchanged selection/training parameters when only evaluation labels change.
The updated combined workflow/preset controls passed58 native assertions; the
annotation renderer passed7 Python tests including separate refined-model
rendering and explicit invalid-projection captions. These remain synthetic or
integration evidence and do not establish real-data accuracy.

## A physical refinement results

Job `refine_A_01` completed in 79.63 s with exit 0 and about 1.2 GB peak sampled
RSS. All original first-fit and annotation inputs remain unchanged. Training-only
equal-block validation selected GCC (163.42 px) over SRP (280.19 px). GCC scores
for omitted blocks 33/35/40/41 were 82.30/445.95/62.07/63.37 px; omitting block35
leaves only nine training points. All four GCC best fold fits converged, but all
folds reached the sound-speed upper bound and camera pose varied by 22.14 degrees.

| GCC physical model | Projected / eligible | In image | RMSE, px | Median, px |
|---|---:|---:|---:|---:|
| Selected model on training |23/23|23|323.99|43.53|
| Selected model on reused block34 |16/16|11|2110.35|414.03|
| Final all39-point refit, descriptive |39/39|34|1677.76|68.90|

Five catastrophic evaluation events remain in all reports. The final fit improves
the median over the same-method uncalibrated D3 baseline (111.35 px), but worsens
RMSE (baseline952.49 px). This is not a reliable general localization calibration.
The selected deployment parameters are c1600 m/s at the bound, horizontal/vertical
FOV77.7691/58.5748 degrees, fx1676.477/fy1355.003 px, centered principal point,
no radial distortion, and relative delays13.5516/2.55757/0 microseconds.
Training-only nuisance-refitted speed objectives change by only about2.75% over
1400–1600 m/s; c is weakly identified and is not a sound-speed measurement.

The refined gallery is
`A/refined_review/ada438af5848_6b04b7370952/index.html`. Saved plots and explicitly
selected diagnostic waveforms are under `A/refinement_v1/diagnostics`.

## B preparation and conservative synchronization review

Job `prepare_B_01` completed with exit0 in6982.87 s and sampled peak RSS
17,630,232,576 bytes. It produced469 candidate events and469 source frames,
initially457 mapped plus12 ambiguous. Seventeen accepted waveform anchors cover
both WAVs and all three chapters; offset is-9.641241734 s, drift65.340131 ppm,
maximum residual146.3 microseconds. The annotation handoff uses the unchanged
469 frame identities and gpt-6.1-terra.

`review_sync_eligibility(DEVICE_OUTPUT)` compares saved waveform anchors after
replacing coarse2 kHz chapter bases with sums of declared camera audio durations.
It retains the larger old/new uncertainty, all original exclusions, and unchanged
source images. It adds an exclusion if the counterfactual changes a source frame
or crosses a frame boundary. It never uses annotation coordinates or fit errors.
The full workflow runs this review before annotation/fitting; `stage=:sync` runs
only this review without reprocessing recordings.

The integrated native review applied no new A exclusions and three additional B
ambiguities, leaving454 synchronization-stable B events before visibility and
block filtering. These three retain their original central images. With zero
additions, canonical preparation remains byte-identical. With additions, exact
original events/manifest are preserved once under `preparation_raw`, canonical
event status cells alone change, and the manifest records explicit postprocessing
lineage and new event digest/count. `sync_eligibility.json` retains every mapping,
uncertainty and reason. This is a conservative metadata sensitivity check, not
new waveform evidence or independent physical clock truth.
Its focused synthetic suite passed40/40 checks. Job
`sync_eligibility_integrated_02` passed24 workflow controls and applied the actual
A0/B3 changes with exit0 in68.08 s. All A canonical and refinement-plan inputs
remain unchanged. B changed only the three status cells and manifest; immutable
originals match their captured SHA256 values, and frames/cache/sync metadata remain
unchanged. The per-file comparison is retained in
`diagnostics/sync_eligibility_integration_check.json` under the run root.

## Algorithm improvements and remaining evidence

The implemented improvement is method-matched physical calibration: corrected
extrema remain the original comparison, while training-block validation compares
GCC and fixed SRP. Invalid rays retain explicit status and coverage; they never
become clamped pixel locations. The refinement does not add optical flexibility
to compensate for unreliable acoustic delays.

A bounded next experiment would compare four predeclared global GCC frequency
masks (current full band, 1–10, 1–20 and 1–30 kHz), retaining the current envelope
guide, window, interpolation, physical limits, event population and fitting
settings. Each candidate needs its own training-block fits; selection must finish
before any evaluation replay. This experiment is recommended, not implemented or
claimed to improve these recordings. A's already inspected evaluation block is
hypothesis-generating evidence and cannot become a fresh validation set.

The five catastrophic A evaluation events have energy centroids above30 kHz,
but accurate training events also occur there. A centroid cutoff is therefore
unjustified. Small pairwise-delay closure only checks additive consistency; it
does not identify the visible sound source. Competing correlation peaks,
fixed-subband disagreement and onset-versus-ringing delays could help distinguish
cycle ambiguity and multipath. Other acoustic sources remain a separate unresolved
possibility. None of these diagnostics currently rejects or relabels an event.

For trustworthy physical parameters, collect clicker-confirmed events across
more independent time blocks and image heights, with pauses separating emissions
from other sources. Independently measured sound speed or additional array/pose
constraints would reduce the speed/delay/camera tradeoff. Keep an independent
recording for validation after estimator changes. The present narrow temporal and
spatial coverage cannot establish a unique sound speed or general camera model.

## Human correction and final deliverables

The user confirmed that A review files00459.jpg through00701.jpg inclusive in
galleryada438af5848 contain the valid underwater clicker sequence; other images
are outside the calibration scope. This is243 selected frame IDs. The old44
underwater-visible annotations are all within this range;199 additional selected
frames previously had uncertain labels. Object presence is user-confirmed;
coordinates still require individual Terra reinspection.

`annotation_revision.py` creates the separate
`calibration_runs/clicker_20260209/20261030_v2/A` revision, with byte-identical
prepared event/frame/cache data, unchanged source-image references, and explicit
parent-manifest lineage. It resolves the inclusive JPG range through the original
gallery HTML frame identities, preserves the old annotations/manifest under
`revision_source`, clears previous coordinates, and records the selection in
`human_review_selection.json`. Every other row is explicitly excluded, rather
than falsely relabeled absent or air. The original A v1 fits remain unchanged.

```sh
python src/clicker_calibration/annotation_revision.py SOURCE_A NEW_A \
  --first-index 459 --last-index 701 --gallery-revision ada438af5848 \
  --reason 'User-confirmed underwater clicker interval'
```

Selected frames must receive independent Terra reinspection before sealing.
Partially covered but localizable housings can receive a point and conservative
uncertainty; uncertain coordinates are shown for human review and do not silently
become fitting truth. Fully unlocalizable objects remain without invented points.
The review defaults to the user-selected frames. Initial overlays show Terra
coordinates/circles; final comparison overlays show Terra in green and calibrated
acoustic estimates in red, with uncertainty/status retained. No baseline marker
obscures that comparison; baseline statistics remain in the reports.

The requested final delivery root is
`/media/spin/anas2/data/calf/Clicktest/clicker_20260209/results`, with independent
A/B results, parameters, individual original frames and comparison overlays,
plus full synchronized annotated videos. Video export is being implemented and
will retain complete chapter frames, explicit timing, and an aligned physical
channel5 review soundtrack. Final A delivery will use the corrected v2 labels;
v1 remains historical evidence. B processing continues in v1 independently.

The revision helper passed12 synthetic tests; annotation validation/rendering
passed8 tests including selection enforcement and stale-selection detection.
Actual revision creation and native annotation handoff completed successfully;
the A timing review again added zero exclusions.

### Revised A uncertainty and explicit fit eligibility

Terra sealed the corrected A revision with233 estimated housing-center points,
10 unlocalizable selected frames, and548 user exclusions. All243 selected frames
carry the original `uncertain` visibility label; typical radii are35–50 px. The
initial strict fit stopped without saving a model because it accepts only visible
annotations. Its failure is recorded in `v2/jobs/finish_A_revision_01`.

The explicit `user_confirmed_uncertain_coordinates` policy permits fitting finite
points in the authenticated user-confirmed interval with their unchanged
uncertainty weights. It does not turn them into new Terra truth. A separate
`20261030_v2_fit/A` directory contains the compatibility input for the unchanged
visible-only fitter, plus `terra_visibility`, eligibility reasons, exact original
truth archives and `annotation_fit_policy.json`. Missing coordinates and all
original exclusions remain excluded. All human galleries and videos continue
to read the unmodified v2 annotations, displaying uncertain estimates in amber.

Use `annotation_fit_policy=:user_confirmed_uncertain_coordinates` with
`fit_input_root=...` in `run_clicker_calibration`; a device-keyed dictionary can
apply it only to A. The default remains strict. The policy requires a sealed
human selection and never infers coordinates or reduces uncertainty. Its14
synthetic tests include233/10/548 eligibility and a strict export-side proof that
every original coordinate/radius is unchanged. This is an exploratory fit to
low-precision visual estimates, not independent ground-truth qualification.

### Synchronized review videos and portable delivery

`video_export.py --annotation-only DEVICE_OUTPUT MEDIA_OUTPUT` produces
`terra_review_synchronized.mp4` and a lossless Float32 mono48 kHz
`synchronized_reference_ch5.wav`, without requiring a calibration model. This
initial review video displays Terra locations/uncertainty circles and synchronized
hydrophone sound, with no acoustic-prediction markers. The calibrated mode takes
`DEVICE_OUTPUT MODEL_OUTPUT MEDIA_OUTPUT` and adds red acoustic estimates.

The exporter retains every source video frame and exact chapter gaps, using
declared audio-duration chapter bases and the saved exact-base timing review.
Output frame count and every PTS are checked. Markers start on their source frame;
a labeled120 ms hold makes them readable without interpolating motion. Invalid,
outside-image, ambiguous-sync and user-excluded estimates are not clamped into
apparently valid locations. The soundtrack uses raw physical channel5, gapless
WAV sample counts and explicit affine sample-time interpolation. It is audible-
band review media; original192 kHz scientific recordings remain unchanged.

Eight short synthetic end-to-end tests passed, covering source-frame timing,
chapter gaps, exact uncertainty-outline pixels, synchronized known impulses
within one48 kHz sample, annotation-only mode, invalid predictions, resumability,
and the real uncertainty-policy validator. Full recording exports are separate
jobs and are not established by these short fixtures.

`package_results.py DEVICE_OUTPUT RESULTS_DEVICE_DIR` copies the initial review,
original frames and scientific artifacts into a portable bundle.
`--model-dir MODEL --video-dir MEDIA --final` additionally requires a completed
matching comparison video. Raw/annotation/comparison images have a frame-index
CSV, and rewritten presentation links resolve inside the bundle. Scientific
tables and manifests remain byte-identical. Unowned or externally edited delivery
files are never silently overwritten. Four focused packaging tests passed.


## Resource policy and immediate review previews

The user requires advance notice of actions and full available CPU/GPU resources
for future runs. Heavy jobs are staggered rather than artificially throttled.
`run_clicker_job.py` now sets Julia/BLAS/OpenMP thread availability from CPU
affinity, uses available RAM at actual launch as a monitoring ceiling, and
serializes heavy commands with a repository queue lock. `--lightweight` is for
short checks, annotation waiting and orchestration; it does not limit CPU access.
The job timeout starts after queue admission. Five launcher tests cover resources,
nonoverlap, failure recording and descendant cleanup.

This machine exposes28 logical CPUs and an RTX4090 with24 GB VRAM. The exporter
defaults to automatic full-CPU selection and qualified NVENC/NVDEC. CUDA decoding
requires one CPU control thread on this driver; the decoding itself uses the GPU.
After download, the original full-range GoPro pixel format is preserved. GPU
qualification matched prepared pixels exactly on all three A chapters; a native
2704x1520 excerpt preserved59/59 frames and timestamps. Ten NTSC/GPU tests passed.
Fresh media/output hashes authenticate saved source PTS without re-decoding all
source frames during preflight; the final output still receives full timestamp
and target-pixel verification. GPU failures are reported with explicit backend
selection and probe logs. Old CPU-limited jobs were intentionally stopped for
replacement; their outputs and interruption receipts remain available.

Review MP4 audio receives one documented constant gain, bounded by512 and a0.6
PCM peak. The aligned Float32 reference WAV remains unscaled. Decoded AAC peaks
are checked for clipping. This supports listening to the quiet hydrophone signal
without altering scientific source recordings.

A and B quick review previews have been delivered under the requested results
root at `A/quick_review/A_terra_quick_review.mp4` and
`B/quick_review/B_terra_quick_review.mp4`. Each contains20 source excerpts,
19.6867 seconds at1352x760, Terra markers/uncertainty and synchronized physical
channel5 audio. All20 target frames matched prepared PNG pixels, with1180 output
frames and zero measured AAC round-trip sample lag. B is explicitly provisional:
its immutable checkpoint contains395 of469 individually inspected frames, and
no final annotation seal or unreviewed labels were invented. The full recordings
remain separate exports. `quick_review_preview.py` now also defaults to automatic
GPU/all-CPU use; its one-event accelerated control passed.

The corrected A uncertainty-weighted run completed in471.57 seconds before this
resource-policy change. It used202 eligible events:163 training and39 reused
evaluation. GCC was selected. Final all-event descriptive median error is44.31 px
and RMSE1755.38 px, with8 outside-image estimates; reused-evaluation median is
176.74 px and RMSE2564.63 px on35/39 projected points. These large failures remain
visible. The final parameters include c1562.177 m/s and FOV75.8906/56.1657 degrees,
but weak identifiability remains. These are exploratory fits to233 low-precision
Terra estimates, not independently measured sound speed or confirmed accuracy.


## Native provisional-observation policy for B

B's sealed469 rows contain2 visible air observations,296 uncertain observations
(234 coordinate-bearing across all media), and171 absent rows. Its calibration
must not silently convert uncertain visual estimates into confirmed positions.
`predeclare_observation_policy(DEVICE_OUTPUT;
policy=:terra_estimated_underwater_coordinates, reason=...)` records an explicit
provisional policy before any fit. It binds the unchanged annotations, completion,
events, frames and fitter source, declares `human_confirmed=false`, and admits
only original underwater visible/uncertain coordinates with finite positive
original radii. All timing, physical, purge and shared-frame exclusions remain.
The shared loader applies this same policy to the first fit, estimator comparison
and refinement. Original visibility/radii/admission are included in evidence;
policy contents and hash are retained in the baseline calibration/report/frozen
evaluation. Without this sidecar, the fitter stays strict visible-only.

The main workflow accepts `observation_policy` as a scalar or device-keyed
mapping. Existing source/input/settings/product mismatches require a new run
version; historical fits are never overwritten. Prepare/SRP/projection code is
unchanged, so A's stored inference contract remains valid. The normal-optimization,
automatic-thread native suite passed139/139 checks. The earlier normal-resource
85-check suite also passed; the old compiler crashes were not reproduced in these
bounded checks. B's provisional-policy fit has not yet been launched at this
checkpoint.
