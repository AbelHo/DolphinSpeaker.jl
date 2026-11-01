# Review: hydrophone sound direction to video pixels

Review date: 2026-10-30. Source snapshot: `d22056c39b949ee5bead6ddd3b9644e4261ab267`.

**Snapshot boundary:** the repository changed concurrently during this review: HEAD advanced to `470776fefc6a7a5bf4ce4fc7e12898be1e2190e2`, restoring a newer `src/clicker_calibration/` workflow, and additional uncommitted integration changes appeared. Findings and line numbers below refer to the original snapshot and the legacy workflow named in the request. The newer path is described separately below; these findings must not be read as a claim that it has no physical calibration implementation. The reproduction script reads definitions from the original Git commit so concurrent edits cannot silently change its subject.

## Assessment and scope

The reviewed legacy Julia workflow produces an approximate acoustic overlay. Its calibration fits two independent pixel-space lines; it does not establish a calibrated physical relationship between hydrophones, camera, and source. The most urgent problems are annotation correspondence, inconsistent frame numbering, and accepting unreliable acoustic estimates. Improving the camera projection before fixing those problems could fit the wrong correspondences more convincingly.

This review covers the parent Julia repository: `src/aspod.jl`, `localization.jl`, `cam_calib.jl`, the relevant device configuration, synchronization/export/rendering helpers, and the experimental camera fitter in `dev2_sync.jl`. It also checks what the newer localization helpers already address. The Python package is outside scope. No production source or recording was changed.

Evidence consists of source inspection and bounded Julia reproductions. The supplied D3 calibration directories under `/media/spin/anas2/data/...` and `/media/spin/anas2/data_res/...` are absent in this environment. Consequently, this review does **not** measure real calibration accuracy, verify the current annotations, or establish which conditional defects occurred in that recording. Historical accuracy/coverage numbers in the improvement READMEs were not rerun.

Supporting artifacts:

- [Runnable source checks](audio_to_video_projection_checks_2026-10-30.jl), which extract selected definitions without executing the workflow scripts.
- [Recorded check output](audio_to_video_projection_checks_2026-10-30.txt).

## What the active algorithm does

| Stage | Current behavior | Source |
|---|---|---|
| Device setup | D3 clicker uses three coplanar receivers with 0.18 m triangle sides; channels 1–2 have `NaN` geometry. Sets image size `(1520,2740)` and overwrites horizontal FOV with 80 degrees. | `src/config.jl:138–163` |
| Acoustic estimate | Selects click-train detections; estimates interchannel arrival times; fits azimuth/elevation by minimizing squared sample-delay residuals. The default method remains `get_tdoa_minmax`. | `src/aspod.jl:165–188`; `src/localization.jl:197–228,684` |
| Pixel estimate | Converts radians to degrees, independently scales azimuth/elevation by horizontal/vertical FOV, and offsets to the image center. | `src/localization.jl:238–245` |
| Video time | Rounds local detection time to a frame and adds 1; subsequently adds a synchronization offset and cumulative file duration in frame units. | `src/aspod.jl:188`; `src/run_example.jl:396–405` |
| Annotation video | Saves detection-frame PNGs, then sorts the filenames into a 4 fps montage. | `src/video.jl:196–247`; `src/pic2vid.jl:54–74`; `src/run_example.jl:493–496` |
| Pixel calibration | Replaces frame identity with detection-row ordinal, joins annotations, removes tag `o`, and fits `x ~ px`, `y ~ py`. | `src/cam_calib.jl:83–150` |

The acoustic model is internally interpretable. With azimuth `a`, elevation `e`, receiver position `r_i`, and sampling frequency `fs`, it uses:

```text
s(a,e) = [sin(a) cos(e), sin(e), cos(a) cos(e)]
predicted_delay_i = -fs/c * (r_i - r_ref) · s
```

The current pixel map, using angles in degrees, is:

```text
u = W/2 + W/FOV_horizontal * a
v = H/2 - H/FOV_vertical   * e
```

There is no explicit camera-to-array rotation, lens model, source range, or housing-refraction model in this active projection.

## Findings, in recommended repair order

### 1. High: annotation pairing depends on an undocumented row-order invariant

**Confirmed logic defect when montage order and detection-row order differ; occurrence in the supplied dataset is unverified.**

`calibrate_from_annotations` filters by `dp.frame`, then creates `dp.frame_index = 0:(length(dp.frame)-1)` and joins on that new column (`src/cam_calib.jl:85–112`). It is joining annotation indices to **detection row numbers**, not to original video frames. This can be intentional for an `_overlaidIMG` montage, but only if every retained detection creates exactly one image in exactly that order.

The montage producer does not guarantee this:

- PNG names contain the rounded source frame and time, without a detection ID (`src/video.jl:241–246`). Two detections in one frame overwrite the same PNG.
- Failed reads are skipped; encountering a frame beyond the duration stops extraction (`src/video.jl:207–220`).
- `pic2vid` sorts filenames numerically (`src/pic2vid.jl:54–59`), while detection CSV rows are grouped by detection type and then audio file (`src/run_example.jl:398–405`). Impulse/tonal rows can therefore be reordered.
- The existing output image directory is reused without an exact image manifest. Old images can enter a later montage.

Example reproduced in the checks: detection frames `[100,100,200]` produce montage frames `[100,200]`; montage index 1 should refer to source frame 200, but the calibration joins it to the second detection at frame 100. Using annotations from the original continuous video also fails because their frame indices are not detection ordinals. A robust regressor cannot repair systematic pairing errors.

**Improve:** write a montage manifest at export time containing `montage_frame_index`, source video identity, source frame/PTS, audio file and sample index, detection ID, and method/type. Join through that manifest. Deduplicate source frames deliberately and associate all detections with the same annotated source point, or create uniquely named images per detection. Validate join cardinality, duplicate annotations, unmatched rows, and annotation-video identity. Preserve both original and montage frame indices. Do not reconstruct old mappings from row order unless equivalence has been checked.

### 2. High: the video renderer and exporter disagree by one frame; time is rounded too early

**Confirmed source inconsistency, with numerical reproduction.**

`process_detections` exports `round(t*fps)+1` (`src/aspod.jl:188`). The image extractor interprets that as one-based and seeks to `(frame-1)/fps` (`src/video.jl:207–208`). In contrast, `overlay_annotations_on_video` takes CSV frame numbers unchanged, then compares them with a loop starting at zero (`src/video.jl:436,596–607,699–710`). An event at `t=0` is exported as frame 1: the image extractor shows the first frame, but the continuous overlay draws on the second frame. At 30 fps the discrepancy is 33.3 ms. The final frame can also fall outside the renderer's loop.

Separately, local time is rounded before adding sync/file offsets, and readers round again. For `t*fps=0.49` and `offset*fps=0.49`, the current two-stage calculation selects one-based frame 1; rounding the complete timestamp selects frame 2. This error is separate from the frame-origin mismatch.

**Improve:** carry audio sample indices and continuous video timestamps through processing. Apply the complete time transform before frame selection. Define one frame convention and explicit conversions at renderer boundaries. Prefer actual video presentation timestamps; for verified constant-rate video, document whether the rule chooses the containing frame or nearest exposure timestamp. Test the first/last frames, fractional offsets, file boundaries, and agreement between montage and continuous overlay. Define pixel origins separately; the current map reaches `W`/`H` at FOV boundaries, so it is not a complete zero-based pixel-center convention either.

### 3. High: the default TDOA estimator's polarity decision is wrong

**Confirmed and reproduced; a corrected helper already exists but is not the default.**

`get_tdoa_minmax` tests `max_val > min_val` (`src/localization.jl:318–328`). For any nonconstant waveform this is true, irrespective of which signed excursion has greater magnitude. It therefore chooses each channel's positive maximum even when the dominant pulse is negative. In the supplied synthetic counterexample it selects indices `[19,21,23]` instead of the dominant negative peaks `[10,11,12]`, changing the inferred interchannel delays.

`get_tdoa_minmax2` fixes the magnitude vote and clips windows (`src/localization.jl:659–666`), but line 684 still sets the old function as the default. The old function also throws a `BoundsError` for an edge detection whose window crosses the recording boundary. The wrapper catches processing failures at audio-file level (`src/run_example.jl:383–388`), so an event-level failure can discard the file's output.

**Improve:** fix or retire the old implementation and make boundary handling consistent. Treat signed-peak selection as a baseline, not a reliable direct-path estimator: polarity correction alone does not prevent different lobes, echoes, or overlapping clicks from being selected across channels. Compare the existing bounded GCC/SRP alternatives on independently paired clicker data before selecting a default. Preserve event IDs across methods so their comparisons use the same examples.

### 4. High: the ordinary overlay converts physically implausible delays into ordinary points

**Confirmed; the alternative multi-TDOA workflow partly addresses it.**

The default solver runs unconstrained Nelder–Mead from `[0,0]`, returns the minimizer, and does not check convergence (`src/localization.jl:208–220`). The active overlay requests residuals but projects the angles without using those residuals to reject or flag points (`src/aspod.jl:179–188`). A synthetic `[100,0,0]` sample-delay vector for the D3 geometry has an unconstrained direction-cosine norm of 5.145, yet the solver returns finite angles, with residual about 6118 samples squared.

For a planar array, solve for `s_x,s_y` first and check `s_x²+s_y² <= 1`, with tolerance derived from timing uncertainty. The D3 0.18 m baseline permits only about 22.44 samples of physical pair delay at 192 kHz and 1540 m/s. The old correlation path searches without that geometry-derived limit (`src/localization.jl:40–58`). Pairwise limits alone are insufficient: all selected pairs must also fit one direction.

**Improve:** return a localization result containing validity, residuals, uncertainty, ambiguity, and rejection reason. Retain invalid detections in the event table without drawing them as confident source positions. Restrict correlation searches to physical lag bounds, retain plausible competing peaks, and compare independent pair delays for closure. The newer `flag_clicks`/GCC/SRP path already supplies useful pieces (`src/run_multitdoa.jl:24–65`); reuse them. Merely switching to `tdoa2dir_ls` is insufficient because that helper clamps infeasible vectors to the unit circle (`src/localization.jl:638–653`).

A small residual or good closure is **not** proof of a correct direct path: a coherent reflection can satisfy the same geometry. Three coplanar receivers provide only two independent relative delays, so a feasible two-dimensional solution can fit even an incorrect arrival selection exactly. Deriving all pair differences from a single fitted arrival vector makes closure zero by construction; quality checks must use independent measurements.

### 5. High accuracy risk: independent linear angle scaling is an unverified camera model

**Confirmed approximation; its actual error depends on the camera's projection and recording mode.**

For a distortion-free pinhole camera aligned with the source coordinate convention, projecting the current direction vector gives:

```text
fx = W / (2 tan(FOV_horizontal/2))
fy = H / (2 tan(FOV_vertical/2))
u = cx + fx tan(a)
v = cy - fy tan(e)/cos(a)
```

The final `1/cos(a)` factor matters: elevation is not an independent vertical image-plane angle. These equations follow by substituting the repository's direction vector into perspective division; they are an analytical comparison with the standard [OpenCV pinhole model](https://docs.opencv.org/4.13.0/d9/d0c/group__calib3d.html).

Using the configured `W=2740`, `H=1520`, horizontal FOV 80 degrees, and vertical FOV 54.4609 degrees:

| Azimuth, elevation | Current pixels | Hypothetical aligned pinhole pixels | Current minus pinhole |
|---|---|---|---|
| 20°, 0° | (2055.000, 760.000) | (1964.255, 760.000) | (+90.745, 0.000) px |
| 30°, 0° | (2397.500, 760.000) | (2312.641, 760.000) | (+84.859, 0.000) px |
| 30°, 20° | (2397.500, 201.801) | (2312.641, 139.306) | (+84.859, +62.495) px |

**These are model discrepancies, not measured GoPro errors.** Do not blindly replace this with a pinhole formula. A fisheye model uses radial angular projection, which also differs from separate linear azimuth/elevation scaling; see the [OpenCV fisheye equations](https://docs.opencv.org/4.13.0/db/d58/group__calib3d__fisheye.html).

**Improve:** identify the exact recording mode, crop/resize, stabilization behavior, and housing. Calibrate the camera and choose an appropriate pinhole, fisheye, or refractive model from held-out reprojection error. Rotate the acoustic unit vector into camera coordinates before projection. Obtain actual image dimensions from metadata, validate them, and transform calibration parameters when images are resized/cropped. Currently `angle2px` uses global `imsize`, and the D3 configuration silently replaces the measured-looking horizontal FOV 71.5412 with 80 degrees (`src/config.jl:149–160`); neither the dimensions nor that replacement is checked against the input video.

### 6. High accuracy risk: a planar array cannot distinguish the two sides of its plane

**Confirmed geometric ambiguity, reproduced with the actual cost function.**

All D3 receiver z coordinates are zero. Directions `(a,e)` and `(π-a,e)` have identical x/y components and therefore identical predicted delays. The source check gives zero cost for both `(30°,10°)` and `(150°,10°)`. An unconstrained optimizer can return either branch; reducing angles modulo `2π` does not select the physical hemisphere.

**Improve:** record front/back ambiguity explicitly. If calibration or deployment guarantees the source lies on the camera side of the array, enforce that as a declared prior and reject cases violating its assumptions. Otherwise additional geometry or independent evidence is needed; a fourth receiver outside the plane is one hardware option. Video association can supply a prior, but selecting the branch from the same annotation later used as the test target would bias evaluation. Do not interpret a camera-side assumption as an acoustic measurement.

### 7. Medium: sound speed and device state are implicit and partly disconnected

**Confirmed configuration defect.**

`src/config.jl:13` defines `c=1540`, but `cost_tdoa2ang` defaults to the separate `default_c` defined in `src/localization.jl:11–12,197`. Changing `c` alone does not change the ordinary solver. The reproduction sets `c=1480`: the default cost stays unchanged, whereas explicitly passing `c=1480` changes it. `tdoa2dir` does not expose a sound-speed keyword, so the intended calibration parameter must otherwise be passed through a custom cost closure.

Other state is also inherited: the clicker setter does not reset every option changed by the dolphin setter, including `PARAM_TONALSPREAD` (`src/config.jl:138–195`). Calling setters in a different sequence can therefore change behavior or the population entering calibration.

**Improve:** use one explicit configuration value for sound speed, units, geometry, channel map, reference channel, sample rate, detector settings, and projection calibration. Pass it through all estimator and residual paths and serialize the effective values with results. Measure geometry and use an appropriate measured/constrained sound speed; do not use FOV to absorb a sound-speed error. If fitting per-channel delays, fix one delay to zero because only relative offsets are observable.

### 8. Medium: the envelope alternative uses the wrong sample rate at 192 kHz

**Confirmed; a corrected helper exists.**

`envelope_lpfiltered` defaults to `fs=500_000` (`src/localization.jl:149`). `get_tdoa_envelope_filtered`, used for Impulse3, invokes it without forwarding the recording rate (`src/localization.jl:152–162`; `src/aspod.jl:237–239`). For a 192 kHz recording, the nominal 10 kHz low-pass consequently corresponds to 3.84 kHz in the actual signal. This can change delay estimates and bias comparisons between methods.

**Improve:** thread `fs` through the complete call chain, using `get_tdoa_envelope_filtered2` or its corrected behavior (`src/localization.jl:669–678`). Test that sample-rate changes preserve the intended physical bandwidth. Express analysis windows and tolerances in seconds where appropriate, then derive sample counts. For context, one sample on a 0.18 m baseline at 192 kHz corresponds to about 2.55° near broadside; sub-sample estimation can matter, but interpolation does not establish accuracy in multipath.

### 9. Medium: two independent regressions cannot express the needed geometric corrections

**Confirmed model limitation.**

`x=a_x+b_x*px`, `y=a_y+b_y*py` can correct independent scale and offset only (`src/cam_calib.jl:122–142`). They cannot represent roll, cross-axis coupling, perspective curvature, radial distortion, or range-dependent parallax. Even a small in-plane rotation needs terms involving both coordinates. Because the upstream map is linear in angles, these fits are also just separate affine functions of azimuth and elevation.

**Improve:** retain this model as a clearly named empirical baseline. A two-dimensional affine map is a useful intermediate comparison, but a physical direction-to-camera model is preferable for portability. Fit roll/yaw/pitch in the unit-vector transform, and camera parameters in the camera model. For noisy TDOAs, compare fitting measured delays against camera-derived directions with fitting already estimated pixels; regressing noisy estimated pixels as if they were exact predictors can distort slope estimates.

For a flat-port underwater housing, refraction can invalidate a single-center camera approximation and make calibration range-dependent. Calibrate in the deployed housing and medium. Escalate to an explicit refractive model if held-out errors require it; the exact housing type is unverified here. This is supported by [Jordt-Sedlazeck and Koch, Refractive Calibration of Underwater Cameras](https://www.geomar.de/fileadmin/personal/fb2/mg/ajordt/finalEccv.pdf).

### 10. Medium: the acoustic model assumes far-field sources and a shared camera/array viewpoint

**Conditional modeling risk; source distances and camera offset were not available.**

The delay formula is a plane-wave model, while `angle2px` also implicitly treats the acoustic bearing as a camera bearing. At finite distance, the exact delay is:

```text
delay_i [samples] = fs/c * (norm(P-r_i) - norm(P-r_ref))
                   + fs * (channel_delay_i - channel_delay_ref)
```

This reduces to the current signed plane-wave expression in the far field. The camera sees `R_CA*P + t_CA`, not just a rotated unit vector. Without source range, translation between camera and array cannot generally be corrected from direction alone. As an illustrative scale, a 0.10 m transverse offset at 1 m range gives a 5.71° bearing difference; this is not a measurement of this device.

**Improve:** measure the camera-array transform and calibration-target distances. Quantify finite-range delay and parallax errors over the intended distance range. Use a spherical-wave model and range/known target positions when necessary; otherwise document the validated far-field range and include the remaining uncertainty. Three receiver arrival differences do not independently determine an arbitrary 3D source position. Do not add range as a free per-click parameter and assume it becomes identifiable.

### 11. Medium: synchronization is one offset, with no drift or continuity validation

**Confirmed implementation assumption; drift/gaps in the actual data are unverified.**

`run_analysis_split_vidau` derives one delay from the first video/audio pair and writes confidence (`src/run_example.jl:329–336`). The writer then adds file durations as though audio files are contiguous (`src/run_example.jl:400–405`). It uses the first video's frame rate; `get_fps` reads `r_frame_rate` (`src/media_info.jl:5–20`). The confidence is not used there to reject a bad synchronization. Recursive discovery can also collect unrelated subrecordings; the localization README already notes this for `not_so_good/` in the clicker folder (`improvements/tdoa_localization_2026-10-21/README.md:247`).

**Improve:** use an explicit recording manifest and check ordering, overlaps, gaps, frame-rate compatibility, and common coverage. Estimate synchronization at several times and fit `t_video = offset + scale*t_audio` if recorder-clock drift is supported. Reject ambiguous matches and inspect residuals. Carry per-segment transforms if recordings are discontinuous. Account for differences in acoustic propagation paths to the synchronization microphones when source distance makes them relevant. A spatial calibration should not compensate for a timing error on a moving target.

### 12. Medium: reported calibration errors are training errors after manual exclusions

**Confirmed evaluation limitation.**

The same joined rows fit the regressions and supply RMS/median/mean errors (`src/cam_calib.jl:122–150`). Tag `o` is excluded beforehand, without a reported rejection denominator or exclusion breakdown. There is no held-out split, raw-error baseline in `stats`, uncertainty, or breakdown by position/time/range. Nearby clicks and repeated frames are correlated, so a random click split could also exaggerate generalization.

**Improve:** report raw and calibrated pixel errors on held-out recordings or independent clicker poses/trials; freeze fitting, method selection, and quality thresholds before evaluating them. Include median, RMS, p90/p95, signed x/y bias, accepted/rejected counts, and matching coverage. Report quality-gated accuracy alongside coverage. Keep manual reasons such as invisible target, non-clicker sound, reflection ambiguity, and poor audio separate; do not remove difficult but correctly labeled cases merely for large residuals. Sample annotation frames independently of overlay closeness, ideally with predictions hidden, and group uncertainty estimates by trial/pose rather than treating adjacent clicks as independent.

To report angular error, unproject annotations using an independently calibrated camera model, then compare unit vectors with a clipped dot product. Simply inverting the current linear FOV mapping assumes the disputed projection is true. Include annotation uncertainty and target identity; pixel accuracy does not require pretending the target lies on a known metric plane.

### 13. Medium, experimental path: the separate camera fitter currently fails even on exact synthetic pairs

**Confirmed and reproduced. This fitter is not called by the ordinary overlay/calibration path.**

`pixel_to_world_correction` is a camera-ray intersection with a world plane, not an acoustic-bearing projector (`src/dev2_sync.jl:297–340`). Its transform algebra is reasonable under that stated assumption. However, the fitting implementation has several distinct problems:

1. `_camera_model(p,xdata)` reverses LsqFit's callback order, which is `model(xdata,p)` (`src/dev2_sync.jl:407,468,483,506`; [LsqFit documentation](https://github.com/JuliaNLSolvers/LsqFit.jl)). With eight point pairs and 15 parameters, the wrong order makes the function return 14 values for 16 observations. The actual `curve_fit` call raises `DimensionMismatch`; the full fitter's fallback fails too. This was checked with the installed LsqFit 0.15.1.
2. Fixing argument order still leaves `Vector{Float64}` output buffers (`src/dev2_sync.jl:346–365`), incompatible with ForwardDiff dual-valued outputs. A direct Jacobian evaluation reproduces `MethodError`.
3. `rodrigues` handles zero rotation with a special branch around `norm(r)` and division by that norm (`src/dev2_sync.jl:379–388`), at the fitter's initial value. Its ForwardDiff Jacobian at zero fails to reproduce the finite, nonzero central-difference Jacobian. Use a differentiable small-angle expansion, including a numerically safe treatment of the norm near zero.
4. RANSAC uses `rand(1:N,sample_size)`, which samples with replacement (`src/dev2_sync.jl:463`). Duplicated points reduce independent constraints. It lacks adequate degeneracy/minimum-inlier checks, parameter bounds, safe ray-plane denominator checks, and explicit handling of an empty inlier set. The default `inlier_thr=0.05` has world-coordinate units; it would mean 0.05 pixels if pixel annotations were supplied as `x,y`.

**Improve:** first extract this into a side-effect-free module and verify callback shape, differentiation, and synthetic parameter recovery. Use unique samples, deterministic seeds, meaningful residual units, bounded parameters, and explicit failure reporting. Do not treat fitted `K,R,t` as a physical calibration when both inputs and targets are pixel coordinates. Even genuine correspondences on one plane do not uniquely determine all free camera intrinsics and pose: without distortion, they constrain a homography. Use independent camera calibration/multiple target poses or constrain the unidentifiable parameters. For this task, a ray-to-pixel projection is the relevant primitive; a pixel-to-world-plane mapper is not a drop-in replacement.

### 14. Lower priority: calibration API and output behavior can mislead callers

**Confirmed from source.**

- `vidfname` and `annotated_csv` are optional positional arguments, but the examples and CLI pass them as keywords (`src/cam_calib.jl:54,68–74,405`). That call does not match the function signature. The standalone file also does not supply its media-info helpers.
- `annotated_frame_col` is unused; `frame` is unconditionally renamed. Defaults permit `nothing` paths even though the body immediately reads them (`src/cam_calib.jl:85–95`).
- `flag_save_csv=true` references undefined `out_prefix` before the writer's `try` block (`src/cam_calib.jl:202–209`). `write_csv=true` alone does not save output because `flag_save_csv` defaults to false.
- Corrections are calculated only for annotated/matched rows after `dropmissing!`. The active `process_detections` path does not load or apply these fitted coefficients to future detections; current calibration is an offline diagnostic, not an automatically deployed correction.

**Improve:** use explicit required inputs, honor column arguments, validate finite values/sample coverage, fix the CLI and output path, and separate fit from apply. Save a versioned calibration artifact with coordinate convention, camera mode, geometry, sound speed, channel mapping, method, training recordings, and held-out metrics. Apply it only to compatible inputs and preserve raw outputs for comparison.

## Proposed calibration design

Use separate acoustic, temporal, and optical transforms with explicit conventions:

```text
audio samples -> candidate arrival times -> valid direction(s) in array coordinates
                                                |
                                      camera-array rotation
                                                |
                                  camera/housing projection -> (u,v)

audio sample index -> synchronized video timestamp -> selected frame/PTS
manual montage index -> export manifest -> same source frame and detection ID
```

For a validated far-field, central-camera approximation, fit `R_CA` and permitted relative delay offsets against calibrated camera rays. Predict a relative arrival delay with:

```text
s_array = transpose(R_CA) * camera_ray(annotation_pixel)
delay_i = -fs/c * (r_i-r_ref) · s_array + fs*(delta_i-delta_ref)
```

Set `delta_ref=0`; state the delay sign and units. Weight residuals by measured timing uncertainty and use a robust loss across observations. A common-reference delay vector has correlated errors, so use its covariance when available rather than assuming all pairs are independent. For near-field or noncentral refractive geometry, replace the direction-only approximation with source positions and the appropriate ray/path model.

Constrain the fit in stages: camera intrinsics/housing first, measured receiver geometry and channel order second, synchronization third, array-to-camera alignment and identifiable delay offsets next. Relax sound speed or other parameters only when independent observations constrain them. Simultaneously freeing FOV, geometry scale, sound speed, rotation, delays, and time offset can create compensating parameter combinations with low training error but poor transfer. Inspect Jacobian rank/conditioning and held-out performance, not just optimizer convergence.

The existing `get_tdoa_minmax2`, `get_tdoa_envelope_filtered2`, bounded GCC/SRP helpers, `tdoa2uv`, and multi-TDOA quality flags provide useful starting points. Their existence does not mean the standard overlay has adopted them or that their field accuracy is established.

## Implementation and validation sequence

| Order | Work | Evidence required before using the result |
|---|---|---|
| 1 | Fix event/frame identity and time conventions; export the montage manifest. | Exact round trips for duplicate/unsorted frames, skipped images, stale output, both renderers, fractional sync offsets, and segment boundaries. |
| 2 | Repair default polarity/window handling and sample-rate propagation; make configuration explicit. | Delayed broadband pulses with known signs, channel permutations, multiple sampling rates, edge events, and nondefault sound speed. |
| 3 | Introduce validity/ambiguity/uncertainty outputs and compare existing delay methods. | Controls with impossible delays, front/back twins, competing echoes, and low-SNR events; accuracy and rejection coverage on the same clicker examples. |
| 4 | Calibrate camera/housing and camera-array alignment; assess near-field and drift effects. | Independently identified points spanning both image axes, edges and center, several distances/poses, and the exact deployment mode. |
| 5 | Fit on selected trials, freeze the model, evaluate on separate trials, then integrate apply-time calibration. | Raw versus corrected held-out pixel/angular errors, coverage, confidence intervals, and residual plots against position, time, range, and method. |

Useful invariants include common-delay-offset cancellation, receiver-coordinate translation invariance, equivalent results after consistent channel/geometry permutation, and projection/unprojection round trips for valid rays. Synthetic success establishes implementation correctness under its assumptions; it does not establish direct-path selection or real underwater accuracy. Likewise, smooth tracks, low training residuals, or a larger fraction of points inside the frame are insufficient accuracy measures.

## Concurrent newer workflow: avoid duplicating implemented work

A bounded inspection of the concurrently restored/edited files found the following relevant changes. This is an implementation inventory, not a completed review or execution of the new workflow.

| Review concern | Newer implementation observed |
|---|---|
| Annotation identity and frame time | `src/clicker_calibration/prepare.jl:5–10,471–478` records event/frame IDs and source PTS; `physical_fit.jl:561–567` checks annotation frame identity/uniqueness. |
| Synchronization drift | `prepare.jl:262–280` fits an affine waveform-derived time mapping and reports validity/residuals. |
| Physical projection and invalid delays | `physical_fit.jl:25–84` applies camera rotation and perspective/radial projection, with explicit invalid states and an assumed positive-z array branch. |
| Held-out evidence | `physical_fit.jl:712–747` separates a training-only evaluation calibration from the final all-eligible deployment refit. Only the former supplies held-out evidence. |
| Optional legacy integration | New A/B presets and `preset_bridge.jl:40–90` load calibration parameters, set both sound-speed globals, select corrected peak extraction, and provide calibrated direction/pixel functions. Modified `localization.jl` routes through them when calibration is active. |

The ordinary `set_device__ophk_acoustic_D3_clicker()` now explicitly clears the active calibration; the newly added `_A`/`_B` presets opt into the calibrated route. Thus the exact legacy call shown in the request still needs to be distinguished from these new entry points. Retain the present review as the baseline explanation and use the newer implementation as the starting point for improvements rather than building a competing calibration subsystem.

Remaining integration questions include whether all output consumers handle invalid/NaN projections, whether old montage/continuous renderers retain their frame-origin mismatch, and whether deployment metadata validates camera mode and timing as well as geometry. The new physical projector explicitly remains far-field and direction-only, so finite-range, housing-refraction, front/back-prior, and field-validation limitations still require evidence. No fitted A/B calibration artifacts, newer-workflow test results, or real-recording accuracy were verified in this review. Concurrent source edits were preserved.

## Reproduction and limits

The bounded checks use Julia 1.12.6, Optim 1.13.3, LsqFit 0.15.1, and ForwardDiff 1.4.6. From the repository root:

```bash
julia --startup-file=no --project=. reviews/audio_to_video_projection_checks_2026-10-30.jl
```

The review environment's `julia` launcher could not create its configuration lock on the read-only home filesystem. The checks therefore used the installed binary directly with `--compiled-modules=existing`; no environment or package installation was changed:

```bash
/home/spin/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia --startup-file=no --compiled-modules=existing --project=. reviews/audio_to_video_projection_checks_2026-10-30.jl
```

The checks reproduce defects in the recorded baseline intentionally; they are review evidence, not proposed permanent regression expectations after those defects are fixed. The camera/projection/acoustic checks evaluate selected actual source definitions read with `git show` at the review commit. Frame/montage examples reproduce the indexing arithmetic and filename identity rules without decoding media. No full workflow, video-rendering, annotation refit, or field-accuracy evaluation was run.
