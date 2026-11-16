# DolphinSpeaker.jl — Clicker Camera–Audio Calibration Integration Status

**Repository:** `hugo-67/DolphinSpeaker.jl`  
**Branch reviewed:** `dev2`  
**Review date:** 2026-11-16  
**Scope:** `src/clicker_calibration/` integration with the main audio localization and video-overlay pipeline.

> This document summarizes the repository inspection and existing implementation. No source code was modified or tests executed during the review.

## 1. Executive Summary

**Status: Implemented and optionally selectable, with important limitations.**

The camera–audio calibration system in `src/clicker_calibration/` has been integrated into the main DolphinSpeaker.jl code through a calibration bridge.

The existing pipeline can load a fitted calibration profile that modifies:

- TDOA estimation (method-matched GCC or SRP).
- TDOA-to-direction conversion using fitted timing offsets, sound speed and geometry.
- Direction-to-camera-pixel conversion using fitted camera orientation and optical parameters.
- Image dimensions and related device configuration.

Two device-specific fitted profiles, **A (GCC)** and **B (SRP)**, are committed.

Calibration is not active by default. The caller explicitly selects the profile using an A/B device setter or `with_clicker_calibration()`.

**Important caveats:**

1. Some alternative localization methods bypass the fitted estimator and/or direction solver.
2. Current fitted profiles are provisional and exhibit significant outliers and weak physical identifiability.
3. Profiles are specific to the original D3 camera/hydrophone mounting and should not be assumed transferable.
4. The review confirms source-code integration and the presence of tests, not successful execution of the full pipeline in this review.

## 2. Repository File Map

All hyperlinks below are relative to the repository root.

### 2.1 Main integration and runtime

| File | Purpose |
|---|---|
| [src/DolphinSpeaker.jl](../src/DolphinSpeaker.jl) | Package entry point; exports calibrated A/B setters, scoped activation and calibration workflows |
| [src/config.jl](../src/config.jl) | Legacy device settings; includes calibration bridge |
| [src/clicker_calibration/preset_bridge.jl](../src/clicker_calibration/preset_bridge.jl) | Primary integration layer: profile loading, activation, restoration and runtime adapters |
| [src/localization.jl](../src/localization.jl) | TDOA estimators, direction solvers, `detection2angle()`, `tdoa2dir()`, `angle2px()` |
| [src/aspod.jl](../src/aspod.jl) | Detection/localization processing and conversion to pixel coordinates |
| [src/run_example.jl](../src/run_example.jl) | Main processing entry points, including `run_contiguous_folders()` |
| [src/video.jl](../src/video.jl) | Annotation rendering and video overlay functions |
| [src/run_multitdoa.jl](../src/run_multitdoa.jl) | Alternative multi-TDOA localization methods; not fully calibration-aware |

### 2.2 Calibration implementation

| File | Purpose |
|---|---|
| [src/clicker_calibration/ClickerCalibration.jl](../src/clicker_calibration/ClickerCalibration.jl) | Calibration module, includes and exports |
| [src/clicker_calibration/workflow.jl](../src/clicker_calibration/workflow.jl) | Workflow integration and `load_calibration()` validation |
| [src/clicker_calibration/physical_fit.jl](../src/clicker_calibration/physical_fit.jl) | Physical fitting, TDOA projection and camera projection |
| [src/clicker_calibration/method_refinement.jl](../src/clicker_calibration/method_refinement.jl) | GCC/SRP estimator specifications and method-matched refinement |
| [src/clicker_calibration/prepare.jl](../src/clicker_calibration/prepare.jl) | Recording/event preparation and calibration input processing |
| [src/clicker_calibration/srp_review.jl](../src/clicker_calibration/srp_review.jl) | SRP method implementation and review |
| [src/clicker_calibration/sync_eligibility.jl](../src/clicker_calibration/sync_eligibility.jl) | Synchronization eligibility and uncertainty checks |
| [src/clicker_calibration/full_workflow.jl](../src/clicker_calibration/full_workflow.jl) | End-to-end calibration workflow |
| [src/clicker_calibration/event_projection.jl](../src/clicker_calibration/event_projection.jl) | Projection of prepared events using a selected fitted model |
| [src/clicker_calibration/support.jl](../src/clicker_calibration/support.jl) | Shared I/O, JSON utilities, hashing and logging |
| [src/clicker_calibration.jl](../src/clicker_calibration.jl) | Standalone calibration entry point |

Additional Python scripts in `src/clicker_calibration/` support annotation review, diagnostics, result packaging and media export.

### 2.3 Committed calibration profiles

```text
src/clicker_calibration/calibrations/
├── clicker_20260209_A.json
└── clicker_20260209_B.json
```

- [Device A calibration](../src/clicker_calibration/calibrations/clicker_20260209_A.json) — GCC.
- [Device B calibration](../src/clicker_calibration/calibrations/clicker_20260209_B.json) — SRP.

These profiles contain fitted geometry, timing offsets, camera intrinsics, camera-array rotation, sound speed, estimator specifications and provenance information.

### 2.4 Documentation and tests

| File | Purpose |
|---|---|
| [docs/clicker_calibration.md](../docs/clicker_calibration.md) | Main calibration documentation, results and limitations |
| [docs/clicker_calibration_quickstart.md](../docs/clicker_calibration_quickstart.md) | Quick-start calibration workflow |
| [reviews/plans/persistent-clicker-AB-setters_2026-11-02.md](../reviews/plans/persistent-clicker-AB-setters_2026-11-02.md) | Design rationale for persistent A/B setters |
| [test/clicker_calibration/test_presets.jl](../test/clicker_calibration/test_presets.jl) | Calibration activation, restoration and projection integration tests |
| [test/clicker_calibration/test_physical_fit.jl](../test/clicker_calibration/test_physical_fit.jl) | Physical-model tests |
| [test/clicker_calibration/test_method_refinement.jl](../test/clicker_calibration/test_method_refinement.jl) | Estimator refinement tests |
| [test/clicker_calibration/test_full_workflow.jl](../test/clicker_calibration/test_full_workflow.jl) | Full workflow tests |

## 3. How to Select a Calibration

### 3.1 Persistent activation

```julia
using DolphinSpeaker

# Device A, GCC
set_device__ophk_acoustic_D3_clicker_A()

run_contiguous_folders(folder; res_dir=res_dir)
```

Select device B using:

```julia
set_device__ophk_acoustic_D3_clicker_B()
```

### 3.2 Custom fitted calibration JSON

```julia
set_device__ophk_acoustic_D3_clicker_A(
    calibration_path="path/to/A/calibration.json"
)
```

The JSON must pass calibration validation and identify device A.

The equivalent B setter requires a device B calibration.

### 3.3 Scoped activation

```julia
using DolphinSpeaker

with_clicker_calibration(:A) do
    run_contiguous_folders(folder; res_dir=res_dir)
end
```

This restores the previous configuration when the block exits, including on an exception.

Avoid leaving asynchronous tasks running beyond the block.

### 3.4 Legacy, uncalibrated configuration

For the comparable D3 clicker configuration:

```julia
set_device__ophk_acoustic_D3_clicker()
```

For the regular D3 detector settings:

```julia
set_device__ophk_acoustic_D3()
```

The second setter also changes detection parameters and is therefore not an otherwise identical calibration comparison.

### 3.5 Check active profile

```julia
cal = DolphinSpeaker.active_clicker_calibration()

if cal === nothing
    println("Calibration inactive")
else
    @show cal["device"]
    @show cal["method"]
    @show cal["fit_status"]
end
```

**Activation lifecycle:** A subsequent `set_device__*()` call deactivates the fitted profile. Set the calibrated device last, immediately before processing.

## 4. Runtime Integration Architecture

```text
User selects device A/B calibration
              |
              v
  preset_bridge.jl
  - validates selected JSON
  - sets D3 device configuration
  - stores active calibration
  - installs fitted TDOA dispatcher
              |
              v
  run_contiguous_folders()
              |
              v
  process_detections()
              |
              v
  detection2angle()
              |
              +--> Default TDOA estimator
              |    - A: GCC
              |    - B: SRP
              |
              v
  tdoa2dir()
  - fitted channel timing offsets
  - fitted sound speed
  - calibrated geometry
  - physically valid direction
              |
              v
  angle2px()
  - camera-array rotation
  - pinhole projection
  - camera intrinsics
              |
              v
  Pixel coordinates
              |
              v
  detection_pixels.csv
              |
              v
  Video annotation overlay
```

### Relevant runtime dispatch

In [src/localization.jl](../src/localization.jl):

- `tdoa2dir()` checks `active_clicker_calibration()`.
- `angle2px()` checks `active_clicker_calibration()`.
- `detection2angle()` uses `default_getTDOA_func`, replaced with a calibration-aware dispatcher during activation.

In [src/clicker_calibration/preset_bridge.jl](../src/clicker_calibration/preset_bridge.jl):

- `active_clicker_calibration()` determines whether the fitted profile remains active.
- `_clicker_tdoas()` selects the fitted estimator or legacy fallback.
- `calibrated_tdoa2dir()` applies the calibrated acoustic inversion.
- `calibrated_angle2px()` applies calibrated camera projection.
- `with_clicker_calibration()` restores configuration after scoped execution.

Calibration corrections are intended to be applied exactly once.

## 5. Older vs Calibrated Localization

| Stage | Legacy method | Calibrated method |
|---|---|---|
| TDOA | Extrema-based default | GCC (A) / SRP (B) |
| Sound speed | Nominal constant | Fitted from calibration |
| Timing offsets | No fitted correction | Per-channel correction |
| Direction | Numerical optimization | Calibrated geometric inversion |
| Array-to-camera alignment | Not fitted | 3D rotation |
| Pixels | Linear FOV mapping | Fitted pinhole projection |
| Camera distortion | Not calibrated | Supported radial parameters |
| Impossible TDOA | May still return a direction | Explicit invalid result |

The calibration is **far-field, direction-only**. It does not estimate source range or near-field position.

Invalid directions and projections may produce `NaN` coordinates. These must not be interpreted as valid video positions.

## 6. Coverage and Known Integration Gaps

| Path | Calibration coverage | Notes |
|---|---|---|
| Standard impulsive detection | Full intended path | Uses fitted GCC/SRP, direction and camera projection |
| Standard video overlay | Integrated | Consumes resulting pixel coordinates |
| `angle2px()` | Calibrated | Uses active profile |
| `tdoa2dir()` | Calibrated | Uses active profile |
| `angle2pxd()` | Legacy | Does not use calibrated camera projection |
| `tdoa2dir_ls()` | Legacy for planar D3 | Does not apply calibrated inversion |
| Alternative impulsive methods 2–4 | Partial | Original delay estimators combined with calibrated direction/camera |
| Tonal localization | Partial | Original TDOA estimation with calibrated direction/camera |
| `run_multitdoa.jl` | Partial/inconsistent | Uses `tdoa2dir_ls()` in some paths |

### Main concern: mixed estimator calibration

The acoustic calibration is method-specific.

Using an alternative TDOA estimator with timing offsets fitted using GCC or SRP can introduce systematic errors.

The standard impulsive path is the most appropriate current deployment path.

Follow-up work should explicitly identify which processing variants are fully calibrated and prevent incompatible method combinations from being silently treated as calibrated.

## 7. Calibration Constraints

Current profiles require:

- D3 equilateral hydrophone array, side length **0.18 m**.
- Physical audio channels **3, 4, 5**.
- Physical reference channel **5**, corresponding to selected reference channel 3.
- Audio sampling rate **192,000 Hz**.
- Method-matched TDOA window **`-300:700`** samples.
- Expected preprocessing: fourth-order zero-phase high-pass filtering at 1 kHz for the fitted GCC/SRP estimators.
- Camera resolution **2704 × 1520** pixels.
- The matching camera–hydrophone physical mounting.

The fitted estimator checks its expected window and sampling rate. Changes to these parameters require compatible recalibration.

The existing `c`, `default_c`, `default_window`, `angle2pxd()` and `tdoa2dir_ls()` are not automatically replaced by fitted equivalents.

Changing other detection globals by direct assignment can preserve activation, but calling another device setter deactivates the profile.

## 8. Current Fitted Profiles and Accuracy

| Property | Device A | Device B |
|---|---:|---:|
| TDOA estimator | GCC | SRP |
| Fitted sound speed | 1562.18 m/s | 1401.00 m/s |
| Horizontal FOV | 75.89° | 78.77° |
| Vertical FOV | 56.17° | 150.00° |
| Reused evaluation median | 176.74 px | 90.97 px |
| Reused evaluation RMSE | 2564.63 px | 996.61 px |
| Deployment descriptive median | 44.31 px | 101.17 px |
| Deployment descriptive RMSE | 1755.38 px | 551.19 px |
| Fit status | `fitted_weak_identifiability` | `fitted_weak_identifiability` |

These are original-resolution pixel errors from the repository's reported results.

**Interpretation:**

- Both calibrations have substantial outliers.
- The reused evaluation datasets were previously inspected; they are not untouched test sets.
- B reaches the vertical-FOV boundary, suggesting parameter instability.
- Small descriptive median errors do not imply reliable physical calibration.
- These models should be considered provisional.

Do not use the existing A/B results as independent validation of localization accuracy.

See [docs/clicker_calibration.md](../docs/clicker_calibration.md) for full methodology, experimental history and diagnostic results.

## 10. Recommended Starting Order for Another Agent

1. Read [docs/clicker_calibration.md](../docs/clicker_calibration.md), especially **Apply a fitted device** and **Current A/B results**.
2. Read [src/clicker_calibration/preset_bridge.jl](../src/clicker_calibration/preset_bridge.jl) to understand the activation and dispatch mechanism.
3. Trace `detection2angle()` → `tdoa2dir()` → `angle2px()` in [src/localization.jl](../src/localization.jl).
4. Trace `process_detections()` in [src/aspod.jl](../src/aspod.jl) and `run_contiguous_folders()` in [src/run_example.jl](../src/run_example.jl).
5. Read [test/clicker_calibration/test_presets.jl](../test/clicker_calibration/test_presets.jl) to understand the verified behavioral contracts.
6. Investigate partial calibration paths in [src/run_multitdoa.jl](../src/run_multitdoa.jl).
7. Run targeted tests and real-recording comparisons before changing algorithms.
