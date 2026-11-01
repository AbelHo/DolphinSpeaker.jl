# Make the clicker-calibration commits less intrusive to existing code (draft v0)

> Earlier draft: kept persistent A/B setters and the per-setter `clear_clicker_calibration!()` calls. It was superseded by the scoped `with_clicker_calibration` design in the approved plan.

## Context
The 29 unpushed `dev2` commits add a clicker calibration subsystem. They also change existing code beyond what the pipeline needs:
- `detection2angle` special-cases `getTDOA_func === _clicker_method_tdoas` to pass `fs`.
- Calibration hooks were added to four localization functions. Two of them, `angle2pxd` and `tdoa2dir_ls`, are not on the pipeline path and quietly override explicit `fov_angle` or geometry arguments.
- `set_device__ophk_acoustic_D3` clears the calibration twice.
- Lazy module loading is written out twice.
- The package exports only `run_clicker_calibration`.

Choices at this stage: hooks only in `tdoa2dir` and `angle2px`; keep the setter clears but drop the duplicate; one entry point with `using DolphinSpeaker` as the main path; `exercise_error_logging!` gets only a comment; new commits on top.

## Commit 1: restrict the legacy localization hooks to `tdoa2dir` and `angle2px`
- `localization.jl`: restore `angle2pxd`, `tdoa2dir_ls` and `detection2angle`. `angle2px` gets the original body back, plus one hook line.
- `preset_bridge.jl`: add `_clicker_sample_rate`. `_clicker_method_tdoas` drops its `fs` kwarg. `calibrated_tdoa2dir` checks the sample rate.
- `test_presets.jl`: `angle2pxd` and `tdoa2dir_ls` stay uncalibrated; the fs=96000 error test is kept.

## Commit 2: one lazy loader, package entry points, duplicate clear
- Add `_ensure_clicker_module()` in `preset_bridge.jl`.
- `DolphinSpeaker.jl`: wrappers for `run_clicker_calibration` and `run_full_clicker_calibration`.
- Remove the redundant clear in `set_device__ophk_acoustic_D3`; the other 10 setter clears stay.
- Add SHA/Random to `Project.toml` if the package load check fails.
- Docs: "Run and resume" shows `using DolphinSpeaker` first.

## Commit 3: comment on `exercise_error_logging!`

## Verification
- `git diff` against the pre-calibration code: only the hook lines in `localization.jl`; in `config.jl`, the include, the A/B setters and 10 clear lines.
- The clicker Julia tests, the package load check, and the CLI usage check.
