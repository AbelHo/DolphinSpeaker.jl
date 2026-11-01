# Make the clicker-calibration commits less intrusive to existing code

## Context
The 29 clicker-calibration commits on `dev2` (`d22056c..d6e9ab3`, now pushed; BASE = `d22056c`, the pre-calibration commit) add a clicker calibration subsystem. They also change existing code more than needed:
- `clear_clicker_calibration!()` is added to all 11 `set_device__*` setters.
- `detection2angle` special-cases `getTDOA_func === _clicker_method_tdoas`.
- Calibration hooks were added to four localization functions. Two of them, `angle2pxd` and `tdoa2dir_ls`, are not on the pipeline path and quietly override explicit `fov_angle` or geometry arguments.
- Lazy module loading is written out twice.
- The package exports only `run_clicker_calibration`.

Your choices:
- Turn calibration on with a **scoped block** instead of persistent setters, so no setter needs to change.
- Keep hooks only in `tdoa2dir` and `angle2px`.
- One loader, with `using DolphinSpeaker` as the main documented path.
- Leave `exercise_error_logging!` in place and only add a comment.
- New commits on a new branch `refactor/click-calibration` off `origin/dev2` (6276235; dev2 is now pushed).

Out of scope: moving the agent/Telegram code, the Python helper module, splitting the docs, test runners.

End state vs BASE:
- `config.jl`: +1 line (the `preset_bridge.jl` include).
- `localization.jl`: the hook in `tdoa2dir` and `angle2px` only.
- `DolphinSpeaker.jl`: exports and two thin wrappers.

## Step 0: before any code edit
Copy this plan unchanged to `reviews/plans/ask-me-relevant-questions-gleaming-teapot.md`, and save the earlier draft (persistent setters + clears) as `..._v0.md`.

## Commit 1: restrict the legacy localization hooks to `tdoa2dir` and `angle2px`
`src/localization.jl`: restore everything to BASE, then re-add only these hooks:
- `tdoa2dir`: keep the `calibration=active_clicker_calibration()` kwarg and the one-line early return.
- `angle2px`: add the kwarg and `calibration === nothing || return calibrated_angle2px(ang,calibration,imsize)` above the **original** `@debug` and `return angle2pxd(...)` lines.
- `angle2pxd`, `tdoa2dir_ls`, `detection2angle`: identical to BASE.

`src/clicker_calibration/preset_bridge.jl`:
- Add `_clicker_sample_rate(cal) = Float64(get(get(cal,"estimator_spec",Dict()),"sample_rate_hz",192000.0))`. The D3/Zoom F6 runs at 192 kHz, the same default as `_mr_estimator_spec` and `prep_gcc`.
- `_clicker_method_tdoas(data,peaks;window,ref_channel=3)`: drop the `fs` kwarg and use `_clicker_sample_rate(calibration)`. It then has the same signature as the other TDOA functions.
- `calibrated_tdoa2dir`: add `Float64(fs)==_clicker_sample_rate(calibration) || error("sample rate differs from active clicker calibration")`, so the fs-mismatch guard now lives here.

`test/clicker_calibration/test_presets.jl`:
- `angle2pxd` and `tdoa2dir_ls` must now stay **uncalibrated** while a calibration is active.
- The `detection2angle(...; fs=96000)` → `ErrorException` test is kept.

Edge case, to note in a comment: `tdoa2dir_ls` on a non-planar array falls back to `tdoa2dir`. While a calibration is active, that hits the D3 geometry check (ArgumentError).

## Commit 2: scoped activation replaces the persistent A/B setters and all setter clears
`src/clicker_calibration/preset_bridge.jl`, which is new code:
- `with_clicker_calibration(f, device::AbstractString; calibration_path=_clicker_default_path(device))`:
  1. `_ensure_clicker_module()` (see commit 3; define it here), then `load_calibration(path; device)`. Validation happens before any global changes.
  2. `before = _clicker_global_snapshot()`: a `Dict` of every **non-const** global in `@__MODULE__`, from `names(mod; all=true)` filtered by `isdefined && !isconst`. It stores references, not copies. Function *definitions* are const and are skipped. Function-valued variables such as `default_getTDOA_func` are included.
  3. `_activate_clicker_calibration!(device, calibration)`: the current body of `_set_clicker_calibration` without the restore record. It calls `set_device__ophk_acoustic_D3_clicker()`, then sets `ACTIVE_CLICKER_CALIBRATION`, `c`, `default_c`, `ref_channel`, `default_window`, `imsize`, `fov_angle`, `default_getTDOA_func` and `FILE_device_ID`.
  4. `changed` = the names in `before` whose binding is now `!==` the old one. That set includes `ACTIVE_CLICKER_CALIBRATION` and the D3-clicker globals.
  5. `try return f() finally` rebind each changed name to its old value with `setglobal!`.
- Effects:
  - Nesting (an A block inside a B block) and exceptions both restore exactly.
  - Only globals the **activation** changed are restored. Globals that pipeline code writes inside `f` are left alone.
- `_clicker_default_path(device)` maps `"A"`/`"B"` to `calibrations/clicker_20260209_{A,B}.json`. It's a function, not a `const`, because `config.jl` and this file are re-included about 10 times.
- Delete `clear_clicker_calibration!`, `CLICKER_CALIBRATION_RESTORE` and `_set_clicker_calibration`. `ACTIVE_CLICKER_CALIBRATION` (default `nothing`) and `active_clicker_calibration()` stay.

`src/config.jl`: reset to BASE plus the single `include("clicker_calibration/preset_bridge.jl")` line. This removes all 11 clears and the `set_device__ophk_acoustic_D3_clicker_A/B` setters.

`src/DolphinSpeaker.jl`: replace the A/B setter export with `export with_clicker_calibration`.

`test/clicker_calibration/test_presets.jl`: rewrite to use blocks. Check:
- A and B both activate and project as before.
- After a block: `active_clicker_calibration()===nothing`, and `c`, `default_c`, `default_getTDOA_func`, `default_window`, `imsize`, `fov_angle`, `FILE_device_ID` and `impulsive_band_pass` are back to their pre-block bindings.
- A nested B-inside-A block restores A afterwards.
- An exception thrown inside a block still restores.
- A bad `calibration_path` changes nothing.
- The method-matched `detection2angle` checks are kept.

`docs/clicker_calibration.md`, "Apply a fitted device" (~line 139): change the examples to `with_clicker_calibration("A") do … end`. Add the usage rules:
- Run calibrated work inside the block.
- Don't leave tasks running after the block ends.
- Don't call `set_device__*` inside it; the D3 geometry check would throw.

## Commit 3: one lazy loader and the package entry points
- `src/DolphinSpeaker.jl`: two thin wrappers, `run_clicker_calibration` and `run_full_clicker_calibration`. Each is `(_ensure_clicker_module(); _clicker_call(:name, args...; kwargs...))`. Export both. The committed inline `include` body goes away.
- Likely real bug: `ClickerCalibration.jl` has `using SHA, Random`, but neither is in `Project.toml` `[deps]`, and including it into the package module needs both. Run the package check below first. If it fails, `Pkg.add(["SHA","Random"])` and commit `Project.toml`.
- `src/clicker_calibration.jl` and `tools/run_clicker_calibration.jl` stay as the lightweight path the CLI uses.
- `docs/clicker_calibration.md` "Run and resume": show `using DolphinSpeaker` + `run_full_clicker_calibration(...)` first; give the standalone include as the CLI's lightweight alternative. Add a one-line pointer in `docs/clicker_calibration_quickstart.md`.

## Commit 4: comment on `exercise_error_logging!`
`src/clicker_calibration/support.jl`: docstring only. It says the function is a manual pre-launch check with no callers in the workflow, that it **rewrites** the log after archiving the marker, and that it should only be run before any workers start.

## Verification
- `git diff d22056c -- src/config.jl`: exactly one added include line. `-- src/localization.jl`: only the `tdoa2dir` and `angle2px` hook lines.
- Tests, each run as `julia --project=. test/clicker_calibration/<file>`: `test_presets.jl`, `test_workflow.jl`, `test_full_workflow.jl`, `test_event_projection.jl`, `test_method_refinement.jl`.
- Package path:
  `julia --project=. -e 'using DolphinSpeaker; c0=DolphinSpeaker.c; with_clicker_calibration("A") do; @assert DolphinSpeaker.active_clicker_calibration()["device"]=="A"; end; @assert DolphinSpeaker.active_clicker_calibration()===nothing && DolphinSpeaker.c===c0; @assert isdefined(Main,:run_full_clicker_calibration)'`
- CLI still loads: `julia --project=. tools/run_clicker_calibration.jl` prints the usage error.
- `julia --version` is at least 1.9, which `getglobal`/`setglobal!` need.
- `tools/replay_clicker_legacy.jl` is not affected: it uses a frozen copy of the old code.

## Execution (herdr, Sonnet implementer, me supervising)
1. Me: `git switch -c refactor/click-calibration origin/dev2`. The local `dev2` is equal to `origin/dev2`; the only working-tree change is the untracked `docs/error.log`, which has since been committed upstream. Copy the plan into `reviews/plans/` (Step 0) as the first commit on the branch.
2. Me: `hagent init clickcal-refactor-lead`, then
   `hagent spawn claude clickcal-refactor-impl --model sonnet --effort high -- "<task>"`.
   The task gives:
   - the plan path;
   - the branch (no switching or pushing);
   - commits 1–4 in order, each one compiling and passing its tests before the next, with the session attribution trailer;
   - the "leave existing functions untouched; add new functions" rule;
   - the Julia command: `julia --project=. --startup-file=no`. If precompile segfaults, as in `docs/error.log`, add `--compiled-modules=no --pkgimages=no`;
   - failures go to `docs/error.log` per the long-task-error-log skill;
   - send blockers and decision questions to the parent, and finish with a `DONE:` summary (commits, test output, deviations).
3. Me: supervise with `hagent wait`/`hagent read` and answer questions. Out-of-plan changes need my approval.
4. Me, final review:
   - read `git diff origin/dev2..refactor/click-calibration`;
   - check `git diff d22056c -- src/config.jl src/localization.jl` is minimal;
   - re-run the verification commands myself;
   - run `/code-review` at medium on the branch, fix or send back findings, then report to you.
   
   No push, no merge, no PR unless you ask.
