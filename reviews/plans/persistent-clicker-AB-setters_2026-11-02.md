# Persistent clicker A/B setters that switch themselves off

## Goal
`set_device__ophk_acoustic_D3_clicker_A()` / `_B()` then `run_contiguous_folders(...)` directly, with no
`with_clicker_calibration` block. No existing function or setter changes. A dolphin-calibrated preset is out of scope.
The user will later build one by calling the clicker A/B setter and then changing detection globals.

## Design (only `src/clicker_calibration/preset_bridge.jl` is new code)
- Activation records `ACTIVE_CLICKER_TOKEN = (rx_vect, fov_angle)`, the array objects it just set.
  `active_clicker_calibration()` returns the calibration only while both are still the same objects (`===`).
  Every existing `set_device__*`, and re-including `config.jl`, rebinds at least one of them, so the calibration
  switches off with no clear calls in the setters.
- `default_getTDOA_func = _clicker_tdoas`, a wrapper. It uses the calibrated estimator while the calibration is
  active and `CLICKER_UNCALIBRATED_TDOA` (the function before activation) otherwise, so no stale estimator is left behind.
- Activation no longer sets `c`, `default_c` or `default_window`. The calibrated `tdoa2dir`/`angle2px` read the sound
  speed and geometry from the JSON, and most setters never reset these globals, so they would be left stale.
  The fitted GCC/SRP estimator still checks its window: callers pass `window_impulsive` (-300:700), as the pipeline does.
- The setters are thin: load and validate the JSON first, then activate. A bad path changes nothing.
- `with_clicker_calibration` is unchanged and still restores everything (including the token) after the block.
- Cost per `tdoa2dir`/`angle2px`/TDOA call: two extra `===` checks.

## Also
- Export both setters from `DolphinSpeaker.jl`; update `docs/clicker_calibration.md` and add one README line.
- Tests (`test_presets.jl`): persistent A/B project as before; `set_device__ophk_acoustic_D3()` switches the calibration off;
  the stale TDOA wrapper falls back; a block inside persistent A restores A; a bad path changes nothing.
- Only `test_presets.jl` is run (niced): the machine is fully loaded.
