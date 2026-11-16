# Python environment (uv)

The repository's Python helpers run from a repo-local uv environment in `.venv/`. These include the calibration review tools, the video export, the checkerboard intrinsics tool and the GoPro metadata reader. The environment is defined by `pyproject.toml` and pinned by `uv.lock`. This replaces the miniconda base interpreter, where `scipy.optimize` failed to import: stale scipy 1.16.2 extension files sat next to scipy 1.18.0, and three OpenCV distributions were installed at once.

```sh
uv sync                      # create or update .venv from uv.lock (add --group dev for pytest)
source .venv/bin/activate    # then `python3` is the repo environment
```

The Julia code calls Python through `CLICKER_PYTHON` (annotation, review, GoPro lens tools) and `CLICKER_PLOT_PYTHON` (diagnostic plots); both default to `python3`. Either activate `.venv` before starting Julia, or set:

```sh
export CLICKER_PYTHON="$PWD/.venv/bin/python"
export CLICKER_PLOT_PYTHON="$PWD/.venv/bin/python"
```

Tests: `uv run --group dev pytest test/clicker_calibration tools/test_gpmf_lens_info.py`.

Optional extras:
- `uv sync --extra synth` adds `scaper`, `pyroomacoustics` and `soundfile` for `synth_bioacoustic.py`. It needs the system library `libsox-dev` first (`sudo apt install libsox-dev`).
- cuML (`test_cuml_direct.py`) is not included; it needs the RAPIDS conda channel.

The Julia `CondaPkg.toml` environment (`.CondaPkg/`, used by PythonCall) is separate and unchanged.
