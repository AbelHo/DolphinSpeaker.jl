# Clicker calibration: quick start

For **underwater D3 recordings**: 0.18 m array, WAV channels **3–5** (reference: 5).
Requires the repo’s Julia dependencies, Python with NumPy/Pillow/Matplotlib,
`/usr/bin/ffmpeg` and `ffprobe`, and configured Herdr/Telegram.

**1. Keep one continuous session per device, including all WAV parts and MP4 chapters.**

```text
new_dataset/
├── A/   # WAV + MP4 files for device A
└── B/   # WAV + MP4 files for device B
```

**2. From the repository root, run in a Herdr pane.**
Initialize the pane once if needed: `hagent init clickcal-new-01` (choose a unique name).
Edit the two paths; use a fresh run name.

```bash
CLICKER_DATA="/path/to/new_dataset"
CLICKER_RUN="calibration_runs/new_dataset_v1"

python tools/run_clicker_job.py \
  --job-dir "$CLICKER_RUN/jobs/01" --timeout-s 86400 --wake -- \
  julia --project=. --startup-file=no --compiled-modules=no --pkgimages=no --threads=auto -O2 \
  tools/run_clicker_calibration.jl "$CLICKER_DATA" "$CLICKER_RUN/work" \
  --results="$CLICKER_DATA/results" \
  --observation-policy=terra_estimated_underwater_coordinates
```

Terra annotates; the pipeline handles everything else. The policy flag includes
uncertain points with their uncertainty weights; omit it for visible-only fitting.

**3. Check the Telegram links, then open `new_dataset/results/A/index.html` and `B/index.html`.**
Find synchronized videos, frames and fitted parameters there. Check fit errors before use.

**Resume:** rerun with the same paths and a new job folder (`jobs/02`).
**One device:** add `--device=A` or `--device=B`.
**Apply:** use the new `results/A/calibration.json` (or B) as
[`calibration_path`](clicker_calibration.md#apply-a-fitted-device).

[Full reference](clicker_calibration.md)
