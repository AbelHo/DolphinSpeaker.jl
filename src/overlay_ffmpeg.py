try:
    import pandas as pd
except Exception:
    pd = None

import csv
import subprocess
import shlex
from typing import Optional


def overlay_with_boxes(csv_path: str,
                       input_video: str,
                       output_video: str,
                       radius: int = 100,
                       dry_run: bool = False,
                       verbose: bool = False,
                       use_gpu: bool = False) -> Optional[str]:
    """Read detection CSV and overlay filled red boxes on matching frames using ffmpeg.

    Parameters
    - csv_path: path to CSV containing columns 'frame', 'px', 'py'
    - input_video: path to input video file
    - output_video: desired output path
    - radius: diameter of the square/box to draw (default 100)
    - dry_run: if True, don't execute ffmpeg; return the constructed command
    - verbose: if True, print progress and the filter string

    Returns the ffmpeg command string when dry_run=True, otherwise returns None.
    """
    filters = []

    if pd is not None:
        df = pd.read_csv(csv_path)
        # detect column names
        cols = set(df.columns.astype(str))
        px_name = 'px' if 'px' in cols else ('x' if 'x' in cols else None)
        py_name = 'py' if 'py' in cols else ('y' if 'y' in cols else None)
        if 'frame' not in cols:
            raise KeyError("CSV must contain a 'frame' column")
        if px_name is None or py_name is None:
            raise KeyError("CSV must contain 'px'/'py' or 'x'/'y' columns")
        rows = (r for _, r in df.iterrows())
        col_names = (px_name, py_name)
    else:
        # fallback to csv.DictReader if pandas isn't installed
        def _gen_rows():
            with open(csv_path, newline='') as fh:
                reader = csv.DictReader(fh)
                cols = set(reader.fieldnames or [])
                px_name = 'px' if 'px' in cols else ('x' if 'x' in cols else None)
                py_name = 'py' if 'py' in cols else ('y' if 'y' in cols else None)
                if 'frame' not in cols:
                    raise KeyError("CSV must contain a 'frame' column")
                if px_name is None or py_name is None:
                    raise KeyError("CSV must contain 'px'/'py' or 'x'/'y' columns")
                for r in reader:
                    yield r, px_name, py_name

        # when using csv fallback we yield tuples (row, px_name, py_name)
        rows = _gen_rows()
        col_names = None

    for item in rows:
        # row may be a pandas Series or a tuple from csv fallback
        if pd is not None:
            row = item
            px_name, py_name = col_names
        else:
            row, px_name, py_name = item

        frame = int(round(float(row['frame'])))
        x = int(round(float(row[px_name])))
        y = int(round(float(row[py_name])))
        # draw a filled box centered at (x, y)
        filters.append(
            f"drawbox=x={x-radius//2}:y={y-radius//2}:w={radius}:h={radius}:color=red@0.5:t=fill:enable='eq(n,{frame})'"
        )

    filters_str = ",".join(filters)

    if verbose:
        print("Constructed filter string:")
        print(filters_str)

    # Safely build the ffmpeg command; wrap the filter string in double quotes for shell use
    # Use subprocess (no shell) by passing args list; ffmpeg expects the filter as a single argument
    # Build base command
    cmd = ["ffmpeg"]

    # If GPU rendering is requested try to auto-detect a usable encoder and add args
    if use_gpu:
        # try to detect common ffmpeg GPU encoders/devices
        encoder, extra_args = _detect_ffmpeg_gpu(verbose=verbose)
        if encoder is None:
            if verbose:
                print("No supported GPU encoder detected; falling back to CPU rendering")
        else:
            # For some encoders (e.g., vaapi) ffmpeg expects a different input setup
            # We'll append any extra args before the input, and specify video codec later
            if extra_args:
                cmd.extend(extra_args)

    cmd.extend(["-i", input_video, "-vf", filters_str])

    # If we detected an encoder, set codec for video output accordingly
    if use_gpu and 'encoder' in locals() and encoder is not None:
        # encoder is the codec name like 'h264_nvenc', which goes with -c:v
        cmd.extend(["-c:v", encoder])
        # Copy audio
        cmd.extend(["-c:a", "copy"])
    else:
        cmd.extend(["-codec:a", "copy"])

    cmd.append(output_video)

    # Return the command string for dry-run
    if dry_run:
        return " ".join(shlex.quote(p) for p in cmd)

    # Execute ffmpeg and raise on error
    proc = subprocess.run(cmd, check=False, capture_output=not verbose, text=True)
    if proc.returncode != 0:
        # when not verbose, include stderr in the exception message
        stderr = proc.stderr if not verbose else ""
        raise RuntimeError(f"ffmpeg failed with exit code {proc.returncode}: {stderr}")

    return None


def _detect_ffmpeg_gpu(verbose: bool = False):
    """Detect a suitable ffmpeg GPU encoder.

    Returns a tuple (encoder_name_or_None, extra_args_list).
    extra_args_list can contain args that need to appear before -i (for VAAPI for example).
    """
    # Prefer order: nvenc (nvidia), qsv (intel), vaapi (intel/amd), h264_amf (amd)
    candidates = [
        ("h264_nvenc", []),
        ("hevc_nvenc", []),
        # macOS hardware encoders: VideoToolbox and newer MPS encoders
        ("h264_videotoolbox", []),
        ("hevc_videotoolbox", []),
        ("mps_h264", []),
        ("mps_hevc", []),
        ("h264_qsv", []),
        ("hevc_qsv", []),
        ("h264_vaapi", ["-vaapi_device", "/dev/dri/renderD128"]),
        ("hevc_vaapi", ["-vaapi_device", "/dev/dri/renderD128"]),
        ("h264_amf", []),
        ("hevc_amf", []),
    ]

    try:
        # ask ffmpeg for its encoders list
        p = subprocess.run(["ffmpeg", "-encoders"], capture_output=True, text=True)
        out = p.stdout + (p.stderr or "")
    except Exception as e:
        if verbose:
            print(f"Failed to run ffmpeg -encoders: {e}")
        return None, []

    for enc, extra in candidates:
        # ffmpeg lists encoders like "V..... h264_nvenc"
        if enc in out:
            if verbose:
                print(f"Detected encoder: {enc}")
            return enc, extra

    if verbose:
        print("No known GPU encoders found in ffmpeg -encoders output")
    return None, []


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Overlay detection boxes on a video using ffmpeg.")
    parser.add_argument("csv", help="CSV file containing detection pixels (columns: frame, px, py)")
    parser.add_argument("input", help="Input video file path")
    parser.add_argument("output", help="Output video file path")
    parser.add_argument("--radius", type=int, default=100, help="Box diameter in pixels (default: 100)")
    parser.add_argument("--dry-run", action="store_true", help="Print the ffmpeg command and exit without running it")
    parser.add_argument("--verbose", action="store_true", help="Print additional debug output")
    parser.add_argument("--use-gpu", action="store_true", help="Attempt to use a GPU encoder if available")

    args = parser.parse_args()

    cmdstr = overlay_with_boxes(args.csv, args.input, args.output, radius=args.radius, dry_run=args.dry_run, verbose=args.verbose, use_gpu=args.use_gpu)
    if args.dry_run:
        print(cmdstr)