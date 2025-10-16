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
                       verbose: bool = False) -> Optional[str]:
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
            f"drawbox=x={x-radius//2}:y={y-radius//2}:w={radius}:h={radius}:color=red@0.7:t=fill:enable='eq(n,{frame})'"
        )

    filters_str = ",".join(filters)

    if verbose:
        print("Constructed filter string:")
        print(filters_str)

    # Safely build the ffmpeg command; wrap the filter string in double quotes for shell use
    # Use subprocess (no shell) by passing args list; ffmpeg expects the filter as a single argument
    cmd = [
        "ffmpeg",
        "-i",
        input_video,
        "-vf",
        filters_str,
        "-codec:a",
        "copy",
        output_video,
    ]

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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Overlay detection boxes on a video using ffmpeg.")
    parser.add_argument("csv", help="CSV file containing detection pixels (columns: frame, px, py)")
    parser.add_argument("input", help="Input video file path")
    parser.add_argument("output", help="Output video file path")
    parser.add_argument("--radius", type=int, default=100, help="Box diameter in pixels (default: 100)")
    parser.add_argument("--dry-run", action="store_true", help="Print the ffmpeg command and exit without running it")
    parser.add_argument("--verbose", action="store_true", help="Print additional debug output")

    args = parser.parse_args()

    cmdstr = overlay_with_boxes(args.csv, args.input, args.output, radius=args.radius, dry_run=args.dry_run, verbose=args.verbose)
    if args.dry_run:
        print(cmdstr)