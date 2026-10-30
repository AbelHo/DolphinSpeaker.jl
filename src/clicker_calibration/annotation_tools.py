#!/usr/bin/env python3
"""Deterministic clicker annotation validation and review images (Pillow only).

Numerical calibration is Julia-owned. Pixel coordinates here are original image
coordinates, zero-origin x right/y down. Imports perform no work.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import sys
from datetime import datetime, timezone
from PIL import Image, ImageDraw, ImageFont

FIELDS = ("frame_id", "x", "y", "visibility", "medium", "uncertainty_px", "note", "annotator")
VISIBILITY = {"visible", "uncertain", "occluded", "absent", "excluded"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, obj: object) -> None:
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("w") as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        return list(reader)


def finite(value: str, name: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"{name} must be finite")
    return x


def validate(directory: Path, *, seal: bool = False) -> tuple[list[dict], dict]:
    frames_path, annotations_path = directory / "frames.csv", directory / "annotations.csv"
    frames, annotations = read_csv(frames_path), read_csv(annotations_path)
    ids = [f["frame_id"] for f in frames]
    if len(set(ids)) != len(ids):
        raise ValueError("frames.csv contains duplicate frame_id")
    annotation_ids = [a.get("frame_id") for a in annotations]
    if len(set(annotation_ids)) != len(annotation_ids):
        raise ValueError("annotations.csv contains duplicate frame_id")
    if set(annotation_ids) != set(ids):
        raise ValueError(f"Annotation identities differ: missing={len(set(ids)-set(annotation_ids))}, extra={len(set(annotation_ids)-set(ids))}")
    by_id = {a["frame_id"]: a for a in annotations}
    selection_path = directory / "human_review_selection.json"
    selection = json.loads(selection_path.read_text()) if selection_path.exists() else None
    accepted = None
    if selection is not None:
        if selection.get("schema_version", selection.get("schema")) != "clicker-human-review-selection-v1":
            raise ValueError("Unsupported human review selection schema")
        if selection.get("frames_sha256") != sha(frames_path):
            raise ValueError("Human review selection does not match frames.csv")
        accepted = set(selection["accepted_frame_ids"])
        if not accepted or not accepted <= set(ids) or len(accepted) != len(selection["accepted_frame_ids"]):
            raise ValueError("Human review selection has invalid frame identities")
    joined, counts = [], {v: 0 for v in sorted(VISIBILITY)}
    for f in frames:
        a = by_id[f["frame_id"]]
        if set(FIELDS) - set(a):
            raise ValueError("Missing required annotation columns")
        status = a["visibility"]
        if status not in VISIBILITY:
            raise ValueError(f"Invalid visibility {status!r}")
        if a["medium"] not in {"underwater", "air", "mixed", "unknown"}:
            raise ValueError(f"Invalid medium: {a['medium']!r}")
        if not a["annotator"].strip():
            raise ValueError("Annotator identity is required")
        if accepted is not None:
            if f["frame_id"] not in accepted:
                if status != "excluded" or a["x"].strip() or a["y"].strip() or a["uncertainty_px"].strip():
                    raise ValueError("Frames outside the human-selected interval must remain excluded without coordinates")
            elif selection.get("coordinate_reannotation_required", False):
                if "gpt-6.1-terra" not in a["annotator"] or status not in {"visible", "uncertain", "occluded"}:
                    raise ValueError("Selected frame still requires Terra coordinate reinspection")
                if a["medium"] != "underwater":
                    raise ValueError("Human-confirmed underwater interval has a conflicting medium label")
        width, height = int(f["width"]), int(f["height"])
        if width <= 0 or height <= 0:
            raise ValueError("Invalid frame dimensions")
        x = y = uncertainty = None
        if status == "visible" or a["x"].strip() or a["y"].strip():
            x, y = finite(a["x"], "x"), finite(a["y"], "y")
            if not (0 <= x < width and 0 <= y < height):
                raise ValueError(f"Coordinate outside original image: {f['frame_id']}")
        if status == "visible" or a["uncertainty_px"].strip():
            uncertainty = finite(a["uncertainty_px"], "uncertainty_px")
            if uncertainty <= 0:
                raise ValueError("uncertainty_px must be positive")
        path = Path(f["image_path"])
        if not path.is_absolute():
            path = directory / path
        if not path.is_file():
            raise ValueError(f"Missing image: {path}")
        with Image.open(path) as im:
            if im.size != (width, height):
                raise ValueError(f"Recorded dimensions differ from image: {path}")
        joined.append({**f, **a, "x_value": x, "y_value": y,
                       "uncertainty_value": uncertainty, "resolved_image_path": path})
        counts[status] += 1
    manifest = {"schema_version": "clicker-annotations-v1",
                "frames_sha256": sha(frames_path),
                "annotations_sha256": sha(annotations_path),
                "n_frames": len(frames), "counts": counts,
                "annotator": sorted({a["annotator"] for a in annotations}),
                "landmark_convention": "visible clicker housing center; uncertainty covers possible emitter position within housing",
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "coordinate_convention": "zero-origin original pixels; x right, y down"}
    if selection is not None:
        manifest["human_review_selection_sha256"] = sha(selection_path)
        manifest["human_selected_frame_count"] = len(accepted)
    if seal:
        atomic_json(directory / "annotations.complete.json", manifest)
    else:
        existing = json.loads((directory / "annotations.complete.json").read_text())
        for key in ("frames_sha256", "annotations_sha256", "n_frames", "human_review_selection_sha256"):
            if existing.get(key) != manifest.get(key):
                raise ValueError(f"Stale annotation completion manifest: {key}")
    return joined, manifest


def font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size)
    except OSError:
        return ImageFont.load_default(size=size)


def cross(draw: ImageDraw.ImageDraw, xy: tuple[float, float], color: str, r: int = 15) -> None:
    x, y = xy
    draw.line((x-r, y, x+r, y), fill=color, width=3)
    draw.line((x, y-r, x, y+r), fill=color, width=3)
    draw.ellipse((x-r, y-r, x+r, y+r), outline=color, width=3)


def gallery(directory: Path, *, fitted: bool = False, model_dir: Path | None = None) -> dict:
    rows, manifest = validate(directory)
    model_dir = directory if model_dir is None else model_dir
    if not fitted and model_dir != directory:
        raise ValueError("A model directory is only valid for a fitted gallery")
    revision = manifest["annotations_sha256"][:12]
    if fitted:
        revision += "_" + sha(model_dir / "calibration.json")[:12]
    name = ("refined_review" if model_dir != directory else "fit_review") if fitted else "annotation_review"
    out = directory / name / revision
    out.mkdir(parents=True, exist_ok=True)
    predictions: dict[str, list[dict]] = {}
    if fitted:
        for p in read_csv(model_dir / "predictions.csv"):
            predictions.setdefault(p["frame_id"], []).append(p)
    cards, thumbs = [], []
    for i, row in enumerate(rows):
        im = Image.open(row["resolved_image_path"]).convert("RGB")
        draw = ImageDraw.Draw(im)
        point = None if row["x_value"] is None else (row["x_value"], row["y_value"])
        if point is not None:
            color = "#00efb0" if row["visibility"] == "visible" else "#ffcb52"
            cross(draw, point, color)
            u = row["uncertainty_value"]
            if u:
                x, y = point
                draw.ellipse((x-u, y-u, x+u, y+u), outline=color, width=2)
        for p in predictions.get(row["frame_id"], []):
            for keys, color in [(("predicted_x", "predicted_y"), "#ff557d")]:
                try:
                    xy = tuple(float(p[k]) for k in keys)
                except (KeyError, TypeError, ValueError):
                    continue
                if all(math.isfinite(v) for v in xy):
                    if 0 <= xy[0] < im.width and 0 <= xy[1] < im.height:
                        cross(draw, xy, color, 11)
        # Preserve original pixel geometry and never obscure a clicker near an
        # image edge with a title banner; identity lives in the gallery caption.
        filename = f"{i:05d}.jpg"
        im.save(out / filename, quality=93)
        im.thumbnail((640, 380))
        thumbnail = f"{i:05d}_thumb.jpg"
        im.save(out / thumbnail, quality=85)
        thumbs.append((im.copy(), row["frame_id"]))
        raw_link = os.path.relpath(row["resolved_image_path"], out)
        position = "" if point is None else f"({point[0]:.1f}, {point[1]:.1f}), ±{row['uncertainty_value']} px"
        fit_notes = []
        if fitted:
            for prediction in predictions.get(row["frame_id"], []):
                status = prediction.get("final_status", prediction.get("predicted_status", "unknown"))
                status = {"impossible_tdoa": "no physically feasible direction",
                          "outside_image": "outside image", "behind_camera": "behind camera"}.get(status, status)
                error = prediction.get("final_error_px", prediction.get("predicted_error_px", ""))
                try:
                    numeric_error = float(error)
                    error_text = f"{numeric_error:.1f} px" if math.isfinite(numeric_error) else "unavailable"
                except (ValueError, TypeError):
                    error_text = "unavailable"
                fit_notes.append(f'Final fit: {status}; error {error_text}; {prediction.get("split", "")}; {prediction.get("eligibility_status", "")}')
            if not predictions.get(row["frame_id"]):
                fit_notes.append("Not part of this fitted event population")
        fit_caption = ''.join(f'<p class="legend">{html.escape(note)}</p>' for note in fit_notes)
        cards.append(f'<article data-status="{row["visibility"]}" data-medium="{row["medium"]}"><a href="{filename}" target="_blank"><img loading="lazy" src="{thumbnail}" alt="{html.escape(row["frame_id"])}"></a><div><b>{html.escape(row["frame_id"])}</b><p>{row["visibility"]} / {row["medium"]} {position}</p><p>{html.escape(row["note"])}</p>{fit_caption}<a href="{html.escape(raw_link, quote=True)}" target="_blank">Original image</a></div></article>')
    counts = " · ".join(f"{k}: {v}" for k, v in manifest["counts"].items())
    selected_button = '<button onclick="filter(\'selected\')">User-selected frames</button>' if "human_selected_frame_count" in manifest else ""
    initial_filter = "filter('selected');" if selected_button else ""
    title = f"Device {directory.name}: {'fit' if fitted else 'annotation'} review"
    body = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{html.escape(title)}</title><style>body{{font:16px system-ui;margin:28px;background:#101822;color:#e8eff5}}header{{max-width:1100px}}a{{color:#74d8ff}}main{{display:grid;grid-template-columns:repeat(auto-fill,minmax(350px,1fr));gap:18px}}article{{background:#1d2a38;border-radius:10px;overflow:hidden}}article img{{width:100%;display:block}}article div{{padding:14px}}button{{padding:8px 12px;margin:5px;border:0;border-radius:5px}}p{{line-height:1.5}}.legend{{color:#a9bccb}}</style></head><body><header><h1>{html.escape(title)}</h1><p>{counts}</p><p class="legend">Green: visible clicker annotation. Amber: uncertain annotation. Red: calibrated acoustic estimate. Baseline comparisons remain in the reports. Coordinates refer to the original image, with x right and y down.</p><p>Click an image for full resolution. Review version: {revision}. <a href="../../annotations.csv">Annotations CSV</a>. Send corrections by frame ID; refitting will preserve this version.</p><p>{selected_button}<button onclick="filter('all')">All</button><button onclick="filter('underwater-visible')">Visible underwater</button><button onclick="filter('visible')">All visible</button><button onclick="filter('uncertain')">Uncertain</button><button onclick="filter('absent')">Absent</button><button onclick="filter('occluded')">Occluded</button></p></header><main>{''.join(cards)}</main><script>function filter(v){{document.querySelectorAll('article').forEach(e=>e.style.display=(v==='all'||(v==='selected'?e.dataset.status!=='excluded':v==='underwater-visible'?e.dataset.status==='visible'&&e.dataset.medium==='underwater':e.dataset.status===v))?'':'none')}}{initial_filter}</script></body></html>'''
    (out / "index.html").write_text(body)
    if thumbs:
        selected = [thumbs[round(i*(len(thumbs)-1)/min(11,len(thumbs)-1))] for i in range(min(12,len(thumbs)))] if len(thumbs)>1 else thumbs
        sheet = Image.new("RGB", (1280, math.ceil(len(selected)/3)*270), "#101822")
        for j, (im, label) in enumerate(selected):
            im.thumbnail((420, 235))
            x, y = (j%3)*426, (j//3)*270
            sheet.paste(im,(x,y))
            ImageDraw.Draw(sheet).text((x+5,y+238),label,fill="white",font=font(15))
        sheet.save(out / "overview.jpg",quality=90)
    result = {"gallery": str((out / "index.html").resolve()), "revision": revision,
              "overview": str((out / "overview.jpg").resolve()), **manifest}
    if fitted:
        result["model_directory"] = str(model_dir.resolve())
        result["model_sha256"] = sha(model_dir / "calibration.json")
    atomic_json(directory / f"{name}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("seal", "validate", "gallery", "fit-gallery"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--model-dir", type=Path)
    args = parser.parse_args()
    if args.action in ("seal", "validate"):
        _, result = validate(args.directory, seal=args.action == "seal")
    else:
        result = gallery(args.directory, fitted=args.action == "fit-gallery", model_dir=args.model_dir)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
