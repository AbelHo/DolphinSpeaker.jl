"""Authenticate optional all-event inference without changing frozen fit products."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


def _sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def resolve_predictions(original_dir, model_dir):
    """Return CSV path and provenance; partial or stale inference fails closed.

    Without sidecars, retain the frozen cohort table. An all-event sidecar must
    bind the current original truth, current frozen model/cohort/report, every
    recorded input, and its producing source. It never replaces fit metrics.
    """
    original = Path(original_dir).resolve()
    model = Path(model_dir).resolve()
    manifest_path = model / "inference_events.json"
    csv_path = model / "inference_events.csv"
    if not manifest_path.exists() and not csv_path.exists():
        return model / "predictions.csv", {"scope": "frozen_fit_cohort", "tracking_paths": []}
    if not manifest_path.is_file() or not csv_path.is_file():
        raise ValueError("Incomplete all-event inference sidecars")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "clicker-all-event-projection-v1" or manifest.get("status") != "complete":
        raise ValueError("Unsupported or incomplete all-event inference manifest")
    if (manifest.get("original_dir") != str(original) or manifest.get("model_dir") != str(model)
            or manifest.get("output_csv") != str(csv_path)):
        raise ValueError("All-event inference directory identity differs")
    if manifest.get("output_csv_sha256") != _sha(csv_path):
        raise ValueError("All-event inference CSV digest differs")
    report = model / ("report.json" if (model / "report.json").exists() else "fit_report.json")
    required = {
        "events_sha256": original / "events.csv",
        "frames_sha256": original / "frames.csv",
        "annotations_sha256": original / "annotations.csv",
        "annotations_complete_sha256": original / "annotations.complete.json",
        "prepare_manifest_sha256": original / "prepare_manifest.json",
        "calibration_sha256": model / "calibration.json",
        "cohort_predictions_sha256": model / "predictions.csv",
        "model_report_sha256": report,
    }
    inputs = manifest.get("inputs", {})
    for key, path in required.items():
        digest = _sha(path)
        if manifest.get("hashes", {}).get(key) != digest or inputs.get(str(path)) != digest:
            raise ValueError(f"All-event inference binding differs: {key}")
    for path, digest in inputs.items():
        if _sha(path) != digest:
            raise ValueError(f"All-event inference input changed: {path}")
    report_data = json.loads(report.read_text())
    for name, alternate in (("calibration.json", "calibration_sha256"), ("predictions.csv", "predictions_sha256")):
        expected = report_data.get("outputs_sha256", {}).get(name, report_data.get(alternate))
        if expected != _sha(model / name):
            raise ValueError(f"Frozen model report binding differs: {name}")
    expected_sources = {"event_projection.jl", "physical_fit.jl", "prepare.jl", "srp_review.jl",
                        "method_refinement.jl", "support.jl", "workflow.jl"}
    if set(manifest.get("source_hashes", {})) != expected_sources:
        raise ValueError("All-event inference lacks producer identities")
    source_paths = []
    for name, digest in manifest.get("source_hashes", {}).items():
        if Path(name).name != name:
            raise ValueError("Unexpected inference source path")
        source = Path(__file__).resolve().with_name(name)
        if _sha(source) != digest:
            raise ValueError(f"All-event inference producer changed: {name}")
        source_paths.append(str(source))
    events, predictions = _rows(original / "events.csv"), _rows(csv_path)
    identities = lambda rows: [(r["event_id"], r.get("frame_id", "")) for r in rows]
    if (identities(events) != identities(predictions)
            or len({r["event_id"] for r in events}) != len(events)
            or manifest.get("counts", {}).get("n_original_events") != len(events)):
        raise ValueError("All-event inference does not retain every event in original order")
    return csv_path, {
        "scope": "all_original_events",
        "manifest_sha256": _sha(manifest_path),
        "csv_sha256": _sha(csv_path),
        "tracking_paths": sorted(set([str(csv_path), str(manifest_path), str(Path(__file__).resolve()), *inputs, *source_paths])),
    }
