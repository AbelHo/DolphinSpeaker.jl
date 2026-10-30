"""Explicit derived fitter inputs; original human annotations remain authoritative.

This stdlib adapter changes eligible uncertain visibility to visible solely for a
frozen fitter. It never infers coordinates, shrinks uncertainty or emits new Terra
truth. Exporters must validate the bridge and render original annotation status.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
from datetime import datetime, timezone
from functools import lru_cache
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import tempfile


POLICY_NAME = "user_confirmed_uncertain_coordinates"
POLICY_SCHEMA = "clicker-annotation-fit-policy-v1"
AUTHORITY = "internal_eligibility_adapter_not_new_terra_annotation"
VISIBILITIES = ("absent", "excluded", "occluded", "uncertain", "visible")
MEDIUMS = {"underwater", "air", "mixed", "unknown"}
EXTRA_FIELDS = ("terra_visibility", "fit_eligibility_reason")
TRUTH_FILES = ("prepare_manifest.json", "annotations.csv", "annotations.complete.json", "human_review_selection.json")


@lru_cache(maxsize=1)
def _helpers():
    # Reuse only byte-copy/atomic-publication helpers; never call the human
    # revision constructor or manufacture a replacement user selection.
    path = Path(__file__).with_name("annotation_revision.py")
    spec = importlib.util.spec_from_file_location("_clicker_fit_policy_clone_helpers", path)
    if spec is None or spec.loader is None:
        raise ImportError("Missing frozen annotation revision helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha(path: Path) -> str:
    return _helpers()._sha(path)


def _json(path: Path) -> dict:
    return _helpers()._read_json(path)


def _annotations(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not fields or len(set(fields)) != len(fields) or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("Malformed annotation CSV schema/rows")
    if not set(_helpers().ANNOTATION_FIELDS) <= set(fields):
        raise ValueError("Missing required human annotation columns")
    return fields, rows


def _number(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _eligibility(row: dict, frame: dict, accepted: set[str]) -> tuple[bool, str]:
    if row["frame_id"] not in accepted:
        return False, "outside_user_confirmed_selection"
    if row["visibility"] not in {"visible", "uncertain"}:
        return False, "original_visibility_" + row["visibility"]
    if row["medium"] != "underwater":
        return False, "not_underwater"
    if not row["x"].strip() or not row["y"].strip():
        return False, "missing_coordinates"
    x, y = _number(row["x"]), _number(row["y"])
    if x is None or y is None:
        return False, "nonfinite_or_invalid_coordinates"
    if not (0 <= x < int(frame["width"]) and 0 <= y < int(frame["height"])):
        return False, "coordinates_out_of_bounds"
    uncertainty = _number(row["uncertainty_px"])
    if uncertainty is None or uncertainty <= 0:
        return False, "missing_or_invalid_positive_uncertainty"
    return True, ("admitted_user_confirmed_uncertain_coordinates" if row["visibility"] == "uncertain"
                  else "admitted_original_visible_coordinates")


def _source(original: Path) -> dict:
    h = _helpers()
    if (original / "annotation_fit_policy.json").exists():
        raise ValueError("Original directory must be human truth, not another fit adapter")
    state = {name: h._snapshot(original / name) for name in (*h.PREPARED_FILES, *TRUTH_FILES)}
    manifest = _json(original / "prepare_manifest.json")
    seal = _json(original / "annotations.complete.json")
    selection = _json(original / "human_review_selection.json")
    frames = h._read_csv(original / "frames.csv")
    fields, rows = _annotations(original / "annotations.csv")
    if set(fields) & set(EXTRA_FIELDS) or seal.get("annotation_authority") == AUTHORITY:
        raise ValueError("Original annotations already carry an internal fitter adapter")
    if seal.get("schema_version") != "clicker-annotations-v1":
        raise ValueError("Unsupported original annotation completion schema")
    hashes = {name: item["sha256"] for name, item in state.items()}
    for key, filename in (("frames_sha256", "frames.csv"), ("annotations_sha256", "annotations.csv"),
                          ("human_review_selection_sha256", "human_review_selection.json")):
        if seal.get(key) != hashes[filename]:
            raise ValueError(f"Stale original annotation completion: {key}")
    if seal.get("n_frames") != len(frames) or not frames:
        raise ValueError("Original annotation completion/frame count differs")
    ids = [f["frame_id"] for f in frames]
    annotation_ids = [r["frame_id"] for r in rows]
    if (not all(isinstance(i, str) and i.strip() for i in ids)
            or len(set(ids)) != len(ids) or len(set(annotation_ids)) != len(annotation_ids)
            or set(ids) != set(annotation_ids)):
        raise ValueError("Original annotation/frame identities differ or duplicate")
    if (selection.get("schema_version", selection.get("schema")) != "clicker-human-review-selection-v1"
            or selection.get("frames_sha256") != hashes["frames.csv"]):
        raise ValueError("Original human selection schema/frames differ")
    chosen = selection.get("accepted_frame_ids")
    if (not isinstance(chosen, list) or not chosen or not all(isinstance(i, str) for i in chosen)
            or len(set(chosen)) != len(chosen) or not set(chosen) <= set(ids)):
        raise ValueError("Original human selection contains invalid frame IDs")
    if "human_selected_frame_count" in seal and seal["human_selected_frame_count"] != len(chosen):
        raise ValueError("Sealed human selection count differs")
    if manifest.get("status") != "complete" or manifest.get("frame_count") != len(frames):
        raise ValueError("Original preparation is incomplete or frame count differs")
    for filename in h.PREPARED_FILES:
        entry = h._output_entry(manifest, original, filename)
        if entry["sha256"] != hashes[filename] or entry["bytes"] != state[filename]["bytes"]:
            raise ValueError(f"Prepared source differs from original manifest: {filename}")
    for filename, key in (("events.csv", "events_csv"), ("frames.csv", "frames_csv")):
        if Path(manifest[key]).resolve() != original / filename:
            raise ValueError(f"Original preparation path differs: {key}")
    cache = manifest["snippet_cache"]
    if (Path(cache["path"]).resolve() != original / "event_snippets.f64"
            or cache["sha256"] != hashes["event_snippets.f64"]
            or cache["bytes"] != state["event_snippets.f64"]["bytes"]):
        raise ValueError("Original cache binding differs")
    frame_map = {f["frame_id"]: f for f in frames}
    for f in frames:
        if int(f["width"]) <= 0 or int(f["height"]) <= 0:
            raise ValueError("Invalid original frame dimensions")
        path = Path(f["image_path"])
        if not path.is_absolute() or not path.is_file():
            raise ValueError("Original image references must be existing absolute paths")
    # Same finite coordinate/radius validation as annotation sealing, without
    # rewriting truth or assigning coordinates to missing/uncertain rows.
    for row in rows:
        if row["visibility"] not in VISIBILITIES or row["medium"] not in MEDIUMS or not row["annotator"].strip():
            raise ValueError("Invalid original annotation visibility/medium/annotator")
        frame = frame_map[row["frame_id"]]
        if row["visibility"] == "visible" or row["x"].strip() or row["y"].strip():
            x, y = _number(row["x"]), _number(row["y"])
            if x is None or y is None or not (0 <= x < int(frame["width"]) and 0 <= y < int(frame["height"])):
                raise ValueError("Invalid/nonfinite/out-of-bounds original coordinates")
        if row["visibility"] == "visible" or row["uncertainty_px"].strip():
            uncertainty = _number(row["uncertainty_px"])
            if uncertainty is None or uncertainty <= 0:
                raise ValueError("Invalid original positive uncertainty")
    counts = {v: sum(r["visibility"] == v for r in rows) for v in VISIBILITIES}
    if "counts" in seal and seal["counts"] != counts:
        raise ValueError("Original completion status counts differ")
    return dict(directory=original, state=state, hashes=hashes, manifest=manifest, seal=seal,
                selection=selection, frames=frames, frame_map=frame_map, fields=fields, rows=rows,
                accepted=set(chosen), counts=counts)


def _adapter(source: dict) -> dict:
    derived, decisions, admitted, converted, missing = [], [], [], [], []
    for row in source["rows"]:
        eligible, reason = _eligibility(row, source["frame_map"][row["frame_id"]], source["accepted"])
        # Frozen fitter ignores our reason column. Fail closed on any preserved
        # visible row that could enter it despite being outside this policy.
        if (not eligible and row["visibility"] == "visible" and row["medium"] == "underwater"
                and _number(row["x"]) is not None and _number(row["y"]) is not None):
            raise ValueError("Unsafe nonadmitted original-visible row would enter frozen fitter")
        adapted = dict(row)
        if eligible:
            admitted.append(row["frame_id"])
            if row["visibility"] == "uncertain":
                adapted["visibility"] = "visible"
                converted.append(row["frame_id"])
        if reason == "missing_coordinates":
            missing.append(row["frame_id"])
        adapted.update(terra_visibility=row["visibility"], fit_eligibility_reason=reason)
        derived.append(adapted)
        decisions.append(dict(frame_id=row["frame_id"], terra_visibility=row["visibility"],
                              fit_visibility=adapted["visibility"], admitted=eligible, reason=reason))
    return dict(rows=derived, decisions=decisions, admitted_frame_ids=admitted,
                converted_uncertain_frame_ids=converted, missing_coordinate_frame_ids=missing,
                original_statuses={r["frame_id"]: r["visibility"] for r in source["rows"]},
                original_visibility_counts=source["counts"],
                fit_visibility_counts={v: sum(r["visibility"] == v for r in derived) for v in VISIBILITIES})


def _assert_source_unchanged(source: dict) -> None:
    for entry in source["state"].values():
        path = Path(entry["path"])
        if _helpers()._fingerprint(path) != tuple(entry["fingerprint"]) or _sha(path) != entry["sha256"]:
            raise ValueError(f"Original human/prepared source changed during adapter operation: {path}")


def create_fit_annotation_input(
    original_dir: str | Path,
    fit_dir: str | Path,
    *,
    policy_name: str,
    reason: str,
) -> dict:
    """Create an explicit absent-directory fitter adapter, never a new annotation.

    Opt-in policy_name must be user_confirmed_uncertain_coordinates. Only the
    authenticated user selection with original visible/uncertain, underwater,
    finite in-bounds coordinates and finite positive radii can be adapter-visible.
    Original coordinate/radius/medium strings, notes, annotator, frame IDs and
    nonadmitted status are retained exactly. Original truth is archived byte-exact
    and never changed. Event timing, physical checks and split purges still apply
    in the frozen fitter: admitted annotation count is not a fitted-event count.
    """
    if policy_name != POLICY_NAME:
        raise ValueError(f"Explicit opt-in policy_name must be {POLICY_NAME!r}")
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("Explicit user-authorization/policy reason is required")
    original = Path(original_dir).resolve(strict=True)
    destination = Path(fit_dir).absolute()
    destination = destination.parent.resolve() / destination.name
    if os.path.lexists(destination):
        raise FileExistsError(f"Derived fit directory already exists: {destination}")
    if destination == original or original in destination.parents:
        raise ValueError("Derived fit directory must be outside original human truth")
    source = _source(original)
    adapted = _adapter(source)
    h = _helpers()
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.fit-adapter-", dir=destination.parent))
    try:
        for name in h.PREPARED_FILES:
            h._copy_exact(original / name, staging / name, source["hashes"][name])
        h._copy_exact(original / "human_review_selection.json", staging / "human_review_selection.json",
                      source["hashes"]["human_review_selection.json"])
        archive = staging / "fit_policy_source"
        archive.mkdir()
        for name in TRUTH_FILES:
            h._copy_exact(original / name, archive / name, source["hashes"][name])
        with (staging / "annotations.csv").open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=source["fields"] + list(EXTRA_FIELDS))
            writer.writeheader()
            writer.writerows(adapted["rows"])
            stream.flush()
            os.fsync(stream.fileno())
        created = datetime.now(timezone.utc).isoformat()
        helper_sha = _sha(Path(__file__))
        clone_helper_sha = _sha(Path(__file__).with_name("annotation_revision.py"))
        fit_annotation_sha = _sha(staging / "annotations.csv")
        policy = {
            "schema_version": POLICY_SCHEMA, "policy_name": POLICY_NAME,
            "annotation_authority": AUTHORITY, "not_new_terra_annotation": True,
            "original_dir": str(original), "fit_dir": str(destination),
            "created_utc": created, "reason": reason, "source_sha256": helper_sha,
            "clone_helper_sha256": clone_helper_sha,
            "original_annotations_sha256": source["hashes"]["annotations.csv"],
            "original_completion_sha256": source["hashes"]["annotations.complete.json"],
            "original_prepare_manifest_sha256": source["hashes"]["prepare_manifest.json"],
            "frames_sha256": source["hashes"]["frames.csv"],
            "selection_sha256": source["hashes"]["human_review_selection.json"],
            "fit_annotations_sha256": fit_annotation_sha,
            "prepared_source_sha256": {name: source["hashes"][name] for name in h.PREPARED_FILES},
            "admitted_frame_ids": adapted["admitted_frame_ids"],
            "admitted_count": len(adapted["admitted_frame_ids"]),
            "converted_uncertain_frame_ids": adapted["converted_uncertain_frame_ids"],
            "missing_coordinate_frame_ids": adapted["missing_coordinate_frame_ids"],
            "missing_coordinate_count": len(adapted["missing_coordinate_frame_ids"]),
            "all_original_statuses": adapted["original_statuses"],
            "original_visibility_counts": adapted["original_visibility_counts"],
            "fit_visibility_counts": adapted["fit_visibility_counts"],
            "per_frame": adapted["decisions"],
            "limitations": [
                "Eligibility adapter for an unchanged visible-only fitter; not a replacement or reseal of Terra truth.",
                "No coordinate inference, interpolation, uncertainty shrinkage, label-driven tuning or admission outside authenticated user selection.",
                "Admitted annotation IDs may still fail event synchronization, physics, time-boundary or blocked-evaluation checks in the frozen fitter.",
                "Final human review must display original annotations and their original uncertain status, using this bridge only for model/input provenance.",
            ],
        }
        h._write_json(staging / "annotation_fit_policy.json", policy)
        policy_sha = _sha(staging / "annotation_fit_policy.json")
        complete = copy.deepcopy(source["seal"])
        complete.update(
            annotations_sha256=fit_annotation_sha, frames_sha256=source["hashes"]["frames.csv"],
            n_frames=len(source["frames"]), counts=adapted["fit_visibility_counts"], created_utc=created,
            human_review_selection_sha256=source["hashes"]["human_review_selection.json"],
            human_selected_frame_count=len(source["accepted"]), annotation_authority=AUTHORITY,
            annotation_fit_policy_sha256=policy_sha, original_truth_directory=str(original),
            original_annotations_sha256=source["hashes"]["annotations.csv"],
            original_completion_sha256=source["hashes"]["annotations.complete.json"],
            original_selection_sha256=source["hashes"]["human_review_selection.json"],
        )
        h._write_json(staging / "annotations.complete.json", complete)
        cloned = copy.deepcopy(source["manifest"])
        cloned.pop("postprocessing", None)
        cloned.pop("producer_manifest_sha256", None)
        copied_paths = {str(original / name): name for name in h.PREPARED_FILES}
        for entry in cloned["outputs"]:
            name = copied_paths.get(str(Path(entry["path"]).resolve()))
            if name is not None:
                stat = (staging / name).stat()
                entry.update(path=str(destination / name), bytes=stat.st_size, sha256=source["hashes"][name],
                             mtime=stat.st_mtime, ctime=stat.st_ctime, inode=stat.st_ino)
        cloned["events_csv"] = str(destination / "events.csv")
        cloned["frames_csv"] = str(destination / "frames.csv")
        cloned["snippet_cache"]["path"] = str(destination / "event_snippets.f64")
        cloned["postprocessing"] = [{
            "stage": "annotation_fit_policy", "policy_name": POLICY_NAME,
            "annotation_authority": AUTHORITY, "not_new_terra_annotation": True,
            "created_utc": created, "source_manifest_sha256": source["hashes"]["prepare_manifest.json"],
            "saved_parent_manifest_path": str(destination / "fit_policy_source/prepare_manifest.json"),
            "saved_parent_postprocessing_locator": "fit_policy_source/prepare_manifest.json#/postprocessing",
            "original_annotations_sha256": source["hashes"]["annotations.csv"],
            "original_completion_sha256": source["hashes"]["annotations.complete.json"],
            "original_selection_sha256": source["hashes"]["human_review_selection.json"],
            "original_selection_copied_byte_exact": True,
            "fit_annotations_sha256": fit_annotation_sha,
            "fit_completion_sha256": _sha(staging / "annotations.complete.json"),
            "policy_sha256": policy_sha, "helper_sha256": helper_sha,
            "policy_path": str(destination / "annotation_fit_policy.json"), "reason": reason,
        }]
        h._write_json(staging / "prepare_manifest.json", cloned)
        _assert_source_unchanged(source)
        # Validate staged bytes against FINAL path identities before atomic publication.
        _validate(original, staging, logical_fit=destination)
        h._directory_fsync(archive)
        h._directory_fsync(staging)
        h._publish_absent(staging, destination)
    finally:
        if staging.is_dir():
            shutil.rmtree(staging)
    return validate_fit_annotation_policy(original, destination)


def _validate(original: Path, storage: Path, *, logical_fit: Path) -> dict:
    h = _helpers()
    source = _source(original)
    expected = _adapter(source)
    policy = _json(storage / "annotation_fit_policy.json")
    helper_sha = _sha(Path(__file__))
    dependency_sha = _sha(Path(__file__).with_name("annotation_revision.py"))
    if (policy.get("schema_version") != POLICY_SCHEMA or policy.get("policy_name") != POLICY_NAME
            or policy.get("annotation_authority") != AUTHORITY or policy.get("not_new_terra_annotation") is not True
            or policy.get("original_dir") != str(original) or policy.get("fit_dir") != str(logical_fit)
            or policy.get("source_sha256") != helper_sha or policy.get("clone_helper_sha256") != dependency_sha
            or not isinstance(policy.get("reason"), str) or not policy["reason"].strip()):
        raise ValueError("Unsupported/unbound fit-annotation policy or helper/source identity")
    for key, filename in (("original_annotations_sha256", "annotations.csv"),
                          ("original_completion_sha256", "annotations.complete.json"),
                          ("original_prepare_manifest_sha256", "prepare_manifest.json"),
                          ("frames_sha256", "frames.csv"), ("selection_sha256", "human_review_selection.json")):
        if policy.get(key) != source["hashes"][filename]:
            raise ValueError(f"Policy is stale against original truth: {key}")
    for filename in TRUTH_FILES:
        if (storage / "fit_policy_source" / filename).read_bytes() != (original / filename).read_bytes():
            raise ValueError(f"Archived original truth differs: {filename}")
    if (storage / "human_review_selection.json").read_bytes() != (original / "human_review_selection.json").read_bytes():
        raise ValueError("Original human selection was replaced or edited")
    prepared_hashes = {name: source["hashes"][name] for name in h.PREPARED_FILES}
    if policy.get("prepared_source_sha256") != prepared_hashes:
        raise ValueError("Policy prepared input bindings differ")
    for name in h.PREPARED_FILES:
        if _sha(storage / name) != prepared_hashes[name]:
            raise ValueError(f"Fit prepared bytes differ from original: {name}")
    fields, actual = _annotations(storage / "annotations.csv")
    if fields != source["fields"] + list(EXTRA_FIELDS) or len(actual) != len(expected["rows"]):
        raise ValueError("Fit adapter annotation schema/row count differs")
    for old, wanted, got in zip(source["rows"], expected["rows"], actual):
        # Exact CSV cell strings preserve numerical spelling/bits including -0,
        # exponents and whitespace; no tolerance permits moved coordinates.
        if got != wanted:
            raise ValueError(f"Fit annotation row/coordinate/status equivalence failed: {old['frame_id']}")
    for key, value in (
        ("admitted_frame_ids", expected["admitted_frame_ids"]), ("admitted_count", len(expected["admitted_frame_ids"])),
        ("converted_uncertain_frame_ids", expected["converted_uncertain_frame_ids"]),
        ("missing_coordinate_frame_ids", expected["missing_coordinate_frame_ids"]),
        ("missing_coordinate_count", len(expected["missing_coordinate_frame_ids"])),
        ("all_original_statuses", expected["original_statuses"]),
        ("original_visibility_counts", expected["original_visibility_counts"]),
        ("fit_visibility_counts", expected["fit_visibility_counts"]), ("per_frame", expected["decisions"]),
    ):
        if policy.get(key) != value:
            raise ValueError(f"Policy eligibility differs from sealed human truth: {key}")
    fit_sha = _sha(storage / "annotations.csv")
    if policy.get("fit_annotations_sha256") != fit_sha:
        raise ValueError("Fit annotation hash differs from policy")
    policy_sha = _sha(storage / "annotation_fit_policy.json")
    completion = _json(storage / "annotations.complete.json")
    required = dict(
        schema_version="clicker-annotations-v1", frames_sha256=source["hashes"]["frames.csv"],
        annotations_sha256=fit_sha, n_frames=len(source["frames"]), counts=expected["fit_visibility_counts"],
        human_review_selection_sha256=source["hashes"]["human_review_selection.json"],
        human_selected_frame_count=len(source["accepted"]), annotation_authority=AUTHORITY,
        annotation_fit_policy_sha256=policy_sha, original_truth_directory=str(original),
        original_annotations_sha256=source["hashes"]["annotations.csv"],
        original_completion_sha256=source["hashes"]["annotations.complete.json"],
        original_selection_sha256=source["hashes"]["human_review_selection.json"],
    )
    expected_completion = copy.deepcopy(source["seal"])
    expected_completion.update(required, created_utc=policy.get("created_utc"))
    if completion != expected_completion:
        raise ValueError("Derived completion is stale or lacks original-truth/policy lineage")
    clone = _json(storage / "prepare_manifest.json")
    for key, value in source["manifest"].items():
        if key in {"outputs", "events_csv", "frames_csv", "snippet_cache", "postprocessing", "producer_manifest_sha256"}:
            continue
        if clone.get(key) != value:
            raise ValueError(f"Prepared producer binding changed: {key}")
    for filename, key in (("events.csv", "events_csv"), ("frames.csv", "frames_csv")):
        if clone.get(key) != str(logical_fit / filename):
            raise ValueError("Fit canonical table path differs")
    if clone.get("snippet_cache") != {**source["manifest"]["snippet_cache"], "path": str(logical_fit / "event_snippets.f64")}:
        raise ValueError("Fit cache identity/path changed")
    original_outputs = source["manifest"]["outputs"]
    if len(clone.get("outputs", [])) != len(original_outputs):
        raise ValueError("Prepared output inventory length changed")
    copied = {str(original / name): name for name in h.PREPARED_FILES}
    for before, after in zip(original_outputs, clone["outputs"]):
        filename = copied.get(str(Path(before["path"]).resolve()))
        if filename is None:
            if after != before:
                raise ValueError("Original image/noncopied producer identity changed")
        elif (after.get("path") != str(logical_fit / filename)
              or after.get("sha256") != before["sha256"] or after.get("bytes") != before["bytes"]):
            raise ValueError("Copied preparation output identity changed")
    history = clone.get("postprocessing", [])
    if len(history) != 1 or history[0].get("stage") != "annotation_fit_policy" or "producer_manifest_sha256" in clone:
        raise ValueError("Fit clone carries incompatible original timing/revision lineage")
    lineage = history[0]
    fit_completion_sha = _sha(storage / "annotations.complete.json")
    checks = {
        "policy_name": POLICY_NAME, "annotation_authority": AUTHORITY, "not_new_terra_annotation": True,
        "source_manifest_sha256": source["hashes"]["prepare_manifest.json"],
        "saved_parent_manifest_path": str(logical_fit / "fit_policy_source/prepare_manifest.json"),
        "saved_parent_postprocessing_locator": "fit_policy_source/prepare_manifest.json#/postprocessing",
        "original_annotations_sha256": source["hashes"]["annotations.csv"],
        "original_completion_sha256": source["hashes"]["annotations.complete.json"],
        "original_selection_sha256": source["hashes"]["human_review_selection.json"],
        "original_selection_copied_byte_exact": True, "fit_annotations_sha256": fit_sha,
        "fit_completion_sha256": fit_completion_sha, "policy_sha256": policy_sha,
        "helper_sha256": helper_sha, "policy_path": str(logical_fit / "annotation_fit_policy.json"), "reason": policy["reason"],
    }
    if any(lineage.get(key) != value for key, value in checks.items()):
        raise ValueError("Fit preparation postprocessing lineage differs")
    _assert_source_unchanged(source)
    return dict(
        valid=True, policy_name=POLICY_NAME, annotation_authority=AUTHORITY,
        original_dir=str(original), fit_dir=str(logical_fit),
        original_annotations_sha256=source["hashes"]["annotations.csv"], fit_annotations_sha256=fit_sha,
        original_completion_sha256=source["hashes"]["annotations.complete.json"],
        fit_completion_sha256=fit_completion_sha, policy_sha256=policy_sha, helper_sha256=helper_sha,
        clone_helper_sha256=dependency_sha, frames_sha256=source["hashes"]["frames.csv"],
        selection_sha256=source["hashes"]["human_review_selection.json"],
        original_prepare_manifest_sha256=source["hashes"]["prepare_manifest.json"],
        fit_prepare_manifest_sha256=_sha(storage / "prepare_manifest.json"),
        policy_path=str(logical_fit / "annotation_fit_policy.json"),
        admitted_frame_ids=expected["admitted_frame_ids"], admitted_count=len(expected["admitted_frame_ids"]),
        missing_coordinate_frame_ids=expected["missing_coordinate_frame_ids"],
        missing_coordinate_count=len(expected["missing_coordinate_frame_ids"]),
        converted_uncertain_frame_ids=expected["converted_uncertain_frame_ids"],
        original_visibility_counts=expected["original_visibility_counts"], fit_visibility_counts=expected["fit_visibility_counts"],
    )


def validate_fit_annotation_policy(original_dir: str | Path, fit_dir: str | Path) -> dict:
    """Return a strict original-truth/derived-input bridge, or raise ValueError.

    Recomputes eligibility from the currently sealed original selection, checks
    all archived and prepared byte identities, all original/derived row strings,
    policy/helper hashes and both completion/manifest lineages. This validates no
    model by itself: exporters must ALSO match model annotation/frames hashes to
    the returned fit-input hashes, and render original truth (uncertain = amber).
    """
    original = Path(original_dir).resolve(strict=True)
    fit = Path(fit_dir).resolve(strict=True)
    if original == fit or original in fit.parents:
        raise ValueError("Derived fit directory is not separate from original human truth")
    return _validate(original, fit, logical_fit=fit)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    create = sub.add_parser("create")
    create.add_argument("original_dir", type=Path)
    create.add_argument("fit_dir", type=Path)
    create.add_argument("--policy-name", required=True)
    create.add_argument("--reason", required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("original_dir", type=Path)
    validate.add_argument("fit_dir", type=Path)
    args = vars(parser.parse_args())
    action = args.pop("action")
    result = create_fit_annotation_input(**args) if action == "create" else validate_fit_annotation_policy(**args)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
