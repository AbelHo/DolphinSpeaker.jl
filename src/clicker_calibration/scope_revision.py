"""Restrict a sealed Terra review to a user-confirmed gallery interval.

Selected CSV records retain their original bytes, coordinates, uncertainty and
status. This is a scope correction, not annotation or coordinate reinspection.
Only new, atomically published directories are supported. Stdlib only.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile


SCHEMA = "clicker-annotation-scope-revision-v1"
AUTHORITY = "original_terra_coordinates_with_user_scope_exclusions"
VISIBILITIES = {"visible", "uncertain", "occluded", "absent", "excluded"}


@lru_cache(maxsize=1)
def _helpers():
    path = Path(__file__).with_name("annotation_revision.py")
    spec = importlib.util.spec_from_file_location("_clicker_scope_clone_helpers", path)
    if spec is None or spec.loader is None:
        raise ImportError("Missing immutable annotation revision helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _annotation_records(path: Path):
    """Return parsed fields/rows plus exact UTF-8 CSV header/record bytes.

    DictReader line_num includes quoted multiline records; StringIO.readlines
    recognizes CSV physical line endings without splitting other Unicode chars.
    A UTF-8 BOM is retained in the header bytes but removed from its field name.
    """
    text = path.read_bytes().decode("utf-8")
    lines = io.StringIO(text, newline="").readlines()
    # Strip the BOM before CSV parsing, otherwise it precedes an opening quote
    # and changes quoted-header parsing. Raw spans below still retain that BOM.
    stream = io.StringIO(text.removeprefix("\ufeff"), newline="")
    reader = csv.DictReader(stream)
    fields = list(reader.fieldnames or [])
    if fields:
        fields[0] = fields[0].removeprefix("\ufeff")
        reader.fieldnames = fields
    if not fields or len(set(fields)) != len(fields) or not set(_helpers().ANNOTATION_FIELDS) <= set(fields):
        raise ValueError("Missing/duplicate original annotation columns")
    end = reader.line_num
    header = "".join(lines[:end]).encode("utf-8")
    records = []
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("Malformed original annotation record")
        raw = "".join(lines[end:reader.line_num]).encode("utf-8")
        end = reader.line_num
        records.append((row, raw))
    ids = [row["frame_id"] for row, _ in records]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError("Empty/duplicate original annotation identities")
    return fields, header, records


def _finite(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid original {name}") from error
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite original {name}")
    return number


def _validate_selected(row, frame):
    if row["medium"] != "underwater":
        raise ValueError(f"Selected frame {row['frame_id']} is not already labeled underwater; refusing to relabel its medium")
    if row["visibility"] not in VISIBILITIES or not row["annotator"].strip():
        raise ValueError("Selected original annotation has invalid visibility/annotator")
    width, height = int(frame["width"]), int(frame["height"])
    if width <= 0 or height <= 0:
        raise ValueError("Invalid original frame dimensions")
    has_coordinates = bool(row["x"].strip() or row["y"].strip())
    if has_coordinates or row["visibility"] == "visible":
        x, y = _finite(row["x"], "x"), _finite(row["y"], "y")
        if not (0 <= x < width and 0 <= y < height):
            raise ValueError("Original selected coordinates are outside the source image")
    if row["uncertainty_px"].strip() or row["visibility"] == "visible":
        if _finite(row["uncertainty_px"], "uncertainty_px") <= 0:
            raise ValueError("Original uncertainty radius must be positive")
    return has_coordinates


def _replace_owned_json(staging: Path, name: str, value: dict):
    """Replace a known constructor-owned file before public publication only."""
    temporary = staging / (name + ".scope-tmp")
    _helpers()._write_json(temporary, value)
    os.replace(temporary, staging / name)


def _assert_snapshots(snapshots):
    h = _helpers()
    for item in snapshots.values():
        path = Path(item["path"])
        if h._fingerprint(path) != tuple(item["fingerprint"]) or h._sha(path) != item["sha256"]:
            raise ValueError(f"Source changed during scope revision: {path}")


def create_scope_revision(source_dir, destination_dir, *, first_index, last_index,
                          gallery_revision, reason):
    """Create a sealed, scope-only revision; never replace an existing directory.

    `first_index` and `last_index` are inclusive numbered gallery JPG filenames.
    The existing clone helper authenticates the gallery HTML/frame mapping and
    all five prepared inputs. Selected rows remain byte-exact original Terra
    records, including uncertain/blank estimates; their medium must already be
    underwater. Outside rows become excluded with blank x/y/radius, retaining
    other original fields. All frames/events remain for provenance; the human
    selection marks publish_selected_only=true for presentation consumers.

    Work is private staging until a final fsynced annotation seal and immutable
    source recheck, followed by atomic no-replace publication. No new annotation,
    image inference, fitting, image decode, or mutation of the source occurs.
    """
    h = _helpers()
    source = Path(source_dir).resolve(strict=True)
    destination = Path(destination_dir).absolute()
    destination = destination.parent.resolve() / destination.name
    if os.path.lexists(destination):
        raise FileExistsError(f"Scope revision destination already exists: {destination}")
    if not source.is_dir() or source == destination or source in destination.parents:
        raise ValueError("Scope destination must be separate from the immutable source")
    if isinstance(first_index, bool) or isinstance(last_index, bool) or not isinstance(first_index, int) or not isinstance(last_index, int):
        raise ValueError("Gallery indices must be integers")
    if not 0 <= first_index <= last_index:
        raise ValueError("Gallery range must be nonnegative, ordered and inclusive")
    if not isinstance(gallery_revision, str) or not re.fullmatch(r"[0-9a-f]{12}", gallery_revision):
        raise ValueError("Expected a sealed 12-character gallery revision")
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("An explicit user-confirmed scope reason is required")
    critical = (*h.PREPARED_FILES, *h.ARCHIVED_FILES, "annotation_review.json")
    snapshots = {name: h._snapshot(source / name) for name in critical}
    original_seal = h._read_json(source / "annotations.complete.json")
    if gallery_revision != original_seal.get("annotations_sha256", "")[:12]:
        raise ValueError("Requested gallery revision differs from the original sealed annotations")
    gallery_path = source / "annotation_review" / gallery_revision / "index.html"
    snapshots["gallery_index"] = h._snapshot(gallery_path)
    snapshots["scope_helper"] = h._snapshot(Path(__file__).resolve())
    snapshots["clone_helper"] = h._snapshot(Path(h.__file__).resolve())
    fields, header, records = _annotation_records(source / "annotations.csv")
    destination.parent.mkdir(parents=True, exist_ok=True)
    outer = Path(tempfile.mkdtemp(prefix=f".{destination.name}.scope-", dir=destination.parent))
    staging = outer / "revision"
    try:
        # The frozen helper authenticates producer bytes and resolves filenames
        # through HTML. Its temporary pending reannotation is never published.
        h.create_annotation_revision(source, staging, first_index=first_index, last_index=last_index,
                                     gallery_revision=gallery_revision, reason=reason)
        _assert_snapshots(snapshots)
        selection = h._read_json(staging / "human_review_selection.json")
        accepted = selection["accepted_frame_ids"]
        accepted_set = set(accepted)
        frame_map = {row["frame_id"]: row for row in h._read_csv(staging / "frames.csv")}
        coordinates = sum(_validate_selected(row, frame_map[row["frame_id"]])
                          for row, _ in records if row["frame_id"] in accepted_set)
        created = datetime.now(timezone.utc).isoformat()
        authority = dict(
            source_annotations_path=str(source / "annotations.csv"),
            source_annotations_sha256=snapshots["annotations.csv"]["sha256"],
            source_completion_path=str(source / "annotations.complete.json"),
            source_completion_sha256=snapshots["annotations.complete.json"]["sha256"],
            archived_annotations_path=str(destination / "revision_source/annotations.csv"),
            archived_completion_path=str(destination / "revision_source/annotations.complete.json"),
            selected_records_byte_exact=True, coordinates_reestimated=False,
            visibility_preserved=True, uncertainty_preserved=True, notes_and_annotator_preserved=True,
        )
        selection.update(coordinate_reannotation_required=False, publish_selected_only=True,
                         scope_authority="user_confirmed_underwater_gallery_interval",
                         coordinate_authority="unchanged_original_sealed_terra_estimates",
                         original_coordinate_source=authority, scope_revision_schema=SCHEMA)
        _replace_owned_json(staging, "human_review_selection.json", selection)
        derived_rows = []
        selected_record_hashes = {}
        with (staging / "annotations.csv").open("xb") as stream:
            stream.write(header)
            for original, raw in records:
                row = dict(original)
                if row["frame_id"] in accepted_set:
                    stream.write(raw)
                    selected_record_hashes[row["frame_id"]] = hashlib.sha256(raw).hexdigest()
                else:
                    row.update(x="", y="", uncertainty_px="", visibility="excluded", annotator="user_scope")
                    buffer = io.StringIO(newline="")
                    csv.DictWriter(buffer, fieldnames=fields, lineterminator="\r\n").writerow(row)
                    stream.write(buffer.getvalue().encode("utf-8"))
                derived_rows.append(row)
            stream.flush()
            os.fsync(stream.fileno())
        (staging / "annotations.pending.csv").unlink()
        counts = {status: sum(row["visibility"] == status for row in derived_rows) for status in sorted(VISIBILITIES)}
        selection_sha = h._sha(staging / "human_review_selection.json")
        annotation_sha = h._sha(staging / "annotations.csv")
        scope = dict(schema_version=SCHEMA, status="complete", created_utc=created,
                     source_directory=str(source), output_directory=str(destination), reason=reason,
                     first_index=first_index, last_index=last_index, gallery_revision=gallery_revision,
                     original_gallery=str(gallery_path), original_gallery_sha256=snapshots["gallery_index"]["sha256"],
                     accepted_frame_ids=accepted, accepted_frame_count=len(accepted),
                     selected_coordinate_count=coordinates, selected_missing_coordinate_count=len(accepted)-coordinates,
                     excluded_frame_count=len(records)-len(accepted), visibility_counts=counts,
                     annotation_authority=AUTHORITY, original_coordinate_source=authority,
                     new_annotation=False, coordinate_reannotation_required=False, publish_selected_only=True,
                     source_manifest_sha256=snapshots["prepare_manifest.json"]["sha256"],
                     original_annotations_sha256=snapshots["annotations.csv"]["sha256"],
                     original_completion_sha256=snapshots["annotations.complete.json"]["sha256"],
                     frames_sha256=snapshots["frames.csv"]["sha256"], annotations_sha256=annotation_sha,
                     human_review_selection_sha256=selection_sha, selected_record_sha256=selected_record_hashes,
                     prepared_files_sha256={name: snapshots[name]["sha256"] for name in h.PREPARED_FILES},
                     clone_helper_sha256=snapshots["clone_helper"]["sha256"],
                     scope_helper_sha256=snapshots["scope_helper"]["sha256"],
                     limits=["User scope confirmation does not upgrade uncertain coordinates to visible or establish acoustic accuracy.",
                             "All original frame/event rows remain authenticated; presentation must honor accepted_frame_ids."])
        h._write_json(staging / "annotation_scope_revision.json", scope)
        scope_sha = h._sha(staging / "annotation_scope_revision.json")
        cloned = h._read_json(staging / "prepare_manifest.json")
        owned = {str(staging / name): name for name in h.PREPARED_FILES}
        for entry in cloned["outputs"]:
            name = owned.get(entry["path"])
            if name is not None:
                entry["path"] = str(destination / name)
        cloned["events_csv"] = str(destination / "events.csv")
        cloned["frames_csv"] = str(destination / "frames.csv")
        cloned["snippet_cache"]["path"] = str(destination / "event_snippets.f64")
        lineage = cloned["postprocessing"]
        if len(lineage) != 1 or lineage[0]["stage"] != "annotation_revision":
            raise ValueError("Unexpected clone helper lineage")
        lineage[0].update(
            saved_parent_manifest_path=str(destination / "revision_source/prepare_manifest.json"),
            human_review_selection_path=str(destination / "human_review_selection.json"),
            human_review_selection_sha256=selection_sha, coordinate_reannotation_required=False,
            publish_selected_only=True, annotation_authority=AUTHORITY,
            scope_revision_manifest_path=str(destination / "annotation_scope_revision.json"),
            scope_revision_manifest_sha256=scope_sha,
            scope_revision_source_sha256=snapshots["scope_helper"]["sha256"],
            original_coordinate_source=authority,
        )
        # review_source_sha256 remains the actual immutable cloning helper SHA.
        if lineage[0]["review_source_sha256"] != snapshots["clone_helper"]["sha256"]:
            raise ValueError("Clone helper identity changed")
        _replace_owned_json(staging, "prepare_manifest.json", cloned)
        # Validate resulting selected records, not only parsed float values.
        _, _, final_records = _annotation_records(staging / "annotations.csv")
        if [row for row, _ in final_records] != derived_rows:
            raise ValueError("Derived annotation serialization changed field values")
        for row, raw in final_records:
            if row["frame_id"] in accepted_set and hashlib.sha256(raw).hexdigest() != selected_record_hashes[row["frame_id"]]:
                raise ValueError("Selected original CSV record bytes changed")
        for name in h.PREPARED_FILES:
            if h._sha(staging / name) != snapshots[name]["sha256"]:
                raise ValueError("Scope revision changed prepared bytes")
        seal = dict(schema_version="clicker-annotations-v1", frames_sha256=snapshots["frames.csv"]["sha256"],
                    annotations_sha256=annotation_sha, n_frames=len(records), counts=counts,
                    annotator=sorted({row["annotator"] for row in derived_rows}), created_utc=created,
                    coordinate_convention=original_seal.get("coordinate_convention", "zero-origin original pixels; x right, y down"),
                    landmark_convention=original_seal.get("landmark_convention", "unchanged source Terra landmark convention"),
                    human_review_selection_sha256=selection_sha, human_selected_frame_count=len(accepted),
                    annotation_scope_revision_sha256=scope_sha, prepare_manifest_sha256=h._sha(staging / "prepare_manifest.json"),
                    annotation_authority=AUTHORITY, original_coordinate_source=authority,
                    coordinate_reannotation_required=False, publish_selected_only=True)
        # Seal last; the destination remains absent until everything is complete.
        h._write_json(staging / "annotations.complete.json", seal)
        _assert_snapshots(snapshots)
        h._directory_fsync(staging / "revision_source")
        h._directory_fsync(staging)
        h._publish_absent(staging, destination)
    finally:
        shutil.rmtree(outer, ignore_errors=True)
    return dict(status="complete", schema_version=SCHEMA, source_dir=str(source), output_dir=str(destination),
                accepted_frame_count=len(accepted), selected_coordinate_count=coordinates,
                selected_missing_coordinate_count=len(accepted)-coordinates, excluded_frame_count=len(records)-len(accepted),
                coordinate_reannotation_required=False, publish_selected_only=True,
                annotations_sha256=annotation_sha, human_review_selection_sha256=selection_sha,
                annotation_scope_revision_sha256=scope_sha,
                annotations_complete_sha256=h._sha(destination / "annotations.complete.json"),
                prepare_manifest_sha256=h._sha(destination / "prepare_manifest.json"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("destination_dir", type=Path)
    parser.add_argument("--first-index", type=int, required=True)
    parser.add_argument("--last-index", type=int, required=True)
    parser.add_argument("--gallery-revision", required=True)
    parser.add_argument("--reason", required=True)
    print(json.dumps(create_scope_revision(**vars(parser.parse_args()))))


if __name__ == "__main__":
    main()
