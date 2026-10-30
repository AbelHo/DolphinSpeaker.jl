"""Fork prepared clicker inputs for a user-selected gallery interval.

Only the five prepared data files are copied; coordinates are never inherited.
Original images remain absolute, authenticated producer-manifest references.
No import-time I/O and no third-party dependencies.
"""
from __future__ import annotations

import argparse
import copy
import csv
import ctypes
from datetime import datetime, timezone
import errno
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile


PREPARED_FILES = (
    "events.csv", "frames.csv", "synchronization.json", "prepare_metadata.json",
    "event_snippets.f64",
)
ANNOTATION_FIELDS = (
    "frame_id", "x", "y", "visibility", "medium", "uncertainty_px", "note", "annotator",
)
ARCHIVED_FILES = ("prepare_manifest.json", "annotations.csv", "annotations.complete.json")


def _sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Missing or duplicate CSV header: {path}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"Malformed CSV row: {path}")
    return rows


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _fingerprint(path: Path) -> tuple[int, int, int, int, int]:
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def _snapshot(path: Path) -> dict:
    before = _fingerprint(path)
    digest = _sha(path)
    if _fingerprint(path) != before:
        raise ValueError(f"Source changed while hashing: {path}")
    return {"path": str(path), "sha256": digest, "bytes": before[2],
            "fingerprint": before}


def _write_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def _copy_exact(source: Path, destination: Path, expected_sha256: str) -> None:
    with source.open("rb") as original, destination.open("xb") as clone:
        shutil.copyfileobj(original, clone, length=1024 * 1024)
        clone.flush()
        os.fsync(clone.fileno())
    if _sha(destination) != expected_sha256:
        raise ValueError(f"Copied bytes differ from authenticated source: {source}")


def _directory_fsync(path: Path) -> None:
    if os.name == "posix":
        descriptor = os.open(path, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def _publish_absent(source: Path, destination: Path) -> None:
    """Publish atomically without replacing even an empty raced-in directory."""
    if sys.platform.startswith("linux"):
        libc = ctypes.CDLL(None, use_errno=True)
        rename = getattr(libc, "renameat2", None)
        if rename is None:
            raise OSError(errno.ENOSYS, "Atomic no-replace directory publication unavailable")
        rename.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint)
        rename.restype = ctypes.c_int
        # Linux AT_FDCWD and RENAME_NOREPLACE; neither path can be overwritten.
        result = rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
        if result:
            code = ctypes.get_errno()
            raise OSError(code, os.strerror(code), str(destination))
    elif os.name == "nt":
        os.rename(source, destination)  # Windows refuses an existing destination.
    else:
        raise OSError(errno.ENOSYS, "Atomic absent-directory publication unsupported on this platform")
    _directory_fsync(destination.parent)


class _GalleryCards(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_article = False
        self.index: int | None = None
        self.frame: str | None = None
        self.cards: list[tuple[int, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        fields = dict(attrs)
        if tag == "article":
            if self.in_article:
                raise ValueError("Nested gallery article")
            self.in_article, self.index, self.frame = True, None, None
        elif self.in_article and tag == "a":
            match = re.fullmatch(r"([0-9]{5,})\.jpg", fields.get("href") or "")
            if match:
                value = int(match.group(1))
                if self.index is not None or fields["href"] != f"{value:05d}.jpg":
                    raise ValueError("Duplicate/noncanonical gallery filename index")
                self.index = value
        elif self.in_article and tag == "img":
            if self.index is None or self.frame is not None:
                raise ValueError("Gallery card image has ambiguous filename identity")
            if fields.get("src") != f"{self.index:05d}_thumb.jpg" or not fields.get("alt"):
                raise ValueError("Gallery thumbnail/frame identity differs")
            self.frame = fields["alt"]

    def handle_endtag(self, tag: str) -> None:
        if tag == "article":
            if not self.in_article or self.index is None or self.frame is None:
                raise ValueError("Incomplete gallery frame card")
            self.cards.append((self.index, self.frame))
            self.in_article = False


def _gallery_mapping(source: Path, revision: str, frames: list[dict], seal: dict,
                     manifest: dict) -> tuple[dict[int, str], Path]:
    expected = source / "annotation_review" / revision / "index.html"
    if manifest.get("revision") != revision or revision != seal["annotations_sha256"][:12]:
        raise ValueError("Requested gallery does not match the current sealed annotation revision")
    if Path(manifest.get("gallery", "")).resolve() != expected.resolve():
        raise ValueError("Requested gallery path does not match annotation_review.json")
    for key in ("frames_sha256", "annotations_sha256", "n_frames"):
        if manifest.get(key) != seal.get(key):
            raise ValueError(f"Gallery is stale against annotation completion: {key}")
    parser = _GalleryCards()
    parser.feed(expected.read_text(encoding="utf-8"))
    parser.close()
    if parser.in_article:
        raise ValueError("Unclosed gallery frame card")
    mapping = dict(parser.cards)
    ids = [row["frame_id"] for row in frames]
    if len(mapping) != len(parser.cards) or set(mapping) != set(range(len(ids))):
        raise ValueError("Gallery filename indices are not a unique zero-origin sequence")
    if len(set(mapping.values())) != len(ids) or set(mapping.values()) != set(ids):
        raise ValueError("Gallery frame identities do not match current frames.csv")
    for index in mapping:
        if not (expected.parent / f"{index:05d}.jpg").is_file():
            raise ValueError(f"Missing reviewed gallery filename: {index:05d}.jpg")
    return mapping, expected


def _output_entry(manifest: dict, source: Path, name: str) -> dict:
    matches = [entry for entry in manifest.get("outputs", [])
               if Path(entry["path"]).resolve() == source / name]
    if len(matches) != 1:
        raise ValueError(f"Prepared manifest lacks a unique output binding: {name}")
    return matches[0]


def create_annotation_revision(
    source_dir: str | Path,
    destination_dir: str | Path,
    *,
    first_index: int,
    last_index: int,
    gallery_revision: str,
    reason: str,
) -> dict:
    """Create an absent revision directory from an authenticated review gallery.

    Indices are the inclusive, zero-origin numbered JPG filenames in the named
    annotation gallery, resolved through its HTML frame IDs (not assumed CSV
    row order). Completion/frame/annotation/gallery and copied producer output
    hashes must agree. Prior image references must already be absolute.

    The five prepared inputs remain byte-exact. All original preparation
    postprocessing stays in revision_source/prepare_manifest.json; the clone has
    a fresh annotation_revision lineage, so native timing review can run anew.
    All pending coordinates are empty. This function neither seals annotations
    nor copies/decodes images nor mutates the source. Producer image identities
    are retained, not rehashed as a new scientific qualification.

    Publication uses Linux renameat2(RENAME_NOREPLACE), or Windows no-replace
    rename, after fsynced staging in the destination parent. Unsupported atomic
    publication fails explicitly; an existing revision is never overwritten.
    """
    source = Path(source_dir).resolve(strict=True)
    destination = Path(destination_dir).absolute()
    if not source.is_dir():
        raise ValueError("Annotation source must be a directory")
    # Resolve the parent without following/reinterpreting a destination symlink.
    destination = destination.parent.resolve() / destination.name
    if os.path.lexists(destination):
        raise FileExistsError(f"Revision destination already exists: {destination}")
    if source == destination or source in destination.parents:
        raise ValueError("Revision destination cannot be inside the immutable source")
    if isinstance(first_index, bool) or isinstance(last_index, bool) or not isinstance(first_index, int) or not isinstance(last_index, int):
        raise ValueError("Gallery indices must be integers")
    if not 0 <= first_index <= last_index:
        raise ValueError("Gallery range must be nonnegative, ordered and inclusive")
    if not isinstance(gallery_revision, str) or not re.fullmatch(r"[0-9a-f]{12}", gallery_revision):
        raise ValueError("Expected the 12-character sealed annotation gallery revision")
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("An explicit user-selection reason is required")

    critical = (*PREPARED_FILES, *ARCHIVED_FILES, "annotation_review.json")
    snapshots = {name: _snapshot(source / name) for name in critical}
    parent = _read_json(source / "prepare_manifest.json")
    seal = _read_json(source / "annotations.complete.json")
    gallery = _read_json(source / "annotation_review.json")
    frames = _read_csv(source / "frames.csv")
    annotations = _read_csv(source / "annotations.csv")
    frame_ids = [row["frame_id"] for row in frames]
    annotation_ids = [row["frame_id"] for row in annotations]
    if not frame_ids or not all(isinstance(identity, str) and identity.strip() for identity in frame_ids) or len(set(frame_ids)) != len(frame_ids):
        raise ValueError("Frame IDs must be nonempty and unique")
    if len(set(annotation_ids)) != len(annotation_ids) or set(annotation_ids) != set(frame_ids):
        raise ValueError("Sealed annotations must cover exactly the current frame IDs")
    if seal.get("schema_version") != "clicker-annotations-v1":
        raise ValueError("Unsupported annotation completion schema")
    if (seal.get("frames_sha256") != snapshots["frames.csv"]["sha256"]
            or seal.get("annotations_sha256") != snapshots["annotations.csv"]["sha256"]
            or seal.get("n_frames") != len(frames)):
        raise ValueError("Source annotation completion hashes/current frames differ")
    if parent.get("status") != "complete":
        raise ValueError("Revision source preparation is incomplete")
    events = _read_csv(source / "events.csv")
    if parent.get("event_count") != len(events) or parent.get("frame_count") != len(frames):
        raise ValueError("Preparation row counts differ from current tables")
    for name in PREPARED_FILES:
        entry = _output_entry(parent, source, name)
        if entry["sha256"] != snapshots[name]["sha256"] or entry["bytes"] != snapshots[name]["bytes"]:
            raise ValueError(f"Prepared source hash differs from its manifest: {name}")
    for name, key in (("events.csv", "events_csv"), ("frames.csv", "frames_csv")):
        if Path(parent.get(key, "")).resolve() != source / name:
            raise ValueError(f"Source manifest has a different canonical {key} path")
    cache = parent["snippet_cache"]
    if (Path(cache["path"]).resolve() != source / "event_snippets.f64"
            or cache["sha256"] != snapshots["event_snippets.f64"]["sha256"]
            or cache["bytes"] != snapshots["event_snippets.f64"]["bytes"]):
        raise ValueError("Source snippet cache binding differs")
    by_id = {row["frame_id"]: row for row in annotations}
    frame_by_id = {row["frame_id"]: row for row in frames}
    for event in events:
        identity = event.get("frame_id", "")
        if identity and (identity not in frame_by_id or event.get("image_path") != frame_by_id[identity]["image_path"]):
            raise ValueError("Prepared event/frame original image pointers differ")
    for frame in frames:
        image = Path(frame["image_path"])
        if not image.is_absolute() or not image.is_file():
            raise ValueError("Source frames must retain existing absolute original image paths")
        if by_id[frame["frame_id"]].get("medium") not in {"underwater", "air", "mixed", "unknown"}:
            raise ValueError("Prior annotation has an invalid medium")
    mapping, gallery_path = _gallery_mapping(source, gallery_revision, frames, seal, gallery)
    snapshots["gallery_index"] = _snapshot(gallery_path)
    if last_index >= len(mapping):
        raise ValueError("Gallery interval is outside the zero-origin filename range")
    accepted = [mapping[index] for index in range(first_index, last_index + 1)]
    accepted_set = set(accepted)

    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent))
    try:
        for name in PREPARED_FILES:
            _copy_exact(source / name, staging / name, snapshots[name]["sha256"])
        archive = staging / "revision_source"
        archive.mkdir()
        for name in ARCHIVED_FILES:
            _copy_exact(source / name, archive / name, snapshots[name]["sha256"])
        created = datetime.now(timezone.utc).isoformat()
        selection = {
            "schema_version": "clicker-human-review-selection-v1",
            "frames_sha256": snapshots["frames.csv"]["sha256"],
            "previous_annotations_sha256": snapshots["annotations.csv"]["sha256"],
            "original_gallery": str(gallery_path),
            "original_gallery_revision": gallery_revision,
            "original_gallery_sha256": snapshots["gallery_index"]["sha256"],
            "original_gallery_manifest_sha256": snapshots["annotation_review.json"]["sha256"],
            "first_index": first_index, "last_index": last_index,
            "accepted_frame_ids": accepted, "reason": reason,
            "coordinate_reannotation_required": True, "created_utc": created,
        }
        _write_json(staging / "human_review_selection.json", selection)
        with (staging / "annotations.pending.csv").open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=ANNOTATION_FIELDS)
            writer.writeheader()
            for frame in frames:
                identity = frame["frame_id"]
                selected = identity in accepted_set
                writer.writerow({
                    "frame_id": identity, "x": "", "y": "", "uncertainty_px": "",
                    "visibility": "uncertain" if selected else "excluded",
                    "medium": "underwater" if selected else by_id[identity]["medium"],
                    "note": ("User-selected interval; manual coordinate reannotation required: " if selected else "User exclusion: ") + reason,
                    "annotator": "pending_manual_reannotation" if selected else "user_scope",
                })
            stream.flush()
            os.fsync(stream.fileno())
        cloned = copy.deepcopy(parent)
        # Original paths in previous timing lineage cannot describe this clone.
        # Preserve all such history in the byte-exact saved parent manifest.
        cloned.pop("postprocessing", None)
        cloned.pop("producer_manifest_sha256", None)
        copied_sources = {str(source / name): name for name in PREPARED_FILES}
        for entry in cloned["outputs"]:
            name = copied_sources.get(str(Path(entry["path"]).resolve()))
            if name is not None:
                info = (staging / name).stat()
                entry.update(path=str(destination / name), bytes=info.st_size,
                             sha256=snapshots[name]["sha256"], mtime=info.st_mtime,
                             ctime=info.st_ctime, inode=info.st_ino)
        cloned["events_csv"] = str(destination / "events.csv")
        cloned["frames_csv"] = str(destination / "frames.csv")
        cloned["snippet_cache"]["path"] = str(destination / "event_snippets.f64")
        cloned["postprocessing"] = [{
            "stage": "annotation_revision", "schema": "clicker-annotation-revision-v1",
            "created_utc": created, "source_directory": str(source),
            "source_manifest_sha256": snapshots["prepare_manifest.json"]["sha256"],
            "saved_parent_manifest_path": str(destination / "revision_source/prepare_manifest.json"),
            "saved_parent_postprocessing_locator": "revision_source/prepare_manifest.json#/postprocessing",
            "previous_annotations_sha256": snapshots["annotations.csv"]["sha256"],
            "previous_annotation_completion_sha256": snapshots["annotations.complete.json"]["sha256"],
            "review_source_sha256": _sha(Path(__file__)),
            "human_review_selection_path": str(destination / "human_review_selection.json"),
            "human_review_selection_sha256": _sha(staging / "human_review_selection.json"),
            "original_gallery": str(gallery_path), "gallery_revision": gallery_revision,
            "first_index": first_index, "last_index": last_index,
            "coordinate_reannotation_required": True, "reason": reason,
        }]
        _write_json(staging / "prepare_manifest.json", cloned)
        # Recheck all authenticated source inputs before publishing. Source edits
        # or resealing during construction abort rather than silently rebase.
        for state in snapshots.values():
            path = Path(state["path"])
            if _fingerprint(path) != tuple(state["fingerprint"]) or _sha(path) != state["sha256"]:
                raise ValueError(f"Source changed during annotation revision: {path}")
        _directory_fsync(archive)
        _directory_fsync(staging)
        _publish_absent(staging, destination)
    finally:
        if staging.is_dir():
            shutil.rmtree(staging)
    return {
        "output_dir": str(destination), "source_dir": str(source),
        "manifest": str(destination / "prepare_manifest.json"),
        "selection": str(destination / "human_review_selection.json"),
        "pending_annotations": str(destination / "annotations.pending.csv"),
        "accepted_frame_count": len(accepted), "excluded_frame_count": len(frames) - len(accepted),
        "source_manifest_sha256": snapshots["prepare_manifest.json"]["sha256"],
        "coordinate_reannotation_required": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("destination_dir", type=Path)
    parser.add_argument("--first-index", required=True, type=int)
    parser.add_argument("--last-index", required=True, type=int)
    parser.add_argument("--gallery-revision", required=True)
    parser.add_argument("--reason", required=True)
    arguments = parser.parse_args()
    print(json.dumps(create_annotation_revision(**vars(arguments))))


if __name__ == "__main__":
    main()
