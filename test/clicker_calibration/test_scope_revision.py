"""Stdlib synthetic scope controls; no real media, annotations, or calibration."""
import csv
import hashlib
import html
import importlib.util
import io
import json
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest import mock
import zlib


spec = importlib.util.spec_from_file_location(
    "scope_revision", Path(__file__).parents[2] / "src/clicker_calibration/scope_revision.py")
scope = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scope)
h = scope._helpers()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def image_bytes(width=320, height=180):
    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data))
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress((b"\0" + bytes(width * 3)) * height)) + chunk(b"IEND", b""))


def snapshot(root):
    return {str(path.relative_to(root)): (h._sha(path), path.stat().st_mtime_ns, path.stat().st_ino)
            for path in root.rglob("*") if path.is_file()}


class ScopeRevisionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.source.mkdir()
        self.destination = self.root / "revision"
        images = self.source / "frames"
        images.mkdir()
        self.frames, self.annotations, events = [], [], []
        for i in range(6):
            image = images / f"B-f{i}.png"
            image.write_bytes(image_bytes())
            identity = f"B-f{i}"
            self.frames.append(dict(frame_id=identity, video_file="synthetic.mp4", frame_index=i,
                                    frame_pts_s=i / 60, video_global_time_s=i / 60,
                                    image_path=str(image), width=320, height=180))
            self.annotations.append(dict(frame_id=identity, x="" if i == 3 else "-0.0" if i == 1 else "+001.2300e2",
                                         y="" if i == 3 else "80.2500", visibility="visible" if i == 2 else "uncertain",
                                         medium="air" if i == 0 else "unknown" if i == 5 else "underwater",
                                         uncertainty_px="" if i == 3 else "3.500e1",
                                         note='Original Terra: "blurred", retain\nsecond line' if i == 3 else "original visual note",
                                         annotator="gpt-6.1-terra", extra="UTF-8 α; keep\vfield"))
            events.append(dict(event_id=f"B-e{i}", frame_id=identity, image_path=str(image), fs=192000,
                               gcc_d1=1.25, gcc_d2=-0.5, gcc_d3=0.0, sync_status="mapped"))
        table(self.source / "frames.csv", self.frames)
        table(self.source / "events.csv", events)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=list(self.annotations[0]), quoting=csv.QUOTE_ALL, lineterminator="\r\n")
        writer.writeheader()
        # Both annotation order and HTML index order differ from frames.csv.
        writer.writerows(self.annotations[i] for i in (5, 2, 0, 4, 1, 3))
        (self.source / "annotations.csv").write_bytes(b"\xef\xbb\xbf" + buffer.getvalue().removesuffix("\r\n").encode())
        dump(self.source / "synchronization.json", {"status": "valid", "scale": 1.00006})
        dump(self.source / "prepare_metadata.json", {"videos": [{"path": "synthetic.mp4", "pts": [0, 0.01]}]})
        (self.source / "event_snippets.f64").write_bytes(bytes(6 * 48))
        (self.source / "calibration.json").write_text("old fit must not be copied")
        (self.source / "rejected_candidates.csv").write_text("event_id\n")
        self.seal = dict(schema_version="clicker-annotations-v1", frames_sha256=h._sha(self.source / "frames.csv"),
                         annotations_sha256=h._sha(self.source / "annotations.csv"), n_frames=6,
                         landmark_convention="source housing center", coordinate_convention="zero-origin original pixels; x right, y down")
        dump(self.source / "annotations.complete.json", self.seal)
        self.gallery_revision = self.seal["annotations_sha256"][:12]
        gallery = self.source / "annotation_review" / self.gallery_revision
        gallery.mkdir(parents=True)
        self.order = [4, 0, 3, 1, 5, 2]
        cards = []
        for index, number in enumerate(self.order):
            (gallery / f"{index:05d}.jpg").write_bytes(b"synthetic gallery file identity")
            identity = f"B-f{number}"
            cards.append(f'<article><a href="{index:05d}.jpg"><img src="{index:05d}_thumb.jpg" alt="{html.escape(identity)}"></a></article>')
        self.gallery_path = gallery / "index.html"
        self.gallery_path.write_text("<main>" + "".join(cards) + "</main>")
        dump(self.source / "annotation_review.json", {**self.seal, "revision": self.gallery_revision, "gallery": str(self.gallery_path)})
        outputs = [dict(path=str(path), sha256=h._sha(path), bytes=path.stat().st_size,
                        mtime=path.stat().st_mtime, inode=path.stat().st_ino)
                   for path in [*(self.source / n for n in h.PREPARED_FILES), self.source / "rejected_candidates.csv", *images.iterdir()]]
        self.prepare = dict(status="complete", event_count=6, frame_count=6, mapped_count=6,
                            events_csv=str(self.source / "events.csv"), frames_csv=str(self.source / "frames.csv"),
                            outputs=outputs, request={"source_sha256": "a" * 64}, completed_utc="2026-01-01T00:00:00Z",
                            snippet_cache=dict(path=str(self.source / "event_snippets.f64"), sha256=h._sha(self.source / "event_snippets.f64"),
                                               bytes=6 * 48, dtype="Float64", endianness="little", layout="column-major"),
                            postprocessing=[{"stage": "sync_eligibility", "prior": "retained only in archive"}], producer_manifest_sha256="c" * 64)
        dump(self.source / "prepare_manifest.json", self.prepare)

    def create(self, **changes):
        kwargs = dict(first_index=2, last_index=3, gallery_revision=self.gallery_revision,
                      reason="User confirmed only numbered gallery images 00002–00003 are valid underwater clicker")
        kwargs.update(changes)
        return scope.create_scope_revision(self.source, self.destination, **kwargs)

    def test_reordered_html_inclusive_selection_preserves_exact_original_records(self):
        before = snapshot(self.source)
        result = self.create()
        self.assertEqual(result["accepted_frame_count"], 2)
        self.assertEqual(result["selected_coordinate_count"], 1)
        self.assertEqual(result["selected_missing_coordinate_count"], 1)
        self.assertEqual(result["excluded_frame_count"], 4)
        selected = h._read_json(self.destination / "human_review_selection.json")
        self.assertEqual(selected["accepted_frame_ids"], ["B-f3", "B-f1"])
        self.assertFalse(selected["coordinate_reannotation_required"])
        self.assertTrue(selected["publish_selected_only"])
        fields, header, original = scope._annotation_records(self.source / "annotations.csv")
        new_fields, new_header, derived = scope._annotation_records(self.destination / "annotations.csv")
        self.assertEqual(fields, new_fields)
        self.assertEqual(header, new_header)
        self.assertEqual([r["frame_id"] for r, _ in original], [r["frame_id"] for r, _ in derived])
        for (old, old_bytes), (new, new_bytes) in zip(original, derived):
            if old["frame_id"] in selected["accepted_frame_ids"]:
                self.assertEqual(new, old)
                self.assertEqual(new_bytes, old_bytes)
                self.assertEqual(new["visibility"], "uncertain")
            else:
                self.assertEqual((new["x"], new["y"], new["uncertainty_px"], new["visibility"], new["annotator"]),
                                 ("", "", "", "excluded", "user_scope"))
                self.assertEqual(new["note"], old["note"])
                self.assertEqual(new["medium"], old["medium"])
                self.assertEqual(new["extra"], old["extra"])
        self.assertEqual(snapshot(self.source), before)
        self.assertFalse((self.destination / "annotations.pending.csv").exists())
        self.assertFalse((self.destination / "calibration.json").exists())

    def test_five_prepared_files_archives_and_only_owned_paths_relocated(self):
        self.create()
        for name in h.PREPARED_FILES:
            self.assertEqual((self.source / name).read_bytes(), (self.destination / name).read_bytes())
        for name in h.ARCHIVED_FILES:
            self.assertEqual((self.source / name).read_bytes(), (self.destination / "revision_source" / name).read_bytes())
        prepared = h._read_json(self.destination / "prepare_manifest.json")
        self.assertEqual(prepared["request"], self.prepare["request"])
        self.assertEqual(prepared["events_csv"], str(self.destination / "events.csv"))
        self.assertEqual(prepared["frames_csv"], str(self.destination / "frames.csv"))
        self.assertEqual(prepared["snippet_cache"]["path"], str(self.destination / "event_snippets.f64"))
        old_outputs = {r["path"]: r for r in self.prepare["outputs"]}
        for entry in prepared["outputs"]:
            path = Path(entry["path"])
            if path.parent == self.destination:
                self.assertIn(path.name, h.PREPARED_FILES)
                self.assertEqual(entry["sha256"], h._sha(path))
            else:
                self.assertEqual(entry, old_outputs[entry["path"]])
        self.assertNotIn(".scope-", json.dumps(prepared))
        lineage = prepared["postprocessing"][0]
        self.assertEqual(lineage["stage"], "annotation_revision")
        self.assertEqual(lineage["review_source_sha256"], h._sha(Path(h.__file__)))
        self.assertEqual(lineage["scope_revision_source_sha256"], h._sha(Path(scope.__file__)))
        self.assertFalse(lineage["coordinate_reannotation_required"])
        self.assertEqual(lineage["saved_parent_manifest_path"], str(self.destination / "revision_source/prepare_manifest.json"))

    def test_seal_and_scope_lineage_hashes_are_consistent(self):
        result = self.create()
        seal = h._read_json(self.destination / "annotations.complete.json")
        audit = h._read_json(self.destination / "annotation_scope_revision.json")
        selection = h._read_json(self.destination / "human_review_selection.json")
        lineage = h._read_json(self.destination / "prepare_manifest.json")["postprocessing"][0]
        for key, name in (("frames_sha256", "frames.csv"), ("annotations_sha256", "annotations.csv"),
                          ("human_review_selection_sha256", "human_review_selection.json"),
                          ("annotation_scope_revision_sha256", "annotation_scope_revision.json"),
                          ("prepare_manifest_sha256", "prepare_manifest.json")):
            self.assertEqual(seal[key], h._sha(self.destination / name))
        self.assertEqual(seal["n_frames"], 6)
        self.assertEqual(seal["human_selected_frame_count"], 2)
        self.assertEqual(seal["counts"], {"absent": 0, "excluded": 4, "occluded": 0, "uncertain": 2, "visible": 0})
        self.assertEqual(seal["annotation_authority"], scope.AUTHORITY)
        self.assertEqual(audit["original_annotations_sha256"], self.seal["annotations_sha256"])
        self.assertEqual(audit["original_completion_sha256"], h._sha(self.source / "annotations.complete.json"))
        self.assertEqual(audit["human_review_selection_sha256"], lineage["human_review_selection_sha256"])
        self.assertEqual(lineage["scope_revision_manifest_sha256"], result["annotation_scope_revision_sha256"])
        self.assertEqual(audit["original_coordinate_source"], selection["original_coordinate_source"])
        self.assertFalse(audit["new_annotation"])

    def test_first_and_last_single_filename_boundaries_are_inclusive(self):
        result = self.create(first_index=0, last_index=0)
        self.assertEqual(result["accepted_frame_count"], 1)
        self.assertEqual(h._read_json(self.destination / "human_review_selection.json")["accepted_frame_ids"], ["B-f4"])
        other = self.root / "last"
        scope.create_scope_revision(self.source, other, first_index=5, last_index=5,
                                    gallery_revision=self.gallery_revision, reason="last numbered image only")
        self.assertEqual(h._read_json(other / "human_review_selection.json")["accepted_frame_ids"], ["B-f2"])
        self.assertEqual({r["frame_id"]: r for r in h._read_csv(other / "annotations.csv")}["B-f2"]["visibility"], "visible")

    def test_selected_nonunderwater_medium_fails_without_relabeling(self):
        before = snapshot(self.source)
        with self.assertRaisesRegex(ValueError, "not already labeled underwater"):
            self.create(first_index=1, last_index=1)
        self.assertFalse(self.destination.exists())
        self.assertEqual(snapshot(self.source), before)
        self.assertFalse(list(self.root.glob(".revision.scope-*")))

    def test_existing_destination_empty_directory_or_broken_symlink_rejected(self):
        self.destination.mkdir()
        with self.assertRaises(FileExistsError):
            self.create()
        self.destination.rmdir()
        self.destination.symlink_to(self.root / "absent")
        with self.assertRaises(FileExistsError):
            self.create()
        self.assertTrue(self.destination.is_symlink())

    def test_stale_annotations_gallery_and_prepared_bytes_rejected(self):
        path = self.source / "annotations.csv"
        old = path.read_bytes()
        path.write_bytes(old + b"\n")
        with self.assertRaisesRegex(ValueError, "completion hashes"):
            self.create()
        path.write_bytes(old)
        with self.assertRaisesRegex(ValueError, "gallery revision"):
            self.create(gallery_revision="0" * 12)
        manifest = h._read_json(self.source / "annotation_review.json")
        manifest["frames_sha256"] = "bad"
        dump(self.source / "annotation_review.json", manifest)
        with self.assertRaisesRegex(ValueError, "Gallery is stale"):
            self.create()
        dump(self.source / "annotation_review.json", {**self.seal, "revision": self.gallery_revision, "gallery": str(self.gallery_path)})
        (self.source / "event_snippets.f64").write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "Prepared source hash"):
            self.create()
        self.assertFalse(self.destination.exists())

    def test_malformed_gallery_mapping_rejected(self):
        self.gallery_path.write_text(self.gallery_path.read_text().replace('alt="B-f3"', 'alt="B-f1"'))
        with self.assertRaisesRegex(ValueError, "Gallery frame identities"):
            self.create()

    def test_source_mutation_after_inner_clone_aborts_publication(self):
        original = scope._replace_owned_json
        def mutate(staging, name, value):
            original(staging, name, value)
            if name == "human_review_selection.json":
                with (self.source / "annotations.csv").open("ab") as stream:
                    stream.write(b"\n")
        with mock.patch.object(scope, "_replace_owned_json", side_effect=mutate):
            with self.assertRaisesRegex(ValueError, "Source changed"):
                self.create()
        self.assertFalse(self.destination.exists())
        self.assertFalse(list(self.root.glob(".revision.scope-*")))

    def test_publication_is_sealed_last_and_destination_race_not_overwritten(self):
        publish = h._publish_absent
        def race(staging, destination):
            if destination == self.destination:
                self.assertFalse(destination.exists())
                seal = h._read_json(staging / "annotations.complete.json")
                self.assertEqual(seal["prepare_manifest_sha256"], h._sha(staging / "prepare_manifest.json"))
                destination.mkdir()
            publish(staging, destination)
        with mock.patch.object(h, "_publish_absent", side_effect=race):
            with self.assertRaises(FileExistsError):
                self.create()
        self.assertEqual(list(self.destination.iterdir()), [])
        self.assertFalse(list(self.root.glob(".revision.scope-*")))

    def test_failure_before_seal_leaves_no_partial_public_destination(self):
        write = h._write_json
        def fail(path, value):
            if path.name == "annotations.complete.json":
                raise OSError("intentional seal failure")
            write(path, value)
        with mock.patch.object(h, "_write_json", side_effect=fail):
            with self.assertRaisesRegex(OSError, "intentional seal failure"):
                self.create()
        self.assertFalse(self.destination.exists())
        self.assertFalse(list(self.root.glob(".revision.scope-*")))

    def test_invalid_ranges_and_source_child_destination_rejected(self):
        for changes in (dict(first_index=True), dict(first_index=-1), dict(last_index=1), dict(last_index=6), dict(reason=" ")):
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    self.create(**changes)
                self.assertFalse(self.destination.exists())
        with self.assertRaisesRegex(ValueError, "immutable source"):
            scope.create_scope_revision(self.source, self.source / "new", first_index=2, last_index=3,
                                        gallery_revision=self.gallery_revision, reason="user interval")


if __name__ == "__main__":
    unittest.main(verbosity=2)
