import copy
import csv
import hashlib
import html
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock


spec = importlib.util.spec_from_file_location(
    "annotation_revision", Path(__file__).parents[2] / "src/clicker_calibration/annotation_revision.py"
)
revision = importlib.util.module_from_spec(spec)
spec.loader.exec_module(revision)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def tree_state(directory):
    return {str(p.relative_to(directory)): (revision._sha(p), p.stat().st_mtime_ns, p.stat().st_ino)
            for p in directory.rglob("*") if p.is_file()}


class AnnotationRevisionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.source.mkdir()
        image_dir = self.source / "frames"
        image_dir.mkdir()
        self.frames, self.annotations, self.events = [], [], []
        media = ["air", "underwater", "mixed", "unknown", "air", "underwater"]
        for i in range(6):
            image = image_dir / f"raw-{i}.png"
            image.write_bytes(b"immutable original image " + str(i).encode())
            self.frames.append(dict(frame_id=f"A-f{i}", video_file="original.MP4", frame_index=i,
                                    frame_pts_s=i / 60, video_global_time_s=i / 60,
                                    image_path=str(image), width=2704, height=1520))
            self.annotations.append(dict(frame_id=f"A-f{i}", x=100 + i, y=200 + i,
                                         visibility="visible", medium=media[i], uncertainty_px=4,
                                         note="Previous estimate, must not carry coordinates", annotator="old_annotator"))
            self.events.append(dict(event_id=f"A-e{i}", frame_id=f"A-f{i}", image_path=str(image),
                                    sample_index=i * 19200, fs=192000, minmax2_d1=-1.25,
                                    minmax2_d2=0.5, minmax2_d3=0.0,
                                    snippet_file="event_snippets.f64", snippet_byte_offset=i * 48,
                                    snippet_samples=2, snippet_channels=3))
        write_csv(self.source / "frames.csv", self.frames)
        write_csv(self.source / "events.csv", self.events)
        write_csv(self.source / "annotations.csv", self.annotations)
        write_json(self.source / "synchronization.json", {"status": "valid", "scale": 1.00006})
        write_json(self.source / "prepare_metadata.json", {"videos": [{"path": "original.MP4", "pts": [0, 0.01]}]})
        (self.source / "event_snippets.f64").write_bytes(bytes(6 * 48))
        (self.source / "rejected_candidates.csv").write_text("event_id\n")
        (self.source / "calibration.json").write_text("previous fit, do not clone\n")
        (self.source / "predictions.csv").write_text("old fitted prediction\n")
        self.seal = dict(schema_version="clicker-annotations-v1",
                         frames_sha256=revision._sha(self.source / "frames.csv"),
                         annotations_sha256=revision._sha(self.source / "annotations.csv"), n_frames=6)
        write_json(self.source / "annotations.complete.json", self.seal)
        self.gallery_revision = self.seal["annotations_sha256"][:12]
        gallery_dir = self.source / "annotation_review" / self.gallery_revision
        gallery_dir.mkdir(parents=True)
        # Deliberately differ from frames.csv order: selection must follow JPG
        # filenames and the HTML frame IDs, never a guessed CSV slice.
        self.gallery_order = [3, 0, 5, 1, 4, 2]
        cards = []
        for index, frame_index in enumerate(self.gallery_order):
            (gallery_dir / f"{index:05d}.jpg").write_bytes(b"review JPEG")
            identity = self.frames[frame_index]["frame_id"]
            cards.append(f'<article><a href="{index:05d}.jpg"><img src="{index:05d}_thumb.jpg" alt="{html.escape(identity)}"></a><b>{identity}</b></article>')
        self.gallery_path = gallery_dir / "index.html"
        self.gallery_path.write_text("<!doctype html><main>" + "".join(cards) + "</main>")
        write_json(self.source / "annotation_review.json", {
            **self.seal, "revision": self.gallery_revision, "gallery": str(self.gallery_path),
        })
        outputs = []
        for path in [*(self.source / name for name in revision.PREPARED_FILES),
                     self.source / "rejected_candidates.csv", *image_dir.iterdir()]:
            outputs.append(dict(path=str(path), sha256=revision._sha(path), bytes=path.stat().st_size,
                                mtime=path.stat().st_mtime, inode=path.stat().st_ino))
        self.manifest = dict(
            status="complete", event_count=6, frame_count=6, mapped_count=5,
            events_csv=str(self.source / "events.csv"), frames_csv=str(self.source / "frames.csv"),
            request={"schema": "clicker-prepare-v1", "source_sha256": "a" * 64,
                     "inputs": [{"path": "original.MP4", "sha256": "b" * 64}]},
            snippet_cache=dict(path=str(self.source / "event_snippets.f64"),
                               sha256=revision._sha(self.source / "event_snippets.f64"), bytes=6 * 48,
                               dtype="Float64", endianness="little", layout="column-major"),
            outputs=outputs, completed_utc="2026-01-01T00:00:00Z",
            postprocessing=[dict(stage="sync_eligibility", raw_events_path=str(self.source / "preparation_raw/events.csv"),
                                 raw_manifest_sha256="c" * 64)], producer_manifest_sha256="c" * 64,
        )
        write_json(self.source / "prepare_manifest.json", self.manifest)
        self.destination = self.root / "revision"

    def create(self, **overrides):
        arguments = dict(first_index=1, last_index=3, gallery_revision=self.gallery_revision,
                         reason="User confirmed only this clicker interval")
        arguments.update(overrides)
        return revision.create_annotation_revision(self.source, self.destination, **arguments)

    def test_inclusive_gallery_filename_mapping_and_empty_reannotation(self):
        before = tree_state(self.source)
        result = self.create()
        self.assertEqual(result["accepted_frame_count"], 3)
        self.assertEqual(result["excluded_frame_count"], 3)
        selected = json.loads((self.destination / "human_review_selection.json").read_text())
        self.assertEqual(selected["schema_version"], "clicker-human-review-selection-v1")
        self.assertEqual(selected["accepted_frame_ids"], ["A-f0", "A-f5", "A-f1"])
        self.assertEqual((selected["first_index"], selected["last_index"]), (1, 3))
        self.assertEqual(selected["frames_sha256"], self.seal["frames_sha256"])
        self.assertEqual(selected["previous_annotations_sha256"], self.seal["annotations_sha256"])
        self.assertEqual(selected["original_gallery"], str(self.gallery_path))
        self.assertTrue(selected["coordinate_reannotation_required"])
        pending = revision._read_csv(self.destination / "annotations.pending.csv")
        self.assertEqual([row["frame_id"] for row in pending], [row["frame_id"] for row in self.frames])
        old = {row["frame_id"]: row for row in self.annotations}
        for row in pending:
            self.assertEqual((row["x"], row["y"], row["uncertainty_px"]), ("", "", ""))
            if row["frame_id"] in selected["accepted_frame_ids"]:
                self.assertEqual((row["visibility"], row["medium"], row["annotator"]),
                                 ("uncertain", "underwater", "pending_manual_reannotation"))
            else:
                self.assertEqual(row["visibility"], "excluded")
                self.assertEqual(row["medium"], old[row["frame_id"]]["medium"])
                self.assertEqual(row["annotator"], "user_scope")
                self.assertIn("User exclusion:", row["note"])
        self.assertFalse((self.destination / "annotations.csv").exists())
        self.assertFalse((self.destination / "annotations.complete.json").exists())
        self.assertEqual(tree_state(self.source), before)

    def test_copy_hash_manifest_lineage_and_original_image_pointers(self):
        original_manifest = (self.source / "prepare_manifest.json").read_bytes()
        self.create()
        for name in revision.PREPARED_FILES:
            self.assertEqual((self.source / name).read_bytes(), (self.destination / name).read_bytes())
        expected_top = set(revision.PREPARED_FILES) | {"prepare_manifest.json", "revision_source", "human_review_selection.json", "annotations.pending.csv"}
        self.assertEqual({p.name for p in self.destination.iterdir()}, expected_top)
        self.assertEqual({p.name for p in (self.destination / "revision_source").iterdir()}, set(revision.ARCHIVED_FILES))
        for name in revision.ARCHIVED_FILES:
            self.assertEqual((self.source / name).read_bytes(), (self.destination / "revision_source" / name).read_bytes())
        clone = json.loads((self.destination / "prepare_manifest.json").read_text())
        self.assertEqual(clone["request"], self.manifest["request"])
        self.assertEqual(clone["mapped_count"], self.manifest["mapped_count"])
        self.assertEqual(clone["completed_utc"], self.manifest["completed_utc"])
        self.assertEqual(clone["events_csv"], str(self.destination / "events.csv"))
        self.assertEqual(clone["frames_csv"], str(self.destination / "frames.csv"))
        self.assertEqual(clone["snippet_cache"], {**self.manifest["snippet_cache"], "path": str(self.destination / "event_snippets.f64")})
        self.assertNotIn("producer_manifest_sha256", clone)
        self.assertEqual(len(clone["postprocessing"]), 1)
        lineage = clone["postprocessing"][0]
        self.assertEqual(lineage["stage"], "annotation_revision")
        self.assertEqual(lineage["source_manifest_sha256"], hashlib.sha256(original_manifest).hexdigest())
        self.assertEqual(lineage["saved_parent_manifest_path"], str(self.destination / "revision_source/prepare_manifest.json"))
        saved = json.loads((self.destination / "revision_source/prepare_manifest.json").read_text())
        self.assertEqual(saved["postprocessing"], self.manifest["postprocessing"])
        self.assertEqual(saved["producer_manifest_sha256"], self.manifest["producer_manifest_sha256"])
        by_original = {entry["path"]: entry for entry in self.manifest["outputs"]}
        for entry in clone["outputs"]:
            name = Path(entry["path"]).name
            if entry["path"].startswith(str(self.destination) + os.sep):
                self.assertIn(name, revision.PREPARED_FILES)
                self.assertEqual(entry["sha256"], revision._sha(self.destination / name))
                self.assertEqual(entry["bytes"], (self.destination / name).stat().st_size)
            else:
                self.assertEqual(entry, by_original[entry["path"]])
        cloned_frames = revision._read_csv(self.destination / "frames.csv")
        self.assertEqual([f["image_path"] for f in cloned_frames], [f["image_path"] for f in self.frames])
        self.assertTrue(all(Path(f["image_path"]).is_absolute() for f in cloned_frames))
        self.assertFalse((self.destination / "frames").exists())
        self.assertFalse((self.destination / "calibration.json").exists())
        self.assertFalse((self.destination / "rejected_candidates.csv").exists())
        self.assertFalse((self.destination / "sync_eligibility.json").exists())

    def test_single_index_inclusive_and_both_endpoints(self):
        self.create(first_index=5, last_index=5)
        selection = json.loads((self.destination / "human_review_selection.json").read_text())
        self.assertEqual(selection["accepted_frame_ids"], ["A-f2"])

    def test_existing_destination_refused_even_empty_or_broken_symlink(self):
        self.destination.mkdir()
        before = tree_state(self.source)
        with self.assertRaises(FileExistsError):
            self.create()
        self.assertTrue(self.destination.is_dir())
        self.assertEqual(tree_state(self.source), before)
        self.destination.rmdir()
        self.destination.symlink_to(self.root / "missing")
        with self.assertRaises(FileExistsError):
            self.create()
        self.assertTrue(self.destination.is_symlink())

    def test_raced_existing_directory_is_not_overwritten(self):
        original_publish = revision._publish_absent
        def raced(staging, destination):
            destination.mkdir()
            original_publish(staging, destination)
        with mock.patch.object(revision, "_publish_absent", side_effect=raced):
            with self.assertRaises(FileExistsError):
                self.create()
        self.assertTrue(self.destination.is_dir())
        self.assertEqual(list(self.destination.iterdir()), [])
        self.assertFalse(list(self.root.glob(".revision.staging-*")))

    def test_reject_source_subdirectory_destination(self):
        with self.assertRaisesRegex(ValueError, "immutable source"):
            revision.create_annotation_revision(self.source, self.source / "new", first_index=0, last_index=0,
                                                gallery_revision=self.gallery_revision, reason="reason")
        self.assertFalse((self.source / "new").exists())

    def test_stale_annotation_completion_rejected(self):
        with (self.source / "annotations.csv").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(ValueError, "completion hashes"):
            self.create()
        self.assertFalse(self.destination.exists())

    def test_wrong_gallery_revision_and_stale_gallery_manifest_rejected(self):
        with self.assertRaisesRegex(ValueError, "gallery.*revision"):
            self.create(gallery_revision="0" * 12)
        gallery = json.loads((self.source / "annotation_review.json").read_text())
        gallery["frames_sha256"] = "0" * 64
        write_json(self.source / "annotation_review.json", gallery)
        with self.assertRaisesRegex(ValueError, "Gallery is stale"):
            self.create()
        self.assertFalse(self.destination.exists())

    def test_prepared_cache_corruption_rejected(self):
        (self.source / "event_snippets.f64").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "Prepared source hash"):
            self.create()
        self.assertFalse(self.destination.exists())

    def test_gallery_mapping_duplicates_rejected(self):
        text = self.gallery_path.read_text()
        text = text.replace('alt="A-f0"', 'alt="A-f3"')
        self.gallery_path.write_text(text)
        with self.assertRaisesRegex(ValueError, "Gallery frame identities"):
            self.create()

    def test_ranges_and_reason_validation(self):
        for changes in [dict(first_index=-1), dict(first_index=4, last_index=3), dict(last_index=6),
                        dict(first_index=True), dict(reason=" ")]:
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    self.create(**changes)
                self.assertFalse(self.destination.exists())

    def test_source_mutation_during_copy_aborts_without_partial_destination(self):
        original_copy = revision._copy_exact
        changed = False
        def changing(source, destination, digest):
            nonlocal changed
            original_copy(source, destination, digest)
            if not changed:
                changed = True
                with (self.source / "annotations.csv").open("a") as stream:
                    stream.write("\n")
        with mock.patch.object(revision, "_copy_exact", side_effect=changing):
            with self.assertRaisesRegex(ValueError, "Copied bytes differ|Source changed"):
                self.create()
        self.assertFalse(self.destination.exists())
        self.assertFalse(list(self.root.glob(".revision.staging-*")))


if __name__ == "__main__":
    unittest.main()
