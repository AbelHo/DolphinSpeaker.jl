"""Synthetic stdlib tests; never opens actual calibration recordings/directories."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


spec = importlib.util.spec_from_file_location(
    "fit_annotation_policy", Path(__file__).parents[2] / "src/clicker_calibration/fit_annotation_policy.py"
)
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def csv_write(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def snapshot(directory):
    return {str(p.relative_to(directory)): (policy._sha(p), p.stat().st_mtime_ns, p.stat().st_ino)
            for p in directory.rglob("*") if p.is_file()}


def make_source(root, *, admitted=3, missing=1, excluded=2):
    """Return a complete synthetic sealed human source for bridge integration."""
    source = Path(root) / "human"
    source.mkdir()
    image_dir = source / "frames"
    image_dir.mkdir()
    n = admitted + missing + excluded
    frames, annotations, events = [], [], []
    fields = list(policy._helpers().ANNOTATION_FIELDS)
    for i in range(n):
        image = image_dir / f"{i}.png"
        image.write_bytes(b"frozen original source image " + str(i).encode())
        identity = f"A-frame{i}"
        frames.append(dict(frame_id=identity, video_file="source.mp4", frame_index=i,
                           frame_pts_s=i / 60, video_global_time_s=i / 60,
                           image_path=str(image), width=320, height=180))
        selected = i < admitted + missing
        coordinates = i < admitted
        annotations.append(dict(frame_id=identity, x=("-0.0" if i == 0 else "1.2300e2") if coordinates else "",
                                y="80.2500" if coordinates else "", visibility="uncertain" if selected else "excluded",
                                medium="underwater" if selected else "air", uncertainty_px="3.500e1" if coordinates else "",
                                note='Original Terra estimate; "uncertain", not new truth', annotator="gpt-6.1-terra"))
        events.append(dict(event_id=f"event{i}", frame_id=identity, image_path=str(image),
                           fs=192000, gcc_d1=1.25, gcc_d2=-0.375, gcc_d3=0.0, sync_status="mapped"))
    csv_write(source / "frames.csv", list(frames[0]), frames)
    csv_write(source / "events.csv", list(events[0]), events)
    csv_write(source / "annotations.csv", fields, annotations)
    dump(source / "synchronization.json", {"status": "valid", "scale": 1.00006})
    dump(source / "prepare_metadata.json", {"videos": [{"pts": [0.0, 0.01]}]})
    (source / "event_snippets.f64").write_bytes(bytes(48 * n))
    selected_ids = [r["frame_id"] for r in annotations[:admitted + missing]]
    dump(source / "human_review_selection.json", dict(
        schema_version="clicker-human-review-selection-v1", frames_sha256=policy._sha(source / "frames.csv"),
        accepted_frame_ids=selected_ids, previous_annotations_sha256="c" * 64,
        original_gallery="frozen-original-review/index.html", first_index=0, last_index=len(selected_ids) - 1,
        coordinate_reannotation_required=True, reason="User confirmed selected clicker sequence"))
    outputs = []
    for path in [*(source / name for name in policy._helpers().PREPARED_FILES), *image_dir.iterdir()]:
        outputs.append(dict(path=str(path), sha256=policy._sha(path), bytes=path.stat().st_size))
    dump(source / "prepare_manifest.json", dict(
        status="complete", event_count=n, frame_count=n, mapped_count=n,
        request=dict(source_sha256="a" * 64, inputs=[dict(path="source.mp4", sha256="b" * 64)]),
        events_csv=str(source / "events.csv"), frames_csv=str(source / "frames.csv"), outputs=outputs,
        snippet_cache=dict(path=str(source / "event_snippets.f64"), sha256=policy._sha(source / "event_snippets.f64"),
                           bytes=48 * n, dtype="Float64", endianness="little", layout="column-major"),
        completed_utc="2026-01-01T00:00:00Z", postprocessing=[dict(stage="annotation_revision", source_manifest_sha256="d" * 64)]))
    reseal(source)
    return source


def reseal(source):
    _, rows = policy._annotations(source / "annotations.csv")
    selection = policy._json(source / "human_review_selection.json")
    dump(source / "annotations.complete.json", dict(
        schema_version="clicker-annotations-v1", frames_sha256=policy._sha(source / "frames.csv"),
        annotations_sha256=policy._sha(source / "annotations.csv"), n_frames=len(rows),
        counts={v: sum(r["visibility"] == v for r in rows) for v in policy.VISIBILITIES},
        annotator=sorted({r["annotator"] for r in rows}), created_utc="2026-01-01T00:00:01Z",
        human_review_selection_sha256=policy._sha(source / "human_review_selection.json"),
        human_selected_frame_count=len(selection["accepted_frame_ids"]),
        landmark_convention="estimated clicker housing center, original uncertainty retained",
        coordinate_convention="zero-origin original pixels; x right, y down"))


class FitAnnotationPolicyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.original = make_source(self.root)
        self.fit = self.root / "fit"

    def create(self, **changes):
        kwargs = dict(policy_name=policy.POLICY_NAME, reason="User explicitly authorized fitting selected uncertain coordinates while checking")
        kwargs.update(changes)
        return policy.create_fit_annotation_input(self.original, self.fit, **kwargs)

    def edit_original(self, callback):
        fields, rows = policy._annotations(self.original / "annotations.csv")
        callback(rows)
        csv_write(self.original / "annotations.csv", fields, rows)
        reseal(self.original)

    def test_real_scale_counts_233_admitted_10_missing_and_no_source_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            original = make_source(directory, admitted=233, missing=10, excluded=548)
            before = snapshot(original)
            fit = Path(directory) / "fit"
            result = policy.create_fit_annotation_input(original, fit, policy_name=policy.POLICY_NAME, reason="Synthetic user authorization")
            self.assertTrue(result["valid"])
            self.assertEqual(result["admitted_count"], 233)
            self.assertEqual(result["missing_coordinate_count"], 10)
            self.assertEqual(result["original_visibility_counts"]["uncertain"], 243)
            self.assertEqual(result["fit_visibility_counts"]["visible"], 233)
            self.assertEqual(result["fit_visibility_counts"]["uncertain"], 10)
            self.assertEqual(result["fit_visibility_counts"]["excluded"], 548)
            self.assertEqual(snapshot(original), before)

    def test_coordinates_medium_radius_and_all_truth_strings_preserved(self):
        before = snapshot(self.original)
        result = self.create()
        old_fields, original = policy._annotations(self.original / "annotations.csv")
        fields, derived = policy._annotations(self.fit / "annotations.csv")
        self.assertEqual(fields, old_fields + list(policy.EXTRA_FIELDS))
        admitted = set(result["admitted_frame_ids"])
        for old, new in zip(original, derived):
            for field in old_fields:
                self.assertEqual(new[field], "visible" if field == "visibility" and old["frame_id"] in admitted else old[field])
            self.assertEqual(new["terra_visibility"], old["visibility"])
            self.assertTrue(new["fit_eligibility_reason"])
        self.assertEqual(derived[0]["x"], "-0.0")
        self.assertEqual(derived[1]["x"], "1.2300e2")
        self.assertEqual(derived[1]["uncertainty_px"], "3.500e1")
        self.assertEqual(snapshot(self.original), before)
        for name in policy.TRUTH_FILES:
            self.assertEqual((self.original / name).read_bytes(), (self.fit / "fit_policy_source" / name).read_bytes())
        self.assertEqual((self.original / "human_review_selection.json").read_bytes(), (self.fit / "human_review_selection.json").read_bytes())

    def test_only_selected_allowed_visibility_underwater_and_valid_radius_enter(self):
        def update(rows):
            rows[0]["visibility"] = "occluded"
            rows[1]["visibility"] = "absent"
            rows[2]["uncertainty_px"] = ""  # uncertain with optional missing radius remains noneligible
            rows[4].update(visibility="uncertain", medium="underwater", x="10", y="11", uncertainty_px="40")
        self.edit_original(update)
        result = self.create()
        self.assertEqual(result["admitted_count"], 0)
        _, derived = policy._annotations(self.fit / "annotations.csv")
        self.assertEqual([r["visibility"] for r in derived], ["occluded", "absent", "uncertain", "uncertain", "uncertain", "excluded"])
        self.assertEqual(derived[4]["fit_eligibility_reason"], "outside_user_confirmed_selection")

    def test_existing_visible_inside_selection_keeps_original_status(self):
        self.edit_original(lambda rows: rows[0].update(visibility="visible"))
        result = self.create()
        self.assertEqual(result["admitted_count"], 3)
        self.assertEqual(len(result["converted_uncertain_frame_ids"]), 2)
        _, derived = policy._annotations(self.fit / "annotations.csv")
        self.assertEqual(derived[0]["visibility"], "visible")
        self.assertEqual(derived[0]["terra_visibility"], "visible")
        self.assertEqual(derived[0]["fit_eligibility_reason"], "admitted_original_visible_coordinates")

    def test_nonunderwater_and_excluded_selected_never_enter(self):
        def update(rows):
            rows[0]["medium"] = "air"
            rows[1]["visibility"] = "excluded"
        self.edit_original(update)
        result = self.create()
        self.assertEqual(result["admitted_frame_ids"], ["A-frame2"])
        _, rows = policy._annotations(self.fit / "annotations.csv")
        self.assertEqual(rows[0]["visibility"], "uncertain")
        self.assertEqual(rows[1]["visibility"], "excluded")

    def test_unsafe_original_visible_outside_selection_fails_closed(self):
        self.edit_original(lambda rows: rows[4].update(visibility="visible", medium="underwater", x="10", y="20", uncertainty_px="40"))
        with self.assertRaisesRegex(ValueError, "Unsafe nonadmitted original-visible"):
            self.create()
        self.assertFalse(self.fit.exists())

    def test_explicit_policy_and_absent_external_destination_required(self):
        with self.assertRaises(ValueError):
            self.create(policy_name="automatic_uncertain")
        with self.assertRaises(ValueError):
            self.create(reason=" ")
        self.fit.mkdir()
        with self.assertRaises(FileExistsError):
            self.create()
        with self.assertRaises(ValueError):
            policy.create_fit_annotation_input(self.original, self.original / "derived", policy_name=policy.POLICY_NAME, reason="authorization")

    def test_stale_or_unauthenticated_original_selection_refused(self):
        selection = policy._json(self.original / "human_review_selection.json")
        selection["accepted_frame_ids"].append("A-frame4")
        dump(self.original / "human_review_selection.json", selection)
        with self.assertRaisesRegex(ValueError, "Stale original annotation completion"):
            self.create()
        self.assertFalse(self.fit.exists())

    def test_invalid_numeric_source_rows_are_not_silently_sealed(self):
        fields, original = policy._annotations(self.original / "annotations.csv")
        for key, value in [("x", "nan"), ("y", "180"), ("uncertainty_px", "-1"), ("uncertainty_px", "inf")]:
            with self.subTest(key=key, value=value):
                rows = [dict(r) for r in original]
                rows[0][key] = value
                csv_write(self.original / "annotations.csv", fields, rows)
                reseal(self.original)
                with self.assertRaises(ValueError):
                    self.create()
                self.assertFalse(self.fit.exists())
        csv_write(self.original / "annotations.csv", fields, original)
        reseal(self.original)

    def test_validator_bridge_hashes_manifest_and_no_images_copied(self):
        result = self.create()
        required = {"original_annotations_sha256", "fit_annotations_sha256", "original_completion_sha256",
                    "fit_completion_sha256", "policy_sha256", "helper_sha256", "frames_sha256", "selection_sha256"}
        self.assertTrue(required <= set(result))
        self.assertNotEqual(result["original_annotations_sha256"], result["fit_annotations_sha256"])
        self.assertEqual(result["helper_sha256"], policy._sha(Path(policy.__file__)))
        self.assertEqual(result["original_annotations_sha256"], policy._sha(self.original / "annotations.csv"))
        self.assertEqual(result["fit_annotations_sha256"], policy._sha(self.fit / "annotations.csv"))
        self.assertEqual(result["policy_sha256"], policy._sha(self.fit / "annotation_fit_policy.json"))
        self.assertEqual(result, policy.validate_fit_annotation_policy(self.original, self.fit))
        for name in policy._helpers().PREPARED_FILES:
            self.assertEqual((self.original / name).read_bytes(), (self.fit / name).read_bytes())
        manifest = policy._json(self.fit / "prepare_manifest.json")
        parent = policy._json(self.original / "prepare_manifest.json")
        self.assertEqual(manifest["request"], parent["request"])
        self.assertEqual(manifest["postprocessing"][0]["stage"], "annotation_fit_policy")
        self.assertEqual(len(manifest["postprocessing"]), 1)
        self.assertEqual(policy._json(self.fit / "fit_policy_source/prepare_manifest.json"), parent)
        self.assertFalse((self.fit / "frames").exists())
        original_images = [r["image_path"] for r in policy._helpers()._read_csv(self.original / "frames.csv")]
        self.assertEqual(original_images, [r["image_path"] for r in policy._helpers()._read_csv(self.fit / "frames.csv")])
        completion = policy._json(self.fit / "annotations.complete.json")
        self.assertEqual(completion["annotation_authority"], policy.AUTHORITY)
        self.assertEqual(completion["human_review_selection_sha256"], result["selection_sha256"])
        self.assertEqual(completion["annotation_fit_policy_sha256"], result["policy_sha256"])

    def test_changed_original_truth_even_resealed_invalidates_export_bridge(self):
        self.create()
        self.edit_original(lambda rows: rows[0].update(x="1.5"))
        with self.assertRaisesRegex(ValueError, "stale against original truth"):
            policy.validate_fit_annotation_policy(self.original, self.fit)

    def test_fit_coordinate_status_radius_or_metadata_tampering_rejected(self):
        self.create()
        path = self.fit / "annotations.csv"
        original_bytes = path.read_bytes()
        fields, rows = policy._annotations(path)
        for key, value in [("x", "1.0"), ("visibility", "uncertain"), ("uncertainty_px", "3"),
                           ("medium", "air"), ("terra_visibility", "visible"), ("fit_eligibility_reason", "manual")]:
            with self.subTest(key=key):
                changed = [dict(r) for r in rows]
                changed[0][key] = value
                csv_write(path, fields, changed)
                with self.assertRaisesRegex(ValueError, "row/coordinate/status equivalence"):
                    policy.validate_fit_annotation_policy(self.original, self.fit)
                path.write_bytes(original_bytes)

    def test_policy_selection_archive_prepared_and_seal_tampering_rejected(self):
        self.create()
        variants = [
            ("annotation_fit_policy.json", lambda p: dump(p, {**policy._json(p), "admitted_frame_ids": []})),
            ("human_review_selection.json", lambda p: p.write_text("{}")),
            ("fit_policy_source/annotations.csv", lambda p: p.write_text("changed")),
            ("event_snippets.f64", lambda p: p.write_bytes(b"changed")),
            ("annotations.complete.json", lambda p: dump(p, {**policy._json(p), "annotation_authority": "terra"})),
            ("prepare_manifest.json", lambda p: dump(p, {**policy._json(p), "mapped_count": 999})),
        ]
        for name, change in variants:
            with self.subTest(name=name):
                path = self.fit / name
                original = path.read_bytes()
                change(path)
                with self.assertRaises(ValueError):
                    policy.validate_fit_annotation_policy(self.original, self.fit)
                path.write_bytes(original)
        self.assertTrue(policy.validate_fit_annotation_policy(self.original, self.fit)["valid"])

    def test_classifier_explicitly_disallows_nonfinite_bounds_and_nonpositive_radius(self):
        _, rows = policy._annotations(self.original / "annotations.csv")
        frame = policy._helpers()._read_csv(self.original / "frames.csv")[0]
        accepted = {r["frame_id"] for r in rows}
        for changes in [dict(x="nan"), dict(x="320"), dict(y="inf"), dict(uncertainty_px="0"),
                        dict(uncertainty_px="nan"), dict(visibility="occluded"), dict(visibility="absent"),
                        dict(visibility="excluded"), dict(medium="mixed")]:
            with self.subTest(changes=changes):
                self.assertFalse(policy._eligibility({**rows[0], **changes}, frame, accepted)[0])


if __name__ == "__main__":
    unittest.main()
