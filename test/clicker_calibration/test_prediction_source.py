"""All-event presentation controls: no fit or estimator changes."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PS = load("prediction_source", ROOT / "src/clicker_calibration/prediction_source.py")
MEDIA = load("video_fixture", Path(__file__).with_name("test_video_export.py"))
ANNOTATION = load("annotation_fixture", ROOT / "src/clicker_calibration/annotation_tools.py")


class PredictionSourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.device, self.model, self.output = MEDIA.fixture(self.root)
        annotations = PS._rows(self.device / "annotations.csv")
        for row in annotations:
            if row["visibility"] == "excluded":
                row.update(x="", y="", uncertainty_px="")
        MEDIA.write_csv(self.device / "annotations.csv", annotations)
        ANNOTATION.validate(self.device, seal=True)
        calibration = json.loads((self.model / "calibration.json").read_text())
        for name in ("annotations.csv", "annotations.complete.json"):
            calibration["hashes"][name] = PS._sha(self.device / name)
        (self.model / "calibration.json").write_text(json.dumps(calibration))
        report = dict(outputs_sha256={name: PS._sha(self.model / name) for name in ("calibration.json", "predictions.csv")})
        (self.model / "report.json").write_text(json.dumps(report))
        self.cohort_before = (self.model / "predictions.csv").read_bytes()
        events = PS._rows(self.device / "events.csv")
        predicted = [dict(event_id=e["event_id"], frame_id=e["frame_id"],
                          predicted_x=105, predicted_y=90, predicted_valid="true",
                          predicted_status="ok", sync_status=e["sync_status"]) for e in events]
        MEDIA.write_csv(self.model / "inference_events.csv", predicted)
        required = {
            "events_sha256": self.device / "events.csv", "frames_sha256": self.device / "frames.csv",
            "annotations_sha256": self.device / "annotations.csv",
            "annotations_complete_sha256": self.device / "annotations.complete.json",
            "prepare_manifest_sha256": self.device / "prepare_manifest.json",
            "calibration_sha256": self.model / "calibration.json",
            "cohort_predictions_sha256": self.model / "predictions.csv",
            "model_report_sha256": self.model / "report.json",
        }
        self.manifest = dict(schema="clicker-all-event-projection-v1", status="complete",
            original_dir=str(self.device), model_dir=str(self.model),
            output_csv=str(self.model / "inference_events.csv"),
            output_csv_sha256=PS._sha(self.model / "inference_events.csv"),
            hashes={k: PS._sha(v) for k, v in required.items()},
            inputs={str(v): PS._sha(v) for v in required.values()},
            source_hashes={name: PS._sha(ROOT / "src/clicker_calibration" / name) for name in
                ("event_projection.jl", "physical_fit.jl", "prepare.jl", "srp_review.jl", "method_refinement.jl", "support.jl", "workflow.jl")},
            counts=dict(n_original_events=len(events)))
        self.save()

    def save(self):
        (self.model / "inference_events.json").write_text(json.dumps(self.manifest))

    def tearDown(self):
        self.assertEqual(self.cohort_before, (self.model / "predictions.csv").read_bytes())
        self.tmp.cleanup()

    def test_media_consumes_new_noncohort_event_and_keeps_original_truth(self):
        loaded = MEDIA.VE._load(self.device, self.model, self.output, "/usr/bin/ffmpeg", "/usr/bin/ffprobe", self.root)
        marks = {row["event_id"]: row for v in loaded["videos"] for row in v["events"]}
        self.assertIsNotNone(marks["event5"]["acoustic"])
        self.assertIsNone(marks["event4"]["acoustic"])  # Explicitly user-excluded.
        self.assertEqual(marks["event3"]["visibility"], "uncertain")
        self.assertEqual(loaded["prediction_provenance"]["scope"], "all_original_events")
        self.assertIn(str(self.model / "inference_events.json"), loaded["hashes"])
        review = ANNOTATION.gallery(self.device, fitted=True, model_dir=self.model)
        self.assertIn(self.manifest["output_csv_sha256"][:12], review["revision"])
        self.assertEqual(review["prediction_provenance"]["scope"], "all_original_events")
        with ANNOTATION.Image.open(Path(review["gallery"]).parent / "00005.jpg") as image:
            r, g, b = image.getpixel((105, 90))
            self.assertGreater(r, g + 50)
        with ANNOTATION.Image.open(Path(review["gallery"]).parent / "00004.jpg") as image:
            r, g, b = image.getpixel((105, 90))
            self.assertLess(abs(r-g), 10)

    def test_tampered_csv_rejected(self):
        with (self.model / "inference_events.csv").open("a") as out: out.write("\n")
        with self.assertRaisesRegex(ValueError, "digest"):
            PS.resolve_predictions(self.device, self.model)

    def test_missing_reordered_events_rejected_even_with_updated_csv_digest(self):
        path = self.model / "inference_events.csv"
        MEDIA.write_csv(path, list(reversed(PS._rows(path))))
        self.manifest["output_csv_sha256"] = PS._sha(path); self.save()
        with self.assertRaisesRegex(ValueError, "original order"):
            PS.resolve_predictions(self.device, self.model)

    def test_changed_source_and_changed_annotations_rejected(self):
        self.manifest["source_hashes"]["workflow.jl"] = "0" * 64; self.save()
        with self.assertRaisesRegex(ValueError, "producer"):
            PS.resolve_predictions(self.device, self.model)
        self.manifest["source_hashes"]["workflow.jl"] = PS._sha(ROOT / "src/clicker_calibration/workflow.jl"); self.save()
        with (self.device / "annotations.csv").open("a") as out: out.write("\n")
        with self.assertRaisesRegex(ValueError, "binding"):
            PS.resolve_predictions(self.device, self.model)

    def test_partial_sidecar_rejected_and_cohort_fallback_requires_both_absent(self):
        (self.model / "inference_events.json").unlink()
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            PS.resolve_predictions(self.device, self.model)
        (self.model / "inference_events.csv").unlink()
        path, provenance = PS.resolve_predictions(self.device, self.model)
        self.assertEqual(path, self.model / "predictions.csv")
        self.assertEqual(provenance["scope"], "frozen_fit_cohort")


if __name__ == "__main__":
    unittest.main()
