import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from PIL import Image


def module(name):
    path = Path(__file__).parents[2] / "src/clicker_calibration" / (name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


pack = module("package_results")
ann = module("annotation_tools")


class PortableResults(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source" / "A"
        self.source.mkdir(parents=True)
        self.destination = self.root / "results" / "A"
        image = self.source / "frame.png"
        Image.new("RGB", (320, 180), "#456789").save(image)
        frame = dict(frame_id="A-frame", video_file="camera.mp4", frame_index=12,
                     frame_pts_s=.2, image_path=str(image), width=320, height=180)
        annotation = dict(frame_id="A-frame", x=100, y=80, visibility="visible",
                          medium="underwater", uncertainty_px=20, note="fixture", annotator="gpt-6.1-terra")
        for name, row in (("frames", frame), ("annotations", annotation)):
            with (self.source / (name + ".csv")).open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(row))
                writer.writeheader()
                writer.writerow(row)
        (self.source / "prepare_manifest.json").write_text(json.dumps(dict(outputs=[dict(path=str(image), sha256=pack.sha(image))])))
        ann.validate(self.source, seal=True)
        ann.gallery(self.source)

    def test_portable_links_original_bytes_and_resume(self):
        result = pack.package_device(self.source, self.destination)
        self.assertEqual(result["status"], "review_available")
        self.assertEqual(result["n_frames"], 1)
        self.assertEqual(pack.sha(self.source / "annotations.csv"), pack.sha(self.destination / "annotations.csv"))
        self.assertEqual(pack.sha(self.source / "frame.png"), pack.sha(self.destination / "frames/raw/00000.png"))
        page = (self.destination / "annotation_review/index.html").read_text()
        self.assertIn('../frames/raw/00000.png', page)
        self.assertIn('../frames/terra/00000.jpg', page)
        self.assertIn('../annotations.csv', page)
        self.assertNotIn(str(self.source), page)
        self.assertTrue((self.destination / "frames/index.csv").is_file())
        again = pack.package_device(self.source, self.destination)
        self.assertEqual(result["files_sha256"], again["files_sha256"])

    def test_external_edits_are_not_overwritten(self):
        pack.package_device(self.source, self.destination)
        edited = self.destination / "annotations.csv"
        edited.write_text("a human edit\n")
        with self.assertRaisesRegex(ValueError, "externally edited"):
            pack.package_device(self.source, self.destination)
        self.assertEqual(edited.read_text(), "a human edit\n")

    def test_frame_tamper_and_stale_annotations_are_rejected(self):
        (self.source / "frame.png").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "checksum"):
            pack.package_device(self.source, self.destination)
        with (self.source / "annotations.csv").open("a") as f:
            f.write("\n")
        with self.assertRaisesRegex(ValueError, "stale"):
            pack.package_device(self.source, self.destination)

    def test_final_requires_model_and_video_and_source_stays_separate(self):
        with self.assertRaisesRegex(ValueError, "Final delivery"):
            pack.package_device(self.source, self.destination, final=True)
        with self.assertRaisesRegex(ValueError, "separate"):
            pack.package_device(self.source, self.source / "delivery")


if __name__ == "__main__":
    unittest.main()
