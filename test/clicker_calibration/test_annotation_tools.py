import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from PIL import Image

spec = importlib.util.spec_from_file_location("annotation_tools", Path(__file__).parents[2] / "src/clicker_calibration/annotation_tools.py")
tools = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tools)


class AnnotationContract(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        image = self.root / "source.png"
        Image.new("RGB", (320, 180), "#456789").save(image)
        self.frames = [dict(frame_id="A-test-f1", video_file="source.mp4", frame_index=1,
                            frame_pts_s=0.04, video_global_time_s=0.04,
                            image_path=str(image), width=320, height=180)]
        self.annotations = [dict(frame_id="A-test-f1", x=160, y=90, visibility="visible",
                                 medium="underwater", uncertainty_px=3,
                                 note="fixture", annotator="test")]
        self.write()

    def tearDown(self):
        self.temp.cleanup()

    def write(self):
        for name, rows in (("frames", self.frames), ("annotations", self.annotations)):
            with (self.root / f"{name}.csv").open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=rows[0])
                writer.writeheader()
                writer.writerows(rows)

    def test_seal_and_overlay_preserve_original_geometry(self):
        _, manifest = tools.validate(self.root, seal=True)
        self.assertEqual(manifest["counts"]["visible"], 1)
        review = tools.gallery(self.root)
        gallery = Path(review["gallery"])
        self.assertTrue(gallery.is_file())
        with Image.open(gallery.parent / "00000.jpg") as im:
            self.assertEqual(im.size, (320, 180))
            r, g, b = im.getpixel((160, 90))
            self.assertGreater(g, r+70)
        self.assertEqual(tools.sha(self.root / "source.png"), tools.sha(Path(self.frames[0]["image_path"])))

    def test_stale_seal_rejects_edits(self):
        tools.validate(self.root, seal=True)
        self.annotations[0]["x"] = 161
        self.write()
        with self.assertRaisesRegex(ValueError, "Stale"):
            tools.validate(self.root)

    def test_duplicate_and_missing_identity(self):
        self.annotations.append(self.annotations[0].copy())
        self.write()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            tools.validate(self.root, seal=True)
        self.annotations = [{**self.annotations[0], "frame_id": "different"}]
        self.write()
        with self.assertRaisesRegex(ValueError, "identities differ"):
            tools.validate(self.root, seal=True)

    def test_invalid_coordinates_and_medium(self):
        for key, value in (("x", 320), ("y", -1), ("x", "nan"),
                           ("uncertainty_px", 0), ("medium", "pool?")):
            with self.subTest(key=key, value=value):
                old = self.annotations[0][key]
                self.annotations[0][key] = value
                self.write()
                with self.assertRaises(ValueError):
                    tools.validate(self.root, seal=True)
                self.annotations[0][key] = old

    def test_absent_rows_remain_reviewable(self):
        self.annotations[0].update(x="", y="", uncertainty_px="", visibility="absent", medium="air")
        self.write()
        _, manifest = tools.validate(self.root, seal=True)
        self.assertEqual(manifest["counts"]["absent"], 1)
        self.assertTrue(Path(tools.gallery(self.root)["gallery"]).exists())

    def test_actual_dimensions_are_verified(self):
        self.frames[0]["width"] = 640
        self.write()
        with self.assertRaisesRegex(ValueError, "dimensions differ"):
            tools.validate(self.root, seal=True)


if __name__ == "__main__":
    unittest.main()
