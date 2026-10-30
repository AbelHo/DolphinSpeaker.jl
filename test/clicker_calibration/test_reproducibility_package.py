"""Synthetic Git/copy controls; never packages actual recordings or runs models."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock


SOURCE = Path(__file__).parents[2] / "tools/package_clicker_reproducibility.py"
spec = importlib.util.spec_from_file_location("reproducibility_package", SOURCE)
package = importlib.util.module_from_spec(spec)
spec.loader.exec_module(package)


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(package.json_bytes(value))


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


class ReproducibilityPackageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = self.root / "repository"
        self.repo.mkdir()
        self.g("init", "-q")
        self.g("config", "user.name", "Synthetic Fixture")
        self.g("config", "user.email", "fixture@example.invalid")
        write(self.repo / "Project.toml", 'name = "SyntheticOnly"\n')
        write(self.repo / "Manifest.toml", 'julia_version = "1.12.6"\n')
        write(self.repo / "README.md", "Synthetic source, no executable research claim\n")
        write(self.repo / "docs/clicker_calibration.md", "Synthetic reproduction documentation\n")
        write(self.repo / "test/fixture_test.py", "# synthetic test source\n")
        write(self.repo / "test/fixtures/delays.csv", "d1,d2,d3\n0,0,0\n")
        write(self.repo / "tools/package_clicker_reproducibility.py", SOURCE.read_text())
        self.module = self.repo / "src/clicker_calibration"
        for name in ("physical_fit.jl", "prepare.jl", "srp_review.jl", "method_refinement.jl", "support.jl"):
            write(self.module / name, f"# frozen synthetic {name}\n")
        write(self.module / "video_export.py", "# original exporter bytes\n")
        self.old_physical = package.sha(self.module / "physical_fit.jl")
        self.old_exporter = package.sha(self.module / "video_export.py")
        self.g("add", ".")
        self.g("commit", "-qm", "original synthetic source")
        self.old_commit = self.g("rev-parse", "HEAD").strip()
        write(self.module / "physical_fit.jl", "# current synthetic physical fitter\n")
        write(self.module / "video_export.py", "# current exporter bytes\n")
        self.g("add", ".")
        self.g("commit", "-qm", "current synthetic source")
        self.head = self.g("rev-parse", "HEAD").strip()
        self.original = self.root / "numerical/A"
        self.model = self.original / "refinement_v1"
        self.model.mkdir(parents=True)
        write(self.original / "events.csv", "event_id,frame_id\ne1,f1\n")
        write(self.original / "frames.csv", "frame_id,width,height\nf1,320,180\n")
        write(self.original / "annotations.csv", "frame_id,x,y,visibility,medium,uncertainty_px,note,annotator\nf1,1.2300e2,80.2500,uncertain,underwater,3.500e1,original point,gpt-6.1-terra\n")
        dump(self.original / "annotations.complete.json", {
            "annotations_sha256": package.sha(self.original / "annotations.csv"),
            "frames_sha256": package.sha(self.original / "frames.csv"), "n_frames": 1})
        dump(self.original / "human_review_selection.json", {"accepted_frame_ids": ["f1"], "publish_selected_only": True})
        write(self.original / "event_snippets.f64", "synthetic derived cache bytes")
        write(self.original / "predictions.csv", "event_id,split\ne1,train\n")
        dump(self.original / "fit_report.json", {"heldout_metrics": {"median_px": 7.0}, "observation_policy": {"provisional": True},
            "hashes": {"physical_fit_source_sha256": self.old_physical, "support_source_sha256": package.sha(self.module / "support.jl")}})
        self.calibration = dict(device="A", method="gcc", fit_status="fitted_weak_identifiability", sound_speed_m_s=1500.0,
                                fov_degrees=[80.0, 55.0], relative_delay_offsets_s=[0.0, 0.0, 0.0],
                                parameter_bound_hits_normalized=[1], software={"Julia": "1.12.6"},
                                hashes={"source/" + name: self.old_physical if name == "physical_fit.jl" else package.sha(self.module / name)
                                        for name in ("physical_fit.jl", "prepare.jl", "srp_review.jl", "method_refinement.jl", "support.jl")})
        dump(self.model / "calibration.json", self.calibration)
        dump(self.model / "evaluation.json", {"interpretation": "reused descriptive evaluation"})
        dump(self.model / "selection.json", {"uses_heldout_labels": False})
        write(self.model / "predictions.csv", "event_id,predicted_x\ne1,125.0\n")
        self.report = dict(selected_method="gcc", fit_status="fitted_weak_identifiability", n_training=1,
                           n_reused_evaluation=1, n_deployment=2, training_metrics={"median_px": 2.0},
                           reused_evaluation_metrics={"median_px": 3.0, "rmse_px": 20.0, "coverage": 0.5},
                           deployment_metrics_descriptive={"median_px": 2.5}, limits=["uncertain visual estimates"],
                           hashes=self.calibration["hashes"], outputs_sha256={})
        self.bind_model()
        self.results = self.root / "results"
        for phase, exporter in (("review", self.old_exporter), ("comparison", package.sha(self.module / "video_export.py"))):
            dump(self.results / "A/videos" / (phase + "_export_manifest.json"),
                 {"status": "complete", "settings": {"exporter_sha256": exporter}})
        write(self.results / "A/videos/full_movie.mp4", "do not duplicate full movies")
        write(self.results / "A/index.html", "existing device bundle untouched")
        self.runs = self.root / "run"
        self.never = self.root / "MUST_NOT_BE_EXECUTED"
        dump(self.runs / "jobs/job01/result.json", {"status": "failed", "command": ["touch", str(self.never)]})
        write(self.runs / "jobs/job01/process.log", "retained failure trace\n")
        self.history = self.root / "old_A1"
        dump(self.history / "calibration.json", {"method": "minmax2", "fit_status": "superseded"})
        write(self.history / "predictions.csv", "event_id,x\nold,20.0\n")
        self.arguments = dict(devices={"A": {"original_dir": self.original, "model_dir": self.model}},
                              run_roots={"run": self.runs}, history_dirs={"A1": self.history})

    def g(self, *arguments):
        return subprocess.check_output(["git", "-C", str(self.repo), *arguments], text=True, stderr=subprocess.STDOUT)

    def bind_model(self):
        self.report["outputs_sha256"] = {name: package.sha(self.model / name) for name in ("calibration.json", "predictions.csv", "evaluation.json", "selection.json")}
        dump(self.model / "report.json", self.report)

    def run_package(self, **changes):
        args = {**self.arguments, **changes}
        return package.package_reproducibility(self.repo, self.results, **args)

    def generation(self, result):
        return Path(result["root"]) / result["generation"]

    def test_canonical_code_numerical_history_and_summary_are_preserved(self):
        original_bytes = (self.original / "annotations.csv").read_bytes()
        bundle_index = (self.results / "A/index.html").read_bytes()
        result = self.run_package()
        gen = self.generation(result)
        manifest = package.read_json(gen / "generation_manifest.json")
        for relative, expected in manifest["files_sha256"].items():
            self.assertEqual(package.sha(gen / relative), expected)
        self.assertEqual((gen / "code/Project.toml").read_bytes(), (self.repo / "Project.toml").read_bytes())
        self.assertEqual((gen / "code/test/fixtures/delays.csv").read_bytes(), (self.repo / "test/fixtures/delays.csv").read_bytes())
        self.assertEqual((gen / "latest/A/original/annotations.csv").read_bytes(), original_bytes)
        self.assertEqual((gen / "superseded_history/A1/predictions.csv").read_bytes(), (self.history / "predictions.csv").read_bytes())
        summary = package.read_json(gen / "summary.json")
        self.assertEqual(summary["git_head"], self.head)
        self.assertEqual(summary["latest"]["A"]["evaluation"]["reused_after_method_selection"], self.report["reused_evaluation_metrics"])
        self.assertIn("not untouched", summary["latest"]["A"]["evaluation"]["reused_label"])
        self.assertIn("In-sample", summary["latest"]["A"]["evaluation"]["deployment_label"])
        self.assertEqual(summary["latest"]["A"]["original_annotations"]["visibility_counts"], {"uncertain": 1})
        self.assertEqual(summary["history"]["A1"]["status"], "superseded_not_current_accuracy")
        self.assertEqual((self.results / "A/index.html").read_bytes(), bundle_index)
        self.assertEqual((self.original / "annotations.csv").read_bytes(), original_bytes)
        self.assertFalse(self.never.exists())

    def test_old_fitter_and_two_media_exporter_versions_recover_exact_commits(self):
        gen = self.generation(self.run_package())
        versions = package.read_json(gen / "source_versions.json")
        physical = next(row for row in versions["latest/A"] if row["path"].endswith("physical_fit.jl"))
        self.assertEqual(physical["sha256"], self.old_physical)
        self.assertEqual(physical["commit"], self.old_commit)
        self.assertEqual(package.sha(gen / physical["snapshot"]), self.old_physical)
        self.assertEqual(versions["latest/A/media/review"][0]["commit"], self.old_commit)
        self.assertEqual(versions["latest/A/media/comparison"][0]["commit"], self.head)
        self.assertNotEqual(versions["latest/A/media/review"][0]["sha256"], versions["latest/A/media/comparison"][0]["sha256"])

    def test_repeat_reuses_verified_generation_new_evidence_preserves_history(self):
        first = self.run_package()
        pointer = Path(first["root"]) / "manifest.json"
        before = pointer.read_bytes()
        second = self.run_package()
        self.assertTrue(second["reused"])
        self.assertEqual(pointer.read_bytes(), before)
        write(self.runs / "jobs/job02/process.log", "new terminal evidence\n")
        third = self.run_package()
        self.assertFalse(third["reused"])
        self.assertNotEqual(third["generation"], first["generation"])
        self.assertTrue(self.generation(first).is_dir())
        self.assertEqual(len(third["generations"]), 2)

    def test_foreign_payload_edits_and_unowned_output_files_are_rejected(self):
        result = self.run_package()
        target = self.generation(result) / "code/README.md"
        target.write_text("foreign edit")
        with self.assertRaisesRegex(ValueError, "externally edited"):
            self.run_package()
        self.assertEqual(target.read_text(), "foreign edit")
        target.write_bytes((self.repo / "README.md").read_bytes())
        write(Path(result["root"]) / "unowned.txt", "leave this alone")
        with self.assertRaisesRegex(ValueError, "Unowned"):
            self.run_package()

    def test_dirty_selected_code_and_uncommitted_packager_are_rejected(self):
        write(self.module / "prepare.jl", "dirty source must not be packaged as HEAD\n")
        with self.assertRaisesRegex(ValueError, "Worktree/source bytes differ"):
            self.run_package()
        self.g("checkout", "--", "src/clicker_calibration/prepare.jl")
        self.g("rm", "--cached", "tools/package_clicker_reproducibility.py")
        self.g("commit", "-qm", "packager is not committed")
        with self.assertRaisesRegex(ValueError, "Commit Project/Manifest"):
            self.run_package()

    def test_missing_frozen_source_and_bounded_history_fail_without_substitution(self):
        with self.assertRaisesRegex(ValueError, "Cannot recover frozen source"):
            self.run_package(max_git_commits=1)
        self.calibration["hashes"]["source/physical_fit.jl"] = "0" * 64
        dump(self.model / "calibration.json", self.calibration)
        self.report["hashes"] = self.calibration["hashes"]
        self.bind_model()
        with self.assertRaisesRegex(ValueError, "Cannot recover frozen source"):
            self.run_package()
        self.assertFalse((self.results / "reproducibility/manifest.json").exists())

    def test_recordings_private_environments_symlinks_and_secret_logs_are_not_copied(self):
        for rel in ("recording.WAV", "video.MP4", "runtime/julia_depot/config.json", "scratch/private.json", "frames/raw.png", ".env"):
            write(self.original / rel, "must not copy " + rel)
        write(self.runs / "jobs/job01/scratch/secret.json", "private child scratch")
        outside = self.root / "private.json"
        write(outside, "private symlink target")
        (self.original / "linked.json").symlink_to(outside)
        secret = "bot" + "123456789" + ":" + "x" * 35
        write(self.runs / "jobs/job01/credential_error.log", secret)
        result = self.run_package()
        gen = self.generation(result)
        manifest = package.read_json(gen / "generation_manifest.json")
        names = list(manifest["files_sha256"])
        self.assertFalse(any(name.lower().endswith((".wav", ".mp4", ".env")) for name in names))
        self.assertFalse(any("scratch" in name or "julia_depot" in name or "linked.json" in name for name in names))
        self.assertEqual(len(manifest["excluded_files"]), 1)
        self.assertTrue(all(secret.encode() not in p.read_bytes() for p in gen.rglob("*") if p.is_file()))
        self.assertEqual((gen / "runs/run/jobs/job01/process.log").read_bytes(), (self.runs / "jobs/job01/process.log").read_bytes())

    def test_latest_model_and_annotation_seal_tampering_are_rejected(self):
        path = self.model / "predictions.csv"
        original = path.read_bytes()
        path.write_bytes(original + b"modified")
        with self.assertRaisesRegex(ValueError, "model digest"):
            self.run_package()
        path.write_bytes(original)
        with (self.original / "annotations.csv").open("ab") as stream:
            stream.write(b"\n")
        with self.assertRaisesRegex(ValueError, "not currently sealed"):
            self.run_package()

    def test_mismatched_explicit_original_model_mapping_is_rejected(self):
        self.calibration["hashes"]["events.csv"] = "0" * 64
        dump(self.model / "calibration.json", self.calibration)
        self.bind_model()
        with self.assertRaisesRegex(ValueError, "model/input mapping differs"):
            self.run_package()

    def test_copy_bounds_and_private_explicit_roots_are_enforced(self):
        with self.assertRaisesRegex(ValueError, "copy bound"):
            self.run_package(max_bytes=100)
        private = self.root / "runtime"
        write(private / "private.json", "do not copy")
        with self.assertRaisesRegex(ValueError, "cannot be explicit"):
            self.run_package(history_dirs={"private": private})
        with self.assertRaises(ValueError):
            package.mappings(["../escape=/tmp"])

    def test_active_jobs_are_not_snapshotted_as_completed_reproducibility_evidence(self):
        dump(self.runs / "jobs/active/result.json", {"status": "running"})
        write(self.runs / "jobs/active/process.log", "mutable active output")
        result = self.run_package(run_roots={"jobs": self.runs / "jobs"})
        gen = self.generation(result)
        manifest = package.read_json(gen / "generation_manifest.json")
        self.assertTrue((gen / "runs/jobs/job01/result.json").is_file())
        self.assertFalse(any("/active/" in name for name in manifest["files_sha256"]))
        self.assertTrue(any(row.get("status") == "running" for row in manifest["excluded_files"]))

    def test_canonical_source_mutation_during_copy_leaves_no_active_pointer(self):
        original_write = package.write_file
        changed = False
        def mutate(path, data):
            nonlocal changed
            original_write(path, data)
            if not changed and "latest" in path.parts:
                changed = True
                with (self.original / "annotations.csv").open("ab") as stream:
                    stream.write(b"\n")
        with mock.patch.object(package, "write_file", side_effect=mutate):
            with self.assertRaisesRegex(ValueError, "changed"):
                self.run_package()
        self.assertFalse((self.results / "reproducibility/manifest.json").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
