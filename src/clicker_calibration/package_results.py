"""Copy calibration results into a portable, hash-checked review bundle.

Canonical scientific files are copied byte-for-byte. Only presentation HTML is
rewritten to link to bundled frames. Partial annotation-only delivery is explicit;
a final delivery additionally requires a model, comparison gallery and video.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import io
import json
import os
from pathlib import Path
import re
import shutil


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


class Bundle:
    def __init__(self, destination, source):
        self.root = Path(destination).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        manifest = self.root / "delivery_manifest.json"
        self.previous = read_json(manifest) if manifest.exists() else {}
        if self.previous and self.previous.get("source_directory") != str(source):
            raise ValueError("This delivery directory belongs to a different annotation revision")
        self.files = dict(self.previous.get("files_sha256", {}))

    def _target(self, relative):
        path = self.root / relative
        if self.root not in path.resolve().parents:
            raise ValueError("Bundle path escapes its destination")
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _check_replace(self, path, digest):
        if not path.exists():
            return False
        actual = sha(path)
        if actual == digest:
            return True
        previous = self.files.get(str(path.relative_to(self.root)))
        if previous is None or previous != actual:
            raise ValueError(f"Refusing to replace an unowned or externally edited result: {path}")
        return False

    def copy(self, source, relative, expected=None):
        source = Path(source)
        before = source.stat()
        digest = sha(source)
        if expected is not None and digest != expected:
            raise ValueError(f"Source checksum differs: {source}")
        target = self._target(relative)
        if not self._check_replace(target, digest):
            temporary = target.with_name(target.name + f".partial-{os.getpid()}")
            try:
                with source.open("rb") as original, temporary.open("wb") as copied:
                    shutil.copyfileobj(original, copied, 1024 * 1024)
                    copied.flush()
                    os.fsync(copied.fileno())
                after = source.stat()
                if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino) or sha(temporary) != digest:
                    raise ValueError(f"Source changed while copying: {source}")
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
        self.files[str(target.relative_to(self.root))] = digest

    def write(self, relative, text):
        data = text.encode()
        digest = hashlib.sha256(data).hexdigest()
        target = self._target(relative)
        if not self._check_replace(target, digest):
            temporary = target.with_name(target.name + f".partial-{os.getpid()}")
            with temporary.open("wb") as stream:
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, target)
        self.files[str(target.relative_to(self.root))] = digest

    def tree(self, source, relative, excluded=()):
        source = Path(source)
        for path in sorted(source.rglob("*")):
            if path.is_file() and not any(part.startswith(".") or part in excluded for part in path.relative_to(source).parts):
                self.copy(path, str(Path(relative) / path.relative_to(source)))


def portable_gallery(bundle, source_meta, annotations_sha, frames, label, model_sha=None, published_ids=None):
    metadata = read_json(source_meta)
    if metadata["annotations_sha256"] != annotations_sha:
        raise ValueError("Review gallery does not bind the current Terra annotations")
    if model_sha is not None and metadata.get("model_sha256") != model_sha:
        raise ValueError("Comparison gallery does not bind the selected calibration")
    source = Path(metadata["gallery"])
    count = 0

    def article(match):
        nonlocal count
        text = match.group(0)
        index_match = re.search(r'href="([0-9]+)\.jpg"', text)
        if index_match is None:
            raise ValueError("Gallery card has no numbered image identity")
        index = int(index_match.group(1))
        if index != count or index >= len(frames):
            raise ValueError("Gallery ordering differs from original frames")
        if f'alt="{html.escape(frames[index]["frame_id"])}"' not in text:
            raise ValueError("Gallery card frame identity differs")
        count += 1
        if published_ids is not None and frames[index]["frame_id"] not in published_ids:
            return ""
        number = f"{index:05d}"
        for suffix in (".jpg", "_thumb.jpg"):
            bundle.copy(source.parent / (number + suffix), f"frames/{label}/{number}{suffix}")
        text = text.replace(f'href="{number}.jpg"', f'href="../frames/{label}/{number}.jpg"')
        text = text.replace(f'src="{number}_thumb.jpg"', f'src="../frames/{label}/{number}_thumb.jpg"')
        raw = number + Path(frames[index]["image_path"]).suffix
        text, substitutions = re.subn(r'<a href="[^"]+"(?: target="_blank")?>Original image</a>', f'<a href="../frames/raw/{raw}" target="_blank">Original image</a>', text)
        if substitutions != 1:
            raise ValueError("Gallery original-image link is ambiguous")
        return text

    text = re.sub(r"<article\b.*?</article>", article, source.read_text(), flags=re.DOTALL)
    if count != len(frames):
        raise ValueError("Gallery does not contain every prepared frame")
    text = text.replace('href="../../annotations.csv"', 'href="../annotations.csv"')
    if published_ids is not None:
        text = text.replace("<main>", f'<p>Published review: only the {len(published_ids)} user-selected images. Original numbering is preserved; excluded rows remain in the audit CSV.</p><main>')
    directory = "annotation_review" if label == "terra" else "calibrated_review"
    bundle.write(f"{directory}/index.html", text)
    return f"{directory}/index.html"


def package_device(device_dir, destination_dir, *, model_dir=None, video_dir=None, final=False):
    source = Path(device_dir).resolve(strict=True)
    model = None if model_dir is None else Path(model_dir).resolve(strict=True)
    video = None if video_dir is None else Path(video_dir).resolve(strict=True)
    if final and (model is None or video is None):
        raise ValueError("Final delivery requires both a calibration model and completed video export")
    destination = Path(destination_dir).resolve()
    if destination == source or source in destination.parents or destination in source.parents:
        raise ValueError("Delivery must be separate from the scientific source directory")
    frames = rows(source / "frames.csv")
    annotations = rows(source / "annotations.csv")
    seal = read_json(source / "annotations.complete.json")
    annotation_sha = sha(source / "annotations.csv")
    if seal["annotations_sha256"] != annotation_sha or seal["frames_sha256"] != sha(source / "frames.csv"):
        raise ValueError("Annotation completion is stale")
    by_id = {row["frame_id"]: row for row in annotations}
    if set(by_id) != {row["frame_id"] for row in frames} or len(by_id) != len(annotations):
        raise ValueError("Annotation/frame identities differ")
    published_ids = None
    selection_path = source / "human_review_selection.json"
    if selection_path.exists():
        selection = read_json(selection_path)
        if selection.get("publish_selected_only", False):
            if seal.get("human_review_selection_sha256") != sha(selection_path) or selection.get("frames_sha256") != sha(source / "frames.csv"):
                raise ValueError("Published selection is not bound to the annotation seal/frames")
            published_ids = set(selection["accepted_frame_ids"])
            if not published_ids or not published_ids <= set(by_id) or len(published_ids) != len(selection["accepted_frame_ids"]):
                raise ValueError("Invalid published frame selection")
            if any(a["visibility"] != "excluded" or a["x"].strip() or a["y"].strip() or a["uncertainty_px"].strip()
                   for fid, a in by_id.items() if fid not in published_ids):
                raise ValueError("Omitted frames must be explicitly excluded without coordinates")
    prepared = read_json(source / "prepare_manifest.json")
    original_hashes = {str(Path(entry["path"]).resolve()): entry["sha256"] for entry in prepared["outputs"]}
    bundle = Bundle(destination, source)
    if bundle.previous and bundle.previous.get("annotations_sha256") != annotation_sha:
        raise ValueError("Changed annotations require a separate delivery revision")
    for path in sorted(source.iterdir()):
        if path.is_file() and path.suffix in {".json", ".csv", ".txt", ".f64"}:
            bundle.copy(path, "artifacts/" + path.name)
    for name in ("diagnostics", "revision_source", "preparation_raw", "annotation_work"):
        if (source / name).is_dir():
            bundle.tree(source / name, "artifacts/" + name)
    for name in ("annotations.csv", "annotations.complete.json", "human_review_selection.json"):
        if (source / name).exists():
            bundle.copy(source / name, name)
    for index, frame in enumerate(frames):
        if published_ids is not None and frame["frame_id"] not in published_ids:
            continue
        original = Path(frame["image_path"]).resolve()
        expected = original_hashes.get(str(original))
        if expected is None:
            raise ValueError("Original frame has no preparation checksum")
        bundle.copy(original, f"frames/raw/{index:05d}{original.suffix}", expected)
    annotation_gallery = portable_gallery(bundle, source / "annotation_review.json", annotation_sha, frames, "terra", published_ids=published_ids)
    comparison_gallery = bundle.previous.get("comparison_gallery")
    video_link = audio_link = None
    review_video_link = bundle.previous.get("review_video")
    comparison_video_link = bundle.previous.get("comparison_video")
    model_sha = bundle.previous.get("calibration_sha256")
    quality = bundle.previous.get("quality_message", "Annotation review ready; calibration and full video are pending.")
    if model is not None:
        calibration = read_json(model / "calibration.json")
        model_sha = sha(model / "calibration.json")
        expected = calibration.get("hashes", {}).get("annotations.csv", calibration.get("hashes", {}).get("annotations_sha256"))
        if expected != annotation_sha:
            from fit_annotation_policy import validate_fit_annotation_policy
            policy = validate_fit_annotation_policy(source, model.parent)
            if policy.get("valid") is not True or policy["fit_annotations_sha256"] != expected:
                raise ValueError("Model and original annotations have no valid explicit fit-policy bridge")
            bundle.tree(model.parent, "fit_input", excluded=("frames", "annotation_review", "fit_review", "refined_review", "annotation_work", "refinement_v1"))
        bundle.tree(model, "model")
        for name in ("calibration.json", "predictions.csv"):
            bundle.copy(model / name, name)
        gallery_meta = source / ("fit_review.json" if model == source else "refined_review.json")
        comparison_gallery = portable_gallery(bundle, gallery_meta, annotation_sha, frames, "comparison", model_sha, published_ids=published_ids)
        quality = "Fit status: " + str(calibration.get("fit_status", "unknown")) + ". Review errors and coverage in the model reports; fitted parameters are not independent physical measurements."
    if video is not None:
        export = read_json(video / "export_manifest.json")
        review_only = export.get("settings", {}).get("annotation_only", False)
        if export.get("status") != "complete" or (model is None and not review_only):
            raise ValueError("Video export is incomplete or has no matching model")
        if final and review_only:
            raise ValueError("Final delivery needs the calibrated comparison video")
        bindings = export["inputs_sha256"]
        if bindings.get(str((source / "annotations.csv").resolve())) != annotation_sha or (not review_only and bindings.get(str((model / "calibration.json").resolve())) != model_sha):
            raise ValueError("Video export belongs to a different annotation/model revision")
        for name, digest in export["outputs_sha256"].items():
            bundle.copy(video / name, "videos/" + name, digest)
        manifest_name = "review_export_manifest.json" if review_only else "comparison_export_manifest.json"
        bundle.copy(video / "export_manifest.json", "videos/" + manifest_name)
        video_link = "videos/" + str(Path(export["video_path"]).relative_to(video))
        audio_link = "videos/" + str(Path(export["reference_audio_path"]).relative_to(video))
        if review_only:
            review_video_link = video_link
        else:
            comparison_video_link = video_link
    index_stream = io.StringIO(newline="")
    index_writer = csv.DictWriter(index_stream, fieldnames=["review_index", "frame_id", "source_video", "source_frame_index", "source_pts_s", "visibility", "medium", "terra_x", "terra_y", "uncertainty_px", "raw_image", "terra_overlay", "comparison_overlay"])
    index_writer.writeheader()
    for index, frame in enumerate(frames):
        if published_ids is not None and frame["frame_id"] not in published_ids:
            continue
        annotation = by_id[frame["frame_id"]]
        index_writer.writerow(dict(review_index=index, frame_id=frame["frame_id"], source_video=frame["video_file"],
                                   source_frame_index=frame["frame_index"], source_pts_s=frame["frame_pts_s"],
                                   visibility=annotation["visibility"], medium=annotation["medium"],
                                   terra_x=annotation["x"], terra_y=annotation["y"], uncertainty_px=annotation["uncertainty_px"],
                                   raw_image=f"raw/{index:05d}{Path(frame['image_path']).suffix}",
                                   terra_overlay=f"terra/{index:05d}.jpg",
                                   comparison_overlay=f"comparison/{index:05d}.jpg" if comparison_gallery else ""))
    bundle.write("frames/index.csv", index_stream.getvalue())
    links = [(annotation_gallery, "Terra annotations and uncertainty circles")]
    if comparison_gallery:
        links.append((comparison_gallery, "Terra versus calibrated acoustic locations"))
    if review_video_link:
        links.append((review_video_link, "Synchronized human-review video: Terra overlay and hydrophone audio"))
    if comparison_video_link:
        links.append((comparison_video_link, "Full synchronized comparison video"))
    if audio_link:
        links.append((audio_link, "Lossless synchronized reference-channel review audio"))
    links.extend([("frames/raw/", "Original synchronized frames"), ("artifacts/", "Scientific inputs and reports")])
    items = "".join(f'<li><a href="{html.escape(path)}">{html.escape(label)}</a></li>' for path, label in links)
    status = "complete" if final or bundle.previous.get("status") == "complete" else "review_available"
    bundle.write("index.html", '<!doctype html><meta charset="utf-8"><title>Clicker calibration review</title>'
                 '<style>body{font:18px system-ui;max-width:1000px;margin:40px;line-height:1.6}li{margin:12px 0}</style>'
                 f'<h1>Device {html.escape(source.name)} — {status.replace("_", " ")}</h1><p>{html.escape(quality)}</p>'
                 '<p>Green: Terra location. Amber: uncertain Terra estimate. Red: calibrated acoustic estimate. '
                 'Missing or invalid locations remain explicitly identified. Scientific input files are preserved; presentation links point to bundled images.</p>'
                 f'<ul>{items}</ul>')
    result = dict(schema_version="clicker-delivery-v1", status=status, source_directory=str(source),
                  destination=str(destination), annotations_sha256=annotation_sha, calibration_sha256=model_sha,
                  source_model_directory=bundle.previous.get("source_model_directory") if model is None else str(model), n_frames=len(frames),
                  n_published_frames=len(frames) if published_ids is None else len(published_ids),
                  publication_scope="all_frames" if published_ids is None else "human_selected_only",
                  annotation_gallery=annotation_gallery, comparison_gallery=comparison_gallery,
                  video=video_link, review_video=review_video_link, comparison_video=comparison_video_link,
                  review_audio=audio_link or bundle.previous.get("review_audio"), files_sha256=bundle.files,
                  quality_message=quality, created_utc=datetime.now(timezone.utc).isoformat())
    path = destination / "delivery_manifest.json"
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    os.replace(temporary, path)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("device_dir")
    parser.add_argument("destination_dir")
    parser.add_argument("--model-dir")
    parser.add_argument("--video-dir")
    parser.add_argument("--final", action="store_true")
    result = package_device(**vars(parser.parse_args()))
    print(json.dumps({key: value for key, value in result.items() if key != "files_sha256"}))


if __name__ == "__main__":
    main()
