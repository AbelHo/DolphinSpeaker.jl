"""Source-frame-preserving review media. No fitting or import-time work.

Chapters are encoded sequentially at a microsecond time base, then stream-copy
joined using DECLARED camera audio durations. Source frames are never resampled
to nominal FPS. Green=Terra housing; amber=uncertain Terra; red=acoustic final.
48kHz float PCM/AAC are audible review media, not scientific waveform outputs.
"""
from __future__ import annotations

import argparse
import bisect
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import traceback
import uuid
import wave

import numpy as np


RATE = 48000
HOLD = 0.12


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path):
    return json.loads(Path(path).read_text())


def _write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("w") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def _rows(path, key):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    result = {row[key]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate {key}: {path}")
    return result


def _number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def _true(value):
    return str(value).lower() in {"true", "1", "yes"}


def _log_failure(output, operation, exception):
    record = dict(timestamp=datetime.now(timezone.utc).isoformat(), pid=os.getpid(),
                  worker="video_export", operation=operation,
                  exception_type=type(exception).__name__, exception=str(exception),
                  traceback=traceback.format_exc())
    targets = [Path(output) / "export_errors.jsonl"]
    configured = os.environ.get("CLICKER_ERROR_LOG")
    if configured:
        targets.append(Path(configured))
    for path in targets:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()
            os.fsync(stream.fileno())


def _run(command, work, label, *, capture=False):
    """Keep complete stderr/stdout durably, including failed subprocesses."""
    log = Path(work) / (label + ".log")
    with log.open("wb") as err:
        err.write((json.dumps(list(map(str, command)))+"\n").encode())
        result = subprocess.run(list(map(str, command)), stdout=subprocess.PIPE if capture else err,
                                stderr=err, check=False)
    if result.returncode:
        raise RuntimeError(f"{label} exit={result.returncode}; command={command!r}; log={log}; "
                           + log.read_text(errors="replace")[-8000:])
    return result.stdout if capture else None


def _probe(path, ffprobe, work, label):
    data = _run([ffprobe, "-v", "error", "-threads", "1", "-select_streams", "v:0", "-show_frames", "-show_streams",
                 "-show_entries", "stream=width,height,duration,time_base,nb_frames:frame=best_effort_timestamp_time,pkt_duration_time",
                 "-of", "json", path], work, label, capture=True)
    value = json.loads(data)
    stream = value["streams"][0]
    pts = [float(row["best_effort_timestamp_time"]) for row in value["frames"]]
    if not pts or any(b <= a for a, b in zip(pts, pts[1:])):
        raise ValueError(f"Missing/nonmonotone source frames: {path}")
    return dict(pts=pts, width=int(stream["width"]), height=int(stream["height"]),
                duration=float(stream["duration"]), time_base=stream["time_base"])


def _quote_file(path):
    text = str(Path(path).resolve())
    if "\n" in text or "\r" in text:
        raise ValueError("Newline in media path")
    return "'" + text.replace("'", "'\\''") + "'"


def _filter_path(path):
    return str(Path(path).resolve()).replace("\\", "\\\\").replace(":", "\\:").replace("'", "\\'")


def _timing(metadata, sync, eligibility=None):
    bases = []
    total = Decimal(0)
    for video in metadata["videos"]:
        bases.append(float(total))
        total += Decimal(str(video["audio_duration_s"]))
    if eligibility is not None:
        review = eligibility["review"]
        for video, base, saved in zip(metadata["videos"], bases, review["chapter_bases"], strict=True):
            if (Path(saved["video_file"]).resolve() != Path(video["path"]).resolve()
                    or abs(float(saved["exact_declared_start_s"]) - base) > 1e-9):
                raise ValueError("Exact-base timing review differs from source chapters")
        fit = review["counterfactual_fit"]
        origin = "saved exact-declared-base counterfactual affine fit"
    else:
        for video, base in zip(metadata["videos"], bases):
            if abs(float(video["global_start_s"]) - base) > 1e-9:
                raise ValueError("Stored sync uses coarse chapter bases; exact-base review is required")
        fit = sync
        origin = "saved affine fit with verified exact declared chapter bases"
    offset, scale = float(fit["offset_s"]), float(fit["scale"])
    if not math.isfinite(offset) or not math.isfinite(scale) or scale <= 0 or sync.get("status") != "valid":
        raise ValueError("Invalid saved affine synchronization")
    return dict(bases=bases, offset_s=offset, scale=scale, origin=origin)


def ass_interval(previous_s, target_s, hold_s=HOLD):
    """libass evaluates integer milliseconds; choose a centisecond strictly
    after the preceding frame and at/before the target. None means no such tick
    exists (e.g. some >100fps VFR intervals), requiring frame-index fallback.
    ASS end is rounded upward; nominal120ms hold can end <10ms later.
    """
    target_ms = math.floor(target_s * 1000 + 1e-7)
    previous_ms = -1 if previous_s is None else math.floor(previous_s * 1000 + 1e-7)
    start = max(0, previous_ms // 10 + 1)
    if start * 10 > target_ms:
        return None
    end = max(start + 1, math.ceil((target_s + hold_s) * 100 - 1e-9))
    return start, end


def _ass_time(centiseconds):
    h, rem = divmod(centiseconds, 360000)
    m, rem = divmod(rem, 6000)
    s, cs = divmod(rem, 100)
    return f"{h}:{m:02}:{s:02}.{cs:02}"


def _clean(text):
    return str(text).replace("\\", "/").replace("{", "(").replace("}", ")").replace("\n", " ")


def _event_marks(event, annotation, prediction, width, height, *, accepted=True, sync_status=None,annotation_only=False):
    visibility = annotation.get("visibility", "absent")
    excluded = (not accepted or visibility == "excluded" or _true(annotation.get("user_excluded", False)))
    medium = annotation.get("medium", "unknown")
    truth = None
    ax, ay = _number(annotation.get("x")), _number(annotation.get("y"))
    if visibility in {"visible", "uncertain"} and ax is not None and ay is not None and 0 <= ax < width and 0 <= ay < height:
        truth = dict(x=ax, y=ay, color="00FF00" if visibility == "visible" else "00BFFF", radius=11)
    status = "missing_prediction"
    acoustic = None
    if annotation_only:
        status = None
    elif excluded:
        status = "user_excluded"
    elif medium != "underwater":
        status = "suppressed_medium_" + medium
    elif (sync_status or event.get("sync_status")) != "mapped":
        status = "suppressed_sync_" + str(sync_status or event.get("sync_status", "unknown"))
    elif prediction is not None:
        x = _number(prediction.get("predicted_x", prediction.get("final_x")))
        y = _number(prediction.get("predicted_y", prediction.get("final_y")))
        valid = _true(prediction.get("predicted_valid", prediction.get("final_valid", False)))
        status = prediction.get("predicted_status", prediction.get("final_status", "invalid_prediction"))
        if valid and x is not None and y is not None:
            if 0 <= x < width and 0 <= y < height:
                status = "ok"
                acoustic = dict(x=x, y=y, color="0000FF", radius=6)
            else:
                status = "outside_image"
        elif valid:
            status = "invalid_nonfinite_prediction"
    label=f"{event['event_id']} | Terra:{visibility} u={annotation.get('uncertainty_px','?')}px {medium}"
    if not annotation_only:
        label+=f" | Acoustic:{status}"
    return dict(terra=truth, acoustic=acoustic, acoustic_status=status, visibility=visibility,
                medium=medium, user_excluded=excluded,uncertainty_px=_number(annotation.get("uncertainty_px")),
                annotator=annotation.get("annotator",""),
                label=_clean(label))


def _drawing(mark,*,outline=False):
    r = mark["radius"]
    vertices = [(r + r * math.cos(2 * math.pi * i / 20), r + r * math.sin(2 * math.pi * i / 20)) for i in range(20)]
    path = "m " + " l ".join(f"{x:.3f} {y:.3f}" for x, y in vertices)
    style=f"\\bord2\\1a&HFF&\\3a&H00&\\3c&H{mark['color']}&" if outline else f"\\bord0\\1c&H{mark['color']}&"
    return (f"{{\\an7\\pos({mark['x']-r:.3f},{mark['y']-r:.3f})\\p1\\shad0{style}}}" + path)


def _ass_header(width, height):
    size = max(12, round(width / 100))
    return f"""[Script Info]
ScriptType: v4.00+
PlayResX: {width}
PlayResY: {height}
ScaledBorderAndShadow: yes
[V4+ Styles]
Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding
Style: Text,DejaVu Sans,{size},&H00FFFFFF,&H00FFFFFF,&H00000000,&H80000000,0,0,0,0,100,100,0,0,1,1,0,7,8,8,8,1
[Events]
Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text
"""


def _write_ass(path, width, height, pts, events,*,annotation_only=False):
    lines = [_ass_header(width, height)]
    finish = max(1, math.ceil((pts[-1] + 1) * 100))
    legend = "GREEN Terra housing | AMBER uncertain | outline=uncertainty_px"
    if not annotation_only:
        legend+=" | RED acoustic"
    legend+=" | nominal120ms hold capped at chapter end; no interpolation"
    lines.append(f"Dialogue: 0,0:00:00.00,{_ass_time(finish)},Text,,0,0,0,,{{\\an7\\pos(8,8)}}{legend}\n")
    fallback = []
    for lane, item in enumerate(events):
        index = item["frame_index"]
        interval = ass_interval(pts[index - 1] if index else None, pts[index])
        item["display_hold_nominal_s"] = HOLD
        item["ass_interval_cs"] = list(interval) if interval else None
        if interval is None:
            last = max(index, bisect.bisect_left(pts, pts[index] + HOLD) - 1)
            item["frame_fallback_last_index"] = last
            fallback.append(item)
            continue
        start, end = map(_ass_time, interval)
        for layer, name in ((10, "terra"), (20, "acoustic")):
            if item[name] is not None:
                lines.append(f"Dialogue: {layer},{start},{end},Text,,0,0,0,,{_drawing(item[name])}\n")
        if item["terra"] is not None and item.get("uncertainty_px") is not None and item["uncertainty_px"]>0:
            outline=dict(item["terra"],radius=item["uncertainty_px"])
            lines.append(f"Dialogue: 11,{start},{end},Text,,0,0,0,,{_drawing(outline,outline=True)}\n")
        y = max(26, round(width / 100) * 2) + (lane % 4) * max(15, round(width / 80))
        lines.append(f"Dialogue: 30,{start},{end},Text,,0,0,0,,{{\\an7\\pos(8,{y})}}{item['label']}\n")
    Path(path).write_text("".join(lines))
    return fallback


def _fallback_filters(events):
    filters = []
    for item in events:
        enable = f"between(n,{item['frame_index']},{item['frame_fallback_last_index']})"
        for name in ("terra", "acoustic"):
            mark = item[name]
            if mark:
                r = mark["radius"]
                rgb = "0x" + mark["color"][4:6] + mark["color"][2:4] + mark["color"][:2]
                filters.append(f"drawbox=x={mark['x']-r}:y={mark['y']-r}:w={2*r+1}:h={2*r+1}:color={rgb}:t=fill:enable='{enable}'")
        if item["terra"] is not None and item.get("uncertainty_px") is not None and item["uncertainty_px"]>0:
            mark=item["terra"];r=item["uncertainty_px"]
            condition=f"between(N,{item['frame_index']},{item['frame_fallback_last_index']})*between((X-{mark['x']})^2+(Y-{mark['y']})^2,{max(0,r-1)**2},{(r+1)**2})"
            colors=[int(mark["color"][i:i+2],16) for i in (4,2,0)]
            filters.append("geq="+":".join(f"{name}='if({condition},{color},{name}(X,Y))'" for name,color in zip("rgb",colors)))
        text = item["label"].replace("\\", "\\\\").replace("'", "").replace(":", "\\:").replace("%", "")
        filters.append(f"drawtext=text='{text}':x=8:y=h-20:fontsize=12:fontcolor=white:borderw=1:enable='{enable}'")
    return filters


def _source_audio(metadata, output, work, ffmpeg, threads, source_hashes):
    wavs = metadata["wav_files"]
    if not wavs:
        raise ValueError("No source WAVs")
    fs = int(wavs[0]["fs"])
    total = 0
    for info in wavs:
        with wave.open(str(info["path"]), "rb") as stream:
            if (stream.getframerate() != fs or int(info["fs"]) != fs or stream.getnchannels() < 5
                    or stream.getnframes() != int(info["samples"]) or int(info["start_sample"]) != total):
                raise ValueError("WAV sample/channel/count metadata is not gapless physical5 PCM")
        total += int(info["samples"])
    key = hashlib.sha256(json.dumps([(info["path"], source_hashes[str(Path(info["path"]).resolve())]) for info in wavs]).encode()).hexdigest()[:20]
    cache = Path(output) / "cache"
    cache.mkdir(exist_ok=True)
    raw = cache / (key + ".ch5_48000.f32")
    seal = raw.with_suffix(".json")
    if raw.exists() and seal.exists():
        saved = _json(seal)
        if saved["sha256"] == _sha(raw):
            return raw, saved
    listing = Path(work) / "recorders.ffconcat"
    listing.write_text("ffconcat version 1.0\n" + "".join(f"file {_quote_file(info['path'])}\n" for info in wavs))
    temporary = raw.with_suffix(".partial")
    _run([ffmpeg, "-v", "error", "-nostdin", "-y", "-threads", threads, "-filter_threads", threads,
          "-f", "concat", "-safe", "0", "-i", listing, "-map", "0:a:0", "-af", "pan=mono|c0=c4,aresample=48000:async=0",
          "-ac", "1", "-ar", RATE, "-c:a", "pcm_f32le", "-f", "f32le", temporary], work, "reference_decode")
    if temporary.stat().st_size % 4:
        raise ValueError("Partial reference float cache")
    count = temporary.stat().st_size // 4
    expected = round(total * RATE / fs)
    if abs(count - expected) > 1:
        raise ValueError(f"48kHz decode count {count} differs from gapless expected {expected}")
    os.replace(temporary, raw)
    saved = dict(sha256=_sha(raw), samples=count, expected_rounded_samples=expected,
                 source_samples=total, source_fs=fs, physical_channel=5,
                 review_fs=RATE, decode="single gapless concatenated-WAV resampling; FFmpeg aresample async=0")
    _write_json(seal, saved)
    return raw, saved


def _validate_policy_bridge(original_dir, fit_dir):
    """Call the explicit validator, never substitute derived labels for Terra.
    Lazy import allows consistent-annotation exports to remain independent of
    the optional user-confirmed uncertainty policy implementation.
    """
    helper = Path(__file__).with_name("fit_annotation_policy.py")
    if not helper.is_file():
        raise ValueError("Model does not bind current annotations; explicit fit annotation policy validator is unavailable")
    specification = importlib.util.spec_from_file_location("_video_export_fit_annotation_policy",helper)
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    try:
        validated = module.validate_fit_annotation_policy(original_dir,fit_dir)
    except Exception as exception:
        raise ValueError(f"Model does not bind current annotations and explicit policy validation failed: {exception}") from exception
    if validated.get("valid") is not True:
        raise ValueError("Explicit annotation policy validator did not accept the bridge")
    return validated, helper


def write_affine_audio(raw, destination, duration_s, offset_s, scale, *, chunk_samples=48000):
    """Bounded Float64 linear interpolation from readonly mono48k float32 cache.
    Write lossless IEEE Float32 WAV; samples outside the recorder are zero.
    x=((video_sample/48000-offset)/scale)*48000. No integer clock approximation.
    """
    count = math.ceil(duration_s * RATE - 1e-7)
    if count * 4 > 0xFFFFFFFF - 48:
        raise ValueError("Review audio exceeds RIFF size; RF64 requires a separate implementation")
    source = np.memmap(raw, dtype="<f4", mode="r")
    with Path(destination).open("wb") as stream:
        stream.write(struct.pack("<4sI4s4sIHHIIHH", b"RIFF", 48+count*4, b"WAVE", b"fmt ", 16,
                                 3, 1, RATE, RATE*4, 4, 32))
        stream.write(struct.pack("<4sII4sI", b"fact", 4, count, b"data", count*4))
        for start in range(0, count, chunk_samples):
            samples = np.arange(start, min(count, start+chunk_samples), dtype=np.float64)
            x = ((samples / RATE - offset_s) / scale) * RATE
            left = np.floor(x).astype(np.int64)
            fraction = x - left
            values = np.zeros(len(x), dtype=np.float64)
            domain = (x >= 0) & (x <= len(source)-1)
            valid = domain & (left >= 0) & (left < len(source))
            values[valid] += (1-fraction[valid]) * source[left[valid]]
            valid = domain & (left+1 >= 0) & (left+1 < len(source))
            values[valid] += fraction[valid] * source[left[valid]+1]
            stream.write(values.astype("<f4").tobytes())
        stream.flush()
        os.fsync(stream.fileno())
    del source
    return dict(samples=count, fs=RATE, encoding="IEEE Float32 mono WAV", duration_s=count/RATE,
                interpolation="chunked NumPy Float64 linear interpolation; Float32 review PCM; zero padding outside recording",
                semantics="audible-band48k review only; scientific source recordings unchanged")


def _load(device, model, output, ffmpeg, ffprobe, work,*,annotation_only=False):
    metadata = _json(device / "prepare_metadata.json")
    sync = _json(device / "synchronization.json")
    events = _rows(device / "events.csv", "event_id")
    frames = _rows(device / "frames.csv", "frame_id")
    annotations = _rows(device / "annotations.csv", "frame_id")
    predictions = {} if annotation_only else _rows(model / "predictions.csv", "event_id")
    calibration = None if annotation_only else _json(model / "calibration.json")
    inputs = [device / name for name in ("prepare_manifest.json", "prepare_metadata.json", "synchronization.json", "events.csv", "frames.csv", "annotations.csv", "annotations.complete.json")]
    if not annotation_only:
        inputs += [model / "calibration.json", model / "predictions.csv"]
    for name in ("sync_eligibility.json", "human_review_selection.json"):
        if (device / name).exists():
            inputs.append(device / name)
    if model is not None and (model / "report.json").exists():
        inputs.append(model / "report.json")
    if model is not None and (model / "fit_report.json").exists():
        inputs.append(model / "fit_report.json")
    sources = [Path(info["path"]).resolve() for info in metadata["wav_files"] + metadata["videos"]]
    hashes = {str(path.resolve()): _sha(path) for path in dict.fromkeys(inputs + sources)}
    prepared = _json(device / "prepare_manifest.json")
    recorded_sources = {str(Path(row["path"]).resolve()):row["sha256"] for row in prepared["request"]["inputs"]}
    if any(recorded_sources.get(str(path)) != hashes[str(path)] for path in sources):
        raise ValueError("Source media bytes differ from authenticated preparation")
    if set(annotations) != set(frames):
        raise ValueError("Annotations do not cover exactly the original frame identities")
    seal = _json(device / "annotations.complete.json")
    if seal["annotations_sha256"] != hashes[str((device / "annotations.csv").resolve())] or seal["frames_sha256"] != hashes[str((device / "frames.csv").resolve())]:
        raise ValueError("Annotation completion is stale")
    bindings = {} if annotation_only else calibration.get("hashes", {})
    model_completion = bindings.get("annotations.complete.json",bindings.get("annotations_complete_sha256"))
    if not annotation_only and model_completion is None:
        raise ValueError("Model lacks annotation-completion binding")
    policy_bridge = None
    checks=() if annotation_only else (("annotations.csv", "annotations_sha256"), ("events.csv", "events_sha256"), ("frames.csv", "frames_sha256"))
    for name, alternate in checks:
        expected = bindings.get(name, bindings.get(alternate))
        if expected is None or expected != hashes[str((device / name).resolve())]:
            if name != "annotations.csv" or expected is None:
                raise ValueError(f"Model does not bind current {name}; refusing stale/revised model")
            validated,helper = _validate_policy_bridge(device,model.parent)
            derived = model.parent
            extra = [derived/"annotations.csv",derived/"annotations.complete.json",derived/"annotation_fit_policy.json",
                     derived/"prepare_manifest.json",derived/"human_review_selection.json",helper,helper.with_name("annotation_revision.py")]
            for path in extra:
                hashes[str(path.resolve())] = _sha(path)
            actual = {"original_annotations_sha256":hashes[str((device/"annotations.csv").resolve())],
                      "fit_annotations_sha256":hashes[str((derived/"annotations.csv").resolve())],
                      "original_completion_sha256":hashes[str((device/"annotations.complete.json").resolve())],
                      "fit_completion_sha256":hashes[str((derived/"annotations.complete.json").resolve())],
                      "policy_sha256":hashes[str((derived/"annotation_fit_policy.json").resolve())],
                      "helper_sha256":hashes[str(helper.resolve())],
                      "clone_helper_sha256":hashes[str(helper.with_name("annotation_revision.py").resolve())],
                      "fit_prepare_manifest_sha256":hashes[str((derived/"prepare_manifest.json").resolve())],
                      "original_prepare_manifest_sha256":hashes[str((device/"prepare_manifest.json").resolve())],
                      "frames_sha256":hashes[str((device/"frames.csv").resolve())]}
            if any(validated.get(key)!=value for key,value in actual.items()) or expected!=actual["fit_annotations_sha256"]:
                raise ValueError("Explicit annotation policy/model/source hashes do not agree")
            if model_completion!=actual["fit_completion_sha256"]:
                raise ValueError("Model does not bind the validated fit annotation completion")
            selection = device/"human_review_selection.json"
            if not selection.is_file() or validated.get("selection_sha256")!=_sha(selection):
                raise ValueError("Explicit annotation policy human selection does not agree")
            hashes[str(selection.resolve())] = _sha(selection)
            policy_bridge = validated
    if not annotation_only and policy_bridge is None and model_completion!=hashes[str((device/"annotations.complete.json").resolve())]:
        raise ValueError("Model does not bind current annotation completion")
    report_path = None if annotation_only else model / ("report.json" if (model / "report.json").exists() else "fit_report.json")
    if report_path is not None and report_path.exists():
        report = _json(report_path)
        expected = report.get("outputs_sha256", {}).get("predictions.csv", report.get("predictions_sha256"))
        if expected != hashes[str((model / "predictions.csv").resolve())]:
            raise ValueError("Model prediction digest mismatch")
        expected_model = report.get("outputs_sha256",{}).get("calibration.json",report.get("calibration_sha256"))
        if expected_model is not None and expected_model != hashes[str((model/"calibration.json").resolve())]:
            raise ValueError("Model calibration digest mismatch")
    eligibility = _json(device / "sync_eligibility.json") if (device / "sync_eligibility.json").exists() else None
    if eligibility:
        if eligibility["events_sha256_after"] != hashes[str((device / "events.csv").resolve())]:
            raise ValueError("Sync eligibility events binding is stale")
        saved = eligibility["input_binding"]["saved_inputs"]
        for name, filename in (("metadata", "prepare_metadata.json"), ("synchronization", "synchronization.json")):
            if saved[name]["sha256"] != hashes[str((device / filename).resolve())]:
                raise ValueError("Exact-base synchronization input is stale")
    timing = _timing(metadata, sync, eligibility)
    videos = []
    expected_pts = []
    for chapter, (info, base) in enumerate(zip(metadata["videos"], timing["bases"])):
        probe = _probe(info["path"], ffprobe, work, f"probe_source_{chapter}")
        if len(probe["pts"]) != len(info["pts"]) or any(abs(a-float(b)) > 1.1e-6 for a,b in zip(probe["pts"],info["pts"])):
            raise ValueError("Source PTS/frame count differs from saved identities")
        if (probe["width"],probe["height"]) != (int(info["width"]),int(info["height"])):
            raise ValueError("Source resolution changed")
        local = [x-float(info.get("audio_start_s",0)) for x in probe["pts"]]
        expected_pts.extend(base+x for x in local)
        videos.append(dict(info=info, probe=probe, local_pts=local, base=base, events=[]))
    if any(b <= a for a,b in zip(expected_pts,expected_pts[1:])):
        raise ValueError("Declared chapter bases overlap source frame timestamps; refusing frame dropping/reordering")
    if len({(v["probe"]["width"],v["probe"]["height"]) for v in videos}) != 1:
        raise ValueError("Original-resolution export requires one resolution per device")
    accepted = None
    if (device / "human_review_selection.json").exists():
        accepted = set(_json(device / "human_review_selection.json")["accepted_frame_ids"])
    qc = {row["event_id"]:row.get("final_sync_status") for row in eligibility["review"]["events"]} if eligibility else {}
    by_video = {str(Path(v["info"]["path"]).resolve()):v for v in videos}
    missing = []
    for identity, event in events.items():
        frame = frames.get(event.get("frame_id"))
        if frame is None:
            missing.append(dict(event_id=identity,status="no_source_frame_identity"))
            continue
        video = by_video[str(Path(frame["video_file"]).resolve())]
        index = int(frame["frame_index"])
        if not 0 <= index < len(video["local_pts"]) or abs(float(frame["frame_pts_s"])-video["probe"]["pts"][index]) > 1.1e-6:
            raise ValueError("Event frame index/PTS identity mismatch")
        if event.get("video_file") != frame["video_file"] or int(event["frame_index"]) != index:
            raise ValueError("Event/frame source identity differs")
        annotation = annotations.get(frame["frame_id"], dict(visibility="absent",medium="unknown"))
        marks = _event_marks(event, annotation, predictions.get(identity), video["probe"]["width"],video["probe"]["height"],
                             accepted=accepted is None or frame["frame_id"] in accepted,sync_status=qc.get(identity),annotation_only=annotation_only)
        marks.update(event_id=identity,frame_id=frame["frame_id"],frame_index=index,
                     source_video=frame["video_file"],source_pts_s=video["probe"]["pts"][index],
                     global_pts_s=video["base"]+video["local_pts"][index])
        video["events"].append(marks)
    duration = max(v["base"] + v["probe"]["duration"]-float(v["info"].get("audio_start_s",0)) for v in videos)
    return dict(metadata=metadata, hashes=hashes, timing=timing, videos=videos,
                expected_pts=expected_pts, duration=duration, missing=missing,policy_bridge=policy_bridge)


def export_video(device_dir, model_dir=None, output_dir=None, *, ffmpeg="/usr/bin/ffmpeg", ffprobe="/usr/bin/ffprobe",
                 resume=True, threads=2, crf=20, preset="fast",annotation_only=False):
    """Export a source-bound atomic generation; never modifies input/model data.
    Return a manifest containing paths to the MP4, lossless float-review WAV,
    raw/overlay per-event PNGs and exact frame/timing plan. Full encodes must be
    authorized by the caller. Imports never launch an export.
    """
    if annotation_only and model_dir is not None:
        raise ValueError("annotation_only uses original sealed annotations, not a model directory")
    if output_dir is None or (not annotation_only and model_dir is None):
        raise ValueError("output_dir is required; calibrated mode also requires model_dir")
    device=Path(device_dir).resolve();model=None if model_dir is None else Path(model_dir).resolve();output=Path(output_dir).resolve()
    if not 1 <= threads <= 2:
        raise ValueError("Review export supports one or two threads")
    if output == device or (model is not None and output == model):
        raise ValueError("Use a separate media output directory")
    output.mkdir(parents=True,exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".export-work-",dir=output))
    try:
        loaded = _load(device,model,output,ffmpeg,ffprobe,work,annotation_only=annotation_only)
        settings = dict(threads=threads,crf=crf,preset=preset,nominal_hold_s=HOLD,review_audio_fs=RATE,
                        exporter_sha256=_sha(__file__),ffmpeg=str(ffmpeg),ffprobe=str(ffprobe),numpy_version=np.__version__,
                        ffmpeg_sha256=_sha(shutil.which(str(ffmpeg)) or ffmpeg),ffprobe_sha256=_sha(shutil.which(str(ffprobe)) or ffprobe),annotation_only=annotation_only)
        manifest_path = output / "export_manifest.json"
        if resume and manifest_path.exists():
            prior = _json(manifest_path)
            if prior.get("status")=="complete" and prior.get("inputs_sha256")==loaded["hashes"] and prior.get("settings")==settings:
                if all((output/name).is_file() and _sha(output/name)==digest for name,digest in prior["outputs_sha256"].items()):
                    shutil.rmtree(work)
                    return prior
        stage = work / "generation"
        stage.mkdir();(stage/"frames"/"raw").mkdir(parents=True);(stage/"frames"/"overlay").mkdir(parents=True)
        raw, source_audio = _source_audio(loaded["metadata"],output,work,ffmpeg,threads,loaded["hashes"])
        audio = write_affine_audio(raw,stage/"synchronized_reference_ch5.wav",loaded["duration"],
                                   loaded["timing"]["offset_s"],loaded["timing"]["scale"])
        chapters = []
        frame_map = []
        for chapter, video in enumerate(loaded["videos"]):
            width,height = video["probe"]["width"],video["probe"]["height"]
            for event in video["events"]:
                event["display_hold_capped_by_chapter_end"] = video["local_pts"][event["frame_index"]]+HOLD > video["probe"]["duration"]-float(video["info"].get("audio_start_s",0))
            ass = stage / f"chapter_{chapter:02}.ass"
            fallback = _write_ass(ass,width,height,video["local_pts"],video["events"],annotation_only=annotation_only)
            decorated = f"ass=filename='{_filter_path(ass)}'"
            extra = _fallback_filters(fallback)
            if extra:
                decorated += "," + ",".join(extra)
            indices = sorted({row["frame_index"] for row in video["events"]})
            part = work / f"chapter_{chapter:02}.mp4"
            command = [ffmpeg,"-v","error","-nostdin","-y","-threads",threads,"-filter_complex_threads",1,
                       "-copyts","-i",video["info"]["path"]]
            normalize = f"settb=1/1000000,setpts=PTS-({float(video['info'].get('audio_start_s',0)):.9f})/TB"
            if indices:
                select = "+".join(f"eq(n,{i})" for i in indices)
                graph = (f"[0:v]{normalize},split=2[base][raw];[base]{decorated},split=2[v][overlay];"
                         f"[raw]select='{select}'[r];[overlay]select='{select}'[o]")
            else:
                graph = f"[0:v]{normalize},{decorated}[v]"
            command += ["-filter_complex",graph,"-map","[v]","-an","-c:v","libx264","-threads",threads,
                        "-crf",crf,"-preset",preset,"-pix_fmt","yuv420p","-fps_mode","passthrough",
                        "-enc_time_base","1:1000000","-video_track_timescale",1000000,part]
            if indices:
                for label in ("r","o"):
                    command += ["-map",f"[{label}]","-c:v","png","-threads",threads,"-fps_mode","passthrough",
                                "-start_number",0,work/f"{chapter}_{label}_%06d.png"]
            _run(command,work,f"encode_chapter_{chapter}")
            chapters.append(part)
            for number,index in enumerate(indices):
                raw_png=work/f"{chapter}_r_{number:06}.png";overlay_png=work/f"{chapter}_o_{number:06}.png"
                for event in (item for item in video["events"] if item["frame_index"]==index):
                    safe=hashlib.sha256(event["event_id"].encode()).hexdigest()[:8]
                    slug="".join(c if c.isalnum() or c in "-_." else "_" for c in event["event_id"])[:96]
                    name=f"ch{chapter:02}_f{index:08}_{slug}_{safe}.png"
                    shutil.copyfile(raw_png,stage/"frames"/"raw"/name)
                    shutil.copyfile(overlay_png,stage/"frames"/"overlay"/name)
                    frame_map.append(dict(event_id=event["event_id"],frame_id=event["frame_id"],source_video=event["source_video"],
                                          frame_index=index,source_pts_s=event["source_pts_s"],global_pts_s=event["global_pts_s"],
                                          raw=f"frames/raw/{name}",overlay=f"frames/overlay/{name}",acoustic_status=event["acoustic_status"],
                                          terra_visibility=event["visibility"],uncertainty_px=event["uncertainty_px"],medium=event["medium"],user_excluded=event["user_excluded"]))
        listing=work/"chapters.ffconcat"
        listing.write_text("ffconcat version 1.0\n"+"".join(f"file {_quote_file(part)}\n"+(f"duration {video['info']['audio_duration_s']}\n" if i<len(chapters)-1 else "") for i,(part,video) in enumerate(zip(chapters,loaded["videos"]))))
        video_path=stage/("terra_review_synchronized.mp4" if annotation_only else "synchronized_annotated.mp4")
        _run([ffmpeg,"-v","error","-nostdin","-y","-threads",threads,"-f","concat","-safe","0","-i",listing,
              "-i",stage/"synchronized_reference_ch5.wav","-map","0:v:0","-map","1:a:0","-c:v","copy",
              "-c:a","aac","-threads",threads,"-b:a","128k","-ar",RATE,"-ac",1,"-video_track_timescale",1000000,
              "-movflags","+faststart",video_path],work,"mux_full_video")
        probe=_probe(video_path,ffprobe,work,"probe_output")
        expected=loaded["expected_pts"]
        if len(probe["pts"])!=len(expected):
            raise ValueError(f"Output frame count {len(probe['pts'])} != source {len(expected)}; partial output rejected")
        maximum_error=max(abs(a-b) for a,b in zip(probe["pts"],expected))
        if maximum_error>2.1e-6:
            raise ValueError(f"Output PTS differs by {maximum_error}s; source timing not preserved")
        for path,digest in loaded["hashes"].items():
            if _sha(path)!=digest:
                raise ValueError(f"Input changed during export: {path}")
        _write_json(stage/"frame_manifest.json",dict(events=frame_map,unmapped=loaded["missing"]))
        _write_json(stage/"marker_plan.json",dict(timing=loaded["timing"],chapters=[dict(source=v["info"]["path"],events=v["events"]) for v in loaded["videos"]]))
        (stage/"logs").mkdir()
        for log in work.glob("*.log"):
            shutil.copyfile(log,stage/"logs"/log.name)
        generation="media_"+uuid.uuid4().hex[:12]
        os.replace(stage,output/generation)
        files=sorted(path for path in (output/generation).rglob("*") if path.is_file())
        manifest=dict(schema_version="clicker-video-export/1.0",status="complete",inputs_sha256=loaded["hashes"],settings=settings,
                      mode="annotation_only" if annotation_only else "calibrated_comparison",
                      generation=generation,video_path=str(output/generation/video_path.name),
                      reference_audio_path=str(output/generation/"synchronized_reference_ch5.wav"),
                      source_frames=len(expected),output_frames=len(probe["pts"]),maximum_pts_error_s=maximum_error,
                      video_duration_s=loaded["duration"],chapter_bases_s=loaded["timing"]["bases"],timing=loaded["timing"],
                      review_audio=audio,source_audio_cache=source_audio,event_frame_pairs=len(frame_map),unmapped_events=loaded["missing"],
                      annotation_fit_policy_bridge=loaded["policy_bridge"],
                      display_semantics="Exact original target frame, never preceding frame; nominal120ms fixed-position hold capped at chapter end, ASS endpoint quantization<10ms; no trajectory interpolation. Rare sub10ms intervals use decoded-index filters.",
                      precision="Source frame index/PTS verified; microsecond encoding and concat timestamps, maximum2.1us acceptance tolerance; no CFR conversion, frame drop, frame duplication, shortest trimming, integer asetrate or atempo.",
                      outputs_sha256={str(path.relative_to(output)):_sha(path) for path in files})
        _write_json(manifest_path,manifest)
        shutil.rmtree(work)
        return manifest
    except BaseException as exception:
        _log_failure(output,"export_video",exception)
        raise


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("device_dir");parser.add_argument("model_dir",nargs="?");parser.add_argument("output_dir",nargs="?")
    parser.add_argument("--annotation-only",action="store_true")
    parser.add_argument("--ffmpeg",default="/usr/bin/ffmpeg");parser.add_argument("--ffprobe",default="/usr/bin/ffprobe")
    parser.add_argument("--threads",type=int,default=2);parser.add_argument("--crf",type=int,default=20)
    parser.add_argument("--preset",default="fast");parser.add_argument("--no-resume",action="store_true")
    args=vars(parser.parse_args(argv));args["resume"]=not args.pop("no_resume")
    if args["annotation_only"] and args["output_dir"] is None:
        args["output_dir"]=args["model_dir"];args["model_dir"]=None
    print(json.dumps(export_video(**args),indent=2))


if __name__=="__main__":
    main()
