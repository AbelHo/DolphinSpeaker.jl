#!/usr/bin/env python3
"""Plot saved calibration results and cached waveforms; never fit or reject events."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
from pathlib import Path


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def number(row, key):
    try:
        return float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return float("nan")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def diagnostics(directory, model_dir=None, waveform_ids=None, output_dir=None):
    directory = directory.resolve()
    model_dir = directory if model_dir is None else model_dir.resolve()
    output_dir = model_dir / "diagnostics" if output_dir is None else output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    calibration = json.loads((model_dir / "calibration.json").read_text())
    report_path = model_dir / ("fit_report.json" if model_dir == directory else "report.json")
    report = json.loads(report_path.read_text())
    predictions = [row for row in read_csv(model_dir / "predictions.csv")
                   if row["split"] in {"train", "heldout"}]
    events = {row["event_id"]: row for row in read_csv(directory / "events.csv")}
    features, snippets = [], {}
    for row in predictions:
        event = events[row["event_id"]]
        ns, nc = int(event["snippet_samples"]), int(event["snippet_channels"])
        with (directory / event["snippet_file"]).open("rb") as f:
            f.seek(int(event["snippet_byte_offset"]))
            raw = f.read(ns * nc * 8)
        if len(raw) != ns * nc * 8 or hashlib.sha256(raw).hexdigest() != event["snippet_sha256"]:
            raise ValueError(f"Cached snippet checksum differs: {row['event_id']}")
        data = np.frombuffer(raw, dtype="<f8").reshape((ns, nc), order="F")
        fs = float(event["fs"])
        freq = np.fft.rfftfreq(ns, 1/fs)
        power = np.abs(np.fft.rfft(data * np.hanning(ns)[:, None], axis=0))**2
        average = power.sum(axis=1)
        total = float(average.sum())
        feature = dict(event_id=row["event_id"], split=row["split"],
                       recording_time_s=float(event["recording_time_s"]),
                       energy_centroid_khz=float(freq @ average / total / 1000) if total else None,
                       energy_above_30khz_fraction=float(average[freq > 30000].sum()/total) if total else None,
                       gcc_closure_samples=number(event, "gcc_closure_samples"),
                       evaluation_error_px=number(row, "evaluation_error_px"),
                       deployment_error_px=number(row, "final_error_px"),
                       snippet_sha256=event["snippet_sha256"])
        features.append(feature)
        snippets[row["event_id"]] = (event, data, freq, power)
    with (output_dir / "acoustic_descriptors.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(features[0]) if features else ["event_id"])
        writer.writeheader()
        writer.writerows(features)

    w, h = calibration["width"], calibration["height"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    colors = {"train": "#1971c2", "heldout": "#e8590c"}
    evaluation_label = "reused evaluation" if model_dir != directory else "held-out evaluation"
    for split in ("train", "heldout"):
        rows = [row for row in predictions if row["split"] == split]
        label = "training" if split == "train" else evaluation_label
        x = np.array([number(row, "annotation_x") for row in rows])
        y = np.array([number(row, "annotation_y") for row in rows])
        px = np.array([number(row, "evaluation_x") for row in rows])
        py = np.array([number(row, "evaluation_y") for row in rows])
        valid = np.isfinite(px) & np.isfinite(py)
        axes[0, 0].scatter(x, y, s=28, facecolors="none", edgecolors=colors[split], label=label+" annotations")
        axes[0, 0].scatter(px[valid], py[valid], s=28, marker="x", color=colors[split], label=label+" predictions")
        for a, b, c, d in zip(x[valid], y[valid], px[valid], py[valid]):
            axes[0, 0].plot([a, c], [b, d], color=colors[split], alpha=.25, lw=.7)
        times = [float(events[row["event_id"]]["recording_time_s"]) for row in rows]
        errors = np.array([number(row, "evaluation_error_px") for row in rows])
        finite = np.isfinite(errors)
        axes[0, 1].scatter(np.array(times)[finite], np.maximum(errors[finite], .1), color=colors[split], s=25, label=f"{label}: {finite.sum()}/{len(rows)} projected")
        desc = [item for item in features if item["split"] == split]
        fe = np.array([item["evaluation_error_px"] for item in desc])
        fc = np.array([item["energy_centroid_khz"] for item in desc], dtype=float)
        valid = np.isfinite(fe) & np.isfinite(fc)
        axes[1, 1].scatter(fc[valid], np.maximum(fe[valid], .1), color=colors[split], s=25, label=label)
    axes[0, 0].set(xlim=(0,w), ylim=(h,0), xlabel="Image x (px)", ylabel="Image y (px)", title="Training-fit projection within image bounds")
    axes[0, 0].set_aspect("equal", adjustable="box")
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].set(xlabel="Recorder time (s)", ylabel="Pixel error (log scale)", title="Training-fit errors; invalid rays remain unscored")
    axes[0, 1].set_yscale("log")
    axes[0, 1].legend(fontsize=8)
    profiles = report.get("diagnostics_training_only", {}).get("training_sound_speed_profiles", [])
    if profiles:
        speed = [p["fixed_sound_speed_m_s"] for p in profiles]
        objective = np.array([p["training_objective"] for p in profiles])
        relative = (objective-objective.min())/max(abs(float(objective.min())),1e-12)*100
        axes[1, 0].plot(speed, relative, "o-", color="#7048e8")
        axes[1, 0].set(xlabel="Fixed sound speed (m/s)", ylabel="Objective above best (%)", title="Training-only nuisance-refitted speed profile")
    else:
        axes[1, 0].text(.5,.5,"No sound-speed profile available",ha="center",transform=axes[1,0].transAxes)
    axes[1, 1].set(xlabel="Windowed acoustic energy centroid (kHz)", ylabel="Pixel error (log scale)", title="Post-fit acoustic diagnostic; no event rejection")
    axes[1, 1].set_yscale("log")
    axes[1, 1].legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.2)
    outside = sum(not (0 <= number(row,"evaluation_x") < w and 0 <= number(row,"evaluation_y") < h)
                  for row in predictions if np.isfinite(number(row,"evaluation_x")))
    invalid = sum(not np.isfinite(number(row,"evaluation_x")) for row in predictions)
    fig.suptitle(f"{calibration['device']} / {calibration['method']} — {len(predictions)} annotated events\n{outside} projections outside displayed image; {invalid} invalid. Sound-speed profile is not a confidence interval.",fontsize=11)
    fig.savefig(output_dir / "calibration_diagnostics.png", dpi=180)
    fig.savefig(output_dir / "calibration_diagnostics.svg")
    plt.close(fig)

    if waveform_ids is None:
        ordered = sorted(predictions, key=lambda r: float(events[r["event_id"]]["recording_time_s"]))
        indices = np.linspace(0,len(ordered)-1,min(4,len(ordered))).round().astype(int) if ordered else []
        waveform_ids = [ordered[i]["event_id"] for i in indices]
        selection = "uniform-time examples among the frozen eligible population"
    else:
        selection = "explicit caller-selected examples; post-fit descriptive inspection only"
    if not 1 <= len(waveform_ids) <= 8:
        raise ValueError("Specify one to eight waveform IDs")
    fig, axes = plt.subplots(len(waveform_ids),2,figsize=(12,2.6*len(waveform_ids)),squeeze=False,constrained_layout=True)
    for i, event_id in enumerate(waveform_ids):
        event,data,freq,power = snippets[event_id]
        t = (np.arange(len(data))-int(event["snippet_peak_index"]))/float(event["fs"])*1000
        for channel in range(3):
            scale = np.max(np.abs(data[:,channel]))
            axes[i,0].plot(t,data[:,channel]/scale if scale else data[:,channel],lw=.7,label=f"physical {channel+3}")
            p = power[:,channel]
            db = 10*np.log10(np.maximum(p/max(float(p.max()),1e-300),1e-10))
            axes[i,1].plot(freq/1000,db,lw=.8)
        axes[i,0].set(title=event_id,xlabel="Time from detected reference peak (ms)",ylabel="Per-channel normalized amplitude")
        axes[i,1].set(xlabel="Frequency (kHz)",ylabel="Relative windowed power (dB)",ylim=(-80,2),xlim=(0,float(event["fs"])/2000))
        axes[i,0].legend(fontsize=7,ncol=3)
        for ax in axes[i]:ax.grid(alpha=.2)
    fig.suptitle("Exact prepared HP1 kHz snippets — no source attribution inferred",fontsize=11)
    fig.savefig(output_dir / "waveform_examples.png",dpi=180)
    fig.savefig(output_dir / "waveform_examples.svg")
    plt.close(fig)
    summary = dict(device=calibration["device"],method=calibration["method"],n_events=len(predictions),
                   output_dir=str(output_dir),waveform_selection=selection,waveform_event_ids=waveform_ids,
                   calibration_sha256=digest(model_dir/"calibration.json"),predictions_sha256=digest(model_dir/"predictions.csv"),
                   annotations_sha256=digest(directory/"annotations.csv"),events_sha256=digest(directory/"events.csv"),
                   report_sha256=digest(report_path),source_sha256=digest(Path(__file__)),
                   python_version=platform.python_version(),
                   numpy_version=np.__version__,matplotlib_version=matplotlib.__version__,
                   limitations="Descriptive plots only. No event rejection, source classification, or model fitting. Log plots floor errors at 0.1 px for display. Spectra are relative PCM power, not calibrated acoustic pressure.")
    (output_dir/"diagnostics.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("--model-dir",type=Path)
    parser.add_argument("--output-dir",type=Path)
    parser.add_argument("--waveform-ids",nargs="+")
    args=parser.parse_args()
    print(json.dumps(diagnostics(args.directory,args.model_dir,args.waveform_ids,args.output_dir)))
