#!/usr/bin/env python3
"""Package committed code, frozen source versions and numerical calibration evidence.

Explicit inputs only; recordings, media, runtime environments and credentials are
excluded. Recorded commands are data and are never executed. Stdlib only.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import ctypes
from datetime import datetime, timezone
import errno
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tempfile


SCHEMA = "clicker-reproducibility-v1"
POINTER_SCHEMA = "clicker-reproducibility-pointer-v1"
MODULE_ROOT = "src/clicker_calibration/"
SKIP_DIRS = {".git", ".julia", ".venv", "venv", "runtime", "julia_depot", "scratch", "test_scratch",
             "cache", "__pycache__", "frames", "annotation_review", "fit_review", "refined_review", "videos"}
TEXT_SUFFIXES = {".json", ".jsonl", ".csv", ".tsv", ".txt", ".md", ".rst", ".log", ".jl", ".py", ".toml", ".sh"}
CODE_SUFFIXES = {".jl", ".py", ".toml", ".md", ".rst", ".txt", ".json", ".csv", ".tsv", ".sh", ".html", ".js", ".css", ".yml", ".yaml", ".bib", ".c", ".h", ".cpp", ".hpp"}
MEDIA_SUFFIXES = {".wav", ".flac", ".mp4", ".mkv", ".mov", ".avi", ".mp3", ".ogg", ".aac", ".pcm"}
SECRET = re.compile(rb"(?:\b(?:bot)?\d{6,}:[A-Za-z0-9_-]{25,}\b|\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{24,}\b|\bAKIA[A-Z0-9]{16}\b|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----)")
SOURCE_KEYS = {name.replace(".jl", "_source_sha256"): name for name in
               ("physical_fit.jl", "prepare.jl", "srp_review.jl", "method_refinement.jl", "support.jl", "workflow.jl")}


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def fingerprint(path):
    s = Path(path).stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def git(repository, *arguments, check=True):
    # Only these internally constructed Git inspection arguments are executed.
    result = subprocess.run(["git", "-C", str(repository), *arguments], capture_output=True, timeout=30)
    if check and result.returncode:
        raise RuntimeError("Git inspection failed: " + result.stderr.decode(errors="replace")[:1000])
    return result


def label(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError("Input labels must be simple directory names")
    return value


def relative(value):
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or not path.parts or "\\" in value:
        raise ValueError("Unsafe bundle-relative path")
    return path


def credential_name(path):
    name = Path(path).name.lower()
    return (name == ".env" or name.startswith(".env.") or name in {"credentials", "credentials.json", "id_rsa", "id_ed25519", "token.json"}
            or name.startswith("client_secret") or Path(name).suffix in {".pem", ".key", ".p12", ".pfx"})


def code_paths(repository, head):
    tracked = git(repository, "ls-tree", "-r", "--name-only", "-z", head).stdout.decode().split("\0")
    selected = []
    for name in filter(None, tracked):
        path = relative(name)
        suffix = path.suffix.lower()
        root = path.parts[0]
        include = (name in {"Project.toml", "Manifest.toml", "LICENSE", "LICENSE.md", "README.md"}
                   or root in {"src", "tools", "test", "tests"}
                   or (root == "docs" and suffix in TEXT_SUFFIXES - {".log", ".jsonl"}))
        code_kind = suffix in CODE_SUFFIXES or name in {"LICENSE", "LICENSE.md"} or (root in {"test", "tests"} and suffix in {".png", ".jpg", ".svg"})
        if include and code_kind and suffix not in MEDIA_SUFFIXES and not credential_name(name):
            selected.append(name)
    if not {"Project.toml", "tools/package_clicker_reproducibility.py"} <= set(selected):
        raise ValueError("Commit Project.toml and the reproducibility packager before packaging")
    return sorted(selected)


def walk_artifacts(root, *, evidence=False, excluded=()):
    """Do not descend symlinks, environments, original frame/media trees or scratch."""
    root = Path(root).resolve(strict=True)
    if root.name.lower() in SKIP_DIRS or credential_name(root):
        raise ValueError("Runtime/private/media directories cannot be explicit artifact roots")
    evidence_root = root.name in {"jobs", "logs", "process_logs", "recovery", "source_snapshot"} or (root / "result.json").is_file()
    for directory, dirs, files in os.walk(root, followlinks=False):
        base = Path(directory)
        dirs[:] = sorted(name for name in dirs if name.lower() not in SKIP_DIRS and not name.startswith(".")
                         and not name.startswith("media_") and not (base / name).is_symlink()
                         and not any((base / name).resolve() == p or p in (base / name).resolve().parents for p in excluded))
        for name in sorted(files):
            path = base / name
            rel = path.relative_to(root)
            if evidence and not evidence_root and len(rel.parts) > 1 and rel.parts[0] not in {"jobs", "logs", "process_logs", "recovery", "source_snapshot"}:
                continue
            suffix = path.suffix.lower()
            if (path.is_symlink() or credential_name(path) or name.startswith(".") or suffix in MEDIA_SUFFIXES
                    or any(path.resolve() == p or p in path.resolve().parents for p in excluded)):
                continue
            if suffix in TEXT_SUFFIXES or name == "event_snippets.f64" or ("diagnostics" in rel.parts and suffix in {".png", ".svg", ".pdf"}):
                yield path, rel.as_posix()


def active_job_roots(root):
    """Avoid archiving mutable launcher receipts/logs, including this packager's own job."""
    root = Path(root).resolve(strict=True)
    candidates = [root] if (root / "result.json").is_file() else []
    jobs = root if root.name == "jobs" else root / "jobs"
    if jobs.is_dir():
        candidates.extend(p for p in jobs.iterdir() if p.is_dir() and not p.is_symlink() and (p / "result.json").is_file())
    active = []
    for path in candidates:
        status = read_json(path / "result.json").get("status")
        if status not in {"complete", "failed", "cancelled", "interrupted", "terminated"}:
            active.append((path, status))
    return active


def source_requirements(document):
    result = {}
    for mapping in (document.get("hashes", {}), document.get("estimator_spec", {})):
        for key, value in mapping.items():
            name = key[7:] if key.startswith("source/") else SOURCE_KEYS.get(key)
            if name is None:
                continue
            if Path(name).name != name or not name.endswith((".jl", ".py")) or not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
                raise ValueError("Malformed frozen source identity")
            path = MODULE_ROOT + name
            if path in result and result[path] != value:
                raise ValueError("Conflicting frozen source identities")
            result[path] = value
    return result


def recover_source(repository, head, path, expected, max_commits):
    commits = [head, *git(repository, "log", f"--max-count={max_commits}", "--format=%H", head, "--", path).stdout.decode().splitlines()]
    for commit in dict.fromkeys(commits):
        result = git(repository, "show", f"{commit}:{path}", check=False)
        if result.returncode == 0 and digest(result.stdout) == expected:
            return result.stdout, commit
    raise ValueError(f"Cannot recover frozen source {path} SHA {expected} within {max_commits} commits; refusing substitution")


def _publish_absent(source, destination):
    if sys.platform.startswith("linux"):
        libc = ctypes.CDLL(None, use_errno=True)
        rename = getattr(libc, "renameat2", None)
        if rename is None:
            raise OSError(errno.ENOSYS, "Atomic no-replace publication unavailable")
        rename.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint)
        if rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1):
            code = ctypes.get_errno()
            raise OSError(code, os.strerror(code), str(destination))
    elif os.name == "nt":
        os.rename(source, destination)
    else:
        raise OSError(errno.ENOSYS, "Atomic no-replace publication unsupported")


def write_file(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def verify_previous(root):
    pointer = root / "manifest.json"
    if not root.exists():
        return None, None
    if root.is_symlink() or not root.is_dir():
        raise ValueError("Reproducibility destination is not an owned regular directory")
    entries = {p.name for p in root.iterdir()}
    if not entries:
        return None, None
    if not pointer.is_file() or pointer.is_symlink():
        raise ValueError("Refusing a nonempty reproducibility directory without its ownership manifest")
    current = read_json(pointer)
    if current.get("schema_version") != POINTER_SCHEMA:
        raise ValueError("Foreign reproducibility manifest")
    generations = current.get("generations", {})
    if entries != {"manifest.json", *generations} or current.get("generation") not in generations or current.get("generation_manifest_sha256") != generations.get(current.get("generation")):
        raise ValueError("Unowned reproducibility files/generations; refusing overwrite")
    for name, expected in generations.items():
        if not re.fullmatch(r"generation_[0-9a-f]{20}", name):
            raise ValueError("Unsafe generation name")
        directory = root / name
        manifest_path = directory / "generation_manifest.json"
        if directory.is_symlink() or manifest_path.is_symlink() or not manifest_path.is_file() or sha(manifest_path) != expected:
            raise ValueError("Modified/foreign generation manifest")
        manifest = read_json(manifest_path)
        files = manifest["files_sha256"]
        actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
        if actual != {*files, "generation_manifest.json"}:
            raise ValueError("Unowned/missing generation files")
        for rel, expected_hash in files.items():
            path = directory / relative(rel)
            if path.is_symlink() or sha(path) != expected_hash:
                raise ValueError(f"Refusing externally edited reproducibility file: {rel}")
    return current, sha(pointer)


def device_summary(name, original, model):
    calibration = read_json(model / "calibration.json")
    report_name = "report.json" if (model / "report.json").is_file() else "fit_report.json"
    report = read_json(model / report_name)
    if str(calibration.get("device")) != name:
        raise ValueError("Latest device/model identity differs")
    if report_name == "report.json":
        outputs = report["outputs_sha256"]
        for rel in ("calibration.json", "predictions.csv", "evaluation.json", "selection.json"):
            if outputs.get(rel) != sha(model / rel):
                raise ValueError(f"Latest frozen model digest differs: {name}/{rel}")
        fit_root = model.parent
    else:
        for rel, key in (("calibration.json", "calibration_sha256"), ("predictions.csv", "predictions_sha256")):
            if report.get(key) != sha(model / rel):
                raise ValueError("Latest first-fit output digest differs")
        fit_root = model
    bindings = calibration.get("hashes", {})
    for filename, alternate, directory in (
            ("events.csv", "events_sha256", original), ("frames.csv", "frames_sha256", original),
            ("annotations.csv", "annotations_sha256", fit_root),
            ("annotations.complete.json", "annotations_complete_sha256", fit_root),
            ("event_snippets.f64", "event_snippets_sha256", fit_root)):
        expected = bindings.get(filename, bindings.get(alternate))
        if expected is not None and sha(directory / filename) != expected:
            raise ValueError(f"Latest model/input mapping differs: {name}/{filename}")
    annotations_path = original / "annotations.csv"
    seal = read_json(original / "annotations.complete.json")
    if seal.get("annotations_sha256") != sha(annotations_path) or seal.get("frames_sha256") != sha(original / "frames.csv"):
        raise ValueError("Latest original annotations are not currently sealed")
    with annotations_path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    selection = read_json(original / "human_review_selection.json") if (original / "human_review_selection.json").is_file() else None
    first_report = read_json(fit_root / "fit_report.json") if (fit_root / "fit_report.json").is_file() else report
    summary = dict(original_directory=str(original), model_directory=str(model), fit_input_directory=str(fit_root),
                   model_sha256=sha(model / "calibration.json"), report_sha256=sha(model / report_name),
                   method=calibration.get("method"), fit_status=calibration.get("fit_status"),
                   physical_parameters={key: calibration.get(key) for key in ("sound_speed_m_s", "relative_delay_offsets_s", "rotation_array_to_camera", "fov_degrees", "camera", "geometry_m", "width", "height")},
                   parameter_bound_hits=calibration.get("parameter_bound_hits_normalized"), software=calibration.get("software"),
                   original_annotations=dict(n_rows=len(rows), visibility_counts=dict(Counter(row["visibility"] for row in rows)),
                                             annotations_sha256=sha(annotations_path), completion_sha256=sha(original / "annotations.complete.json"),
                                             uncertainty="Original per-frame uncertainty_px strings are retained unchanged in the copied annotations.csv; they are not asserted to be calibrated confidence intervals."),
                   user_scope=selection, observation_policy=first_report.get("observation_policy"),
                   evaluation=dict(first_fit_heldout=first_report.get("heldout_metrics"),
                                   first_fit_label="Frozen first-fit outer-block evaluation under its original observation/split policy; not independent visual ground truth.",
                                   reused_after_method_selection=report.get("reused_evaluation_metrics"),
                                   reused_label="Reused descriptive evaluation, not untouched confirmatory validation.",
                                   deployment_descriptive=report.get("deployment_metrics_descriptive", report.get("final_all_eligible_metrics_descriptive")),
                                   deployment_label="In-sample deployment errors after all-eligible refit; not heldout accuracy."),
                   training_metrics=report.get("training_metrics"),
                   counts={key: report.get(key) for key in ("n_training", "n_reused_evaluation", "n_deployment")},
                   limits=report.get("limits", []))
    return summary, calibration, report, first_report, fit_root


def package_reproducibility(repository, results_dir, *, devices, run_roots=None, history_dirs=None,
                            max_git_commits=256, max_files=20000, max_bytes=2 * 1024**3):
    """Publish immutable reproducibility generations and an atomic active pointer.

    devices maps A/B to {original_dir, model_dir}; run_roots/history_dirs map
    explicit labels to directories. History is always marked superseded. Only
    tracked HEAD code/tests/docs are snapshotted; selected worktree differences
    fail. The required actual runtime Manifest.toml is separately preserved under
    environment/, including when intentionally ignored/untracked; it is never
    force-added to Git or represented as committed source. Frozen model and media-exporter sources are recovered by exact SHA-256
    from bounded Git history. Missing source versions fail, never substitute.
    Repeat inputs reuse the verified generation; foreign output edits fail.
    """
    repository = Path(repository).resolve(strict=True)
    results = Path(results_dir).resolve()
    root = results / "reproducibility"
    if not isinstance(devices, dict) or not devices or not set(devices) <= {"A", "B"}:
        raise ValueError("devices must explicitly map A/B original and model directories")
    if min(max_git_commits, max_files, max_bytes) <= 0:
        raise ValueError("Inspection/copy bounds must be positive")
    run_roots, history_dirs = run_roots or {}, history_dirs or {}
    if not isinstance(run_roots, dict) or not isinstance(history_dirs, dict):
        raise ValueError("run_roots/history_dirs require labeled explicit directory mappings")
    previous, pointer_sha = verify_previous(root)
    head = git(repository, "rev-parse", "HEAD").stdout.decode().strip()
    entries, snapshots, excluded, versions, requirements = {}, {}, [], {}, []
    total = 0

    def add(path, target, category, *, expected=None):
        nonlocal total
        target = relative(target).as_posix()
        path = Path(path)
        if path.is_symlink():
            raise ValueError("Source symlinks are not copied")
        before = fingerprint(path)
        if before[2] > max_bytes - total or len(entries) >= max_files:
            raise ValueError("Explicit reproducibility copy bound exceeded; narrow inputs or increase the declared bound")
        data = path.read_bytes()
        if fingerprint(path) != before:
            raise ValueError("Input changed while reading")
        actual = digest(data)
        if expected is not None and actual != expected:
            raise ValueError(f"Worktree/source bytes differ from committed or bound content: {path}")
        if SECRET.search(data):
            if category in {"committed_code", "runtime_environment"}:
                raise ValueError(f"Credential-like content in required code/environment; refusing publication: {path}")
            excluded.append(dict(path=str(path), reason="credential-like content; bytes not copied"))
            return
        if target in entries:
            raise ValueError("Duplicate output artifact path")
        entries[target] = dict(source=str(path), sha256=actual, bytes=len(data), category=category)
        snapshots[str(path)] = (before, actual)
        total += len(data)

    committed_paths = code_paths(repository, head)
    for name in committed_paths:
        blob = git(repository, "show", f"{head}:{name}").stdout
        add(repository / name, "code/" + name, "committed_code", expected=digest(blob))
    tool = repository / "tools/package_clicker_reproducibility.py"
    if sha(tool) != sha(Path(__file__)):
        raise ValueError("Executing packager differs from the committed repository copy")
    runtime_manifest = repository / "Manifest.toml"
    if not runtime_manifest.is_file():
        raise FileNotFoundError(f"Actual runtime Manifest.toml is required: {runtime_manifest}")
    add(runtime_manifest, "environment/Manifest.toml", "runtime_environment")
    tracked_manifest = "Manifest.toml" in committed_paths
    environment = dict(
        schema_version="clicker-runtime-environment-v1", source_path=str(runtime_manifest),
        snapshot="environment/Manifest.toml", manifest_sha256=entries["environment/Manifest.toml"]["sha256"],
        project_sha256=entries["code/Project.toml"]["sha256"], git_head=head,
        git_tracked_at_head=tracked_manifest,
        git_ignored=git(repository, "check-ignore", "-q", "--", "Manifest.toml", check=False).returncode == 0,
        provenance="tracked_current_runtime_environment" if tracked_manifest else "untracked_runtime_environment",
        committed_manifest_sha256=entries["code/Manifest.toml"]["sha256"] if tracked_manifest else None,
        historical_execution_environment_asserted=False,
        interpretation="Actual current runtime environment captured at packaging, not proof of historical fit environments. Compare saved calibration software versions and historical execution receipts.")
    latest = {}
    for name, locations in sorted(devices.items()):
        original = Path(locations["original_dir"]).resolve(strict=True)
        model = Path(locations["model_dir"]).resolve(strict=True)
        if any(root == p or p in root.parents or root in p.parents for p in (original, model)):
            raise ValueError("Reproducibility output must be separate from current scientific inputs")
        info, calibration, report, first_report, fit_root = device_summary(name, original, model)
        latest[name] = info
        roots = [("original", original), ("model", model)]
        if fit_root not in (original, model):
            roots.append(("fit_inputs", fit_root))
        for kind, directory in roots:
            skip = (root, model) if directory != model else (root,)
            for path, rel in walk_artifacts(directory, excluded=skip):
                add(path, f"latest/{name}/{kind}/{rel}", "latest_numerical")
        for document in (calibration, report, first_report):
            for path, expected in source_requirements(document).items():
                requirements.append((f"latest/{name}", path, expected))
        for phase in ("review", "comparison"):
            path = results / name / "videos" / f"{phase}_export_manifest.json"
            if path.is_file():
                media = read_json(path)
                if media.get("status") != "complete":
                    raise ValueError("Incomplete historical media manifest")
                expected = media["settings"]["exporter_sha256"]
                if not re.fullmatch(r"[0-9a-f]{64}", expected):
                    raise ValueError("Malformed media exporter source hash")
                add(path, f"latest/{name}/media_manifests/{path.name}", "media_metadata_only")
                requirements.append((f"latest/{name}/media/{phase}", MODULE_ROOT + "video_export.py", expected))
    for kind, mapping in (("runs", run_roots), ("superseded_history", history_dirs)):
        for name, directory in sorted(mapping.items()):
            label(name)
            active = active_job_roots(directory) if kind == "runs" else []
            excluded.extend(dict(path=str(path), reason="nonterminal job evidence is still mutable", status=status) for path, status in active)
            for path, rel in walk_artifacts(Path(directory), evidence=kind == "runs", excluded=(root, *(p for p, _ in active))):
                add(path, f"{kind}/{name}/{rel}", "execution_evidence" if kind == "runs" else "superseded_numerical_history")
    recovered = {}
    for owner, path, expected in sorted(set(requirements)):
        key = (path, expected)
        if key not in recovered:
            data, commit = recover_source(repository, head, path, expected, max_git_commits)
            if SECRET.search(data):
                raise ValueError("Credential-like content in frozen source; refusing publication")
            target = f"frozen_sources/{expected}/{path}"
            if len(entries) >= max_files or len(data) > max_bytes - total:
                raise ValueError("Frozen source copy bound exceeded")
            entries[target] = dict(source=f"git:{commit}:{path}", sha256=expected, bytes=len(data), category="frozen_source", commit=commit)
            recovered[key] = (data, commit, target)
            total += len(data)
        _, commit, target = recovered[key]
        versions.setdefault(owner, []).append(dict(path=path, sha256=expected, commit=commit, snapshot=target))
    configuration = dict(repository=str(repository), git_head=head, devices={k: {n: str(Path(v).resolve()) for n, v in d.items()} for k, d in sorted(devices.items())},
                         run_roots={k: str(Path(v).resolve()) for k, v in sorted(run_roots.items())},
                         superseded_history={k: str(Path(v).resolve()) for k, v in sorted(history_dirs.items())},
                         runtime_environment=environment)
    identity = digest(json_bytes(dict(configuration=configuration, files=entries, excluded=excluded, source_versions=versions)))
    generation = "generation_" + identity[:20]
    if previous and previous["generation"] == generation:
        return {**previous, "root": str(root), "reused": True}
    summary = dict(schema_version=SCHEMA, git_head=head, runtime_environment=environment, latest=latest,
                   history={name: {"source_directory": str(Path(path).resolve()), "status": "superseded_not_current_accuracy"} for name, path in sorted(history_dirs.items())},
                   source_versions=versions, quality="Provisional physical fits to visual estimates; weak identifiability, bounds, coverage and original uncertainty remain material.",
                   exclusions="No input recordings, full videos, private depots, scratch, credential files, or detected credential-bearing logs. Copied canonical bytes are never redacted or transformed.")
    instructions = f"""# Clicker calibration reproducibility

Current `code/` bytes match committed HEAD `{head}`. `summary.json` is extracted
from the supplied frozen reports/calibrations; it is not a new evaluation.
`source_versions.json` identifies exact SHA-256 recovered source bytes and the
Git commit containing each version, including historical video exporters.

Original Terra annotations, uncertainty radii, user scope and explicit fit-input
adapters remain distinct under `latest/`. First-fit outer-block results retain
their original policy. Refinement evaluation is reused and descriptive; final
deployment errors are in-sample. `superseded_history/` is preserved evidence,
not the current calibration or independent accuracy evidence.

To reproduce in a NEW writable workspace:

1. Copy `code/` to that workspace, then place `environment/Manifest.toml`
   alongside `code/`'s Project.toml in the new workspace. The environment copy
   is the actual current runtime Manifest, even when ignored/untracked in Git;
   `environment/provenance.json` records that distinction and its exact SHA.
   `code/` remains strictly HEAD, including its Manifest only if committed.
   This current environment is not asserted to reproduce every historical fit:
   compare the recorded Julia/package versions and execution receipts for each
   saved calibration before reproducing older numerical results.
2. For an existing scientific model, overlay the exact files listed for its
   device in `source_versions.json` from `frozen_sources/`. Do not silently use
   current fitter bytes for older A results. Historical exporter versions are
   identified separately by device and review/comparison phase.
3. Use the archived canonical events, exact filtered snippet cache, original
   annotations/seal and explicit observation policy. A's fit adapter is not
   replacement Terra truth. Keep evaluation splits and plans unchanged.
4. The full fresh-workflow entry point is `run_full_clicker_calibration`; see
   `code/docs/clicker_calibration.md`. Original recordings are NOT bundled:
   obtain their exact bytes using the retained preparation input manifests.
   Choose new output/results directories; do not overwrite preserved runs.
5. Job receipts and process logs under `runs/` record actual arguments and
   failures. This packager never executes those recorded commands. Review a
   command and its paths before explicitly rerunning it. Historical source or
   artifact mismatches must fail rather than rebase a reference.

The current pointer is `../manifest.json`; each generation contains its own
`generation_manifest.json` mapping every payload file to SHA-256. Original
device delivery bundles and videos remain outside this reproducibility folder.
"""
    generated = {"summary.json": json_bytes(summary), "source_versions.json": json_bytes(versions),
                 "environment/provenance.json": json_bytes(environment),
                 "REPRODUCE.md": instructions.encode()}
    root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".staging-", dir=root))
    try:
        frozen_data = {target: data for data, _, target in recovered.values()}
        for target, info in entries.items():
            data = frozen_data[target] if target in frozen_data else Path(info["source"]).read_bytes()
            if digest(data) != info["sha256"]:
                raise ValueError("Canonical input changed before copying")
            write_file(staging / target, data)
        for target, data in generated.items():
            write_file(staging / target, data)
        for path, (before, expected) in snapshots.items():
            if fingerprint(path) != before or sha(path) != expected:
                raise ValueError("Canonical input changed during packaging")
        if git(repository, "rev-parse", "HEAD").stdout.decode().strip() != head:
            raise ValueError("Committed HEAD changed during packaging")
        pointer = root / "manifest.json"
        if (pointer_sha is None and pointer.exists()) or (pointer_sha is not None and sha(pointer) != pointer_sha):
            raise ValueError("Reproducibility ownership pointer changed during packaging")
        files = {target: info["sha256"] for target, info in entries.items()}
        files.update({target: digest(data) for target, data in generated.items()})
        manifest = dict(schema_version=SCHEMA, generation=generation, inventory_sha256=identity,
                        created_utc=datetime.now(timezone.utc).isoformat(), configuration=configuration,
                        files_sha256=files, canonical_artifacts=entries, excluded_files=excluded,
                        payload_bytes=sum((staging / p).stat().st_size for p in files))
        write_file(staging / "generation_manifest.json", json_bytes(manifest))
        destination = root / generation
        if destination.exists():
            # Returning to an earlier already-owned identical generation is safe.
            if not previous or generation not in previous["generations"]:
                raise ValueError("Foreign generation destination exists")
            manifest_sha = previous["generations"][generation]
        else:
            _publish_absent(staging, destination)
            manifest_sha = sha(destination / "generation_manifest.json")
        generations = {} if previous is None else dict(previous["generations"])
        generations[generation] = manifest_sha
        active = dict(schema_version=POINTER_SCHEMA, generation=generation, generation_manifest_sha256=manifest_sha,
                      generations=generations, git_head=head, inventory_sha256=identity)
        temp_pointer = root / (".manifest-" + str(os.getpid()))
        write_file(temp_pointer, json_bytes(active))
        os.replace(temp_pointer, pointer)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return {**active, "root": str(root), "reused": False}


def mappings(values):
    result = {}
    for value in values:
        name, separator, path = value.partition("=")
        if not separator or not path or label(name) in result:
            raise ValueError("Use unique LABEL=PATH mappings")
        result[name] = path
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--devices-json", required=True, help="JSON mapping A/B to original_dir/model_dir")
    parser.add_argument("--run-root", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--history", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--max-git-commits", type=int, default=256)
    args = parser.parse_args()
    print(json.dumps(package_reproducibility(args.repository, args.results_dir,
        devices=read_json(args.devices_json), run_roots=mappings(args.run_root),
        history_dirs=mappings(args.history), max_git_commits=args.max_git_commits)))


if __name__ == "__main__":
    main()
