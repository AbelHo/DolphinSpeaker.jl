#!/usr/bin/env python3
"""Run an explicit calibration command with persistent logs and optional wake-up.

Linux execution helper; scientific work stays in the native Julia entry point.
No command is started on import. Hard resets leave an explicit running receipt.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import re
import signal
import subprocess
import sys
import time
import traceback


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: dict) -> None:
    temp = path.with_suffix(".tmp")
    with temp.open("w") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    temp.replace(path)


def log_error(path: Path, operation: str, error: BaseException, **context) -> None:
    record = dict(timestamp=now(), worker="clicker-job-launcher", pid=os.getpid(),
                  operation=operation, exception_type=type(error).__name__,
                  exception=str(error), traceback=traceback.format_exc(), context=context)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a") as f:
            f.write(json.dumps(record) + "\n")
            f.flush()
            os.fsync(f.fileno())
    except OSError as log_failure:
        print(f"Could not persist failure: {log_failure}; original: {record}", file=sys.stderr)


def child_limits() -> None:
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    requested = 64 * 1024**2
    resource.setrlimit(resource.RLIMIT_STACK, (min(requested, hard) if hard >= 0 else requested, hard))


def sampled_tree_rss(pid: int) -> int:
    """Sample the known process and its child tree; exited processes disappear."""
    queue, seen, total = [pid], set(), 0
    while queue:
        current = queue.pop()
        if current in seen:
            continue
        seen.add(current)
        base = Path("/proc") / str(current)
        try:
            for line in (base / "status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    total += int(line.split()[1]) * 1024
            for children in (base / "task").glob("*/children"):
                queue.extend(map(int, children.read_text().split()))
        except (FileNotFoundError, ProcessLookupError):
            continue
    return total


def group_live(pgid: int) -> bool:
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            fields = (proc / "stat").read_text().rsplit(") ", 1)[1].split()
            if int(fields[2]) == pgid and fields[0] not in {"Z", "X"}:
                return True
        except (FileNotFoundError, ProcessLookupError):
            continue
    return False


def stop_owned(child: subprocess.Popen, grace_s: float = 5) -> None:
    # The session/group is ours even after its leader has exited. Descendants
    # that ignore TERM must not survive merely because wait() reaped the leader.
    if group_live(child.pid):
        try:
            os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + grace_s
    while group_live(child.pid) and time.monotonic() < deadline:
        time.sleep(0.05)
    if group_live(child.pid):
        try:
            os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    child.wait(timeout=10)


class JobInterrupted(KeyboardInterrupt):
    pass


def run_job(command: list[str], job_dir: Path, *, timeout_s: float,
            rss_limit_gib: float, wake: bool, error_log: Path) -> dict:
    job_dir = job_dir.resolve()
    job_dir.mkdir(parents=True, exist_ok=False)
    scratch = job_dir / "scratch"
    scratch.mkdir()
    repository = Path(__file__).resolve().parents[1]
    depot = repository / "calibration_runs/runtime/julia_depot"
    depot.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(TMPDIR=str(scratch), JULIA_DEPOT_PATH=str(depot) + ":" + str(Path.home()/".julia"),
               JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1", CLICKER_ERROR_LOG=str(error_log.resolve()))
    python = repository / "dolphin-whisper/.venv/bin/python"
    if python.is_file():
        env.setdefault("CLICKER_PYTHON", str(python))
    receipt = dict(status="starting", command=command, started_utc=now(),
                   job_dir=str(job_dir), timeout_s=timeout_s, rss_limit_gib=rss_limit_gib,
                   peak_sampled_tree_rss_bytes=0, wake_requested=wake)
    receipt_path = job_dir / "result.json"
    write_json(receipt_path, receipt)
    child = None
    wake_status = None
    started = time.monotonic()
    prior_handlers = {}

    def interrupted(signum, _frame):
        raise JobInterrupted(f"Launcher received signal {signum}")

    for signum in (signal.SIGTERM, signal.SIGINT):
        prior_handlers[signum] = signal.signal(signum, interrupted)
    try:
        with (job_dir / "process.log").open("wb", buffering=0) as log:
            child = subprocess.Popen(command, cwd=repository, env=env, stdout=log,
                                     stderr=subprocess.STDOUT, start_new_session=True,
                                     preexec_fn=child_limits)
            receipt.update(status="running", child_pid=child.pid)
            if wake:
                state_root = job_dir / "watcher_state"
                state_root.mkdir()
                watch_env = {**env, "XDG_STATE_HOME": str(state_root)}
                stat = Path(f"/proc/{child.pid}/stat").read_text()
                ticks = stat.rsplit(") ", 1)[1].split()[19]
                wake_status = state_root / ".hagent/wake-on-exit" / f"{child.pid}-{ticks}.status"
                message = f"Clicker calibration process ended. Inspect {receipt_path} and {job_dir/'process.log'}, then continue the authorized A/B task."
                registrar = subprocess.run(["wake-on-exit", str(child.pid), message], env=watch_env,
                                           text=True, capture_output=True, timeout=15)
                (job_dir / "watcher_registration.log").write_text(registrar.stdout + registrar.stderr)
                if registrar.returncode:
                    raise RuntimeError(f"Completion watcher registration failed ({registrar.returncode})")
                receipt.update(wake_registered=True, wake_status_path=str(wake_status))
            write_json(receipt_path, receipt)
            print(json.dumps(dict(status="running", child_pid=child.pid, receipt=str(receipt_path))), flush=True)
            next_checkpoint = 0.0
            while child.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed > timeout_s:
                    raise TimeoutError(f"Calibration command exceeded {timeout_s}s")
                rss = sampled_tree_rss(child.pid)
                receipt["peak_sampled_tree_rss_bytes"] = max(rss, receipt["peak_sampled_tree_rss_bytes"])
                if rss > rss_limit_gib * 1024**3:
                    raise MemoryError(f"Owned sampled RSS {rss} exceeds {rss_limit_gib}GiB")
                if elapsed >= next_checkpoint:
                    receipt["elapsed_s"] = elapsed
                    write_json(receipt_path, receipt)
                    next_checkpoint = elapsed + 15
                time.sleep(0.5)
            receipt["exit_code"] = child.returncode
            if group_live(child.pid):
                raise RuntimeError("Command exited while owned process-group descendants were still running")
            if child.returncode:
                raise subprocess.CalledProcessError(child.returncode, command)
        receipt["status"] = "complete"
    except BaseException as error:
        for signum in prior_handlers:
            signal.signal(signum, signal.SIG_IGN)
        if child is not None:
            try:
                stop_owned(child)
            except BaseException as cleanup_error:
                receipt["cleanup_error"] = str(cleanup_error)
                log_error(error_log, "calibration subprocess cleanup", cleanup_error, child_pid=child.pid)
            receipt["exit_code"] = child.returncode
        receipt.update(status="failed", exception_type=type(error).__name__, exception=str(error))
        process_log=job_dir/"process.log"
        child_trace=""
        if process_log.exists():
            with process_log.open("rb") as f:
                f.seek(max(0,process_log.stat().st_size-131072))
                child_trace=f.read().decode(errors="replace")
            child_trace=re.sub(r"bot[0-9]+:[A-Za-z0-9_-]+","bot[REDACTED]",child_trace)
        log_error(error_log, "calibration subprocess", error, command=command,
                  process_log=str(process_log), child_output_tail=child_trace,
                  exit_code=receipt.get("exit_code"))
    finally:
        receipt.update(completed_utc=now(), elapsed_s=time.monotonic()-started)
        write_json(receipt_path, receipt)
        # Keep the launching process alive through watcher delivery. This avoids
        # losing a fallback watcher when an enclosing execution session closes.
        if wake_status is not None:
            deadline = time.monotonic() + 60
            while time.monotonic() < deadline:
                value = wake_status.read_text().strip() if wake_status.exists() else "missing"
                if value in {"delivered", "delivery-failed", "monitor-killed"}:
                    break
                time.sleep(0.5)
            receipt["wake_final_status"] = value
            if value != "delivered":
                error = RuntimeError(f"Completion notification ended in {value}")
                log_error(error_log, "calibration completion notification", error, status_path=str(wake_status))
            write_json(receipt_path, receipt)
        for signum, handler in prior_handlers.items():
            signal.signal(signum, handler)
    print(json.dumps(receipt), flush=True)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-dir", type=Path, required=True)
    parser.add_argument("--timeout-s", type=float, default=7200)
    parser.add_argument("--rss-limit-gib", type=float, default=28)
    parser.add_argument("--wake", action="store_true")
    parser.add_argument("--error-log", type=Path, default=Path("docs/error.log"))
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.timeout_s <= 0 or args.rss_limit_gib <= 0:
        parser.error("An explicit command and positive limits are required")
    result = run_job(command, args.job_dir, timeout_s=args.timeout_s,
                     rss_limit_gib=args.rss_limit_gib, wake=args.wake, error_log=args.error_log)
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
