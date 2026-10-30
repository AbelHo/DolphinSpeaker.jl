import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "tools/run_clicker_job.py"
spec = importlib.util.spec_from_file_location("clicker_job", SCRIPT)
job = importlib.util.module_from_spec(spec)
spec.loader.exec_module(job)


class JobLifecycle(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_nonzero_exit_is_durable(self):
        result = job.run_job([sys.executable, "-c", "raise RuntimeError('controlled failure')"],
                             self.root / "failed", timeout_s=5, rss_limit_gib=1,
                             wake=False, error_log=self.root / "errors.log", serialize=False)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["exit_code"], 1)
        record = json.loads((self.root / "errors.log").read_text().splitlines()[0])
        self.assertEqual(record["operation"], "calibration subprocess")
        self.assertIn("controlled failure", (self.root / "failed/process.log").read_text())

    def test_descendant_cleanup_after_leader_exit(self):
        pidfile = self.root / "grandchild.pid"
        grandchild = "import signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);time.sleep(60)"
        parent = f"import subprocess,sys,pathlib,time;p=subprocess.Popen([sys.executable,'-c',{grandchild!r}]);pathlib.Path({str(pidfile)!r}).write_text(str(p.pid));time.sleep(.2)"
        leader = subprocess.Popen([sys.executable, "-c", parent], start_new_session=True)
        leader.wait(timeout=5)
        self.assertTrue(job.group_live(leader.pid))
        job.stop_owned(leader, grace_s=0.2)
        deadline = time.monotonic() + 3
        while job.group_live(leader.pid) and time.monotonic() < deadline:
            time.sleep(0.05)
        self.assertFalse(job.group_live(leader.pid))

    def test_wrapper_sigterm_stops_owned_session(self):
        target = self.root / "signal_job"
        wrapper = subprocess.Popen([sys.executable, str(SCRIPT), "--job-dir", str(target),
                                    "--lightweight",
                                    "--timeout-s", "60", "--error-log", str(self.root/"signal_errors.log"),
                                    "--", sys.executable, "-c", "import time;time.sleep(60)"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        receipt = target / "result.json"
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            if receipt.exists():
                data = json.loads(receipt.read_text())
                if data.get("status") == "running":
                    break
            time.sleep(0.05)
        else:
            wrapper.kill()
            self.fail("wrapper did not start")
        child_pid = data["child_pid"]
        wrapper.send_signal(signal.SIGTERM)
        _, stderr = wrapper.communicate(timeout=15)
        self.assertEqual(wrapper.returncode, 1, stderr.decode())
        result = json.loads(receipt.read_text())
        self.assertEqual(result["exception_type"], "JobInterrupted")
        self.assertFalse(job.group_live(child_pid))

    def test_full_cpu_environment_and_memory_default(self):
        command = "import os,json;print(json.dumps({k:os.environ[k] for k in ['JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}))"
        result = job.run_job([sys.executable, "-c", command], self.root / "resources",
                             timeout_s=5, rss_limit_gib=None, wake=False,
                             error_log=self.root / "errors.log", serialize=False)
        self.assertEqual(result["status"], "complete")
        actual = json.loads((self.root / "resources/process.log").read_text())
        self.assertEqual(set(actual.values()), {str(job.available_resources()[0])})
        self.assertGreater(result["rss_limit_gib"], 0)

    def test_heavy_commands_are_staggered(self):
        # Two concurrently submitted commands must run one at a time, not each
        # consume a fraction of the available machine under artificial caps.
        workers = []
        for i in range(2):
            directory = self.root / f"queued{i}"
            command = "import time,json;start=time.time();time.sleep(.6);print(json.dumps([start,time.time()]))"
            workers.append((directory, subprocess.Popen([sys.executable, str(SCRIPT),
                "--job-dir", str(directory), "--timeout-s", "5", "--error-log", str(self.root / "queue_errors.log"),
                "--", sys.executable, "-c", command], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)))
        intervals = []
        for directory, worker in workers:
            _, stderr = worker.communicate(timeout=15)
            self.assertEqual(worker.returncode, 0, stderr.decode())
            intervals.append(json.loads((directory / "process.log").read_text()))
            self.assertTrue(json.loads((directory / "result.json").read_text())["serialized_heavy_job"])
        intervals.sort()
        self.assertLessEqual(intervals[0][1], intervals[1][0])


if __name__ == "__main__":
    unittest.main()
