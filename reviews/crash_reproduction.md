# Crash reproduction and resource notes

## Question tested

A previous computer/session crash occurred while this repository and unrelated
workloads were active. The most suspicious command actually run by this review
was replayed under a watchdog to determine whether it independently reproduced
the crash:

```bash
JULIA_DEPOT_PATH=/tmp/<isolated-depot>:/home/spin/.julia \
JULIA_PKG_PRECOMPILE_AUTO=0 \
julia --startup-file=no --project=. -e 'using DolphinSpeaker'
```

The temporary depot isolated writable compilation products while allowing the
existing package cache to be read. Automatic package-manager precompilation was
disabled, but `using` still caused Julia to compile/load the dependency graph as
needed.

## Result

The crash did **not** reproduce.

- A 90-second default-parallel run reached 19 processes and about 8.59 GiB
  aggregate resident memory before the watchdog stopped it.
- Continuing with the same temporary depot briefly reached about 15.4 GiB while
  Reactant extensions compiled. Roughly 74 GiB system memory remained available.
- The machine remained responsive and no Julia/Python process survived the
  completed probe.
- Normal loading eventually printed `LOAD_OK` and exited successfully after
  Julia's fallback behavior.

The same load exposed a repository defect:

```text
WARNING: Method definition showall(Any) in src/utils.jl overwritten ...
ERROR: Method overwriting is not permitted during Module precompilation.
```

It also produced many repeated documentation-replacement warnings. Reactant
warned that it does not support Julia 1.12+ and recommended Julia 1.11/LTS.

## Other action replayed

A full ignored-file name enumeration was repeated with output reduced to counts:

```bash
rg --files -uu | wc -l -c
```

It found 99,406 paths and about 7.64 MB of path text in 0.02 seconds at roughly
28 MiB resident memory. This did not crash the computer. Sending the full
unsuppressed path list into a chat/tool payload could still overload that client
or session, so future scans should continue to exclude `.CondaPkg/`, `temp/`, and
`results/` or aggregate their output.

## System evidence

The prior boot ended abruptly around 07:52 and the current boot began at
07:54:56; `last -x` labeled the graphical tty session as a crash. Available logs
did not show an OOM kill, kernel panic, or GPU Xid immediately before the end.
`coredumpctl` was unavailable. Absence of those records is not proof that no
hardware/driver/resource event occurred.

## Conclusion

The actions actually run by this repository review were unable to reproduce the
machine crash. The Julia load is unnecessarily heavy and the package's include
graph is broken for precompilation, but this machine had ample memory during the
replay. The evidence is therefore more consistent with an unrelated concurrent
workload or an unlogged external/system event than with the bounded audit itself.

## Static hazards discovered but not executed

Two files are more dangerous than their names suggest, but they were not run by
this review and therefore are not candidates for “something the review ran”:

- Untracked `src/pinn_ai.jl` automatically starts two very large training
  examples if included without first defining `Main.DISABLE_PINN_EXAMPLES`.
- Untracked `test_cuml_direct.py` runs `pip install --upgrade cuml-cu11`, then
  launches GPU cuML UMAP, KMeans, and PCA; its shebang incorrectly says Julia.

Do not reproduce either hazard unguarded merely to provoke a crash. Convert them
to explicit, bounded entry points first. The two debug PINN scripts already show
the safer pattern by setting the disable flag, and one reduces the workload to
200 points and 10 epochs.

## Operational recommendations

- Default Julia precompile/task concurrency to a bounded value in local run
  instructions until the dependency split is complete.
- Keep large ignored/untracked trees out of recursive tool payloads.
- Split Reactant/PINN/cuML from the core environment.
- Add timeouts and peak-RSS/process-count logging to any future reproduction.
- If another whole-machine crash occurs, immediately preserve the previous boot
  journal, kernel ring buffer if available, GPU driver log, and workload process
  list before reboot history rotates.
