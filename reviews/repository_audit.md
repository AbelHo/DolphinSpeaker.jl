# Repository audit

## Scope and method

This review combined bounded static inspection, Git history, file identity and
reference checks, directory sizing, and one isolated package-load reproduction.
It did not run the existing test suite, because `test/runtests.jl` is still the
PkgTemplates test suite rather than a DolphinSpeaker suite. It did not execute
the untracked PINN or cuML experiments.

The conclusions below use three confidence labels:

- **Confirmed**: identity, history, or direct runtime evidence establishes the
  conclusion.
- **High confidence**: replacement and usage evidence is strong, but a focused
  regression test should precede deletion.
- **Needs characterization**: old or misplaced, but not safe to delete on static
  evidence alone.

## What is good and worth preserving

The repository contains a substantial working domain pipeline, not just
prototypes:

- `src/audio.jl` covers audio ingest/conversion.
- `src/synchronization.jl` and `src/media_info.jl` coordinate audio/video timing
  and metadata.
- `src/detector_impulsive.jl`, `src/detector_tonal.jl`, and
  `src/detector_band.jl` implement complementary detectors.
- `src/dsp.jl` and `src/localization.jl` contain reusable signal-processing,
  TDOA, geometry, and localization work.
- `src/audacity.jl` and `src/raven.jl` preserve useful label interchange.
- `src/video.jl`, `src/pic2vid.jl`, and `overlay_ffmpeg.py` provide an annotation
  and visualization path.
- `src/cam_calib.jl` and `src/cam_projection.jl` preserve camera calibration
  knowledge.
- `run_contiguous_folders` in `src/run_example.jl` is a valuable end-to-end
  workflow: synchronization/concatenation, threaded per-file detection,
  combined labels, video overlay, and optional spectrogram creation.
- The device presets in `src/config.jl` capture hard-won deployment knowledge
  for several hydrophone/camera arrangements. Their contents should survive the
  migration even though their current global-mutation API should not.
- The newer overlay documentation and the synthetic integration draft in
  `test/test_overlay_annotations.jl` are the beginning of a real project test
  corpus.
- The large local data/output trees are not committed as Git blobs. The tracked
  tree is small (about 1.55 MB), so the problem is organization and local
  environment weight rather than repository-object bloat.

Recent history agrees with this assessment: the most actively developed files
include `run_example.jl`, `video.jl`, `media_info.jl`, `config.jl`, `audio.jl`,
`plotting.jl`, `dsp.jl`, and `synchronization.jl`.

## Principal findings

### 1. The majority of tracked file count is foreign PkgTemplates residue

**Confirmed.** The root commit, `31bd58f` (2023-12-01), imported the complete
PkgTemplates repository: 182 files and 10,807 lines. Commit `a382130`, titled
`chore: remove template files`, removed only four files:

- `src/PkgTemplates.jl`
- `src/deprecated.jl`
- `src/interactive.jl`
- `src/template.jl`

At least 173 unambiguous PkgTemplates files, about 381 KB, remain:

| Area | Count | Direct evidence |
|---|---:|---|
| `.github/**` | 5 | Original automation; CI still doctests `PkgTemplates` |
| PkgTemplates docs setup/guides | 6 | `docs/make.jl` builds PkgTemplates; index documents it |
| `templates/**` | 29 | PkgTemplates generator templates |
| `src/plugin.jl`, `src/show.jl`, `src/plugins/**` | 20 | Foreign implementation files |
| Original tests and fixtures | 113 | `test/runtests.jl` imports and tests PkgTemplates |

These files have no repository-specific evolution after the import. They distort
search results, make the test and documentation entry points false, and hide the
actual package surface.

Recommended disposition: remove the foreign implementation, templates, tests,
fixtures, and docs in a provenance-labelled change. Replace—not merely delete—
the test runner, docs build, and CI workflow.

The current `LICENSE` is also the imported PkgTemplates/Invenia MIT copyright.
That is a provenance issue, not just clutter. Preserve it until authorship and
licensing of copied and new DolphinSpeaker code have been reviewed; do not
silently substitute a new license.

### 2. The default branch does not represent the maintained project

**Confirmed.** `origin/HEAD` points to `master` at `9306007` (2024-04-10), while
the working branch `temp_commit` is at `621fe73` (2026-06-29). `temp_commit` is a
linear descendant, 195 commits ahead and zero behind. The README explicitly
installs `#temp_commit`.

This is confusing for users, CI, releases, and code browsers. It is not a hard
merge problem: after preserving and verifying the current dirty working tree,
the maintained lineage can be fast-forwarded into a properly named default
branch. The current CI only reacts to `master`, further reducing the chance that
active code has meaningful automated coverage.

### 3. The include graph repeatedly evaluates the same source files

**Confirmed by static and runtime evidence.** `src/DolphinSpeaker.jl` includes a
small set of files, but those files recursively include siblings, which include
more siblings and sometimes the original files again. Common files have many
include sites: `audio.jl` has 13, `config.jl` and `plotting.jl` have 10 each, and
`utils.jl`, `dsp.jl`, and `run_example.jl` have 7 each.

An isolated `using DolphinSpeaker` reproduction produced a precompile failure:

```text
WARNING: Method definition showall(Any) in src/utils.jl overwritten ...
ERROR: Method overwriting is not permitted during Module precompilation.
```

Julia then fell back and eventually loaded the package, but printed extensive
repeated documentation-replacement warnings. This means “it loads in a running
session” masks a broken package-precompilation contract.

Required rule: `src/DolphinSpeaker.jl` must be the only production include
coordinator. Leaf files should import names through the module and never include
siblings.

### 4. Core and experimental dependency scopes are conflated

**Confirmed.** `Project.toml` contains 73 direct dependencies, only 25 compat
entries, and no Julia compatibility declaration. Heavy differentiable-programming
and accelerator packages—such as Lux, LuxCUDA, Reactant, Zygote,
ComponentArrays, ForwardDiff, and Optimisers—are direct dependencies even though
the apparent PINN implementation is currently an untracked, non-included
experiment.

This makes ordinary package loading resolve and precompile an experimental
machine-learning stack. During reproduction, the process tree briefly reached
about 15.4 GiB resident memory while Reactant extensions were compiling. Reactant
also warned that Julia 1.12+ is unsupported and recommended Julia 1.11/LTS.

Python/Conda and GPU clustering material has the same boundary problem. The
untracked `CondaPkg.toml` declares pandas, matplotlib, and hdbscan—but not cuML—
while `CUML_INSTALLATION_COMPLETE.md` claims cuML is configured and installed.
That claim is not supported by the configuration present in the tree.

Recommended disposition: make the core package CPU-capable and relatively
small; place PINN, clustering, and Python/GPU work in explicit experiment
environments or Julia package extensions.

### 5. Production, workflow, experiment, and scratch code share `src/`

**Confirmed.** Files named `test_run.jl` and `test_make_video.jl` contain
production functions and are included by the live pipeline. Conversely,
`init.jl`, `dev2_sync.jl`, `run_analysisDatai.jl`, scratch/temp files, and other
research scripts have hard-coded local paths or top-level execution.

This makes the `src/` directory an unreliable description of what constitutes
the library. Code should be classified by execution contract:

- reusable, side-effect-free package code in `src/`;
- user-facing recipes in `examples/`;
- parameterized batch entry points in `scripts/`;
- exploratory work with independent environments in `experiments/`;
- test-only material in `test/`.

### 6. Configuration is valuable but hidden in mutable global state

**Confirmed structurally.** Device functions in `src/config.jl` mutate global
values, and the file is repeatedly included. The file also selects `calf_hk`
and invokes its setter during inclusion. A run's behavior therefore depends on
include order and prior calls in the Julia process. For example,
`run_example.jl` assigns a tonal threshold before later nested includes lead
back to `config.jl` and reset configuration.

Preserve the presets as constructors for explicit immutable `DeviceConfig` and
`ProcessingConfig` values. Pass those values to pipeline functions. Temporary
compatibility wrappers can preserve old entry points while emitting a
deprecation warning.

### 7. Several workflows lack bounded failure behavior

Examples found during inspection:

- `run_contiguous_folders` starts an ffmpeg concat process asynchronously and
  polls `ffprobe` without a timeout. A failed producer can leave a consumer
  waiting forever.
- Linux CSV append logic in `src/detector.jl` reads and rewrites the entire
  existing file, increasing latency and data-loss exposure.
- Many analysis files contain machine-specific `/media/spin/...` paths.
- `src/init.jl` executes a `process_folder(...)` call against an absolute path at
  file load.
- The active include chain reaches `src/pinger_calibration.jl`, which creates
  `~/Desktop/results` during inclusion. Importing a library must not create a
  user directory.

These are reliability issues to address after a characterization harness exists.

## Confirmed redundant or obsolete material

| Material | Classification | Evidence | Action |
|---|---|---|---|
| PkgTemplates source/plugins/templates | Confirmed foreign residue | Byte-identical to root import; no relevant history | Delete in one provenance-labelled change |
| PkgTemplates tests/fixtures | Confirmed foreign residue | `test/runtests.jl` imports PkgTemplates | Replace with DolphinSpeaker harness, then delete |
| PkgTemplates docs/build | Confirmed foreign residue | Docs build and index name PkgTemplates | Replace with DolphinSpeaker docs |
| Existing CI | Confirmed obsolete/misaligned | Doctests PkgTemplates; old Julia matrix; watches stale branch | Replace, do not incrementally patch |
| `src/tonal_detector.jl` | High-confidence superseded | Old and new files define `process_audioVideo_tonal1`; `detector_tonal.jl` is larger and actively maintained | Remove old include after regression coverage, then delete |
| `src/init.jl` | High-confidence obsolete/dangerous | Older entry-point copy plus top-level machine-specific execution | Move any recipe to `examples/`, then delete |
| `src/temp_rx_tx.jl` | High-confidence superseded | Unreferenced temporary file; newer expanded `analysis_rx_tx.jl` exists | Preserve useful notes/results, then archive/delete |
| `src/bearing_video.jl` | High-confidence duplicate | Its main implementation is effectively duplicated in `test_run.jl` | Select one implementation during pipeline extraction, then retire the duplicate |
| Root environment snapshots | Superseded snapshots | `Project_old.toml`, `Project_20250728.toml`, `Manifest_20250728.toml`; current project is newer | Move to named legacy environment or rely on a tag |
| `temp_test_detections.csv` | Misplaced artifact | Root-level test/generated naming | Move to asserted fixture or delete |
| Root `test_overlay_annotations.jl` | Misplaced demo | Writes CSVs in current directory; newer test exists under `test/` | Move useful recipe to `examples/` or retire |

## Likely legacy, but characterize before deletion

| Material | Why it looks old/duplicated | Why deletion must wait |
|---|---|---|
| `src/calibration_aspod4_2.jl` | Old monolith duplicates localization, plotting, and media helpers | Still enters the runtime graph through `test_run.jl` |
| `src/test_run.jl` | Duplicates newer pipeline/plot functions and is badly named | Still provides production `process_audioVideo` behavior |
| `src/ambient_noise.jl`, `multimedia.jl`, `map.jl`, `plot_map.jl` | Old ARL import, little/no later history, largely unreferenced | Could preserve niche analysis APIs or field recipes |
| `src/bandlimited_impulse.jl` | Mixed library/experiment structure and duplicate includes | Domain behavior may not exist elsewhere |
| `src/dev2_sync.jl`, `run_analysisDatai.jl`, scratch/temp files | Hard-coded paths and exploratory top-level code | Research provenance may matter; move rather than delete |
| `src/pinn_ai.jl` and clustering/cuML material | Experimental and currently untracked | User work; preserve in dedicated environments |

`src/pinn_ai.jl` deserves special care: if directly included without defining
`Main.DISABLE_PINN_EXAMPLES`, it launches two large training examples (roughly
5,000 points × 3,000 epochs and 8,000 points × 4,000 epochs). It was not run in
this review. Its examples must be opt-in before the file is used or moved.

Likewise, `test_cuml_direct.py` performs a `pip install --upgrade cuml-cu11` and
then launches GPU cuML algorithms. Its shebang says Julia although its contents
are Python. It is not a safe test and was not executed.

Within otherwise current files, `detect_tonal2` and the empty
`quietest_segment` stub in `detector_tonal.jl` have no discovered callers. Treat
them as cleanup candidates, not confirmed deletions, until detector fixtures
cover their intended cases. The older blip/trigger synchronization API and the
newer segment-correlation synchronization are different algorithms; they must
not be silently collapsed merely because the newer contiguous workflow uses the
latter.

## Documentation and test state

The apparent volume of tests and docs is misleading because most belongs to
PkgTemplates. The useful Dolphin-specific test material is small and not wired
into `test/runtests.jl`.

Priorities for a real baseline suite:

1. package import and precompile;
2. configuration preset construction and isolation;
3. deterministic synthetic audio and synchronization fixtures;
4. impulsive, tonal, and band detector contracts;
5. TDOA/localization geometry;
6. Audacity and Raven round trips;
7. overlay generation using tiny generated fixtures;
8. a bounded, mocked workflow test for external ffmpeg/ffprobe behavior.

`test/test_overlay_annotations.jl` is a promising seed, but it includes
`src/video.jl` directly and is not run by the current runner. Another file,
`test/run_overlay_examples.jl`, uses hard-coded local paths and enables a
non-dry-run path; it belongs under `examples/` after parameterization.

## Other hygiene findings

- The package version in `Project.toml` is `0.0.7`, while
  `src/DolphinSpeaker.jl` exposes a hard-coded timestamp-like version string.
  One canonical version source is needed.
- The ignored current `Manifest.toml` was generated with Julia 1.12.6. Ignoring
  a root manifest is normal for a reusable Julia library, but a research
  application needs named, committed experiment manifests. Choose and document
  the policy instead of keeping dated files at the root.
- `.CondaPkg/` is untracked, about 3.7 GB and 54,681 files, but is not ignored.
  It should be excluded from Git once current work is protected.
- `temp/` is about 7.1 GB and 43,768 files; `results/` is about 759 MB. Both are
  ignored. Their contents may be valuable scientific output and were not
  classified as disposable.
- Some source files import undeclared or apparently optional packages
  (`ARLToolkit`, GeoMakie/CairoMakie, Distances, FilePaths, ColorTypes). Moving
  experiments first will clarify which dependencies truly belong in core.
- There are no conventional release tags; a tag with the unusual name
  `bugfix(flac_write)-channel(5,6,7)` is not a replacement for a versioned
  release history.

## Overall recommendation

Treat this as a recovery and boundary-setting migration, not a formatting
exercise. First remove the objectively foreign shell around the project and
establish trustworthy tests. Then linearize the include graph. Only after those
two gates should duplicate algorithms or dependencies be removed. This order
turns “looks obsolete” judgments into evidence-backed decisions and retains the
domain knowledge that makes the repository valuable.
