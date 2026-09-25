# Staged reorganization plan

## Safety principles

- Protect the current dirty working tree before branch or mass-move operations.
- Do not clean untracked files. Several appear to be active experiments.
- Use one concern per commit so removal, movement, and behavioral changes do not
  become inseparable.
- Preserve values and outputs before refactoring APIs.
- Never classify ignored `temp/` or `results/` contents as disposable merely
  because Git ignores them.
- Treat licensing/provenance as a gate, not a cosmetic edit.

## Phase 0 — freeze and inventory the recoverable state

Actions:

1. Make a user-approved backup, WIP branch/commit, or patch bundle containing the
   current tracked edits and untracked experiment files.
2. Record the current head (`621fe73` at review time), Julia version, package
   status, relevant external tool versions, and dataset locations.
3. Inventory which files under ignored `results/`, `temp/`, and `tried_old/` are
   scientific evidence versus reproducible cache.
4. Determine code authorship and licenses for the initial ARL/Dolphin import and
   all retained PkgTemplates-derived code.

Gate: every current modification and scientifically relevant output has a known,
recoverable location. Licensing questions are recorded with an owner.

## Phase 1 — replace the false project shell

Actions:

1. Replace `test/runtests.jl` with a minimal DolphinSpeaker test runner.
2. Replace the docs build and index with DolphinSpeaker-specific material.
3. Replace CI with a small supported-Julia matrix on the maintained branch.
4. Remove the 173 confirmed PkgTemplates implementation/template/test/docs files.
5. Keep or revise `LICENSE` only according to the Phase 0 provenance decision.
6. Add ignores for local environment products such as `.CondaPkg/` and clearly
   generated root outputs, after protecting any current untracked work.

Gate:

- `test/runtests.jl`, docs, and CI mention DolphinSpeaker, not PkgTemplates.
- A clean clone does not attempt to test or document a different project.
- No user experiment or scientific output has been deleted.

## Phase 2 — create characterization coverage

Start with small deterministic fixtures and current behavior, including awkward
behavior that may later change:

1. package import/precompile smoke test;
2. device preset values;
3. short generated waveforms for DSP and each detector family;
4. synchronization edge cases and timestamp conversions;
5. TDOA/localization geometry with synthetic positions;
6. Audacity and Raven read/write round trips;
7. overlay generation using tiny generated media;
8. external process failure, timeout, and missing-tool behavior.

Capture representative outputs from the current pipeline before rearranging the
include graph. For floating-point algorithms, document tolerances and compare
structured values rather than byte-identical plots.

Gate: the active, intended pipeline is named and covered well enough to detect a
structural regression. Expected current failures are explicit, not silently
skipped.

## Phase 3 — linearize package loading

Actions:

1. Choose a dependency order and include every production file once from
   `src/DolphinSpeaker.jl`.
2. Remove all sibling `include(...)` statements from internal files.
3. Move production functions out of `test_run.jl` and `test_make_video.jl` into
   properly named pipeline/visualization files.
4. Separate top-level executable statements from definitions.
5. Use one canonical package version from `Project.toml`.

Gate:

- `Pkg.precompile()` succeeds.
- `using DolphinSpeaker` has no method-overwrite or documentation-replacement
  warnings.
- Including the package performs no data processing, training, package install,
  network operation, directory creation, or external media command.
- Characterization tests remain within their specified tolerances.

## Phase 4 — make configuration explicit

Actions:

1. Define immutable `DeviceConfig` and `ProcessingConfig` structures.
2. Turn current device setters into named constructors that return values.
3. Pass configuration through detector, localization, video, and pipeline APIs.
4. Retain thin compatibility wrappers for existing callers during one migration
   cycle.

Gate:

- Two different device configurations can run in one Julia process without
  changing each other's behavior.
- Preset values match the captured baseline.
- No result depends on source include order.

## Phase 5 — separate experiments and dependency environments

Actions:

1. Move hard-coded analyses, scratch files, and development scripts out of
   `src/` into `experiments/`, `scripts/`, or `examples/` according to contract.
2. Give PINN, clustering/cuML, and other heavy workflows their own
   `Project.toml` and, for reproducibility, committed `Manifest.toml`.
3. Make all example/training calls explicit; in particular, remove automatic
   training from `pinn_ai.jl`.
4. Replace `test_cuml_direct.py` with a correctly named, opt-in environment setup
   and a bounded smoke test. Never install/upgrade packages from a test.
5. Keep Python/GPU integration as an optional extension only if the core API
   truly needs it.
6. Prune core dependencies only after static search plus actual load/tests show
   they are not part of the core path.

Gate:

- The core package resolves, precompiles, and runs its CPU test suite without
  PINN, Reactant, cuML, or a local Conda environment.
- Every experiment documents its Julia/Python/CUDA assumptions and has an
  explicit command to start it.
- No experiment runs simply because its definitions are loaded.

## Phase 6 — retire replacements one pair at a time

Suggested order:

1. Remove the `tonal_detector.jl` include and compare its fixtures with
   `detector_tonal.jl`; then delete it.
2. Move any recipe from `init.jl`, verify no API depends on it, then delete it.
3. Compare `temp_rx_tx.jl` with `analysis_rx_tx.jl`; retain unique provenance,
   then retire the temporary file.
4. Decompose `calibration_aspod4_2.jl`; migrate calibration-only functions and
   verify duplicated helpers against localization/media/plotting modules.
5. Characterize `ambient_noise.jl`, `multimedia.jl`, `map.jl`, `plot_map.jl`, and
   `bandlimited_impulse.jl`; either give each an owner and target module or
   archive it.
6. Consolidate root overlay demos/docs/test artifacts.

Gate for every retirement: no remaining call/reference, replacement coverage is
green, and any changed output has an explained review artifact. Never use age or
filename alone as the deletion criterion.

## Phase 7 — branch, release, and operating policy

Actions:

1. After the current state is protected and the new baseline passes, fast-forward
   the maintained lineage into the real default branch (prefer `main` or a
   clearly named supported branch).
2. Update installation instructions to use the default branch or a release tag.
3. Establish one version source, changelog practice, supported Julia range, and
   semantic release tags.
4. Document whether this is a reusable library, research application, or both;
   use named environments for the latter.
5. Add contribution guidance stating where core code, experiments, examples,
   fixtures, and generated outputs belong.

Gate: a new user can clone the default branch, instantiate the core environment,
run the real tests, build the real docs, and execute a small example without
machine-specific paths.

## Disposition matrix

| Category | Files/areas | Proposed disposition | When |
|---|---|---|---|
| Keep and refactor | active audio, sync, DSP, detection, localization, annotation, video, camera, workflow | Move into clear modules without algorithm changes | Phases 3–4 |
| Preserve domain data | device presets and format conventions | Convert to typed explicit configuration | Phase 4 |
| Delete and replace | PkgTemplates source, templates, tests, docs, CI | Remove as confirmed foreign residue; provide real equivalents | Phase 1 |
| Retire after test | `tonal_detector.jl`, `init.jl`, `temp_rx_tx.jl` | Remove once replacement/recipe is verified | Phase 6 |
| Decompose first | `calibration_aspod4_2.jl`, `test_run.jl`, `test_make_video.jl` | Extract live production functions, then retire shells | Phases 3 and 6 |
| Move, do not discard | analysis/dev/scratch/temp files, PINN, clustering, field scripts | Independent experiments/scripts with explicit entry points | Phase 5 |
| Archive or tag | dated Project/Manifest snapshots | Named legacy environment if still executable; otherwise Git tag/history | Phases 5–6 |
| Preserve pending decision | `results/`, `temp/`, `tried_old/` | Inventory provenance before any cleanup | Phase 0 |
| Verify legally | `LICENSE`, copied origins | Authorship/license audit | Phase 0 |

## Suggested first implementation slice

Keep the first change deliberately small:

1. protect the dirty working state;
2. add a real `test/runtests.jl` containing only an import/precompile smoke test
   plus one tiny Dolphin-specific unit test;
3. replace the CI workflow to run that suite on the maintained branch;
4. delete only the confirmed PkgTemplates tests/templates/plugins/docs;
5. demonstrate that all removed files were foreign by attaching the identity
   inventory to the commit or pull request.

Do not combine that slice with algorithm refactors, dependency pruning, branch
renaming, or global-configuration changes. A clean first boundary makes every
later migration step easier to review and reverse.
