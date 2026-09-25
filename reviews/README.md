# DolphinSpeaker.jl repository review

Review snapshot: 2026-09-21, branch `temp_commit` at `621fe73`, with pre-existing
tracked and untracked working-tree changes left untouched.

This directory contains a read-only architectural and repository-hygiene review.
No package source, tests, documentation, generated output, branch, or environment
was changed as part of the review.

## Documents

- [Repository audit](repository_audit.md) — findings, strengths, confirmed
  residue, likely superseded implementations, and risk assessment.
- [Current repository map](current_repository_map.md) — the repository as it is
  now, including its nested include graph and misplaced material.
- [Proposed repository map](proposed_repository_map.md) — a target layout and
  dependency direction.
- [Migration plan](migration_plan.md) — staged changes, validation gates, and a
  conservative disposition table.
- [Crash reproduction notes](crash_reproduction.md) — what was replayed, the
  measured resource peak, and what the replay did and did not establish.

## Executive conclusion

DolphinSpeaker has a useful bioacoustics pipeline buried inside three different
kinds of repository debt:

1. an incomplete removal of the PkgTemplates source repository;
2. production, research, demo, and scratch code mixed together under `src/`;
3. a cyclic/repeated include structure that prevents normal package
   precompilation and needlessly inflates load cost.

The safest path is not a rewrite. Preserve the active audio, detection,
synchronization, localization, annotation, and video pipeline; first replace the
foreign test/docs/CI shell, then characterize current behavior, make the package
entry point the sole include coordinator, and only then retire duplicate code.

## Important boundary

The working tree was already dirty. In particular, current edits and untracked
experiments under `src/`, `scripts/`, and the root were treated as user work.
They were inventoried but not altered, cleaned, staged, or executed.
