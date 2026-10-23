#!/usr/bin/env bash
# All multi-arrival methods (rounds 1-4) on one recording, as for New_device: viewer CSVs
# detection_pixels__color255_webui__Impulse7..11*.csv next to the combined video, per-click tables and summaries in
# improvements/tdoa_localization_2026-10-22/results/, and optionally Audacity labels.
# Prerequisite: the multitdoa result folder of the recording (run_contiguous_folders_multitdoa; see how_to_run_on_recordings.md).
# usage (from the repository root):
#   bash improvements/tdoa_localization_2026-10-22/scripts/run_on_recording.sh <res_folder> [--audacity] [--n N] [--threads T]
set -euo pipefail
RES=""; AUD=0; N=2; T=12
while [ $# -gt 0 ]; do
  case "$1" in
    --audacity) AUD=1 ;;
    --n) N="$2"; shift ;;
    --threads) T="$2"; shift ;;
    *) RES="$1" ;;
  esac
  shift
done
[ -d "$RES" ] || { echo "usage: $0 <res_folder> [--audacity] [--n N] [--threads T]"; exit 1; }
S=improvements/tdoa_localization_2026-10-22/scripts
J="julia --project=. --startup-file=no -t $T"
echo "== round 1: SRP candidates + tracks (Impulse7_srpTracks, Impulse8_srpCandidates)"
$J $S/run_multiarrival_newdevice.jl "$RES"
echo "== round 2: SIC candidates, windows -50:50 (sic) and -110:50 (sicw) (Impulse7/8/9 sic, sicw)"
$J $S/run_multiarrival2_newdevice.jl sic "$RES"
$J $S/run_multiarrival2_newdevice.jl sicw "$RES"
echo "== round 2: track-rule selection (Impulse9sict / Impulse9sicwt _trackSelected)"
$J $S/write_track_selection_overlay.jl sic "$RES"
$J $S/write_track_selection_overlay.jl sicw "$RES"
echo "== round 3: longest track + followers inherit (Impulse10sicw_default_longesttrack-inherit)"
$J $S/apply_rules_real.jl sicw "default" "longest track + followers inherit" "$RES"
echo "== round 4: dolphin separation, n = $N, tuned parameters (Impulse11sicw_sources_n${N}_tuned, _auto_tuned)"
$J $S/apply_sources_real.jl "$N" tuned "$RES"
if [ "$AUD" = 1 ]; then
  echo "== Audacity labels"
  $J $S/export_audacity.jl "$RES" "$N"
fi
echo "done: viewer CSVs in $RES, tables in improvements/tdoa_localization_2026-10-22/results/"
