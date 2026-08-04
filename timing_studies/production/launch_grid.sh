#!/usr/bin/env bash
# Produce every grid point defined in config.sh, NPAR points at a time,
# and record a manifest. Usage: launch_grid.sh [particles...]
#   launch_grid.sh gamma        # phase 1
#   launch_grid.sh pi-          # phase 2
#   launch_grid.sh              # everything in GRID_PARTICLES
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

PARTICLES="${*:-$GRID_PARTICLES}"
NPAR="${NPAR:-3}"

mkdir -p "$DATA_DIR"
MANIFEST="$DATA_DIR/production_manifest.tsv"
[ -f "$MANIFEST" ] || echo -e "particle\tenergy_gev\tnevents\tseed\tstem" > "$MANIFEST"

pids=()
for p in $PARTICLES; do
    for e in $GRID_ENERGIES; do
        stem="$(file_stem "$p" "$e")"
        grep -q "	${stem}$" "$MANIFEST" || \
            echo -e "$p\t$e\t$GRID_NEVENTS\t$(point_seed "$p" "$e")\t$stem" >> "$MANIFEST"
        "$HERE/run_chain.sh" "$p" "$e" &
        pids+=($!)
        # throttle to NPAR concurrent chains
        while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n; done
    done
done
rc=0
for pid in "${pids[@]}"; do wait "$pid" || rc=1; done
exit $rc
