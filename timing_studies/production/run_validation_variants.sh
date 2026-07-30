#!/usr/bin/env bash
# Run the real-chain validation points (PLAN.md section 2) on all existing sim
# files matching a pattern, then summarize each variant and delete the bulky
# variant output to save disk.
# Usage: run_validation_variants.sh [sim-file-glob]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

SEL_ALL="EcalBarrelSelector EcalEndcapSelector HcalBarrelSelector HcalEndcapSelector"
DIGI_ALL="ECalBarrelDigi ECalEndcapDigi HCalBarrelDigi HCalEndcapDigi"

sel_overrides() { # args: min max
    local out=""
    for s in $SEL_ALL; do out+=" --$s.TimeWindowMin=$1 --$s.TimeWindowMax=$2"; done
    echo "$out"
}
digi_overrides() { # args: max
    local out=""
    for d in $DIGI_ALL; do out+=" --$d.timingWindowMax=$1"; done
    echo "$out"
}

mkdir -p "$DATA_DIR/summaries"
for SIM in ${1:-$DATA_DIR/sim/*_sim.edm4hep.root}; do
    STEM="$(basename "$SIM" _sim.edm4hep.root)"
    while read -r CUTTAG OVERRIDES; do
        SUMMARY="$DATA_DIR/summaries/${STEM}_${CUTTAG}.npz"
        [ -s "$SUMMARY" ] && { echo "have $SUMMARY"; continue; }
        "$HERE/run_digireco_variant.sh" "$SIM" "$CUTTAG" $OVERRIDES
        RECO="$DATA_DIR/reco_variants/${STEM}_reco_${CUTTAG}.edm4hep.root"
        "$HERE/run_container.sh" "python3 $STUDY_DIR/analysis/summarize_reco.py $RECO $SUMMARY"
        rm -f "$DATA_DIR/reco_variants/${STEM}"_*"${CUTTAG}"*.root
    done <<EOF
selwin0p15 $(sel_overrides -0.15 0.15)
selwin1p0 $(sel_overrides -1.0 1.0)
digimax2 $(digi_overrides 2.0)
EOF
done
echo "validation variants done"
