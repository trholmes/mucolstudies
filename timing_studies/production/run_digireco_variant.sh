#!/usr/bin/env bash
# Re-run digitization + reconstruction on an EXISTING sim file with modified
# timing-cut settings (the k4run CLI exposes every algorithm property).
# Usage: run_digireco_variant.sh <sim_file> <cuttag> [extra k4run digi overrides...]
#   e.g.: run_digireco_variant.sh data/sim/photonGun_E10GeV_..._sim.edm4hep.root \
#           selwin0p15 \
#           --EcalBarrelSelector.TimeWindowMin=-0.15 --EcalBarrelSelector.TimeWindowMax=0.15 \
#           --EcalEndcapSelector.TimeWindowMin=-0.15 --EcalEndcapSelector.TimeWindowMax=0.15 \
#           --HcalBarrelSelector.TimeWindowMin=-0.15 --HcalBarrelSelector.TimeWindowMax=0.15 \
#           --HcalEndcapSelector.TimeWindowMin=-0.15 --HcalEndcapSelector.TimeWindowMax=0.15
# All overrides are applied to the DIGI step (both timing cuts live there);
# the reco step runs with defaults on the variant digi output.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

SIM="$(readlink -f "$1")"
CUTTAG="$2"
shift 2
OVERRIDES="$*"

STEM="$(basename "$SIM" _sim.edm4hep.root)"
mkdir -p "$DATA_DIR/reco_variants" "$DATA_DIR/logs"
DIGI="$DATA_DIR/reco_variants/${STEM}_digi_${CUTTAG}.edm4hep.root"
RECO="$DATA_DIR/reco_variants/${STEM}_reco_${CUTTAG}.edm4hep.root"
LOG="$DATA_DIR/logs/${STEM}_digireco_${CUTTAG}.log"

# Number of events: read from the stem (n<NEVENTS>), fall back to config
NEVENTS=$(echo "$STEM" | sed -n 's/.*_n\([0-9]*\)_.*/\1/p')
NEVENTS="${NEVENTS:-$GRID_NEVENTS}"

SETUP="source $BENCHMARKS_DIR/setup_config.sh $BENCHMARKS_DIR MAIA_v0 >/dev/null"

"$HERE/run_container.sh" "$SETUP && cd \$MUCOLL_CONFIG/\$MUCOLL_CONFIG_NAME && \
    k4run digi_steer.py -n $NEVENTS \
        --inputFiles $SIM --outputFile $DIGI \
        --histoFile $DATA_DIR/reco_variants/${STEM}_digi_${CUTTAG}_hists.root \
        $OVERRIDES \
    && k4run reco_steer.py -n $NEVENTS \
        --inputFiles $DIGI --outputFile $RECO \
        --histoFile $DATA_DIR/reco_variants/${STEM}_reco_${CUTTAG}_hists.root" \
    > "$LOG" 2>&1
echo "variant done: $RECO"
