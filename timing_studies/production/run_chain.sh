#!/usr/bin/env bash
# Produce one grid point with the default reconstruction settings:
#   particle gun -> ddsim (MAIA_v0) -> k4run digi -> k4run reco
# Usage: run_chain.sh <particle> <energy_GeV> [nevents] [seed]
#   e.g.: run_chain.sh gamma 10
# Output: data/sim/<stem>_sim.edm4hep.root, data/reco/<stem>_reco_default.edm4hep.root
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

PARTICLE="$1"
ENERGY="$2"
NEVENTS="${3:-$GRID_NEVENTS}"
SEED="${4:-$(point_seed "$PARTICLE" "$ENERGY")}"
STEM="$(particle_tag "$PARTICLE")_E${ENERGY}GeV_th${GRID_THETA}_n${NEVENTS}_s${SEED}"

mkdir -p "$DATA_DIR/sim" "$DATA_DIR/digi" "$DATA_DIR/reco" "$DATA_DIR/logs"
SIM="$DATA_DIR/sim/${STEM}_sim.edm4hep.root"
DIGI="$DATA_DIR/digi/${STEM}_digi_default.edm4hep.root"
RECO="$DATA_DIR/reco/${STEM}_reco_default.edm4hep.root"
LOG="$DATA_DIR/logs/${STEM}"

# All steps run inside the container; setup_config.sh exports MUCOLL_GEO and
# makes the MAIAConfig steering importable. k4run must run from the MAIAConfig
# package dir so the relative PandoraSettings path resolves.
SETUP="source $BENCHMARKS_DIR/setup_config.sh $BENCHMARKS_DIR MAIA_v0 >/dev/null"

if [ ! -s "$SIM" ]; then
    "$HERE/run_container.sh" "$SETUP && export BENCHMARKS_DIR=$BENCHMARKS_DIR && ddsim \
        --steeringFile $HERE/sim_gun_steer.py \
        --compactFile \$MUCOLL_GEO \
        --enableGun --gun.particle '$PARTICLE' \
        --gun.energy '${ENERGY}*GeV' --gun.distribution uniform \
        --gun.thetaMin '${GRID_THETA}*degree' --gun.thetaMax '${GRID_THETA}*degree' \
        --numberOfEvents $NEVENTS --random.seed $SEED --random.enableEventSeed \
        --outputFile $SIM" > "${LOG}_sim.log" 2>&1
    echo "sim done: $SIM"
else
    echo "sim exists, skipping: $SIM"
fi

if [ ! -s "$RECO" ]; then
    "$HERE/run_container.sh" "$SETUP && cd \$MUCOLL_CONFIG/\$MUCOLL_CONFIG_NAME && \
        k4run digi_steer.py -n $NEVENTS --RandSeed $SEED \
            --inputFiles $SIM --outputFile $DIGI --histoFile $DATA_DIR/digi/${STEM}_digi_hists.root \
        && k4run reco_steer.py -n $NEVENTS \
            --inputFiles $DIGI --outputFile $RECO --histoFile $DATA_DIR/reco/${STEM}_reco_hists.root" \
        > "${LOG}_digireco.log" 2>&1
    echo "digi+reco done: $RECO"
else
    echo "reco exists, skipping: $RECO"
fi
