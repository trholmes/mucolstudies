# Shared configuration for the MAIA calorimeter timing-cut study.
# Source this from the production scripts; every knob lives here.

# --- software -----------------------------------------------------------
export MUCOLL_IMAGE="ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v3.0"
# mucoll-benchmarks checkout (cloned by setup_benchmarks.sh, gitignored).
# MAIAConfig comes in as a submodule; commits are recorded in data/benchmarks_commits.txt
export BENCHMARKS_BRANCH="main"
# MAIAConfig@main has tracking properties newer than the v3.0 container's
# k4ActsTracking (CKF_Chi2CutOffOutlier etc.), so pin the newest submodule
# commit contemporary with the image (2026-06-23). The calo digi/selector
# settings are identical to main at this commit.
export MAIACONFIG_COMMIT="d290d02"

# --- paths --------------------------------------------------------------
# Absolute path of the timing_studies directory (works when sourced from anywhere)
export STUDY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA_DIR="${STUDY_DIR}/data"
export BENCHMARKS_DIR="${DATA_DIR}/mucoll-benchmarks"
# Everything is mounted at the same path inside the container, so host and
# container commands can share file names.

# --- sample grid (per study spec: 100 events, theta = 90 deg only) -------
export GRID_NEVENTS=100
export GRID_THETA=90
export GRID_ENERGIES="10 50 200"        # GeV
export GRID_PARTICLES="gamma pi-"       # phase 1: gamma ; phase 2: pi-
export SEED_BASE=42

# Human-readable particle names for file naming
particle_tag() {
    case "$1" in
        gamma) echo "photonGun" ;;
        pi-)   echo "pionGun" ;;
        *)     echo "$(echo "$1" | tr -cd '[:alnum:]')Gun" ;;
    esac
}

# Deterministic per-point seed: base + index of (particle, energy) in the grid
point_seed() { # args: particle energy
    local p="$1" e="$2" idx=0 pp ee
    for pp in $GRID_PARTICLES; do
        for ee in $GRID_ENERGIES; do
            idx=$((idx + 1))
            if [ "$pp" = "$p" ] && [ "$ee" = "$e" ]; then
                echo $((SEED_BASE + 10000 * idx)); return
            fi
        done
    done
    echo "$SEED_BASE"
}

# Canonical file stem, e.g. photonGun_E10GeV_th90_n100_s10042
file_stem() { # args: particle energy
    echo "$(particle_tag "$1")_E${2}GeV_th${GRID_THETA}_n${GRID_NEVENTS}_s$(point_seed "$1" "$2")"
}
