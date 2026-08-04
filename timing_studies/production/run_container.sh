#!/usr/bin/env bash
# Run a command inside the mucoll-sim container with the study area mounted.
# Usage: run_container.sh <command ...>
# Works with apptainer/singularity (clusters) or docker (falls back).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"

# The whole study directory is mounted at the same absolute path inside the
# container, so file names in arguments work unchanged on both sides.
# The image entrypoint execs an interactive shell and would swallow the
# command, so bypass it; /opt/setup_mucoll.sh loads the mucoll stack.
# In this image the CaloHitSelector threshold maps ship with MyBIBUtils, not
# k4Reco, so point MAIAConfig's locator (MUCOLL_CALO_THRESHOLDS_DIR) at them.
PRELUDE='source /opt/setup_mucoll.sh >/dev/null 2>&1;
for d in /opt/spack/opt/spack/*/*/*/*/linux-x86_64/mybibutils-*/share/MyBIBUtils/data; do
    [ -d "$d" ] && export MUCOLL_CALO_THRESHOLDS_DIR="$d";
done;'
if command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1; then
    APP=$(command -v apptainer || command -v singularity)
    exec "$APP" exec --cleanenv -B "$STUDY_DIR" "docker://${MUCOLL_IMAGE}" \
        /bin/bash -c "$PRELUDE $*"
else
    exec docker run --rm --entrypoint /bin/bash --user "$(id -u):$(id -g)" \
        -v "$STUDY_DIR":"$STUDY_DIR" -w "$STUDY_DIR" "${MUCOLL_IMAGE}" \
        -c "$PRELUDE $*"
fi
