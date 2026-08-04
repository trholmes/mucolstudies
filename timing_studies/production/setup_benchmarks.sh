#!/usr/bin/env bash
# Clone mucoll-benchmarks (with the MAIAConfig submodule) into data/ and record
# the exact commits used, so the study is reproducible.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"

mkdir -p "$DATA_DIR"
if [ ! -d "$BENCHMARKS_DIR/.git" ]; then
    git clone --branch "$BENCHMARKS_BRANCH" --recurse-submodules \
        https://github.com/MuonColliderSoft/mucoll-benchmarks.git "$BENCHMARKS_DIR"
fi
# pin MAIAConfig to the container-compatible commit (see config.sh)
git -C "$BENCHMARKS_DIR/configs/MAIAConfig" checkout -q "$MAIACONFIG_COMMIT"

{
    echo "# repositories used for this study (recorded $(date -u +%Y-%m-%dT%H:%M:%SZ))"
    echo "image: ${MUCOLL_IMAGE}"
    git -C "$BENCHMARKS_DIR" log -1 --format="mucoll-benchmarks: %H (%ci)"
    git -C "$BENCHMARKS_DIR/configs/MAIAConfig" log -1 --format="MAIAConfig: %H (%ci)"
} | tee "$DATA_DIR/benchmarks_commits.txt"
