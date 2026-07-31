# Deployment pack: k4Reco calorimeter fixes → upstream PR

Two verified fixes for [MuonColliderSoft/k4Reco](https://github.com/MuonColliderSoft/k4Reco),
found during the timing-cut study (see `../PLAN.md` §1, notes 2 and 6):

1. **`CaloHitSelector.cpp`** — remove the pseudo time-of-flight correction
   (`r/TMath::C()` with `r` in mm and `C` in m/s ≈ 7×10⁻⁶ ns): the digi hit time
   is already TOF-corrected, so a *correct* second correction would double-count.
2. **`RealisticCaloReco.cc`** — `int cellid = hit->getCellID()` truncates the
   64-bit cellID; every Rec/Coned/Sel hit loses the fields above bit 31
   (`x`/`y` in the MAIA encoding). Fixed with `std::uint64_t`.

Prepared against upstream `main` @ `6282f2a884948417c66885d4a8e6de7148040bf3`
(2026-07-10). Validation ran in `ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v3.0`
with MAIAConfig @ `d290d02`; results in `validation/` and summarized in
`PR_BODY.md` (photon sample bit-identical; pion sample differs by exactly two
hits within 7×10⁻⁶ ns of the +0.3 ns selector edge, now handled correctly).

## Contents

- `k4reco-calo-fixes.patch` — the two-file diff (`git apply`-able)
- `PR_BODY.md` — ready-to-use PR title and body (includes the validation summary)
- `validation/*_compare.log` — full stock-vs-patched comparison output
  (produced by `../analysis/compare_reco_files.py`)
- `validation/pion10_edge_hits.log` — the two edge-window hits, explained

## Instructions for a fresh Claude Code session (or a human)

This session could not push to k4Reco (repo-scope limitation), so hand this
directory to a session that can:

1. Fork `MuonColliderSoft/k4Reco` to your account (GitHub UI or
   `gh repo fork MuonColliderSoft/k4Reco`), and start the session with that fork
   as its repository (so it has push access).
2. In the fork (synced with upstream `main`):
   ```bash
   git checkout -b fix/calo-selector-tof-and-cellid origin/main
   git apply k4reco-calo-fixes.patch     # from this directory
   git commit -am "Fix CaloHitSelector pseudo TOF correction and RealisticCaloReco cellID truncation"
   git push -u origin fix/calo-selector-tof-and-cellid
   ```
3. Open a PR against `MuonColliderSoft/k4Reco` `main` with the title and body
   from `PR_BODY.md` (keep the validation section and the attribution footer).
4. If upstream `main` has moved past `6282f2a` and the patch no longer applies
   cleanly, rebase the two hunks by hand (they are small) and mention in the PR
   that validation was performed at `6282f2a`.

## Re-running the validation (optional, ~15 min on 4 cores)

From `timing_studies/` with docker or apptainer available:

```bash
./production/setup_benchmarks.sh
# build the patched k4Reco inside the container
./production/run_container.sh "cd data/k4Reco && cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build -j4"
# re-run digi+reco on an existing sim file with the patched plugins preloaded
# (prepend build dir to GAUDI_PLUGIN_PATH and LD_LIBRARY_PATH; same --RandSeed),
# then compare:
./production/run_container.sh "python3 analysis/compare_reco_files.py \
    data/reco/<sample>_reco_default.edm4hep.root data/patched/<sample>_reco.edm4hep.root"
```
