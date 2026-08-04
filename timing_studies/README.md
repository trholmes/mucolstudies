# timing_studies

Study of the calorimeter timing cuts in the MAIA reconstruction chain
(what they remove from ECAL/HCAL clusters, and what cut values avoid shaping
the cluster energy). **See `PLAN.md` for the study plan, the verified
timing-cut inventory, and results.**

Software: `ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v3.0` +
`mucoll-benchmarks`/`MAIAConfig` (Gaudi chain from the
[mcd-wiki tutorial](https://mcd-wiki.web.cern.ch/software/howto/gaudi/)).
Needs docker or apptainer; analysis needs only numpy + matplotlib.

## Quickstart

```bash
cd timing_studies
./production/setup_benchmarks.sh          # clone steering (pinned in data/benchmarks_commits.txt)
./production/extract_thresholds.sh        # ECAL selector threshold map -> data/thresholds_ecal.npz
./production/run_chain.sh gamma 10 10     # 10-event smoke test
./production/launch_grid.sh gamma         # phase 1: photons E=10/50/200, 100 evts, theta=90

# flatten (inside the container) and check closure
./production/run_container.sh python3 analysis/flatten.py \
    data/reco/photonGun_E10GeV_th90_n100_s10042_reco_default.edm4hep.root \
    data/flat/photonGun_E10GeV_th90_n100_s10042.npz
python3 analysis/sanity_check.py data/flat/photonGun_E10GeV_*.npz \
    --thresholds data/thresholds_ecal.npz

# scan + plots
python3 analysis/emulate.py data/flat/photonGun_*.npz \
    --thresholds data/thresholds_ecal.npz --out data/results/scan_photons.npz
python3 analysis/plot_distributions.py data/results/scan_photons.npz --outdir plots/ecal_photons
python3 analysis/plot_scans.py data/results/scan_photons.npz --outdir plots/ecal_photons

# real-chain validation at one variant point (overrides both windows everywhere)
./production/run_digireco_variant.sh data/sim/photonGun_E50GeV_*_sim.edm4hep.root \
    selwin0p15 \
    --EcalBarrelSelector.TimeWindowMin=-0.15 --EcalBarrelSelector.TimeWindowMax=0.15 \
    --EcalEndcapSelector.TimeWindowMin=-0.15 --EcalEndcapSelector.TimeWindowMax=0.15 \
    --HcalBarrelSelector.TimeWindowMin=-0.15 --HcalBarrelSelector.TimeWindowMax=0.15 \
    --HcalEndcapSelector.TimeWindowMin=-0.15 --HcalEndcapSelector.TimeWindowMax=0.15
./production/run_container.sh python3 analysis/summarize_reco.py \
    data/reco_variants/photonGun_E50GeV_*_reco_selwin0p15.edm4hep.root \
    data/summaries/photonGun_E50GeV_selwin0p15.npz
python3 analysis/validate.py --flat data/flat/photonGun_E50GeV_*.npz \
    --summary data/summaries/photonGun_E50GeV_selwin0p15.npz \
    --thresholds data/thresholds_ecal.npz --sel-min -0.15 --sel-max 0.15
```

Everything under `data/` is gitignored and reproducible from the seeds in
`data/production_manifest.tsv`. Analysis after `flatten.py` is plain
numpy/matplotlib (deliberately not PyROOT, so it runs outside the container).
