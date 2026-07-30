#!/usr/bin/env bash
# Extract the ECAL CaloHitSelector threshold map (th_2dmode_sym) from the ROOT
# file shipped inside the container, into data/thresholds_ecal.npz for the
# python emulation. Also copies the HCAL map for reference.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"
mkdir -p "$DATA_DIR"

"$HERE/run_container.sh" "python3 - <<'EOF'
import glob, os
import numpy as np
import ROOT

out = {}
for det in ('ECAL', 'HCAL'):
    fname = f'{det}_Thresholds_10TeV.root'
    # same search logic as MAIAConfig Common/calo_thresholds.py
    cands = []
    for inc in os.environ.get('ROOT_INCLUDE_PATH', '').split(os.pathsep):
        if 'k4reco' in inc.lower():
            cands.append(os.path.join(os.path.dirname(inc), 'share', 'k4Reco', 'data', fname))
    stack = os.environ.get('MUCOLL_STACK', '')
    if stack:
        arch = os.path.dirname(os.path.dirname(stack))
        cands += glob.glob(os.path.join(arch, 'k4reco-*', 'share', 'k4Reco', 'data', fname))
    cands += glob.glob(f'/opt/**/{fname}', recursive=True)
    path = next((c for c in cands if os.path.isfile(c)), None)
    if path is None:
        print(f'WARNING: {fname} not found, skipping')
        continue
    print(f'{det}: {path}')
    f = ROOT.TFile.Open(path)
    for hname in ('th_2dmode_sym', 'stddev_sym'):
        h = f.Get(hname)
        if not h:
            continue
        nx, ny = h.GetNbinsX(), h.GetNbinsY()
        xedges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
        yedges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
        # include under/overflow so FindBin semantics can be replicated exactly
        vals = np.array([[h.GetBinContent(ix, iy) for iy in range(0, ny + 2)]
                         for ix in range(0, nx + 2)])
        out[f'{det}_{hname}_xedges'] = xedges
        out[f'{det}_{hname}_yedges'] = yedges
        out[f'{det}_{hname}_values'] = vals
    f.Close()

np.savez('$DATA_DIR/thresholds_ecal.npz', **out)
print('wrote $DATA_DIR/thresholds_ecal.npz with keys:', sorted(out))
EOF"
