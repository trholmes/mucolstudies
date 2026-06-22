# Overlay histograms from one or more plotMET output files, e.g. BIB vs. no BIB.
# The heavy lifting (the event loop) happens in plotMET.cpp; this just plots.
#
# Usage: python compareMET.py <label1>:<file1.root> [<label2>:<file2.root> ...]
# e.g.:  python compareMET.py "no BIB":met_noBIB.root "BIB":met_BIB.root

import sys
import ROOT

sys.path.append("..")
import plotHelper

ROOT.gROOT.SetBatch()

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

# Parse label:file pairs from the command line
inputs = []
for arg in sys.argv[1:]:
    label, fname = arg.split(":", 1)
    inputs.append((label, ROOT.TFile.Open(fname)))

hist_names = ["met", "met_low", "met_phi", "met_x", "met_y", "sumet", "npfo",
              "pmiss", "pmiss_low", "pmiss_z", "sume"]
for hname in hist_names:
    h_map = {}
    for label, f in inputs:
        h = f.Get(hname)
        if not h: continue
        h = h.Clone(f"{hname}_{label}")
        h.SetDirectory(0)
        if h.Integral() > 0: h.Scale(1.0/h.Integral())   # Normalize for shape comparison
        h_map[label] = h
    if not h_map:
        print(f"No histograms found for {hname}, skipping.")
        continue
    xlabel = list(h_map.values())[0].GetXaxis().GetTitle()
    plotHelper.plotHistograms(h_map, f"plots/compare_{hname}.png", xlabel, "Fraction of events", logy=True)
