#!/usr/bin/env python3
"""Distributions of the per-event retained-energy fraction f at a few selector
windows, one panel per gun energy — shows WHY the 68% spread vs window is
non-monotonic at low energy (bimodal population: f = 0 "missed shower" events
vs prompt-rich events), and how the effect fades with energy.

Usage:
  python3 plot_fraction_dists.py data/results/scan_pions.npz \
      --stems pionGun_E10GeV_th90_n100_s40042 pionGun_E50GeV_th90_n100_s50042 \
              pionGun_E200GeV_th90_n100_s60042 \
      --region HcalBarrel --out plots/hcal_pions/fraction_dist_HcalBarrel.png
"""

import argparse
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WINDOWS = [0.1, 0.3, 1.0, 5.0]
COLORS = ["#1f77b4", "#C0392B", "#2ca02c", "#7f7f7f"]


def label_of(stem):
    m = re.search(r"E(\d+)GeV", stem)
    return f"{m.group(1)} GeV" if m else stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--stems", nargs="+", required=True)
    ap.add_argument("--region", default="HcalBarrel")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    d = dict(np.load(args.results))
    n = len(args.stems)
    fig, axes = plt.subplots(1, n, figsize=(4.1 * n, 3.6), sharey=False)
    if n == 1:
        axes = [axes]
    bins = np.linspace(0, 1.0001, 21)

    for ax, stem in zip(axes, args.stems):
        ref = d[f"{stem}__sel_{args.region}__dmin-0.5__dmax100.0__syminf__e"]
        ok = ref > 0
        for w, c in zip(WINDOWS, COLORS):
            f = d[f"{stem}__sel_{args.region}__dmin-0.5__dmax10.0__sym{w}__e"][ok] / ref[ok]
            n0 = int(np.sum(f == 0))
            cur = "  (current)" if abs(w - 0.3) < 1e-9 else ""
            ax.hist(np.clip(f, 0, 1), bins=bins, histtype="step", lw=1.8, color=c,
                    label=f"±{w} ns{cur}: med {np.median(f):.2f}, f=0: {n0}")
        ax.set_xlabel("retained fraction f per event")
        ax.set_title(f"{label_of(stem)}  ({int(ok.sum())} events)", fontsize=11)
        ax.legend(fontsize=7.2, loc="upper center", title="selector window", title_fontsize=7.5)
        ax.set_xlim(0, 1)
    axes[0].set_ylabel("events / bin")
    fig.suptitle(f"{args.region}, π⁻: per-event retained energy fraction "
                 "f = E(window) / E(no timing cuts);  digi window fixed at default (−0.5, +10) ns",
                 fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(args.out, dpi=150)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
