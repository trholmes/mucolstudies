#!/usr/bin/env python3
"""Plot the energy-weighted delta_t spectra from an emulate.py results file.

One figure per calorimeter region: contribution-level and cell-earliest-level
spectra overlaid per sample (gun energy), with the current cut values marked.

Usage:
  python3 plot_distributions.py data/results/scan_photons.npz --outdir plots/ecal_photons
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import cutmodel as cm  # noqa: E402


def sample_label(stem):
    m = re.search(r"E(\d+)GeV", stem)
    return f"{m.group(1)} GeV" if m else stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    d = dict(np.load(args.results, allow_pickle=False))
    bins = d["dt_bins"]
    centers = 0.5 * (bins[:-1] + bins[1:])
    widths = np.diff(bins)

    stems = sorted({k.split("__")[0] for k in d if "__dthist_" in k},
                   key=lambda s: (len(s), s))

    LEVELS = (("contrib", "per contribution (no cuts)"),
              ("cellearliest_nodigi", "per cell, earliest contribution (no digi cut)"),
              ("cellearliest", "per cell (earliest accepted time, digi window (-0.5,+10) ns)"))

    # identical y scale across all levels of a region -> directly comparable
    ylims = {}
    for region in cm.REGIONS:
        lo, hi = np.inf, 0.0
        for level, _ in LEVELS:
            for stem in stems:
                key = f"{stem}__dthist_{level}_{region}"
                if key in d and d[key].sum() > 0:
                    h = d[key] / d[key].sum() / widths
                    pos = h[h > 0]
                    lo, hi = min(lo, pos.min()), max(hi, pos.max())
        if hi > 0:
            ylims[region] = (lo * 0.5, hi * 3)

    for level, title in LEVELS:
        for region in cm.REGIONS:
            fig, ax = plt.subplots(figsize=(7, 5))
            found = False
            for stem in stems:
                key = f"{stem}__dthist_{level}_{region}"
                if key not in d or d[key].sum() <= 0:
                    continue
                found = True
                total = d[key].sum()
                h = d[key] / total / widths  # normalized density
                line, = ax.step(centers, h, where="mid", label=sample_label(stem))
                key2 = f"{stem}__dthist2_{level}_{region}"
                if key2 in d:  # statistical error: sqrt(sum of squared weights)
                    err = np.sqrt(d[key2]) / total / widths
                    ax.errorbar(centers, h, yerr=err, fmt="none",
                                ecolor=line.get_color(), elinewidth=0.8,
                                capsize=0, alpha=0.7)
            if not found:
                plt.close(fig)
                continue
            for x, ls, lab in ((cm.SEL_WINDOW_MAX, "--", "selector +/-0.3 ns"),
                               (cm.DIGI_WINDOW_MAX, ":", "digi max 10 ns")):
                ax.axvline(x, color="gray", linestyle=ls, lw=1, label=lab)
            ax.axvline(cm.SEL_WINDOW_MIN, color="gray", linestyle="--", lw=1)
            ax.set_yscale("log")
            ax.set_xlabel(r"$\Delta t = t - r/c$ [ns]")
            ax.set_ylabel("energy-weighted density [1/ns]")
            ax.set_title(f"{region}, {title}", fontsize=11)
            ax.legend(fontsize=8)
            ax.set_xlim(bins[0], 25)
            if region in ylims:
                ax.set_ylim(*ylims[region])
            fig.tight_layout()
            out = os.path.join(args.outdir, f"dt_{level}_{region}.png")
            fig.savefig(out, dpi=150)
            plt.close(fig)
            print("wrote", out)


if __name__ == "__main__":
    main()
