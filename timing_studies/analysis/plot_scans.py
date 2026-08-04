#!/usr/bin/env python3
"""Scan plots from an emulate.py results file: retained energy fraction,
response and resolution vs the timing-cut values.

The denominator ("no timing shaping" reference) is the same emulation with
the widest windows in the scan (digi max = last scan point, selector window
infinite), so thresholds are applied identically in numerator and denominator
and only the timing cuts differ.

Usage:
  python3 plot_scans.py data/results/scan_photons.npz --outdir plots/ecal_photons \
      [--regions EcalBarrel HcalBarrel] [--combined]
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


def sample_sort_key(stem):
    m = re.search(r"E(\d+)GeV", stem)
    return int(m.group(1)) if m else 0


def get(d, stem, region, dmin, dmax, tag):
    key = f"{stem}__sel_{region}__dmin{dmin}__dmax{dmax}__{tag}__e"
    return d.get(key)


def collect(d, stem, regions, dmin, dmax, tag):
    """Summed selected energy across regions, per event; None if missing."""
    tot = None
    for r in regions:
        arr = get(d, stem, r, dmin, dmax, tag)
        if arr is None:
            continue
        tot = arr.copy() if tot is None else tot + arr
    return tot


def robust_stats(x):
    med = np.median(x)
    lo, hi = np.percentile(x, [16, 84])
    return med, 0.5 * (hi - lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--regions", nargs="+", default=["EcalBarrel"],
                    help="regions summed for the retained-energy metric")
    ap.add_argument("--tag", default="", help="suffix for output file names")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    d = dict(np.load(args.results))
    dmins = d["digi_min_scan"]
    dmaxs = d["digi_max_scan"]
    stems = sorted({k.split("__")[0] for k in d if "__sel_" in k},
                   key=sample_sort_key)
    dmin0 = -0.5 if -0.5 in dmins else dmins[0]
    dmax_ref = dmaxs[-1]
    suffix = f"_{args.tag}" if args.tag else ""
    rlabel = "+".join(args.regions)

    def halfwidths():
        out = []
        for k in d:
            m = re.search(r"__sym([0-9.inf]+)__e$", k)
            if m:
                out.append(float(m.group(1)))
        return sorted(set(out))

    def selmaxes():
        out = []
        for k in d:
            m = re.search(r"__asym([0-9.inf]+)__e$", k)
            if m:
                out.append(float(m.group(1)))
        return sorted(set(out))

    # ---- 1) retained fraction vs selector half-width (digi window default)
    for scan, tagfmt, xlabel, fname in (
        (halfwidths(), "sym{}", "selector half-width w [ns]  (cut: |t| < w)",
         f"retained_vs_selwidth{suffix}.png"),
        (selmaxes(), "asym{}", "selector max [ns]  (cut: -0.3 < t < max)",
         f"retained_vs_selmax{suffix}.png"),
    ):
        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(7, 7), sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1]})
        for stem in stems:
            ref = collect(d, stem, args.regions, dmin0, dmax_ref, "syminf")
            if ref is None:
                continue
            ok = ref > 0
            xs, med, sig = [], [], []
            for w in scan:
                tag = tagfmt.format(w if w != np.inf else "inf")
                num = collect(d, stem, args.regions, dmin0, 10.0, tag)
                if num is None:
                    continue
                frac = num[ok] / ref[ok]
                m, s = robust_stats(frac)
                xs.append(w if np.isfinite(w) else max(
                    [v for v in scan if np.isfinite(v)]) * 2)
                med.append(m)
                sig.append(s)
            line, = ax.plot(xs, med, "o-", ms=4, label=sample_label(stem))
            axr.plot(xs, sig, "o-", ms=4, color=line.get_color())
        for a in (ax, axr):
            a.axvline(0.3, color="gray", ls="--", lw=1)
            a.set_xscale("log")
        ax.axhline(1.0, color="gray", lw=0.5)
        ax.axhline(0.99, color="gray", lw=0.5, ls=":")
        ax.set_ylabel(f"median retained E fraction ({rlabel})")
        ax.legend(fontsize=8, title="gun energy")
        ax.set_title("dashed line: current selector cut (0.3 ns); "
                     "denominator: no timing cuts")
        axr.set_ylabel("68% spread of fraction")
        axr.set_xlabel(xlabel)
        fig.tight_layout()
        out = os.path.join(args.outdir, fname)
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print("wrote", out)

    # ---- 2) retained fraction vs digi window max (selector at default and open)
    fig, ax = plt.subplots(figsize=(7, 5))
    for stem in stems:
        ref = collect(d, stem, args.regions, dmin0, dmax_ref, "syminf")
        if ref is None:
            continue
        ok = ref > 0
        for tag, ls, lab in (("syminf", "-", "no selector cut"),
                             ("sym0.3", "--", "selector +/-0.3 ns")):
            xs, med = [], []
            for dm in dmaxs:
                num = collect(d, stem, args.regions, dmin0, float(dm), tag)
                if num is None:
                    continue
                xs.append(dm)
                med.append(np.median(num[ok] / ref[ok]))
            line, = ax.plot(xs, med, ls, marker="o", ms=4,
                            label=f"{sample_label(stem)}, {lab}")
    ax.axvline(10.0, color="gray", ls=":", lw=1)
    ax.axhline(1.0, color="gray", lw=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("digi timingWindowMax [ns] (min = -0.5 ns)")
    ax.set_ylabel(f"median retained E fraction ({rlabel})")
    ax.legend(fontsize=7)
    ax.set_title("dotted line: current digi window max (10 ns)")
    fig.tight_layout()
    out = os.path.join(args.outdir, f"retained_vs_digimax{suffix}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)

    # ---- 3) 2D map: median retained fraction over (digi max, sel half-width)
    ws = [w for w in halfwidths() if np.isfinite(w)]
    for stem in stems:
        ref = collect(d, stem, args.regions, dmin0, dmax_ref, "syminf")
        if ref is None:
            continue
        ok = ref > 0
        grid = np.full((len(ws), len(dmaxs)), np.nan)
        for i, w in enumerate(ws):
            for j, dm in enumerate(dmaxs):
                num = collect(d, stem, args.regions, dmin0, float(dm), f"sym{w}")
                if num is not None:
                    grid[i, j] = np.median(num[ok] / ref[ok])
        fig, ax = plt.subplots(figsize=(7.5, 5))
        im = ax.pcolormesh(np.arange(len(dmaxs) + 1), np.arange(len(ws) + 1),
                           grid, vmin=max(np.nanmin(grid), 0.8), vmax=1.0,
                           cmap="viridis")
        ax.set_xticks(np.arange(len(dmaxs)) + 0.5,
                      [str(v) for v in dmaxs], fontsize=8)
        ax.set_yticks(np.arange(len(ws)) + 0.5, [str(v) for v in ws], fontsize=8)
        ax.set_xlabel("digi timingWindowMax [ns]")
        ax.set_ylabel("selector half-width [ns]")
        ax.set_title(f"median retained E fraction, {sample_label(stem)} ({rlabel})")
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        out = os.path.join(
            args.outdir, f"retained_2d_{sample_label(stem).replace(' ', '')}{suffix}.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
