#!/usr/bin/env python3
"""
Make momentum-distribution plots per PDGID for BIB particles from FLUKA binary input.

Updates in this version:
  - 1D plots are made PER |PDGID| and OVERLAY +PDG and -PDG contributions
    (e.g. +13 vs -13 on the same |p| histogram).
  - For neutrals (where + and - are the same, e.g. 22, 2112), you just get one curve.

Outputs:
  1) For each selected |PDGID|: 1D overlay histogram of |p| for +pid and -pid
  2) For each selected |PDGID|: 2D histogram of |p| vs parent muon z (z_mu)
     (this 2D plot uses BOTH charges together)

Selections:
  - radial cut: keep only particles with r = sqrt(x^2+y^2) >= --r-min
  - PDG selection:
      * if --pdg-list given: interpreted as a list of PDGIDs; internally plots by |PDGID|
      * else: auto-select top-K by COUNT after cuts, grouping by |PDGID|

Assumptions:
  - Record layout matches LINE_DT below
  - E is kinetic energy [GeV]
  - Momentum magnitude |p| computed from kinetic energy T and mass m:
        |p| = sqrt(T^2 + 2 T m)  (c=1)
  - fid maps to PDGID via bib_pdgs.FLUKA_PIDS
  - mass is bib_pdgs.PDG_PROPS[pdgid] = (charge, mass)

Examples:
  python plotPdgMomentum_overlay.py /data/.../summary*.dat --outdir pdg_p --r-min 50
  python plotPdgMomentum_overlay.py file.dat --pdg-list 11 -11 13 -13 22 2112 --r-min 50
  python plotPdgMomentum_overlay.py file.dat --topk 15 --r-min 100 --z-round 3 --logz
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from bib_pdgs import FLUKA_PIDS, PDG_PROPS


LINE_DT = np.dtype([
    ("fid",   np.int32),
    ("fid_mo", np.int32),
    ("E",     np.float64),  # kinetic energy [GeV]
    ("x",     np.float64),
    ("y",     np.float64),
    ("z",     np.float64),
    ("cx",    np.float64),
    ("cy",    np.float64),
    ("cz",    np.float64),
    ("time",  np.float64),
    ("x_mu",  np.float64),
    ("y_mu",  np.float64),
    ("z_mu",  np.float64),  # parent muon z
])


def iter_records(filename: str, chunk_records: int = 1_000_000):
    with open(filename, "rb") as f:
        while True:
            arr = np.fromfile(f, dtype=LINE_DT, count=chunk_records)
            if arr.size == 0:
                break
            yield arr


def pdg_abs_tag(pid_abs: int) -> str:
    return f"pdgAbs_{pid_abs}"


def plot_hist_p_overlay(p_pos, p_neg, outpath: Path, title: str, bins: int, logy: bool,
                        label_pos: str, label_neg: str):
    """
    Overlay +pid and -pid 1D |p| histograms. If one side is empty, it simply won't appear.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    plotted_any = False
    if p_pos is not None and p_pos.size > 0:
        ax.hist(p_pos, bins=bins, histtype="step", linewidth=1.5, label=label_pos)
        plotted_any = True
    if p_neg is not None and p_neg.size > 0:
        ax.hist(p_neg, bins=bins, histtype="step", linewidth=1.5, label=label_neg)
        plotted_any = True

    if not plotted_any:
        # Still write an empty figure with a message
        ax.text(0.5, 0.5, "No entries after cuts", ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("|p| [GeV/c]")
    ax.set_ylabel("BIB particles (count)")
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_hist2_p_vs_zmu(zmu, p, outpath: Path, title: str,
                        bins_z: int, bins_p: int, logz: bool):
    fig, ax = plt.subplots(figsize=(8, 6))
    if logz:
        h = ax.hist2d(zmu, p, bins=[bins_z, bins_p], norm=LogNorm(vmin=1))
        fig.colorbar(h[3], ax=ax, label="BIB particles / bin (log scale)")
    else:
        h = ax.hist2d(zmu, p, bins=[bins_z, bins_p])
        fig.colorbar(h[3], ax=ax, label="BIB particles / bin")
    ax.set_xlabel("Parent muon z_mu [cm]")
    ax.set_ylabel("|p| [GeV/c]")
    ax.set_title(title)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files_in", nargs="+", help="Input FLUKA binary file(s)")
    ap.add_argument("--outdir", default="pdg_momentum_plots", help="Output directory")
    ap.add_argument("--r-min", type=float, default=0.0,
                    help="Only include particles with r = sqrt(x^2+y^2) >= r-min [cm]")
    ap.add_argument("--z-round", type=int, default=None,
                    help="Round parent z_mu to this many decimals (helps if z is discrete)")
    ap.add_argument("--chunk-records", type=int, default=1_000_000, help="Chunk size for reading")
    ap.add_argument("--max-lines", type=int, default=None, help="Optional cap on total records read (debug)")

    # PDG selection
    ap.add_argument("--pdg-list", type=int, nargs="+", default=None,
                    help="Explicit list of PDGIDs to plot; internally grouped by |PDGID| "
                         "(e.g. --pdg-list 13 or --pdg-list 13 -13 are equivalent)")
    ap.add_argument("--topk", type=int, default=20,
                    help="If --pdg-list not provided, plot only top-K |PDGID| by count after cuts (default: 20)")
    ap.add_argument("--min-count", type=int, default=1,
                    help="If auto-selecting, require at least this many entries (default: 1)")

    # Plot controls
    ap.add_argument("--bins-p", type=int, default=160, help="Bins for |p| in 1D plots")
    ap.add_argument("--bins-z", type=int, default=200, help="Bins in z_mu for 2D plots")
    ap.add_argument("--bins-p2d", type=int, default=200, help="Bins in |p| for 2D plots")
    ap.add_argument("--logy", action="store_true", help="Log-y for 1D histograms")
    ap.add_argument("--logz", action="store_true", help="Log color scale for 2D histograms")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Pass 1: count PDGs after cuts (grouping by abs(PDG)) unless user gave explicit list ----
    counts_abs = {}  # abs(pdg) -> count after cuts
    n_read = 0

    def process_chunk_for_counts(chunk):
        nonlocal counts_abs

        # r cut
        x = chunk["x"].astype(np.float64)
        y = chunk["y"].astype(np.float64)
        r = np.sqrt(x * x + y * y)
        mask = (r >= args.r_min)
        if not np.any(mask):
            return

        sub = chunk[mask]

        # fid -> pdg
        fids = sub["fid"].astype(np.int32)
        pdg = np.zeros_like(fids, dtype=np.int32)
        for i, fid in enumerate(fids):
            pdg[i] = FLUKA_PIDS.get(int(fid), 0)

        # keep only PDGs with known mass (so momentum is well-defined)
        known = np.array([int(pid) in PDG_PROPS for pid in pdg], dtype=bool)
        if not np.any(known):
            return
        pdg = pdg[known]

        pdg_abs = np.abs(pdg)
        uniq, cnt = np.unique(pdg_abs, return_counts=True)
        for u, c in zip(uniq, cnt):
            u = int(u)
            counts_abs[u] = counts_abs.get(u, 0) + int(c)

    for fin in args.files_in:
        for chunk in iter_records(fin, chunk_records=args.chunk_records):
            if args.max_lines is not None:
                remaining = args.max_lines - n_read
                if remaining <= 0:
                    break
                if chunk.size > remaining:
                    chunk = chunk[:remaining]

            process_chunk_for_counts(chunk)
            n_read += chunk.size

        if args.max_lines is not None and n_read >= args.max_lines:
            break

    if not counts_abs:
        raise RuntimeError("No particles survived selection (or no known-mass PDGs).")

    # Decide which |PDG| to plot
    if args.pdg_list is not None:
        pdg_abs_to_plot = sorted(set(int(abs(x)) for x in args.pdg_list))
    else:
        items = [(pid_abs, c) for pid_abs, c in counts_abs.items() if c >= args.min_count]
        items.sort(key=lambda t: t[1], reverse=True)
        pdg_abs_to_plot = [pid_abs for pid_abs, _ in items[:args.topk]]

    # Write a CSV of abs(PDG) counts after cuts
    csv_path = outdir / "pdgAbs_counts_after_cuts.csv"
    with open(csv_path, "w") as f:
        f.write("pdg_abs,count_after_cuts\n")
        for pid_abs, c in sorted(counts_abs.items(), key=lambda t: t[1], reverse=True):
            f.write(f"{pid_abs},{c}\n")

    print(f"Wrote: {csv_path}")
    print(f"Selected {len(pdg_abs_to_plot)} |PDGID| to plot:", pdg_abs_to_plot)

    # ---- Pass 2: collect data per abs(PDG), separately for + and - charges (for 1D overlay) ----
    # For 2D, we’ll combine both charges together in one plot per abs(PDG).
    data_p_pos = {pid_abs: [] for pid_abs in pdg_abs_to_plot}
    data_p_neg = {pid_abs: [] for pid_abs in pdg_abs_to_plot}
    data_zmu_all = {pid_abs: [] for pid_abs in pdg_abs_to_plot}
    data_p_all = {pid_abs: [] for pid_abs in pdg_abs_to_plot}

    n_read2 = 0
    abs_set = set(pdg_abs_to_plot)

    for fin in args.files_in:
        for chunk in iter_records(fin, chunk_records=args.chunk_records):
            if args.max_lines is not None:
                remaining = args.max_lines - n_read2
                if remaining <= 0:
                    break
                if chunk.size > remaining:
                    chunk = chunk[:remaining]

            # r cut
            x = chunk["x"].astype(np.float64)
            y = chunk["y"].astype(np.float64)
            r = np.sqrt(x * x + y * y)
            mask = (r >= args.r_min)
            if not np.any(mask):
                n_read2 += chunk.size
                continue

            sub = chunk[mask]

            # parent z_mu (optionally rounded)
            zmu = sub["z_mu"].astype(np.float64)
            if args.z_round is not None:
                zmu = np.round(zmu, args.z_round)

            # fid -> pdg
            fids = sub["fid"].astype(np.int32)
            pdg = np.zeros_like(fids, dtype=np.int32)
            for i, fid in enumerate(fids):
                pdg[i] = FLUKA_PIDS.get(int(fid), 0)

            # keep only PDGs with known mass
            known = np.array([int(pid) in PDG_PROPS for pid in pdg], dtype=bool)
            if not np.any(known):
                n_read2 += chunk.size
                continue

            sub = sub[known]
            zmu = zmu[known]
            pdg = pdg[known]

            pdg_abs = np.abs(pdg)
            sel_abs = np.isin(pdg_abs, pdg_abs_to_plot)
            if not np.any(sel_abs):
                n_read2 += chunk.size
                continue

            sub = sub[sel_abs]
            zmu = zmu[sel_abs]
            pdg = pdg[sel_abs]
            pdg_abs = pdg_abs[sel_abs]

            e_kin = sub["E"].astype(np.float64)
            masses = np.array([float(PDG_PROPS[int(pid)][1]) for pid in pdg], dtype=np.float64)
            p_mag = np.sqrt(e_kin * e_kin + 2.0 * e_kin * masses)

            # Split per abs(PDG)
            for pid_abs in np.unique(pdg_abs):
                pid_abs = int(pid_abs)
                if pid_abs not in abs_set:
                    continue
                m_abs = (pdg_abs == pid_abs)
                if not np.any(m_abs):
                    continue

                # For 2D combined:
                data_p_all[pid_abs].append(p_mag[m_abs])
                data_zmu_all[pid_abs].append(zmu[m_abs])

                # For 1D overlay split by sign:
                pdg_here = pdg[m_abs]
                p_here = p_mag[m_abs]

                m_pos = (pdg_here == pid_abs)      # +pid
                m_neg = (pdg_here == -pid_abs)     # -pid

                if np.any(m_pos):
                    data_p_pos[pid_abs].append(p_here[m_pos])
                if np.any(m_neg):
                    data_p_neg[pid_abs].append(p_here[m_neg])

            n_read2 += chunk.size

        if args.max_lines is not None and n_read2 >= args.max_lines:
            break

    # ---- Plot per abs(PDG) ----
    for pid_abs in pdg_abs_to_plot:
        p_pos = np.concatenate(data_p_pos[pid_abs]) if data_p_pos[pid_abs] else np.array([], dtype=np.float64)
        p_neg = np.concatenate(data_p_neg[pid_abs]) if data_p_neg[pid_abs] else np.array([], dtype=np.float64)

        p_all = np.concatenate(data_p_all[pid_abs]) if data_p_all[pid_abs] else np.array([], dtype=np.float64)
        zmu_all = np.concatenate(data_zmu_all[pid_abs]) if data_zmu_all[pid_abs] else np.array([], dtype=np.float64)

        if p_all.size == 0:
            print(f"Skipping |PDG|={pid_abs}: no entries after cuts.")
            continue

        base = outdir / pdg_abs_tag(pid_abs)

        # 1D overlay (+pid_abs vs -pid_abs)
        # For neutrals: pid_abs == -pid_abs is not a thing in PDG, but neutrals just populate one side typically.
        title_1d = f"|p| for PDGID ±{pid_abs} (N={p_all.size}, r>={args.r_min} cm)"
        label_pos = f"+{pid_abs} (N={p_pos.size})"
        label_neg = f"-{pid_abs} (N={p_neg.size})"

        plot_hist_p_overlay(
            p_pos=p_pos,
            p_neg=p_neg,
            outpath=Path(str(base) + "_p_overlay.png"),
            title=title_1d,
            bins=args.bins_p,
            logy=args.logy,
            label_pos=label_pos,
            label_neg=label_neg,
        )

        # 2D combined (both charges together)
        title_2d = f"|p| vs parent z_mu for PDGID ±{pid_abs} (N={p_all.size}, r>={args.r_min} cm)"
        plot_hist2_p_vs_zmu(
            zmu_all,
            p_all,
            outpath=Path(str(base) + "_p_vs_zmu.png"),
            title=title_2d,
            bins_z=args.bins_z,
            bins_p=args.bins_p2d,
            logz=args.logz,
        )

        print(f"Wrote: {base}_p_overlay.png")
        print(f"Wrote: {base}_p_vs_zmu.png")

    # Summary text file
    summary = outdir / "summary.txt"
    with open(summary, "w") as f:
        f.write(f"r_min_cm={args.r_min}\n")
        if args.z_round is not None:
            f.write(f"z_round={args.z_round}\n")
        f.write(f"pdg_abs_plotted={pdg_abs_to_plot}\n")
        f.write("pdg_abs,N_total,N_pos,N_neg\n")
        for pid_abs in pdg_abs_to_plot:
            n_tot = sum(arr.size for arr in data_p_all[pid_abs]) if data_p_all[pid_abs] else 0
            n_pos = sum(arr.size for arr in data_p_pos[pid_abs]) if data_p_pos[pid_abs] else 0
            n_neg = sum(arr.size for arr in data_p_neg[pid_abs]) if data_p_neg[pid_abs] else 0
            f.write(f"{pid_abs},{n_tot},{n_pos},{n_neg}\n")
    print(f"Wrote: {summary}")


if __name__ == "__main__":
    main()
