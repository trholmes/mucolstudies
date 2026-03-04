#!/usr/bin/env python3
"""
Make momentum-distribution plots per PDGID for BIB particles from FLUKA binary input.

Outputs:
  1) For each selected PDGID: 1D histogram of |p|
  2) For each selected PDGID: 2D histogram of |p| vs parent muon z (z_mu)

Selections:
  - radial cut: keep only particles with r = sqrt(x^2+y^2) >= --r-min
  - PDG selection: either explicitly via --pdg-list, or automatically choose top-K by rate after cuts

Assumptions (same as earlier scripts):
  - Record layout matches LINE_DT below
  - E is kinetic energy [GeV]
  - Momentum magnitude |p| computed from kinetic energy T and mass m:
        |p| = sqrt(T^2 + 2 T m)  (c=1)
  - fid maps to PDGID via bib_pdgs.FLUKA_PIDS
  - mass is bib_pdgs.PDG_PROPS[pdgid] = (charge, mass)

Examples:
  python plotPdgMomentum.py /data/.../summary*.dat --outdir pdg_p --r-min 50
  python plotPdgMomentum.py file.dat --pdg-list 11 -11 13 -13 22 2112 --r-min 50
  python plotPdgMomentum.py file.dat --topk 15 --r-min 100 --z-round 3 --logz
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


def safe_tag(x: float) -> str:
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    s = s.replace("-", "m").replace(".", "p")
    return s


def pdg_tag(pid: int) -> str:
    return f"pdg_{pid}"


def plot_hist_p(p, outpath: Path, title: str, bins: int, logy: bool):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(p, bins=bins, histtype="step", linewidth=1.5)
    ax.set_xlabel("|p| [GeV/c]")
    ax.set_ylabel("BIB particles (count)")
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
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
                    help="Explicit list of PDGIDs to plot (e.g. --pdg-list 11 -11 13 -13 22)")
    ap.add_argument("--topk", type=int, default=20,
                    help="If --pdg-list not provided, plot only top-K PDGIDs by count after cuts (default: 20)")
    ap.add_argument("--min-count", type=int, default=1,
                    help="If auto-selecting PDGs, require at least this many entries (default: 1)")

    # Plot controls
    ap.add_argument("--bins-p", type=int, default=160, help="Bins for |p| in 1D plots")
    ap.add_argument("--bins-z", type=int, default=200, help="Bins in z_mu for 2D plots")
    ap.add_argument("--bins-p2d", type=int, default=200, help="Bins in |p| for 2D plots")
    ap.add_argument("--logy", action="store_true", help="Log-y for 1D histograms")
    ap.add_argument("--logz", action="store_true", help="Log color scale for 2D histograms")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # First pass: determine which PDGs to plot (unless user supplies list)
    pdg_counts = {}
    n_read = 0

    def process_chunk_for_counts(chunk):
        nonlocal pdg_counts

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

        # keep only PDGs with known mass (so momentum is well-defined here)
        known = np.array([int(pid) in PDG_PROPS for pid in pdg], dtype=bool)
        if not np.any(known):
            return
        pdg = pdg[known]

        uniq, cnt = np.unique(pdg, return_counts=True)
        for u, c in zip(uniq, cnt):
            pdg_counts[int(u)] = pdg_counts.get(int(u), 0) + int(c)

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

    if not pdg_counts:
        raise RuntimeError("No particles survived selection (or no known-mass PDGs).")

    # Decide which PDGs to plot
    if args.pdg_list is not None:
        pdgs_to_plot = [int(x) for x in args.pdg_list]
    else:
        items = [(pid, c) for pid, c in pdg_counts.items() if c >= args.min_count]
        items.sort(key=lambda t: t[1], reverse=True)
        pdgs_to_plot = [pid for pid, _ in items[:args.topk]]

    # Write a CSV of PDG counts after cuts for transparency
    csv_path = outdir / "pdg_counts_after_cuts.csv"
    with open(csv_path, "w") as f:
        f.write("pdg_id,count_after_cuts\n")
        for pid, c in sorted(pdg_counts.items(), key=lambda t: t[1], reverse=True):
            f.write(f"{pid},{c}\n")

    print(f"Wrote: {csv_path}")
    print(f"Selected {len(pdgs_to_plot)} PDGIDs to plot:", pdgs_to_plot)

    # Second pass: collect arrays per PDG for plotting
    # Store as lists of arrays to avoid huge reallocs; concatenate at the end.
    data_p = {pid: [] for pid in pdgs_to_plot}
    data_zmu = {pid: [] for pid in pdgs_to_plot}

    n_read2 = 0
    pdg_set = set(pdgs_to_plot)

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

            # Keep only selected PDGs
            sel = np.isin(pdg, pdgs_to_plot)
            if not np.any(sel):
                n_read2 += chunk.size
                continue

            sub = sub[sel]
            zmu = zmu[sel]
            pdg = pdg[sel]

            # mass known filter (should always be true for PDGs we counted, but keep safe)
            known = np.array([int(pid) in PDG_PROPS for pid in pdg], dtype=bool)
            if not np.any(known):
                n_read2 += chunk.size
                continue
            sub = sub[known]
            zmu = zmu[known]
            pdg = pdg[known]

            e_kin = sub["E"].astype(np.float64)
            masses = np.array([float(PDG_PROPS[int(pid)][1]) for pid in pdg], dtype=np.float64)

            p_mag = np.sqrt(e_kin * e_kin + 2.0 * e_kin * masses)

            # Split into per-PDG buffers
            for pid in np.unique(pdg):
                pid = int(pid)
                if pid not in pdg_set:
                    continue
                m = (pdg == pid)
                data_p[pid].append(p_mag[m])
                data_zmu[pid].append(zmu[m])

            n_read2 += chunk.size

        if args.max_lines is not None and n_read2 >= args.max_lines:
            break

    # Plot per PDG
    for pid in pdgs_to_plot:
        if not data_p[pid]:
            print(f"Skipping PDG {pid}: no entries after cuts.")
            continue

        p = np.concatenate(data_p[pid])
        zmu = np.concatenate(data_zmu[pid])

        # Output filenames
        base = outdir / pdg_tag(pid)
        title_1d = f"|p| distribution for PDGID {pid} (N={p.size}, r>={args.r_min} cm)"
        title_2d = f"|p| vs parent z_mu for PDGID {pid} (N={p.size}, r>={args.r_min} cm)"

        plot_hist_p(
            p,
            outpath=Path(str(base) + "_p.png"),
            title=title_1d,
            bins=args.bins_p,
            logy=args.logy,
        )

        plot_hist2_p_vs_zmu(
            zmu,
            p,
            outpath=Path(str(base) + "_p_vs_zmu.png"),
            title=title_2d,
            bins_z=args.bins_z,
            bins_p=args.bins_p2d,
            logz=args.logz,
        )

        print(f"Wrote: {base}_p.png")
        print(f"Wrote: {base}_p_vs_zmu.png")

    # Summary text file (handy for quick scanning)
    summary = outdir / "summary.txt"
    with open(summary, "w") as f:
        f.write(f"r_min_cm={args.r_min}\n")
        if args.z_round is not None:
            f.write(f"z_round={args.z_round}\n")
        f.write(f"pdgs_plotted={pdgs_to_plot}\n")
        f.write("pdg_id,N\n")
        for pid in pdgs_to_plot:
            n = sum(arr.size for arr in data_p[pid]) if data_p[pid] else 0
            f.write(f"{pid},{n}\n")
    print(f"Wrote: {summary}")


if __name__ == "__main__":
    main()
