#!/usr/bin/env python3
"""
Find parent-muon (PosZmu / z_mu) outliers with unusually many BIB particles,
then print summary and make diagnostic plots per outlier parent.

Added:
  --log   -> log-y for 1D hists, log color scale (log-z) for 2D hist2d.

Usage:
  python find_bib_parent_outliers.py file1.bin [file2.bin ...] --threshold 20000 --outdir outliers --log
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from bib_pdgs import FLUKA_PIDS, PDG_PROPS


LINE_DT = np.dtype([
    ('fid',   np.int32),
    ('fid_mo', np.int32),
    ('E',     np.float64),
    ('x',     np.float64),
    ('y',     np.float64),
    ('z',     np.float64),
    ('cx',    np.float64),
    ('cy',    np.float64),
    ('cz',    np.float64),
    ('time',  np.float64),
    ('x_mu',  np.float64),
    ('y_mu',  np.float64),
    ('z_mu',  np.float64),
])


def iter_records(filename: str, chunk_records: int = 1_000_000):
    with open(filename, "rb") as f:
        while True:
            arr = np.fromfile(f, dtype=LINE_DT, count=chunk_records)
            if arr.size == 0:
                break
            yield arr


def safe_tag(z_mu: float) -> str:
    s = f"{z_mu:.6f}".rstrip("0").rstrip(".")
    s = s.replace("-", "m").replace(".", "p")
    return s


def plot_momentum(p_mag, outpath: Path, title: str, bins: int = 120, logy: bool = False):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(p_mag, bins=bins, histtype="step", linewidth=1.5)
    ax.set_xlabel("Momentum magnitude |p| [GeV/c]")
    ax.set_ylabel("BIB particles (count)")
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_pdg_bar(pdg_ids, outpath: Path, title: str, topk: int = 25, logy: bool = False):
    uniq, counts = np.unique(pdg_ids, return_counts=True)
    order = np.argsort(counts)[::-1]
    uniq = uniq[order]
    counts = counts[order]

    if uniq.size > topk:
        uniq_show = uniq[:topk]
        counts_show = counts[:topk]
        other = counts[topk:].sum()
        labels = [str(int(u)) for u in uniq_show] + ["Other"]
        counts_show = np.append(counts_show, other)
    else:
        labels = [str(int(u)) for u in uniq]
        counts_show = counts

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x, counts_show)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_xlabel("PDGID")
    ax.set_ylabel("BIB particles (count)")
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_r_vs_z(r_vals, z_vals, outpath: Path, title: str,
                bins_r: int = 200, bins_z: int = 200, logz: bool = False):
    fig, ax = plt.subplots(figsize=(8, 6))

    if logz:
        # Use LogNorm for the color scale ("log z"). vmin=1 avoids log(0) issues from empty bins.
        h = ax.hist2d(
            z_vals, r_vals,
            bins=[bins_z, bins_r],
            norm=LogNorm(vmin=1)
        )
        cbar = fig.colorbar(h[3], ax=ax, label="BIB particles / bin (log scale)")
    else:
        h = ax.hist2d(z_vals, r_vals, bins=[bins_z, bins_r])
        cbar = fig.colorbar(h[3], ax=ax, label="BIB particles / bin")

    ax.set_xlabel("Particle z [cm]")
    ax.set_ylabel("Particle r = sqrt(x^2 + y^2) [cm]")
    ax.set_title(title)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files_in", nargs="+", help="Input FLUKA binary file(s)")
    ap.add_argument("--threshold", type=int, default=20000,
                    help="Flag parents with N(BIB) > threshold (default: 20000)")
    ap.add_argument("--outdir", default="bib_parent_outliers",
                    help="Output directory for plots (default: bib_parent_outliers)")
    ap.add_argument("--z-round", type=int, default=None,
                    help="Round z_mu to this many decimals when grouping")
    ap.add_argument("--chunk-records", type=int, default=1_000_000,
                    help="Chunk size for reading (default: 1,000,000)")
    ap.add_argument("--max-lines", type=int, default=None,
                    help="Optional maximum number of records to read per file (debug)")
    ap.add_argument("--log", action="store_true",
                    help="Log-y for 1D plots and log color scale (log-z) for 2D plots")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    any_outliers = False

    for fin in args.files_in:
        fin_path = Path(fin)
        #print(f"\n=== Scanning file: {fin_path} ===")

        # Pass 1: count BIB per parent z_mu in THIS FILE
        counts = {}
        n_read = 0

        for chunk in iter_records(fin, chunk_records=args.chunk_records):
            if args.max_lines is not None:
                remaining = args.max_lines - n_read
                if remaining <= 0:
                    break
                if chunk.size > remaining:
                    chunk = chunk[:remaining]

            z_mu = chunk["z_mu"].astype(np.float64)
            if args.z_round is not None:
                z_mu = np.round(z_mu, args.z_round)

            uz, c = np.unique(z_mu, return_counts=True)
            for z, cc in zip(uz, c):
                counts[z] = counts.get(z, 0) + int(cc)

            n_read += chunk.size

        if not counts:
            print("No records found.")
            continue

        outliers = sorted([z for z, c in counts.items() if c > args.threshold])
        if not outliers:
            #print(f"No parents above threshold {args.threshold}. Max count was {max(counts.values())}.")
            continue

        any_outliers = True
        print(f"Found {len(outliers)} outlier parent(s) with N(BIB) > {args.threshold}:")
        for z in outliers:
            print(f"  z_mu={z} : N={counts[z]}")

        # Pass 2: collect per-outlier particle data
        data = {z: {"p": [], "pdg": [], "r": [], "z": []} for z in outliers}
        outlier_set = set(outliers)

        n_read_2 = 0
        for chunk in iter_records(fin, chunk_records=args.chunk_records):
            if args.max_lines is not None:
                remaining = args.max_lines - n_read_2
                if remaining <= 0:
                    break
                if chunk.size > remaining:
                    chunk = chunk[:remaining]

            z_mu = chunk["z_mu"].astype(np.float64)
            if args.z_round is not None:
                z_mu = np.round(z_mu, args.z_round)

            mask_any = np.isin(z_mu, outliers)
            if not np.any(mask_any):
                n_read_2 += chunk.size
                continue

            sub = chunk[mask_any]
            sub_zmu = z_mu[mask_any]

            # FLUKA fid -> PDG
            fids = sub["fid"].astype(np.int32)
            pdg = np.zeros_like(fids, dtype=np.int32)
            for i, fid in enumerate(fids):
                pdg[i] = FLUKA_PIDS.get(int(fid), 0)

            # momentum magnitude from Ekin and mass
            e_kin = sub["E"].astype(np.float64)
            masses = np.zeros_like(e_kin)
            known = np.ones_like(e_kin, dtype=bool)
            for i, pid in enumerate(pdg):
                props = PDG_PROPS.get(int(pid))
                if props is None:
                    known[i] = False
                else:
                    masses[i] = float(props[1])

            if not np.any(known):
                n_read_2 += chunk.size
                continue

            sub = sub[known]
            sub_zmu = sub_zmu[known]
            pdg = pdg[known]
            e_kin = e_kin[known]
            masses = masses[known]

            p_mag = np.sqrt(e_kin * e_kin + 2.0 * e_kin * masses)

            x = sub["x"].astype(np.float64)
            y = sub["y"].astype(np.float64)
            zpos = sub["z"].astype(np.float64)
            r = np.sqrt(x * x + y * y)

            for z_parent in np.unique(sub_zmu):
                if z_parent not in outlier_set:
                    continue
                m = (sub_zmu == z_parent)
                data[z_parent]["p"].append(p_mag[m])
                data[z_parent]["pdg"].append(pdg[m])
                data[z_parent]["r"].append(r[m])
                data[z_parent]["z"].append(zpos[m])

            n_read_2 += chunk.size

        # Report + plots per outlier parent
        for z_parent in outliers:
            p = np.concatenate(data[z_parent]["p"]) if data[z_parent]["p"] else np.array([], dtype=np.float64)
            pdg = np.concatenate(data[z_parent]["pdg"]) if data[z_parent]["pdg"] else np.array([], dtype=np.int32)
            r = np.concatenate(data[z_parent]["r"]) if data[z_parent]["r"] else np.array([], dtype=np.float64)
            zpos = np.concatenate(data[z_parent]["z"]) if data[z_parent]["z"] else np.array([], dtype=np.float64)

            n_out = int(p.size)
            sum_p = float(np.sum(p)) if n_out else 0.0

            print("\n--- Outlier parent ---")
            print(f"File: {fin_path}")
            print(f"Parent z_mu (PosZmu): {z_parent}")
            print(f"Number of outgoing BIB particles: {n_out}")
            print(f"Sum of |p| over outgoing BIB particles: {sum_p:.6g} GeV/c")

            tag = safe_tag(float(z_parent))
            base = outdir / f"{fin_path.stem}_zmu_{tag}"

            plot_momentum(
                p, Path(str(base) + "_p_mag.png"),
                title=f"|p| distribution (parent z_mu={z_parent} cm; N={n_out})",
                logy=args.log
            )

            plot_pdg_bar(
                pdg, Path(str(base) + "_pdg.png"),
                title=f"PDGID distribution (parent z_mu={z_parent} cm; N={n_out})",
                topk=25,
                logy=args.log
            )

            plot_r_vs_z(
                r_vals=r, z_vals=zpos,
                outpath=Path(str(base) + "_r_vs_z.png"),
                title=f"r vs z (parent z_mu={z_parent} cm; N={n_out})",
                logz=args.log
            )

            print("Plots:")
            print(f"  {base}_p_mag.png")
            print(f"  {base}_pdg.png")
            print(f"  {base}_r_vs_z.png")

    if not any_outliers:
        print("\nNo outliers found in any input files.")


if __name__ == "__main__":
    main()
