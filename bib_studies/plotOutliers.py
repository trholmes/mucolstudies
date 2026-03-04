#!/usr/bin/env python3
"""
Find parent-muon (PosZmu / z_mu) outliers with unusually many BIB particles,
then print summary and make diagnostic plots per outlier parent.

Parent definition: unique z_mu values *within each file* (i.e. (file, z_mu)).

Cuts supported (applied consistently in both counting pass and plotting pass):
  - radial cut: keep particles with r = sqrt(x^2+y^2) >= r-min
  - momentum cut: keep particles with |p| in [p-min, p-max] (either can be omitted)

Plots per outlier parent:
  - |p| distribution (of selected particles)
  - PDGID distribution (topK + Other) (of selected particles)
  - 2D histogram of r vs z (of selected particles)

Options:
  --log     : log-y for 1D plots and log color scale (log-z) for 2D plots
  --r-min R : only count/plot particles with r >= R [cm]
  --p-min P : only count/plot particles with |p| >= P [GeV/c]
  --p-max P : only count/plot particles with |p| <= P [GeV/c]
  --z-round : round z_mu to N decimals before grouping

Assumptions:
  - Record layout matches LINE_DT
  - E is kinetic energy [GeV]
  - |p| = sqrt(T^2 + 2 T m) with m from PDG_PROPS, c=1

Usage:
  python plotOutliers.py *.bin --threshold 20000 --outdir out --log --r-min 200 --p-min 1.0
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
        h = ax.hist2d(
            z_vals, r_vals,
            bins=[bins_z, bins_r],
            norm=LogNorm(vmin=1)
        )
        fig.colorbar(h[3], ax=ax, label="BIB particles / bin (log scale)")
    else:
        h = ax.hist2d(z_vals, r_vals, bins=[bins_z, bins_r])
        fig.colorbar(h[3], ax=ax, label="BIB particles / bin")

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

    # Geometry cut
    ap.add_argument("--r-min", type=float, default=0.0,
                    help="Only count/plot particles with r = sqrt(x^2+y^2) >= r-min [cm] (default: 0)")

    # Momentum cuts (incoming BIB particle momentum magnitude)
    ap.add_argument("--p-min", type=float, default=None,
                    help="Only count/plot particles with |p| >= p-min [GeV/c]")
    ap.add_argument("--p-max", type=float, default=None,
                    help="Only count/plot particles with |p| <= p-max [GeV/c]")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    any_outliers = False

    # Helper to apply all cuts to a chunk, returning indices mask and computed |p|
    def compute_masks_and_p(chunk, z_mu_arr):
        # r cut
        x = chunk["x"].astype(np.float64)
        y = chunk["y"].astype(np.float64)
        r = np.sqrt(x * x + y * y)
        mask = (r >= args.r_min)

        # Momentum: need PDG mass, so map fid -> pdg and get mass
        fids = chunk["fid"].astype(np.int32)
        pdg = np.zeros_like(fids, dtype=np.int32)
        for i, fid in enumerate(fids):
            pdg[i] = FLUKA_PIDS.get(int(fid), 0)

        e_kin = chunk["E"].astype(np.float64)
        masses = np.zeros_like(e_kin)
        known = np.ones_like(e_kin, dtype=bool)
        for i, pid in enumerate(pdg):
            props = PDG_PROPS.get(int(pid))
            if props is None:
                known[i] = False
            else:
                masses[i] = float(props[1])

        # require known mass if any momentum cuts are specified (or if we need p for plots)
        # We always compute p for selected subsets anyway, so enforce known here.
        mask &= known

        if not np.any(mask):
            return mask, None, None, None, None, None

        # compute p only for masked entries
        p_mag = np.zeros_like(e_kin)
        idx = np.where(mask)[0]
        p_mag[idx] = np.sqrt(e_kin[idx] * e_kin[idx] + 2.0 * e_kin[idx] * masses[idx])

        # p cuts
        if args.p_min is not None:
            mask &= (p_mag >= args.p_min)
        if args.p_max is not None:
            mask &= (p_mag <= args.p_max)

        if not np.any(mask):
            return mask, None, None, None, None, None

        # Recompute p_mag on final mask (keeps alignment and avoids stale values)
        idx = np.where(mask)[0]
        p_mag = np.sqrt(e_kin[idx] * e_kin[idx] + 2.0 * e_kin[idx] * masses[idx])

        # return arrays aligned to final mask
        pdg_sel = pdg[idx]
        r_sel = r[idx]
        zpos_sel = chunk["z"][idx].astype(np.float64)
        zmu_sel = z_mu_arr[idx]

        return mask, p_mag, pdg_sel, r_sel, zpos_sel, zmu_sel

    for fin in args.files_in:
        fin_path = Path(fin)
        #print(f"\n=== Scanning file: {fin_path} ===")

        # -------- Pass 1: count BIB per parent z_mu (THIS FILE), applying cuts ----------
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

            # Apply r and p cuts
            mask, p_mag, pdg_sel, r_sel, zpos_sel, zmu_sel = compute_masks_and_p(chunk, z_mu)
            if zmu_sel is None or zmu_sel.size == 0:
                n_read += chunk.size
                continue

            uz, c = np.unique(zmu_sel, return_counts=True)
            for z, cc in zip(uz, c):
                counts[z] = counts.get(z, 0) + int(cc)

            n_read += chunk.size

        if not counts:
            #print("No records found (or all records failed cuts).")
            continue

        outliers = sorted([z for z, c in counts.items() if c > args.threshold])
        if not outliers:
            #print(f"No parents above threshold {args.threshold}. Max count was {max(counts.values())}.")
            continue

        any_outliers = True
        cut_desc = f"r>={args.r_min} cm"
        if args.p_min is not None:
            cut_desc += f", |p|>={args.p_min} GeV/c"
        if args.p_max is not None:
            cut_desc += f", |p|<={args.p_max} GeV/c"

        print(f"Found {len(outliers)} outlier parent(s) with N(BIB) > {args.threshold} (after cuts: {cut_desc}):")
        for z in outliers:
            print(f"  z_mu={z} : N={counts[z]}")

        # -------- Pass 2: collect per-outlier particle data (apply same cuts) ----------
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

            # quick z_mu prefilter
            mask_any = np.isin(z_mu, outliers)
            if not np.any(mask_any):
                n_read_2 += chunk.size
                continue

            sub = chunk[mask_any]
            sub_zmu = z_mu[mask_any]

            # Apply r and p cuts on subset
            mask, p_mag, pdg_sel, r_sel, zpos_sel, zmu_sel = compute_masks_and_p(sub, sub_zmu)
            if zmu_sel is None or zmu_sel.size == 0:
                n_read_2 += chunk.size
                continue

            # Split by parent z_mu
            for z_parent in np.unique(zmu_sel):
                if z_parent not in outlier_set:
                    continue
                m = (zmu_sel == z_parent)
                data[z_parent]["p"].append(p_mag[m])
                data[z_parent]["pdg"].append(pdg_sel[m])
                data[z_parent]["r"].append(r_sel[m])
                data[z_parent]["z"].append(zpos_sel[m])

            n_read_2 += chunk.size

        # -------- Report + plots per outlier parent ----------
        for z_parent in outliers:
            p = np.concatenate(data[z_parent]["p"]) if data[z_parent]["p"] else np.array([], dtype=np.float64)
            pdg = np.concatenate(data[z_parent]["pdg"]) if data[z_parent]["pdg"] else np.array([], dtype=np.int32)
            rvals = np.concatenate(data[z_parent]["r"]) if data[z_parent]["r"] else np.array([], dtype=np.float64)
            zpos = np.concatenate(data[z_parent]["z"]) if data[z_parent]["z"] else np.array([], dtype=np.float64)

            n_out = int(p.size)
            sum_p = float(np.sum(p)) if n_out else 0.0

            print("\n--- Outlier parent ---")
            print(f"File: {fin_path}")
            print(f"Parent z_mu (PosZmu): {z_parent}")
            print(f"Number of outgoing BIB particles (after cuts): {n_out}")
            print(f"Sum of |p| over outgoing BIB particles: {sum_p:.6g} GeV/c")

            tag = safe_tag(float(z_parent))
            cut_tag = f"rmin_{args.r_min:g}"
            if args.p_min is not None:
                cut_tag += f"_pmin_{args.p_min:g}"
            if args.p_max is not None:
                cut_tag += f"_pmax_{args.p_max:g}"

            base = outdir / f"{fin_path.stem}_zmu_{tag}_{cut_tag}"

            plot_momentum(
                p, Path(str(base) + "_p_mag.png"),
                title=f"|p| distribution (z_mu={z_parent} cm; N={n_out}; {cut_desc})",
                logy=args.log
            )

            plot_pdg_bar(
                pdg, Path(str(base) + "_pdg.png"),
                title=f"PDGID distribution (z_mu={z_parent} cm; N={n_out}; {cut_desc})",
                topk=25,
                logy=args.log
            )

            plot_r_vs_z(
                r_vals=rvals, z_vals=zpos,
                outpath=Path(str(base) + "_r_vs_z.png"),
                title=f"r vs z (z_mu={z_parent} cm; N={n_out}; {cut_desc})",
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
