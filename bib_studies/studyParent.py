#!/usr/bin/env python3
"""
Muon-collider BIB plots from FLUKA binary (line_dt with z_mu):

(1) # BIB particles vs parent muon z (PosZmu / z_mu)        [particle-weighted]
(2) # unique parent muons vs parent muon z                  [parent-weighted]
(3) average BIB particles per parent vs parent muon z:
        <N_BIB>(z) = N_BIB_particles_in_bin(z) / N_parent_muons_in_bin(z)
(4) distribution of N_BIB per parent muon                   [group sizes]

Unique parent definition matches fluka_remix.py:
  unique parents are defined by unique values of PosZmu (z_mu).
Use --z-round to avoid float-equality issues if z_mu should be discrete.
"""


# Example usage
# python studyParent.py /data/fmeloni/FLUKA/summary*_DET_IP.dat --max-lines 10000


import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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


def hist_step(data, bins, xlabel, ylabel, title, out, logy=False, x_range=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(data, bins=bins, range=x_range, histtype="step", linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def line_plot(x_centers, y, xlabel, ylabel, title, out, logy=False):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x_centers, y, marker="o", linestyle="-", linewidth=1.2, markersize=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files_in", nargs="+", help="Input FLUKA binary file(s)")
    ap.add_argument("--out-prefix", default="poszmu", help="Prefix for output filenames")
    ap.add_argument("--bins", type=int, default=200, help="Number of z bins (if --bin-width not set)")
    ap.add_argument("--bin-width", type=float, default=None, help="z bin width in cm (overrides --bins)")
    ap.add_argument("--range", nargs=2, type=float, default=None, metavar=("ZMIN", "ZMAX"),
                    help="z range in cm for binnings (e.g. --range -600 0). Default: autoscale from data.")
    ap.add_argument("--z-round", type=int, default=None,
                    help="Round z_mu to this many decimal places before grouping/unique (recommended if z is discrete)")
    ap.add_argument("--max-lines", type=int, default=None, help="Max number of records total (across all files)")
    ap.add_argument("--logy", action="store_true", help="Log y-axis for the histogram plots")
    args = ap.parse_args()

    # Read z_mu for all particles (particle-weighted distribution)
    z_particles = []
    n_read = 0

    for fin in args.files_in:
        for chunk in iter_records(fin):
            if args.max_lines is not None:
                remaining = args.max_lines - n_read
                if remaining <= 0:
                    break
                if chunk.size > remaining:
                    chunk = chunk[:remaining]

            z = chunk["z_mu"].astype(np.float64)
            if args.z_round is not None:
                z = np.round(z, args.z_round)

            z_particles.append(z)
            n_read += chunk.size

        if args.max_lines is not None and n_read >= args.max_lines:
            break

    if not z_particles:
        raise RuntimeError("No data read (empty input?).")

    z_particles = np.concatenate(z_particles)
    if z_particles.size == 0:
        raise RuntimeError("No z_mu values found.")

    # Determine z range
    if args.range is not None:
        zmin, zmax = map(float, args.range)
    else:
        zmin = float(np.min(z_particles))
        zmax = float(np.max(z_particles))

    # Define common z bin edges for numerator/denominator ratio
    if args.bin_width is not None:
        nbins = int(np.ceil((zmax - zmin) / args.bin_width))
        nbins = max(nbins, 1)
        z_edges = np.linspace(zmin, zmax, nbins + 1)
    else:
        z_edges = np.linspace(zmin, zmax, args.bins + 1)

    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    out_prefix = Path(args.out_prefix)

    # (1) Particle-weighted: # BIB particles vs z_mu
    hist_step(
        z_particles,
        bins=z_edges,
        xlabel="PosZmu (z_mu) [cm]",
        ylabel="BIB particles entering detector volume (count)",
        title="BIB particle count vs parent muon z (particle-weighted)",
        out=str(out_prefix) + "_bib_particles_vs_poszmu.png",
        logy=args.logy
    )

    # Unique parents defined by unique z_mu values
    # Also compute particles-per-parent = counts per unique z
    unique_z, counts_per_parent = np.unique(z_particles, return_counts=True)

    # (2) Parent-weighted: # unique parent muons vs z_mu (one per unique z)
    hist_step(
        unique_z,
        bins=z_edges,
        xlabel="PosZmu (z_mu) [cm]",
        ylabel="Unique parent muons (count)",
        title="Unique parent muon count vs z (parent-weighted)",
        out=str(out_prefix) + "_parent_muons_vs_poszmu.png",
        logy=args.logy
    )

    # (3) Average BIB per parent vs z (ratio of histograms with same binning)
    # Numerator: particles per z-bin (counts of z_particles)
    n_bib_per_bin, _ = np.histogram(z_particles, bins=z_edges)
    # Denominator: parents per z-bin (counts of unique_z)
    n_parent_per_bin, _ = np.histogram(unique_z, bins=z_edges)

    avg_bib_per_parent = np.full_like(n_bib_per_bin, fill_value=np.nan, dtype=np.float64)
    mask = n_parent_per_bin > 0
    avg_bib_per_parent[mask] = n_bib_per_bin[mask] / n_parent_per_bin[mask]

    line_plot(
        z_centers,
        avg_bib_per_parent,
        xlabel="PosZmu (z_mu) [cm] (bin centers)",
        ylabel="Average BIB particles per parent muon in bin",
        title="Average BIB per parent vs parent muon z",
        out=str(out_prefix) + "_avg_bib_per_parent_vs_poszmu.png",
        logy=False
    )

    # (4) Distribution of # BIB particles per parent muon
    # (counts_per_parent is already N_BIB for each unique parent z)
    # Use a reasonable binning: integers from 0..max
    max_n = int(np.max(counts_per_parent))
    # Cap number of bins if huge
    if max_n <= 2000:
        bins_pp = np.arange(0, max_n + 2)  # integer bins
    else:
        bins_pp = 200  # fallback
    hist_step(
        counts_per_parent,
        bins=bins_pp,
        xlabel="BIB particles per parent muon (group size)",
        ylabel="Parent muons (count)",
        title="Distribution of BIB particles per parent muon",
        out=str(out_prefix) + "_bib_per_parent.png",
        logy=args.logy
    )

    # Quick summary to stdout
    print(f"Read records: {n_read}")
    print(f"Particles histogrammed: {z_particles.size}")
    print(f"Unique parents (unique z_mu): {unique_z.size}")
    print(f"z range used: [{zmin:.6g}, {zmax:.6g}] cm")
    print("Wrote:")
    print(f"  {out_prefix}_bib_particles_vs_poszmu.png")
    print(f"  {out_prefix}_parent_muons_vs_poszmu.png")
    print(f"  {out_prefix}_avg_bib_per_parent_vs_poszmu.png")
    print(f"  {out_prefix}_bib_per_parent.png")
    empty_bins = int(np.sum(~mask))
    if empty_bins:
        print(f"Note: {empty_bins} z-bins had zero parents; avg plot shows NaN there.")


if __name__ == "__main__":
    main()
