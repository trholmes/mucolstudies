#!/usr/bin/env python3
"""
Make a histogram of number of BIB particles entering the detector region
as a function of parent muon z position (PosZmu / z_mu) from FLUKA binary files.

This matches the binary record structure used in fluka_remix.py:
  ... ('z_mu', np.float64)

Usage examples:
    python studyParent.py /data/fmeloni/FLUKA/summary1_DET_IP.dat -o plots/poszmu.png
python bib_poszmu_hist.py input1.bin input2.bin -o poszmu.png
  python bib_poszmu_hist.py input.bin --bin-width 1.0 --range -600 0
  python bib_poszmu_hist.py input.bin --logy
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt


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
    """Yield numpy structured arrays of records from a FLUKA binary file."""
    with open(filename, "rb") as f:
        while True:
            arr = np.fromfile(f, dtype=LINE_DT, count=chunk_records)
            if arr.size == 0:
                break
            yield arr


def main():
    ap = argparse.ArgumentParser(description="Histogram BIB particles vs PosZmu (z_mu).")
    ap.add_argument("files_in", nargs="+", help="Input FLUKA binary file(s)")
    ap.add_argument("-o", "--out", default="poszmu_hist.png", help="Output plot filename (png/pdf/etc.)")

    # Binning controls
    ap.add_argument("--bins", type=int, default=200, help="Number of bins (used if --bin-width is not set)")
    ap.add_argument("--bin-width", type=float, default=None, help="Bin width in cm (overrides --bins)")
    ap.add_argument("--range", nargs=2, type=float, default=None, metavar=("ZMIN", "ZMAX"),
                    help="Histogram x-range in cm (e.g. --range -600 0). Default: autoscale from data.")

    # Optional filters similar in spirit to your example
    ap.add_argument("--t-max", type=float, default=None,
                    help="Max particle time [ns]. Uses t = time*1e9 - bx_time.")
    ap.add_argument("--bx-time", type=float, default=0.0, help="Bunch crossing time [ns] (subtracted)")
    ap.add_argument("--ne-min", type=float, default=None, help="Min kinetic energy for neutrons [GeV] (PDG=2112)")
    ap.add_argument("--pdg", type=int, nargs="+", default=None,
                    help="Keep only these PDG IDs (requires --use-fluka-pid-map)")
    ap.add_argument("--nopdg", type=int, nargs="+", default=None,
                    help="Drop these PDG IDs (requires --use-fluka-pid-map)")
    ap.add_argument("--use-fluka-pid-map", action="store_true",
                    help="Enable FLUKA->PDG mapping via bib_pdgs.py (like fluka_remix.py).")
    ap.add_argument("--max-lines", type=int, default=None, help="Max number of records total (across all files)")
    ap.add_argument("--logy", action="store_true", help="Log-scale y-axis")
    args = ap.parse_args()

    fluka_to_pdg = None
    if args.use_fluka_pid_map:
        # Same mapping module your example uses
        from bib_pdgs import FLUKA_PIDS
        fluka_to_pdg = FLUKA_PIDS

    z_all = []
    n_read = 0

    for fin in args.files_in:
        for chunk in iter_records(fin):
            if args.max_lines is not None:
                remaining = args.max_lines - n_read
                if remaining <= 0:
                    break
                if chunk.size > remaining:
                    chunk = chunk[:remaining]

            # Basic masks (start with all True)
            mask = np.ones(chunk.size, dtype=bool)

            # Time cut (t = time*1e9 - bx_time)
            if args.t_max is not None:
                t_ns = chunk["time"] * 1e9 - args.bx_time
                mask &= (t_ns <= args.t_max)

            # Optional PDG-based cuts need FLUKA->PDG mapping
            if (args.ne_min is not None) or (args.pdg is not None) or (args.nopdg is not None):
                if fluka_to_pdg is None:
                    raise RuntimeError("PDG/neutron filters require --use-fluka-pid-map (needs bib_pdgs.py).")

                # Map FLUKA fid -> PDG (unknown fids become 0 and can be filtered out if desired)
                fids = chunk["fid"]
                pdg = np.zeros_like(fids, dtype=np.int32)
                # vectorized-ish mapping
                for i, fid in enumerate(fids):
                    pdg[i] = fluka_to_pdg.get(int(fid), 0)

                if args.pdg is not None:
                    keep = np.isin(pdg, np.array(args.pdg, dtype=np.int32))
                    mask &= keep
                if args.nopdg is not None:
                    drop = np.isin(pdg, np.array(args.nopdg, dtype=np.int32))
                    mask &= ~drop

                if args.ne_min is not None:
                    # Apply only to neutrons (PDG 2112); others pass
                    is_neutron = (np.abs(pdg) == 2112)
                    mask &= (~is_neutron) | (chunk["E"] >= args.ne_min)

            z_all.append(chunk["z_mu"][mask])
            n_read += chunk.size

        if args.max_lines is not None and n_read >= args.max_lines:
            break

    if not z_all:
        raise RuntimeError("No data read (empty input or all events filtered).")

    z = np.concatenate(z_all)
    if z.size == 0:
        raise RuntimeError("All entries were filtered out; nothing to histogram.")

    # Choose bins
    if args.range is not None:
        zmin, zmax = args.range
    else:
        zmin, zmax = float(np.min(z)), float(np.max(z))

    if args.bin_width is not None:
        nbins = int(np.ceil((zmax - zmin) / args.bin_width))
        nbins = max(nbins, 1)
        bins = np.linspace(zmin, zmax, nbins + 1)
    else:
        bins = args.bins

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(z, bins=bins, range=(zmin, zmax) if isinstance(bins, int) else None, histtype="step", linewidth=1.5)

    ax.set_xlabel("PosZmu (z_mu) [cm]")
    ax.set_ylabel("BIB particles entering detector volume (count)")
    ax.set_title("BIB particle count vs parent muon z position")

    if args.logy:
        ax.set_yscale("log")

    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"Wrote: {args.out}")
    print(f"Records read: {n_read}, histogrammed entries: {z.size}, z range: [{zmin:.3g}, {zmax:.3g}] cm")


if __name__ == "__main__":
    main()
