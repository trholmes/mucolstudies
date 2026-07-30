#!/usr/bin/env python3
"""Emulate the digi + CaloHitSelector chain on flattened simhits, scanning the
timing windows. Pure numpy; runs anywhere.

For every (digi window, selector window) point of the scan grid it computes,
per event and per calorimeter region, the retained selected energy (calibrated
GeV), the number of surviving cells, and (at default digi settings) the
energy-weighted delta_t distributions.

Key structural fact (see PLAN.md): with the digi window MINIMUM fixed, the
stored cell time (earliest accepted contribution) does not depend on the digi
window MAXIMUM; the max edge only changes the windowed energy and therefore
which cells survive the thresholds.

Usage:
  python3 emulate.py flat1.npz [flat2.npz ...] --thresholds data/thresholds_ecal.npz \
      --out data/results/scan_photonGun.npz
"""

import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import cutmodel as cm  # noqa: E402

# ---- scan grid (defaults marked in PLAN.md) --------------------------------
DIGI_MAX_SCAN = [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 25.0, 50.0, 100.0]
DIGI_MIN_SCAN = [-0.25, -0.5, -1.0]
SEL_HALFWIDTH_SCAN = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0,
                      1.5, 2.0, 5.0, np.inf]
SEL_MAX_SCAN = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, np.inf]
DT_BINS = np.concatenate([np.arange(-1.0, 3.0, 0.02), np.arange(3.0, 25.0, 0.25)])


class RegionData:
    """Per-region cell/contribution arrays with contributions sorted by dt
    within each cell, plus cumulative energy for fast window sums."""

    def __init__(self, flat, ireg):
        sel = flat["cell_region"] == ireg
        self.region = cm.REGIONS[ireg]
        self.event = flat["cell_event"][sel]
        self.cellid = flat["cell_cellid"][sel]
        self.layer = flat["cell_layer"][sel]
        self.theta = flat["cell_theta"][sel]
        self.nev = int(flat["ev_gun_e"].shape[0])

        ncon = flat["cell_n_contrib"][sel]
        alloff = np.concatenate([[0], np.cumsum(flat["cell_n_contrib"])]).astype(np.int64)
        # gather this region's contributions, sorted by dt within each cell
        idx = np.concatenate([np.arange(alloff[i], alloff[i + 1])
                              for i in np.flatnonzero(sel)]) if sel.any() else \
            np.array([], dtype=np.int64)
        dt = flat["contrib_dt"][idx]
        e = flat["contrib_e"][idx]
        self.off = np.concatenate([[0], np.cumsum(ncon)]).astype(np.int64)
        cellindex = np.repeat(np.arange(len(ncon)), ncon)
        order = np.lexsort((dt, cellindex))
        self.dt = dt[order]
        self.e = e[order]
        self.cum_e = np.cumsum(self.e)

def windowed_vectorized(rd, wmin, wmax):
    """Per cell: energy in (wmin, wmax) and earliest contribution dt > wmin.
    Within-cell searchsorted done via one global searchsorted on (cell, dt)
    keys: each cell gets a slot of width SPAN, dt values (clipped to +-CLIP)
    are offset into their slot, so per-cell window edges become global keys."""
    n = len(rd.off) - 1
    if n == 0:
        return (np.zeros(0), np.zeros(0), np.zeros(0, dtype=bool))
    SPAN, CLIP = 1e6, 1e5  # |dt| < 1e5 ns always; window edges clipped to 2*CLIP
    cellindex = np.repeat(np.arange(n), np.diff(rd.off))
    key = cellindex.astype(np.float64) * SPAN + np.clip(rd.dt, -CLIP, CLIP)
    slots = np.arange(n) * SPAN
    lo = np.searchsorted(key, slots + np.clip(wmin, -2 * CLIP, 2 * CLIP), side="right")
    hi = np.searchsorted(key, slots + np.clip(wmax, -2 * CLIP, 2 * CLIP), side="left")
    csum = np.concatenate([[0.0], rd.cum_e])
    e_win = csum[hi] - csum[lo]
    has = hi > lo
    # first dt > wmin (index clamped for empty cells; masked by `has`)
    t_cell = np.where(has, rd.dt[np.minimum(lo, max(len(rd.dt) - 1, 0))], np.nan)
    return e_win, t_cell, has


def run_scan(flat, ecal_map, digi_min_scan, digi_max_scan):
    """Full scan for one flat file. Returns dict of result arrays."""
    nev = int(flat["ev_gun_e"].shape[0])
    res = {"gun_e": flat["ev_gun_e"], "gun_theta": flat["ev_gun_theta"]}
    grids = {}

    for ireg in range(len(cm.REGIONS)):
        region = cm.REGIONS[ireg]
        rd = RegionData(flat, ireg)
        if len(rd.event) == 0:
            continue

        for dmin in digi_min_scan:
            e_win_ref, t_cell, _ = windowed_vectorized(rd, dmin, np.inf)
            for dmax in digi_max_scan:
                e_win, _, _ = windowed_vectorized(rd, dmin, dmax)
                # digi: mean response + threshold (on windowed energy)
                e_mip = cm.digi_cell_response_mip(e_win, region)
                pass_digi = cm.digi_threshold_pass(e_mip, region)
                # cell time = earliest contribution > dmin, must itself be in window
                # (a cell whose first accepted dt >= dmax has no in-window
                # contribution -> e_win = 0 -> fails threshold anyway)
                e_rec = cm.rec_energy_gev(e_mip, region)

                for tag, smin, smax in sel_windows():
                    ok = pass_digi & cm.selector_pass(
                        e_rec, t_cell, rd.theta, rd.layer, region,
                        smin, smax, ecal_map)
                    esel = np.bincount(rd.event[ok],
                                       weights=e_rec[ok], minlength=nev)
                    nsel = np.bincount(rd.event[ok], minlength=nev)
                    grids[(region, dmin, dmax, tag)] = (esel, nsel)

        # delta_t histograms at default digi min (energy weighted):
        # contribution level and cell-earliest level
        if abs(cm.DIGI_WINDOW_MIN - (-0.5)) < 1e-9:
            hw, _ = np.histogram(rd.dt, bins=DT_BINS, weights=rd.e)
            res[f"dthist_contrib_{region}"] = hw
            e_win, t_cell, has = windowed_vectorized(
                rd, cm.DIGI_WINDOW_MIN, cm.DIGI_WINDOW_MAX)
            hw2, _ = np.histogram(t_cell[has], bins=DT_BINS,
                                  weights=e_win[has])
            res[f"dthist_cellearliest_{region}"] = hw2

    return res, grids


def sel_windows():
    """(tag, min, max) triplets for the selector scans."""
    out = []
    for w in SEL_HALFWIDTH_SCAN:
        out.append((f"sym{w}", -w, w))
    for m in SEL_MAX_SCAN:
        out.append((f"asym{m}", cm.SEL_WINDOW_MIN, m))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("flats", nargs="+")
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--digi-min", type=float, nargs="+", default=[-0.5])
    ap.add_argument("--digi-max", type=float, nargs="+", default=DIGI_MAX_SCAN)
    args = ap.parse_args()

    ecal_map = cm.SelectorThresholdMap(args.thresholds, "ECAL")

    out = {"dt_bins": DT_BINS,
           "digi_min_scan": np.array(args.digi_min),
           "digi_max_scan": np.array(args.digi_max),
           "sel_tags": np.array([t for t, _, _ in sel_windows()])}
    for fname in args.flats:
        flat = dict(np.load(fname))
        stem = fname.rsplit("/", 1)[-1].replace(".npz", "")
        res, grids = run_scan(flat, ecal_map, args.digi_min, args.digi_max)
        for k, v in res.items():
            out[f"{stem}__{k}"] = v
        for (region, dmin, dmax, tag), (esel, nsel) in grids.items():
            key = f"{stem}__sel_{region}__dmin{dmin}__dmax{dmax}__{tag}"
            out[key + "__e"] = esel
            out[key + "__n"] = nsel
        print(f"scanned {fname}")

    np.savez_compressed(args.out, **out)
    print(f"wrote {args.out} ({len(out)} arrays)")


if __name__ == "__main__":
    main()
