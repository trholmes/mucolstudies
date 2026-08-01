#!/usr/bin/env python3
"""Shower-propagation-corrected time variable dt' for calorimeter timing cuts.

delta_t = t - r/c treats every deposit as if it came promptly from the IP at
the speed of light. Hadronic showers instead develop late energy as they
penetrate: the energy-weighted late quantiles of delta_t grow ~linearly with
the 3D distance r of the cell from the IP. This script:

  1. makes the 2D (energy-weighted) distribution of contribution delta_t vs r
     for pions, with energy-weighted quantile profiles overlaid, to show how
     linear the growth is;
  2. fits the q75/q90 profiles per region and globally ->
        dt' = delta_t - k * max(0, r - R0)
     (equivalently: allow propagation at an effective speed
      v_eff = 1/(1/c + k) beyond the calorimeter front at R0);
  3. quantifies at the CELL level (earliest accepted time, default digi
     window, all thresholds applied) how much more pion energy a selector
     window on dt' retains compared to the same window on delta_t.

Usage:
  python3 shower_speed.py --flats data/flat/pionGun_*.npz \
      --thresholds data/thresholds_ecal.npz --outdir plots/hcal_pions
"""

import argparse
import glob
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import cutmodel as cm  # noqa: E402
from emulate import RegionData, windowed_vectorized  # noqa: E402

R0 = 1857.0  # ECAL barrel front face [mm] (min cell r in the samples)


def wquant(x, w, qs):
    o = np.argsort(x)
    x, w = x[o], w[o]
    cw = np.cumsum(w)
    cw /= cw[-1]
    return np.interp(qs, cw, x)


def gather_contribs(flats):
    DT, E, R = [], [], []
    for f in flats:
        d = dict(np.load(f))
        ncon = d["cell_n_contrib"]
        alloff = np.concatenate([[0], np.cumsum(ncon)]).astype(np.int64)
        for reg in (0, 2):  # EcalBarrel, HcalBarrel
            sel = d["cell_region"] == reg
            if not sel.any():
                continue
            idx = np.concatenate([np.arange(alloff[i], alloff[i + 1])
                                  for i in np.flatnonzero(sel)])
            DT.append(d["contrib_dt"][idx])
            E.append(d["contrib_e"][idx])
            R.append(np.repeat(d["cell_r"][sel], ncon[sel]))
    return np.concatenate(DT), np.concatenate(E), np.concatenate(R)


def profile(dt, e, r, edges, qs=(0.5, 0.75, 0.9)):
    cen, out = [], []
    for i in range(len(edges) - 1):
        s = (r >= edges[i]) & (r < edges[i + 1]) & (dt < 25)
        if e[s].sum() < 1.0:  # require >= 1 GeV in the slice
            continue
        cen.append(0.5 * (edges[i] + edges[i + 1]))
        out.append(wquant(dt[s], e[s], qs))
    return np.array(cen), np.array(out)


def wfit(x, y, w):
    """Weighted linear fit y = a + b x."""
    W = np.sum(w)
    xm, ym = np.sum(w * x) / W, np.sum(w * y) / W
    b = np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2)
    return ym - b * xm, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flats", nargs="+", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    dt, e, r = gather_contribs(args.flats)

    # ---- 2D map + quantile profiles ------------------------------------
    redges = np.arange(1850, 4300, 50)
    tedges = np.concatenate([np.arange(-0.5, 1.0, 0.05),
                             np.geomspace(1.0, 25.0, 40)])
    H, _, _ = np.histogram2d(r, np.clip(dt, -0.49, 24.9),
                             bins=[redges, tedges], weights=e)
    cen, prof = profile(dt, e, r, redges)

    # weights for the fits: energy per slice
    slice_e = np.array([e[(r >= redges[i]) & (r < redges[i + 1]) & (dt < 25)].sum()
                        for i in range(len(redges) - 1)])
    slice_e = np.array([se for i, se in enumerate(slice_e)
                        if se >= 1.0])

    fits = {}
    for iq, qname in ((1, "q75"), (2, "q90")):
        a, b = wfit(cen, prof[:, iq], slice_e)
        fits[qname] = (a, b)
        v_eff = 1.0 / (1.0 / cm.C_MM_NS + b)
        print(f"global {qname} fit: dt = {a:+.3f} + {b*1e3:.3f} ns/m * r; "
              f"slope-> v_eff = {v_eff:.0f} mm/ns = {v_eff/cm.C_MM_NS:.2f} c")

    fig, ax = plt.subplots(figsize=(9, 6))
    X, Y = np.meshgrid(redges, tedges)
    pc = ax.pcolormesh(X, Y, H.T, norm=LogNorm(vmin=max(H[H > 0].min(), 1e-4),
                                               vmax=H.max()), cmap="viridis")
    fig.colorbar(pc, ax=ax, label="energy / bin [GeV]")
    for iq, (qname, ls) in enumerate((("q50", ":"), ("q75", "--"), ("q90", "-"))):
        ax.plot(cen, prof[:, iq], "w" + ls, lw=2, label=f"energy-weighted {qname}")
    for qname, color in (("q75", "orange"), ("q90", "red")):
        a, b = fits[qname]
        rr = np.array([R0, redges[-1]])
        k = b  # ns/mm
        ax.plot(rr, (a + b * R0) + k * (rr - R0), color=color, lw=2,
                label=f"linear fit to {qname}: k = {b*1e3:.2f} ns/m")
    ax.set_xlabel("cell distance from IP  r [mm]")
    ax.set_ylabel(r"$\Delta t = t - r/c$ [ns]")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_ylim(-0.5, 25)
    ax.set_title("pions (10+50+200 GeV), ECAL+HCAL barrel: contribution "
                 r"$\Delta t$ vs $r$" + "\n(energy weighted; ECAL/HCAL boundary ~2.1 m)",
                 fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    out = f"{args.outdir}/dt_vs_r_pion.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)

    # ---- retention: late edge on dt vs on dt' (cell level) -------------
    # The correction must only relax the LATE edge: prompt cores at depth get
    # dt' < 0, so a symmetric dt' window would cut them. Cut definition:
    #   standard:  -0.3 < t_cell < +w
    #   dt'     :  -0.3 < t_cell  and  t_cell - k*(r-R0) < +w
    ecal_map = cm.SelectorThresholdMap(args.thresholds, "ECAL")
    print("\ncell-level retention (default digi window + all thresholds), per-event median,")
    print("ECAL+HCAL barrel summed; early edge fixed at t > -0.3 ns in both variants:")
    print(f"{'sample':>10s} {'late edge':>10s} {'on dt':>8s} {'on dtp(q75)':>12s} {'on dtp(q90)':>12s}")
    for f in args.flats:
        flat = dict(np.load(f))
        nev = len(flat["ev_gun_e"])
        acc = {}
        for ireg, region in ((0, "EcalBarrel"), (2, "HcalBarrel")):
            rd = RegionData(flat, ireg)
            if not len(rd.event):
                continue
            e_win, t_cell, _ = windowed_vectorized(
                rd, cm.DIGI_WINDOW_MIN, cm.DIGI_WINDOW_MAX)
            e_mip = cm.digi_cell_response_mip(e_win, region)
            ok = cm.digi_threshold_pass(e_mip, region)
            e_rec = cm.rec_energy_gev(e_mip, region)
            base = ok & cm.selector_pass(
                e_rec, np.zeros_like(t_cell), rd.theta,
                rd.layer, region, -1, 1, ecal_map)
            cellr = flat["cell_r"][flat["cell_region"] == ireg]
            acc[region] = (rd.event, e_rec, t_cell, cellr, base)
        for w in (0.3, 0.5, 1.0):
            sums = {}
            for tag, kk in (("dt", 0.0), ("q75", fits["q75"][1]), ("q90", fits["q90"][1])):
                tot = np.zeros(nev)
                ref = np.zeros(nev)
                for region, (ev, e_rec, t_cell, cellr, base) in acc.items():
                    tp = t_cell - kk * np.maximum(0.0, cellr - R0)
                    keep = base & (t_cell > -0.3) & (tp < w)
                    tot += np.bincount(ev[keep], weights=e_rec[keep], minlength=nev)
                    ref += np.bincount(ev[base], weights=e_rec[base], minlength=nev)
                okev = ref > 0
                sums[tag] = np.median(tot[okev] / ref[okev])
            label = f.rsplit("/", 1)[-1].split("_")[1]
            print(f"{label:>10s} {'+' + str(w) + ' ns':>10s} {sums['dt']:8.3f} "
                  f"{sums['q75']:12.3f} {sums['q90']:12.3f}")


if __name__ == "__main__":
    main()
