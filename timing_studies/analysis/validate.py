#!/usr/bin/env python3
"""Compare a real-chain variant run (summarize_reco.py output) against the
emulation at the matching grid point (gate 6), and report how the leading
Pandora cluster energy tracks the Sel-level energy.

Usage:
  python3 validate.py --flat flat.npz --summary summary_variant.npz \
      --thresholds data/thresholds_ecal.npz \
      --digi-min -0.5 --digi-max 10 --sel-min -0.15 --sel-max 0.15
"""

import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import cutmodel as cm  # noqa: E402
from emulate import RegionData, windowed_vectorized  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flat", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--digi-min", type=float, default=cm.DIGI_WINDOW_MIN)
    ap.add_argument("--digi-max", type=float, default=cm.DIGI_WINDOW_MAX)
    ap.add_argument("--sel-min", type=float, default=cm.SEL_WINDOW_MIN)
    ap.add_argument("--sel-max", type=float, default=cm.SEL_WINDOW_MAX)
    args = ap.parse_args()

    flat = dict(np.load(args.flat))
    summ = dict(np.load(args.summary))
    ecal_map = cm.SelectorThresholdMap(args.thresholds, "ECAL")
    nev = int(flat["ev_gun_e"].shape[0])
    assert len(summ["gun_e"]) == nev, "event counts differ - wrong file pair?"

    print(f"emulation grid point: digi ({args.digi_min}, {args.digi_max}) ns, "
          f"selector ({args.sel_min}, {args.sel_max}) ns")
    total_emu = np.zeros(nev)
    total_act = np.zeros(nev)
    for ireg, region in enumerate(cm.REGIONS):
        rd = RegionData(flat, ireg)
        if len(rd.event) == 0:
            continue
        e_win, t_cell, _ = windowed_vectorized(rd, args.digi_min, args.digi_max)
        e_mip = cm.digi_cell_response_mip(e_win, region)
        ok = cm.digi_threshold_pass(e_mip, region)
        e_rec = cm.rec_energy_gev(e_mip, region)
        ok &= cm.selector_pass(e_rec, t_cell, rd.theta, rd.layer, region,
                               args.sel_min, args.sel_max, ecal_map)
        emu = np.bincount(rd.event[ok], weights=e_rec[ok], minlength=nev)
        act = summ[f"sel_sum_{region}"]
        total_emu += emu
        total_act += act
        both = act > 1e-6
        if both.any():
            rel = (emu[both] - act[both]) / act[both]
            print(f"[{region}] emul vs actual Sel sum: median rel dev "
                  f"{np.median(rel):+.4f}, 68% |dev| < "
                  f"{np.percentile(np.abs(rel), 68):.4f}")

    both = total_act > 1e-6
    lead = summ["lead_cluster_e"]
    print(f"[all]   emul vs actual Sel sum: median rel dev "
          f"{np.median((total_emu[both] - total_act[both]) / total_act[both]):+.4f}")
    ok2 = both & (lead > 0)
    print(f"[all]   leading cluster / Sel sum: median "
          f"{np.median(lead[ok2] / total_act[ok2]):.4f} "
          f"(clustering leakage; 1.0 = cluster captures all selected energy)")
    print(f"[all]   median leading cluster E = {np.median(lead[lead > 0]):.3f} GeV, "
          f"median gun E = {np.median(summ['gun_e']):.3f} GeV")


if __name__ == "__main__":
    main()
