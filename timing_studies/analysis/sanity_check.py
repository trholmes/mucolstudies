#!/usr/bin/env python3
"""Closure gates: emulated digi/selector at DEFAULT settings vs the hits the
real chain actually produced (stored in the flat file). Run before trusting
any scan. Pure numpy.

Gates (see PLAN.md):
  2. prompt delta_t peaks near 0
  3. per-cell join emulated vs actual Digi hits (cells / energy / time)
  4. per-cell join emulated vs actual Sel hits
  5. event-level emulated vs actual Sel sums

Expected residuals: ECAL energies agree to <~1% per cell (Poisson e-h smear in
the real chain); HCAL cell energies fluctuate by the binomial SiPM smear
(~26%/sqrt(E_mip)); times must agree to float precision for matched cells.

Usage: python3 sanity_check.py flat.npz --thresholds data/thresholds_ecal.npz
"""

import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import cutmodel as cm  # noqa: E402
from emulate import RegionData, windowed_vectorized  # noqa: E402


def join(ev_a, id_a, ev_b, id_b):
    """Index pairs matching (event, cellid) between two hit lists."""
    dt = np.dtype([("ev", np.int64), ("id", np.uint64)])
    ka = np.empty(len(ev_a), dt); ka["ev"], ka["id"] = ev_a, id_a
    kb = np.empty(len(ev_b), dt); kb["ev"], kb["id"] = ev_b, id_b
    common, ia, ib = np.intersect1d(ka, kb, return_indices=True)
    return ia, ib, len(common)



# layer field of the MAIA calo cellID encoding (bits 19-27); survives the
# 32-bit truncation that RealisticCaloReco applies to Rec/Coned/Sel cellIDs
def coder_layer(flat, cellid):
    return ((cellid.astype(np.uint64) >> np.uint64(19))
            & np.uint64((1 << 9) - 1)).astype(np.int64)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("flat")
    ap.add_argument("--thresholds", required=True)
    args = ap.parse_args()

    flat = dict(np.load(args.flat))
    ecal_map = cm.SelectorThresholdMap(args.thresholds, "ECAL")
    nev = int(flat["ev_gun_e"].shape[0])
    print(f"{args.flat}: {nev} events")
    fails = 0

    for ireg, region in enumerate(cm.REGIONS):
        rd = RegionData(flat, ireg)
        if len(rd.event) == 0:
            continue

        # ---- gate 2: prompt peak
        if len(rd.dt):
            med = np.median(rd.dt[rd.e > 0])
            wsum = np.sum(rd.e[np.abs(rd.dt) < 0.5]) / max(np.sum(rd.e), 1e-12)
            print(f"[{region}] gate2 median contribution dt = {med:+.3f} ns; "
                  f"|dt|<0.5ns energy fraction = {wsum:.3f}")

        # ---- emulate defaults
        e_win, t_cell, _ = windowed_vectorized(rd, cm.DIGI_WINDOW_MIN, cm.DIGI_WINDOW_MAX)
        e_mip = cm.digi_cell_response_mip(e_win, region)
        pass_digi = cm.digi_threshold_pass(e_mip, region)
        e_rec = cm.rec_energy_gev(e_mip, region)
        pass_sel = pass_digi & cm.selector_pass(
            e_rec, t_cell, rd.theta, rd.layer, region, ecal_map=ecal_map)

        # ---- gate 3: per-cell join vs actual Digi hits (full 64-bit cellIDs)
        regsel = flat["digi_region"] == ireg
        dev, did = flat["digi_event"][regsel], flat["digi_cellid"][regsel]
        de, dtm = flat["digi_e"][regsel], flat["digi_t"][regsel]
        # HCAL digi stores NPE, emulation e_mip is in MIP -> convert
        e_cmp = e_mip * (cm.PPD_PE_PER_MIP if region.startswith("Hcal") else 1.0)
        ia, ib, ncommon = join(rd.event[pass_digi], rd.cellid[pass_digi], dev, did)
        n_emul, n_act = int(pass_digi.sum()), len(dev)
        only_e, only_a = n_emul - ncommon, n_act - ncommon
        rel = np.abs(e_cmp[pass_digi][ia] - de[ib]) / np.maximum(de[ib], 1e-12)
        trel = np.abs(t_cell[pass_digi][ia] - dtm[ib])
        print(f"[{region}] gate3/Digi: emul {n_emul} / actual {n_act} cells, "
              f"matched {ncommon} (emul-only {only_e}, actual-only {only_a}); "
              f"median |dE|/E = {np.median(rel) if len(rel) else 0:.4f}, "
              f"max |dt| = {np.max(trel) if len(trel) else 0:.2e} ns")
        frac_unmatched = (only_e + only_a) / max(n_act, 1)
        if frac_unmatched > (0.25 if region.startswith("Hcal") else 0.05) or \
                (len(trel) and np.max(trel) > 1e-4):
            fails += 1
            print(f"[{region}] gate3 FAIL")

        # ---- gate 4: vs actual Sel hits. RealisticCaloReco truncates the
        # cellID to 32 bits (k4Reco bug, see PLAN.md) so Rec/Coned/Sel hit IDs
        # lose the x/y fields and a per-cell join is impossible. The layer
        # bits survive, so compare per-(event, layer) counts and energy sums.
        regsel = flat["sel_region"] == ireg
        sev, sid = flat["sel_event"][regsel], flat["sel_cellid"][regsel]
        se = flat["sel_e"][regsel]
        slay = coder_layer(flat, sid)
        n_emul, n_act = int(pass_sel.sum()), len(sev)
        key_e = rd.event[pass_sel].astype(np.int64) * 1000 + rd.layer[pass_sel]
        key_a = sev.astype(np.int64) * 1000 + slay
        allk = np.union1d(key_e, key_a)
        cnt_e = np.array([np.sum(key_e == k) for k in allk])
        cnt_a = np.array([np.sum(key_a == k) for k in allk])
        sum_e = np.array([e_rec[pass_sel][key_e == k].sum() for k in allk])
        sum_a = np.array([se[key_a == k].sum() for k in allk])
        dn = np.abs(cnt_e - cnt_a).sum() / max(cnt_a.sum(), 1)
        big = sum_a > 0.01  # only compare bins with >10 MeV
        rel = (np.abs(sum_e[big] - sum_a[big]) / sum_a[big]) if big.any() else np.array([0.0])
        print(f"[{region}] gate4/Sel (per event x layer): emul {n_emul} / actual "
              f"{n_act} cells, count mismatch fraction = {dn:.4f}, "
              f"median |dE|/E = {np.median(rel):.4f}")
        tol = 0.25 if region.startswith("Hcal") else 0.05
        # too few hits -> single-cell stochastic fluctuations dominate
        if n_act >= 20 and (dn > tol or np.median(rel) > tol):
            fails += 1
            print(f"[{region}] gate4 FAIL")

        # ---- gate 5: event sums (Sel, GeV)
        act = flat[f"ev_sel_sum_{region}"]
        emu = np.bincount(rd.event[pass_sel], weights=e_rec[pass_sel], minlength=nev)
        nsel = flat[f"ev_sel_n_{region}"]
        both = (act > 1e-6) & (nsel >= 5)
        if both.any():
            rel = np.abs(emu[both] - act[both]) / act[both]
            print(f"[{region}] gate5 event Sel sums: median rel dev = "
                  f"{np.median(rel):.4f}, worst = {np.max(rel):.4f}")
            if np.median(rel) > (0.05 if region.startswith("Hcal") else 0.01):
                fails += 1
                print(f"[{region}] gate5 FAIL")

    print("SANITY", "FAIL" if fails else "PASS", f"({fails} gate failures)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
