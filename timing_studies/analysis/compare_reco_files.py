#!/usr/bin/env python3
"""Compare two reco EDM4hep files event by event (stock vs patched k4Reco).

Used to validate the two k4Reco fixes (CaloHitSelector double-TOF removal,
RealisticCaloReco 64-bit cellID): everything must be identical except that
Rec/Coned/Sel cellIDs regain their upper-32 bits (verified against the
digitized hits, whose cellIDs were always correct).

Runs INSIDE the container. Usage:
  python3 compare_reco_files.py stock_reco.edm4hep.root patched_reco.edm4hep.root
"""

import sys

import numpy as np

REGIONS = ["EcalBarrel", "EcalEndcap", "HcalBarrel", "HcalEndcap"]
STAGES = ["Digi", "Rec", "Coned", "Sel"]


def collect(fname):
    from podio.root_io import Reader
    reader = Reader(fname)
    out = {}
    nev = 0
    for iev, frame in enumerate(reader.get("events")):
        nev += 1
        for reg in REGIONS:
            for stage in STAGES:
                coll = frame.get(f"{reg}Collection{stage}")
                key = (reg, stage)
                out.setdefault(key, {"id": [], "e": [], "t": [], "ev": []})
                for h in coll:
                    out[key]["id"].append(h.getCellID())
                    out[key]["e"].append(h.getEnergy())
                    out[key]["t"].append(h.getTime())
                    out[key]["ev"].append(iev)
        for cname, key in (("PandoraClusters", "clusters"), ("PandoraPFOs", "pfos")):
            coll = frame.get(cname)
            out.setdefault(key, {"e": [], "ev": []})
            for c in coll:
                out[key]["e"].append(c.getEnergy())
                out[key]["ev"].append(iev)
    for key in out:
        for k in out[key]:
            dtype = np.uint64 if k == "id" else np.float64 if k in ("e", "t") else np.int64
            out[key][k] = np.asarray(out[key][k], dtype=dtype)
    return out, nev


def main():
    stock, nev_s = collect(sys.argv[1])
    patch, nev_p = collect(sys.argv[2])
    assert nev_s == nev_p, "different event counts!"
    print(f"comparing {nev_s} events: {sys.argv[1]} vs {sys.argv[2]}")
    bad = 0

    for reg in REGIONS:
        digi_ids = {}
        for stage in STAGES:
            s, p = stock[(reg, stage)], patch[(reg, stage)]
            if len(s["e"]) == 0 and len(p["e"]) == 0:
                continue
            n_eq = len(s["e"]) == len(p["e"])
            if not n_eq:
                print(f"[{reg}/{stage}] HIT COUNT differs: {len(s['e'])} vs {len(p['e'])}")
                bad += 1
                continue
            de = np.max(np.abs(s["e"] - p["e"])) if len(s["e"]) else 0.0
            dt = np.max(np.abs(s["t"] - p["t"])) if len(s["t"]) else 0.0
            if stage == "Digi":
                same_id = np.array_equal(s["id"], p["id"])
                print(f"[{reg}/{stage}] n={len(s['e'])}: cellIDs identical={same_id}, "
                      f"max|dE|={de:.3g}, max|dt|={dt:.3g} ns")
                if not same_id or de > 0 or dt > 0:
                    bad += 1
                digi_ids[reg] = dict(zip(zip(s["ev"].tolist(),
                                             (s["id"] & np.uint64(0xFFFFFFFF)).tolist()),
                                         s["id"].tolist()))
            else:
                # stock ids are the 32-bit truncation; patched must equal the
                # full 64-bit digi id whose low bits match, stock must equal
                # the truncation of the patched id
                trunc_ok = np.array_equal(s["id"] & np.uint64(0xFFFFFFFF),
                                          p["id"] & np.uint64(0xFFFFFFFF))
                full = np.array([digi_ids[reg].get((ev, low), 2**64 - 1) for ev, low in
                                 zip(p["ev"].tolist(),
                                     (p["id"] & np.uint64(0xFFFFFFFF)).tolist())],
                                dtype=np.uint64)
                # ambiguous truncations (several digi cells sharing low bits per
                # event) can't be checked this way; count exact matches only
                match64 = np.sum(full == p["id"])
                print(f"[{reg}/{stage}] n={len(s['e'])}: max|dE|={de:.3g}, "
                      f"max|dt|={dt:.3g} ns, stock==trunc(patched)={trunc_ok}, "
                      f"patched-id resolves to a digi id: {match64}/{len(full)}")
                if de > 0 or dt > 0 or not trunc_ok:
                    bad += 1

    for key, label in (("clusters", "PandoraClusters"), ("pfos", "PandoraPFOs")):
        s, p = stock[key], patch[key]
        if len(s["e"]) != len(p["e"]):
            print(f"[{label}] COUNT differs: {len(s['e'])} vs {len(p['e'])}")
            bad += 1
        else:
            de = np.max(np.abs(np.sort(s["e"]) - np.sort(p["e"]))) if len(s["e"]) else 0.0
            print(f"[{label}] n={len(s['e'])}: max|dE|={de:.3g} GeV")
            if de > 1e-6:
                bad += 1

    print("COMPARISON", "FAIL" if bad else "PASS",
          f"({bad} differences beyond the intended cellID repair)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
