#!/usr/bin/env python3
"""Summarize a (variant) reco file into a tiny per-event npz: Sel sums per
region, leading Pandora cluster, gun kinematics. Runs INSIDE the container.

Usage: python3 summarize_reco.py input_reco.edm4hep.root output.npz
"""

import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from cutmodel import REGIONS  # noqa: E402

SEL_COLLECTIONS = {r: r + "CollectionSel" for r in REGIONS}


def main():
    infile, outfile = sys.argv[1], sys.argv[2]
    from podio.root_io import Reader

    reader = Reader(infile)
    events = reader.get("events")

    ev = {k: [] for k in ["gun_e", "gun_theta", "gun_phi",
                          "n_clusters", "lead_cluster_e", "sum_cluster_e"]}
    for r in REGIONS:
        ev[f"sel_sum_{r}"] = []
        ev[f"sel_n_{r}"] = []

    for frame in events:
        gun = next((p for p in frame.get("MCParticles")
                    if p.getGeneratorStatus() == 1), None)
        if gun is not None:
            mom = gun.getMomentum()
            pmag = math.sqrt(mom.x**2 + mom.y**2 + mom.z**2)
            ev["gun_e"].append(gun.getEnergy())
            ev["gun_theta"].append(math.acos(mom.z / pmag) if pmag else 0.0)
            ev["gun_phi"].append(math.atan2(mom.y, mom.x))
        else:
            for k in ("gun_e", "gun_theta", "gun_phi"):
                ev[k].append(0.0)

        try:
            clusters = list(frame.get("PandoraClusters"))
        except Exception:
            clusters = []
        ev["n_clusters"].append(len(clusters))
        ev["sum_cluster_e"].append(sum(c.getEnergy() for c in clusters))
        ev["lead_cluster_e"].append(max((c.getEnergy() for c in clusters),
                                        default=0.0))

        for r in REGIONS:
            try:
                coll = frame.get(SEL_COLLECTIONS[r])
                ev[f"sel_sum_{r}"].append(sum(h.getEnergy() for h in coll))
                ev[f"sel_n_{r}"].append(len(coll))
            except Exception:
                ev[f"sel_sum_{r}"].append(0.0)
                ev[f"sel_n_{r}"].append(0)

    np.savez_compressed(outfile, **{k: np.asarray(v) for k, v in ev.items()})
    print(f"wrote {outfile}: {len(ev['gun_e'])} events")


if __name__ == "__main__":
    main()
