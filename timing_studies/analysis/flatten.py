#!/usr/bin/env python3
"""Flatten an EDM4hep reco file into numpy arrays for the timing-cut study.

Runs INSIDE the mucoll-sim container (needs podio); everything downstream
(emulate.py, plots) needs only numpy.

For each calorimeter region it stores, per event:
  - all SimCalorimeterHit cells (position, r, theta, layer, summed sim energy)
  - all their contributions, with delta_t = t_contrib - r_cell/c
  - the actual Digi and Sel hits (cellID, energy, time) for closure checks
  - gun kinematics from MCParticles, and Pandora cluster summaries

Usage (inside container):
  python3 flatten.py input_reco.edm4hep.root output.npz
"""

import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from cutmodel import C_MM_NS, REGIONS, SIM_COLLECTIONS, BitFieldCoder  # noqa: E402

# NB: sim collections are ECal*/HCal*, digi/sel are Ecal*/Hcal* (REGIONS spelling)
DIGI_COLLECTIONS = {r: r + "CollectionDigi" for r in REGIONS}
SEL_COLLECTIONS = {r: r + "CollectionSel" for r in REGIONS}


def get_encoding(reader, collname):
    """CellIDEncoding for a collection, from the metadata frame."""
    try:
        meta = reader.get("metadata")[0]
        enc = meta.get_parameter(f"{collname}__CellIDEncoding")
        if enc:
            return str(enc)
    except Exception as exc:
        print(f"WARNING: no CellIDEncoding for {collname}: {exc}")
    return None


def main():
    infile, outfile = sys.argv[1], sys.argv[2]
    from podio.root_io import Reader

    reader = Reader(infile)
    events = reader.get("events")
    nevents = len(events)
    print(f"{infile}: {nevents} events")

    # per-event scalars
    ev = {k: [] for k in ["gun_e", "gun_theta", "gun_phi", "gun_pdg",
                          "n_clusters", "lead_cluster_e", "lead_cluster_theta",
                          "lead_cluster_phi", "sum_cluster_e"]}
    for r in REGIONS:
        ev[f"digi_sum_{r}"] = []   # sum of Digi hit energies (digi units)
        ev[f"sel_sum_{r}"] = []    # sum of Sel hit energies (GeV, calibrated)
        ev[f"sel_n_{r}"] = []

    # per-cell / per-contribution jagged arrays
    cells = {k: [] for k in ["event", "region", "cellid", "layer", "r", "theta",
                             "sim_e", "n_contrib"]}
    contribs = {k: [] for k in ["dt", "e"]}
    # actual digi/sel hits for closure
    digihits = {k: [] for k in ["event", "region", "cellid", "e", "t"]}
    selhits = {k: [] for k in ["event", "region", "cellid", "e", "t"]}

    coders = {}
    first_frame_report = True

    for iev, frame in enumerate(events):
        if first_frame_report:
            print("collections:", sorted(frame.getAvailableCollections()))
            first_frame_report = False

        # --- gun particle: first generator-status-1 MCParticle
        mcps = frame.get("MCParticles")
        gun = None
        for p in mcps:
            if p.getGeneratorStatus() == 1:
                gun = p
                break
        if gun is not None:
            mom = gun.getMomentum()
            pmag = math.sqrt(mom.x**2 + mom.y**2 + mom.z**2)
            ev["gun_e"].append(gun.getEnergy())
            ev["gun_theta"].append(math.acos(mom.z / pmag) if pmag else 0.0)
            ev["gun_phi"].append(math.atan2(mom.y, mom.x))
            ev["gun_pdg"].append(gun.getPDG())
        else:
            for k in ["gun_e", "gun_theta", "gun_phi", "gun_pdg"]:
                ev[k].append(0)

        # --- Pandora clusters
        try:
            clusters = frame.get("PandoraClusters")
            cl_e = sorted((c.getEnergy() for c in clusters), reverse=True)
            ev["n_clusters"].append(len(cl_e))
            ev["sum_cluster_e"].append(sum(cl_e))
            if len(cl_e):
                lead = max(clusters, key=lambda c: c.getEnergy())
                ev["lead_cluster_e"].append(lead.getEnergy())
                ev["lead_cluster_theta"].append(lead.getITheta())
                ev["lead_cluster_phi"].append(lead.getPhi())
            else:
                ev["lead_cluster_e"].append(0.0)
                ev["lead_cluster_theta"].append(0.0)
                ev["lead_cluster_phi"].append(0.0)
        except Exception:
            for k in ["n_clusters", "sum_cluster_e", "lead_cluster_e",
                      "lead_cluster_theta", "lead_cluster_phi"]:
                ev[k].append(0)

        # --- calorimeter regions
        for ireg, region in enumerate(REGIONS):
            # sim hits + contributions
            try:
                simcoll = frame.get(SIM_COLLECTIONS[region])
            except Exception:
                simcoll = []
            if region not in coders:
                enc = get_encoding(reader, SIM_COLLECTIONS[region])
                coders[region] = BitFieldCoder(enc) if enc else None
                print(f"{region}: encoding = {enc}")
            coder = coders[region]

            for hit in simcoll:
                pos = hit.getPosition()
                rr = math.sqrt(pos.x**2 + pos.y**2 + pos.z**2)
                tof = rr / C_MM_NS
                nc = 0
                for c in hit.getContributions():
                    contribs["dt"].append(c.getTime() - tof)
                    contribs["e"].append(c.getEnergy())
                    nc += 1
                cells["event"].append(iev)
                cells["region"].append(ireg)
                cells["cellid"].append(hit.getCellID())
                lay = int(coder.get(hit.getCellID(), "layer")) if coder else -1
                cells["layer"].append(lay)
                cells["r"].append(rr)
                cells["theta"].append(math.acos(pos.z / rr) if rr else 0.0)
                cells["sim_e"].append(hit.getEnergy())
                cells["n_contrib"].append(nc)

            # actual digi / sel hits
            for colls, store in ((DIGI_COLLECTIONS, digihits),
                                 (SEL_COLLECTIONS, selhits)):
                esum, n = 0.0, 0
                try:
                    coll = frame.get(colls[region])
                except Exception:
                    coll = []
                for hit in coll:
                    store["event"].append(iev)
                    store["region"].append(ireg)
                    store["cellid"].append(hit.getCellID())
                    store["e"].append(hit.getEnergy())
                    store["t"].append(hit.getTime())
                    esum += hit.getEnergy()
                    n += 1
                if store is digihits:
                    ev[f"digi_sum_{region}"].append(esum)
                else:
                    ev[f"sel_sum_{region}"].append(esum)
                    ev[f"sel_n_{region}"].append(n)

    out = {}
    for prefix, d in (("ev_", ev), ("cell_", cells), ("contrib_", contribs),
                      ("digi_", digihits), ("sel_", selhits)):
        for k, v in d.items():
            dtype = np.uint64 if k == "cellid" else (
                np.int32 if k in ("event", "region", "layer", "n_contrib",
                                  "gun_pdg", "n_clusters") else np.float64)
            out[prefix + k] = np.asarray(v, dtype=dtype)
    np.savez_compressed(outfile, **out)
    ncells = len(out["cell_event"])
    ncontrib = len(out["contrib_dt"])
    print(f"wrote {outfile}: {nevents} events, {ncells} cells, {ncontrib} contributions")

    # flatten self-consistency (verification gate 1): sum of contribution
    # energies must reproduce the sim hit energy
    off = np.concatenate([[0], np.cumsum(out["cell_n_contrib"])])
    sums = np.add.reduceat(out["contrib_e"], off[:-1])
    sums[out["cell_n_contrib"] == 0] = 0.0
    bad = np.abs(sums - out["cell_sim_e"]) > 1e-6 + 1e-4 * out["cell_sim_e"]
    print(f"gate1 contribution-sum closure: {bad.sum()} / {ncells} cells off")


if __name__ == "__main__":
    main()
