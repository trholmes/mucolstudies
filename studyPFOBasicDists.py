import pyLCIO
import glob
import math
import os

# This helper is expected to define:
#  - ROOT
#  - getTLV(obj): returns ROOT.TLorentzVector for MCParticle/ReconstructedParticle/etc.
#  - plotHistograms(hmap, outpath, xtitle, ytitle)
exec(open("./plotHelper.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running
ROOT.gROOT.SetBatch()

# ---- user options ----
max_events = 1000         # total events to process (across all files). Set <=0 for all.
obj_type = "ne"            # choose: "ph", "mu", "el", "ne" (extendable below)
max_E = 12000                 # x-axis max for E/pt plots
pick_event = True
# ----------------------

settings = {
    "fnames": {
        "ph": "/data/fmeloni/DataMuC_MuColl10_v0A/v2/recoBIB/photonGun_E_0*",
        "mu": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun*",
        "el": "/data/fmeloni/DataMuC_MAIA_v0/v2/reco/electronGun_pT_0_50",
        "ne": "/data/fmeloni/DataMuC_MAIA_v0/v8/recoBIB/neutronGun_E_0_*",
    },
    "labelname": {
        "ph": "Photon",
        "mu": "Muon",
        "el": "Electron",
        "ne": "Neutron",
    },
    "plotdir": {
        "ph": "photons_basic",
        "mu": "muons_basic",
        "el": "electrons_basic",
        "ne": "neutrons_basic",
    },
    "pdgid": {
        "ph": 22,
        "mu": 13,
        "el": 11,
        "ne": 2112,
    },
}

if obj_type not in settings["pdgid"]:
    raise ValueError(f"Unknown obj_type='{obj_type}'. Add it to settings['pdgid']/['fnames']/['labelname']/['plotdir'].")


def isGoodTLV(tlv):
    if abs(tlv.Eta()) > 2.4: return False
    if tlv.E() < 5: return False
    return True

print("Running basic PFO distributions on", settings["labelname"][obj_type])

# Gather input files
samples = glob.glob(settings["fnames"][obj_type])
fnames = []
for s in samples:
    fnames += glob.glob(f"{s}/*.slcio")
print("Found %i files." % len(fnames))

# ############## HISTOGRAM DEFINITIONS #############################

variables = {
    "E":   {"nbins": 30, "xmin": 0,    "xmax": max_E,  "title": "E [GeV]"},
    "pt":  {"nbins": 30, "xmin": 0,    "xmax": max_E,  "title": "p_{T} [GeV]"},
    "eta": {"nbins": 30, "xmin": -3.0, "xmax": 3.0,    "title": "#eta"},
    "phi": {"nbins": 30, "xmin": -3.5, "xmax": 3.5,    "title": "#phi"},
    "n":   {"nbins": 40, "xmin": 0,    "xmax": 200,     "title": f"n({settings['labelname'][obj_type]} PFO) per event"},
}

hists = {}
for var, cfg in variables.items():
    hists[var] = ROOT.TH1F(f"pfo_{var}", f"{settings['labelname'][obj_type]} PFOs", cfg["nbins"], cfg["xmin"], cfg["xmax"])

# ############## EVENT LOOP #############################

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["PandoraPFOs"])

select_event = False
iev = 0
for f in fnames:
    print("Opening", f)
    reader.open(f)
    for event in reader:
        if max_events > 0 and iev >= max_events:
            break
        if iev % 100 == 0:
            print(f"Processing event {iev}")

        pfoCollection = event.getCollection("PandoraPFOs")

        n_selected = 0
        if pick_event:
            for pfo in pfoCollection:
                if abs(pfo.getType()) != settings["pdgid"][obj_type]: continue
                tlv = getTLV(pfo)
                if not isGoodTLV(tlv): continue
                if tlv.E()>2000:
                    select_event=True
                    print(iev, tlv.E())

        for pfo in pfoCollection:
            if abs(pfo.getType()) != settings["pdgid"][obj_type]: continue
            tlv = getTLV(pfo)
            if not isGoodTLV(tlv): continue
            n_selected += 1

            if not pick_event or select_event:
                hists["E"].Fill(tlv.E())
                hists["pt"].Fill(tlv.Perp())
                hists["eta"].Fill(tlv.Eta())
                hists["phi"].Fill(tlv.Phi())

        if not pick_event or select_event: hists["n"].Fill(n_selected)

        if pick_event and select_event: break
        iev += 1

    reader.close()
    if pick_event and select_event: break
    if max_events > 0 and iev >= max_events:
        break

print(f"Processed {iev} events.")

# ############## PLOTTING #############################

outdir = f"plots/{settings['plotdir'][obj_type]}"
os.makedirs(outdir, exist_ok=True)

# Use the same plotting helper as the original script to keep style consistent.
for var, cfg in variables.items():
    plotHistograms(
        {settings["labelname"][obj_type]: hists[var]},
        f"{outdir}/pfo_{var}.png",
        cfg["title"],
        "Count",
        logy=True
    )

print("Wrote plots to", outdir)
