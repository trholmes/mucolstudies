import pyLCIO
import glob
import ctypes

exec(open("./plotHelper.py").read())

# ############## SETUP #############################
ROOT.gROOT.SetBatch()

# Set up some options
max_events = 10000
magnetic_field = 5.00
max_E = 100
filenames = "/data/fmeloni/DataMuC_MuColl10_v0A/v2/recoBIB/photonGun_*"
min_energy = 20
max_eta = 2.6

settings = {
        "labelname": {  "ph": "Photon",
                        "mu": "Muon",
                        "el": "Electron",
                        "ne": "Neutron",},
        "plotdir":{ "ph": "photons",
                    "mu": "muons",
                    "el": "electrons",
                    "ne": "neutrons",},
        "pdgid":  { "ph": 22,
                    "mu": 13,
                    "el": 11,
                    "ne": 2112,},
        "mass":   { "ph": 0,
                    "mu": 0.106,
                    "el": 0.000511,
                    "ne": 0.9396,}
}

# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
variables = {}
variables["E"] =  {"nbins": 30, "xmin": 0, "xmax": max_E,   "title": "E [GeV]"}
variables["pt"] =  {"nbins": 30, "xmin": 0, "xmax": max_E,   "title": "p_{T} [GeV]"}
variables["eta"] = {"nbins": 30, "xmin": -3, "xmax": 3,     "title": "#eta"}
variables["phi"] = {"nbins": 30, "xmin": -3.5, "xmax": 3.5, "title": "#phi"}
variables["n"] =   {"nbins": 200, "xmin": 0, "xmax": 1000,     "title": "n"}
variables["nzoom"] =   {"nbins": 50, "xmin": 0, "xmax": 50,     "title": "n"}

hists = {}
for var in variables:
    hists[var] = {}
    for obj in settings["pdgid"]:
        hists[var][obj] = ROOT.TH1F(obj+"_"+var, obj+"_"+var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
samples = glob.glob(filenames)
fnames = []
for s in samples:
    fnames += glob.glob(f"{s}/*.slcio")
print("Found %i files."%len(fnames))


# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks_Refitted"])
i = 0
for f in fnames:
    reader.open(f)
    if max_events > 0 and i >= max_events: break

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        #mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")

        # Initialize counts
        counts = {}
        for obj in settings["pdgid"]: counts[obj] = 0

        # Loop over PFOs
        for pfo in pfoCollection:
            pfo_tlv = getTLV(pfo)
            if pfo_tlv.E() < min_energy: continue
            if abs(pfo_tlv.Eta()) > max_eta: continue

            for obj in settings["pdgid"]:
                if abs(pfo.getType()) == settings["pdgid"][obj]:
                    hists["eta"][obj].Fill(pfo_tlv.Eta())
                    hists["pt"][obj].Fill(pfo_tlv.Perp())
                    hists["E"][obj].Fill(pfo_tlv.E())
                    hists["phi"][obj].Fill(pfo_tlv.Phi())

                    counts[obj] += 1

        for obj in settings["pdgid"]:
            hists["n"][obj].Fill(counts[obj])
            hists["nzoom"][obj].Fill(counts[obj])

        i += 1

# Make plots
for var in hists:
    h_to_plot = {}
    for obj in hists[var]:
        h_to_plot[obj] = hists[var][obj]
    plotHistograms(h_to_plot, f"plots/fakes/comp_{var}.png", variables[var]["title"], "Count", logy=True)


