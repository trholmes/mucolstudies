import pyLCIO
import ROOT
import glob

exec(open("./plotHelper.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1
pdgid = 11
append = "_electrons"
require_particle_id = True
min_true_E = 10
min_reco_E = 5

# Gather input files
#fnames = glob.glob("../v2.11/reco/photonGun_E_0_50/photonGun_E_0_50_reco_0*")
#fnames = glob.glob("../v2.11/reco/neutronGun_E_0_50/neutronGun_E_0_50_reco_0*")
fnames = glob.glob("../v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0*")
#fnames = glob.glob("../v2.11/reco/pionGun_pT_0_50/pionGun_pT_0_50_reco_0*")

# Set up histograms
#keys = ["calib", "uncalib", "calibtr"]
#keys = ["uncalib", "calib", "calibtrfix", "calibpid", "calibtrpid"]
#keys = ["uncalib", "calibtr", "hadreweightuncalib", "hadreweightcalibtrpid"]
keys = ["uncalib", "calibtrpid", "proxim", "assoc"]
comphists = {}
hists = {}
compvariables = {
        "eres": {"nbins": 100, "xmin": -2, "xmax": 2, "xlabel": "(E_{reco}-E_{true})/E_{true}"},
        "cluseres": {"nbins": 100, "xmin": -2, "xmax": 2, "xlabel": "(E_{reco}-E_{true})/E_{true}"},
        }
variables = {
        "pt":  {"nbins": 20, "xmin": 0, "xmax": 50, "xlabel": "p_T [GeV]"},
        "E":  {"nbins": 20, "xmin": 0, "xmax": 200, "xlabel": "E [GeV]"},
        "eta": {"nbins": 20, "xmin": -3, "xmax": 3, "xlabel": "#eta"},
        "phi": {"nbins": 20, "xmin": -3.5, "xmax": 3.5, "xlabel": "#phi"},
        "n": {"nbins": 20, "xmin": 0, "xmax": 20, "xlabel": "n"},
        }
objects = ["pfo", "mcp", "mcpmatch"]

for key in keys:
    comphists[key] = {}
    for var in compvariables:
        comphists[key][var] = ROOT.TH1F(f"{var}{key}", var, compvariables[var]["nbins"], compvariables[var]["xmin"], compvariables[var]["xmax"])

for key in keys:
    hists[key] = {}
    for obj in objects:
        hists[key][obj] = {}
        for var in variables:
            hists[key][obj][var] = ROOT.TH1F(f"{var}{key}{obj}", var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# Loop over events
for f in fnames:
    print(f)
    i = 0
    key = f.split("_")[-1].split(".")[0]
    if key not in keys: continue
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)


    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        pfoCollection = event.getCollection("PandoraPFOs")
        mcpCollection = event.getCollection("MCParticle")

        # Loop over the reconstructed objects and fill histograms
        true_E = -99
        n_mcp = 0
        mymcp = None
        for mcp in mcpCollection:
            if mcp.getGeneratorStatus()==1 and mcp.getEnergy() > min_true_E:
                true_E = mcp.getEnergy()
                mymcp = mcp
                mcptlv = getTLV(mymcp)
                n_mcp += 1
                hists[key]["mcp"]["pt"].Fill(mcptlv.Perp())
                hists[key]["mcp"]["E"].Fill(mcptlv.E())
                hists[key]["mcp"]["eta"].Fill(mcptlv.Eta())
                hists[key]["mcp"]["phi"].Fill(mcptlv.Phi())
        hists[key]["mcp"]["n"].Fill(n_mcp)

        # Pick highest energy particle of the right type in the event
        pfo_E = -99
        n_pfo = 0
        mypfo = None
        for pfo in pfoCollection:
            if not require_particle_id or abs(pfo.getType())==pdgid:
                if pfo.getEnergy() > min_reco_E:
                    n_pfo += 1
                    if pfo_E < pfo.getEnergy():
                        pfo_E = pfo.getEnergy()
                        mypfo = pfo
        hists[key]["pfo"]["n"].Fill(n_pfo)

        # Fill histograms with reco in them when both are available
        if n_mcp > 0 and n_pfo > 0:
            #myclu = mypfo.getClusters()[0]
            try: myclu = mypfo.getClusters()[0]
            except: continue
            comphists[key]["eres"].Fill((pfo_E - true_E)/true_E)
            comphists[key]["cluseres"].Fill((myclu.getEnergy() - true_E)/true_E)

            pfotlv = getTLV(mypfo)
            hists[key]["pfo"]["pt"].Fill(pfotlv.Perp())
            hists[key]["pfo"]["E"].Fill(pfotlv.E())
            hists[key]["pfo"]["eta"].Fill(pfotlv.Eta())
            hists[key]["pfo"]["phi"].Fill(pfotlv.Phi())
            hists[key]["mcpmatch"]["pt"].Fill(mcptlv.Perp())
            hists[key]["mcpmatch"]["E"].Fill(mcptlv.E())
            hists[key]["mcpmatch"]["eta"].Fill(mcptlv.Eta())
            hists[key]["mcpmatch"]["phi"].Fill(mcptlv.Phi())

        i+=1

# Make your plots
for h in compvariables:
    hist_dict = {}
    for key in keys: hist_dict[key] = comphists[key][h]
    #hist_dict = {"calibrated": comphists["calib"][h], "uncalibrated": comphists["uncalib"][h], "calibrated, track comp.": comphists["calibtr"][h]}
    #hist_dict = {"uncalibrated": comphists["uncalib"][h], "calibrated, track comp.": comphists["calibtr"][h]}
    #hist_dict = {"uncalibrated": comphists["uncalib"][h], "calibrated": comphists["calib"][h]}
    plotHistograms(hist_dict, f"plots/calib/calib_{h}{append}.png", compvariables[h]["xlabel"])

for h in variables:
    for obj in objects:
        hist_dict = {}
        for key in keys: hist_dict[key] = hists[key][obj][h]
        #hist_dict = {"uncalibrated": hists["uncalib"][obj][h], "calibrated, track comp.": hists["calibtr"][obj][h]}
        #hist_dict = {"uncalibrated": hists["uncalib"][obj][h], "calibrated": hists["calib"][obj][h]}
        plotHistograms(hist_dict, f"plots/calib/calib_{h}_{obj}{append}.png", variables[h]["xlabel"])

# Make efficiency plots
for h in variables:
    if h=="n": continue
    eff_dict = {}
    for key in keys: eff_dict[key] = ROOT.TEfficiency(hists[key]["mcpmatch"][h], hists[key]["mcp"][h])
    #        "uncalibrated": ROOT.TEfficiency(hists["uncalib"]["mcpmatch"][h], hists["uncalib"]["mcp"][h]),
    #        "calibrated, track comp.": ROOT.TEfficiency(hists["calibtr"]["mcpmatch"][h], hists["calibtr"]["mcp"][h])
    #        #"calibrated": ROOT.TEfficiency(hists["calib"]["mcpmatch"][h], hists["calib"]["mcp"][h])
    #        }
    plotEfficiencies(eff_dict, f"plots/calib/calib_eff_{h}{append}.png", variables[h]["xlabel"])
