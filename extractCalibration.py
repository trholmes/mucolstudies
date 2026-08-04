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
pdgid = 211
append = "_pions"
require_particle_id = False
min_true_E = 5
min_reco_E = 1

central_hadECAL = 1.24223718397
central_hadHCAL = 1.01799349172
fit_slope = -0.0965212
fit_const = 0.884234

#proposed_hadHCAL = central_hadHCAL
#proposed_hadECAL = central_hadECAL
#proposed_hadECAL = central_hadECAL*(1+fit_slope)
#proposed_hadHCAL /= central_hadHCAL
#proposed_hadECAL /= central_hadHCAL

proposed_hadECAL = 1.38
proposed_hadHCAL = 1.25


# Gather input files
fnames = glob.glob("../v2.11/reco/pionGun_pT_0_50/pionGun_pT_0_50_reco_0_uncalib.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v6/reco/neutronGun_E_50_250/*")

# Set up histograms
e_slices = ["0_50", "50_100", "100_150", "150_200", "200_250"]
hists = {}
for sl in e_slices:
    hists[sl] = ROOT.TH2F(sl, sl, 20, 0, 1, 40, 0, 1.5)

# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "PandoraClusters"])
i = 0
for f in fnames:
    reader.open(f)

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        pfoCollection = event.getCollection("PandoraPFOs")
        mcpCollection = event.getCollection("MCParticle")
        i+=1

        # Loop over the reconstructed objects and fill histograms
        true_E = -99
        n_mcp = 0
        mymcp = None
        for mcp in mcpCollection:
            if mcp.getGeneratorStatus()==1 and mcp.getEnergy() > min_true_E:
                mcptlv = getTLV(mcp)
                if abs(mcptlv.Eta()) < 1.25 or abs(mcptlv.Eta())>2.4: continue
                true_E = mcp.getEnergy()
                mymcp = mcp
                n_mcp += 1
        if n_mcp==0: continue

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

        # Fill histograms with reco in them when both are available
        if n_mcp > 0 and n_pfo > 0:
            evectors = mypfo.getClusters()[0].getSubdetectorEnergies()
            hcal_frac = evectors[1]/(evectors[0]+evectors[1])
            corrected_E = evectors[0]*proposed_hadECAL + evectors[1]*proposed_hadHCAL
            #print("PFO E:", pfo_E, "Scaled Raw E:", corrected_E)
            response = corrected_E / true_E
            for sl in e_slices:
                if true_E < int(sl.split("_")[-1]):
                    hists[sl].Fill(hcal_frac, response)

# Make your plots
profiles = {}
for sl in e_slices:
    print(sl)
    can = ROOT.TCanvas("c"+sl)
    hists[sl].Draw()
    hists[sl].GetXaxis().SetTitle("HCAL Fraction")
    hists[sl].GetYaxis().SetTitle("Response")
    profiles[sl] = hists[sl].ProfileX("_pfx")
    f = ROOT.TF1(f"fit{sl}", "pol1", 0.1, 1)
    #f = ROOT.TF1(f"fit{sl}", "pol1")
    f.SetLineColor(ROOT.kRed)
    hists[sl].Fit(f"fit{sl}", "R")
    p = f.GetParameters()
    print(p[0], p[1])
    can.SaveAs(f"plots/calib/response_{sl}.png")

plotHistograms(profiles, "plots/calib/response_profiles.png")

print("Effective values:")
print("Hadronic ECAL:", proposed_hadECAL)
print("Hadronic HCAL:", proposed_hadHCAL)
