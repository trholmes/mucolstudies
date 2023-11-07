import numpy as np
import pyLCIO
import glob
import ctypes

exec(open("helpers.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/reco/*.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/gen_muonGun/recoBIB/*.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/muonGun_1000/recoBIB/*.slcio")
print("Found %i files."%len(fnames))

# Define good particle
def isGood(tlv):
    if abs(tlv.Eta()) < 1.88:
        return True
    return False

# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1 and abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp() < 0.2:
        return True
    return False

#pfo_list = [-211, 13, 22, 2112]
pfo_list = [-14, -13, -11, 11, 12, 13, 14, 22, 211, 2112, 2212, 1000010020, 1000020040, 1000030080, 1000040080, 1000040090]
h_pfoType = ROOT.TH1F("pfo_types", "pfo_types", len(pfo_list), 0, len(pfo_list))
h_pfoType.SetCanExtend(ROOT.TH1.kAllAxes)

pt_hists = {}
for pdg in pfo_list:
    pt_hists["h_%s_pt"%(str(abs(pdg)))] = ROOT.TH1F("h_%s_pt"%(str(abs(pdg))), "", 20, 0, 1)

# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
pfo_types = []
i = 0
for f in fnames:
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")
        trkCollection = event.getCollection("SiTracks_Refitted")

        # Make counter variables
        n_mcp_mu = 0
        n_pfo_mu = 0
        has_mcp_mu = False
        has_pfo_mu = False
        my_pfo_mu = 0
        my_mcp_mu = 0


        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())

            pt_hists["h_%s_pt"%(str(abs(mcp.getPDG())))].Fill(mcp_tlv.Perp())
            h_pfoType.Fill(pfo_list.index(mcp.getPDG()))
            pfo_types.append(mcp.getPDG())

            if abs(mcp.getPDG())==13 and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_mu = True
                n_mcp_mu += 1
                my_mcp_mu = mcp_tlv

        # Loop over the reconstructed objects and fill histograms
        for pfo in pfoCollection:
            pfo_p = pfo.getMomentum()
            pfo_tlv = ROOT.TLorentzVector()
            pfo_tlv.SetPxPyPzE(pfo_p[0], pfo_p[1], pfo_p[2], pfo.getEnergy())

            #pfo_types.append(pfo.getType())
            #pt_hists["h_%s_pt"%(str(abs(pfo.getType())))].Fill(pfo_tlv.Perp())
            #h_pfoType.Fill(pfo_list.index(pfo.getType()))

            if abs(pfo.getType())==13:
                n_pfo_mu += 1
                has_pfo_mu = True
                my_pfo_mu = pfo_tlv     # Storing this to use for matching in the next loop

        #if not has_pfo_mu:
        #    for pfo in pfoCollection:
        #        h_pfoType.Fill(pfo_list.index(pfo.getType()))

        i+=1

        # Fill the histograms comparing properties

print(sorted(np.unique(pfo_types)))

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################

c = ROOT.TCanvas("c", "c")
h_pfoType.Draw("e hist")
for i, val in enumerate(pfo_list):
    h_pfoType.GetXaxis().SetBinLabel(i+1, str(val))
h_pfoType.GetXaxis().SetTitle("PDG ID of MCP Object")
c.SetLogy()
c.SaveAs("plots/pfo_types.png")

for h in pt_hists:
    c = ROOT.TCanvas("c"+h, "")
    pt_hists[h].Draw()
    pt_hists[h].GetXaxis().SetTitle("p_{T} [GeV]")
    latex = ROOT.TLatex()
    latex.DrawLatexNDC(.64, .85, "PDGID: %s"%h.split("_")[1])
    c.SetLogy()
    c.SaveAs("plots/"+h+".png")
