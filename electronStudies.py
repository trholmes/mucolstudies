# Note: This script was to test out opening k4hep files with ROOT. An actual example can be found here:
# https://gitlab.desy.de/ftx-sft-key4hep/tutorials/-/blob/main/edm4hep_analysis/edm4hep_python.ipynb

# Before running this script be sure to run the commands needed to access the software:
# apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest
# apptainer run --no-home -B /collab/project/snowmass21/data/muonc:/data -B /home/$USER k4toroid.sif
# source /setup.sh

import glob
import ROOT

#from ROOT import edm4hep
#from podio.root_io import Reader
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Open the edm4hep files with ROOT
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/k4reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/electronGun*")
samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco_highrange/electronGun*")
files = {}
#slices = ["250_1000"]
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices: files[f"electronGun_pT_{s}"] = []
for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: continue
    #files[sname] = glob.glob(f"{s}/*.root")
    files[sname] = glob.glob(f"{s}/*.slcio")

# Set up histograms
hists = {}
for s in files:
    hists[s] = {}
    for obj in ["pfo", "pfo_el", "mcp", "mcp_el", "mcp_el_match", "trk", "trk_el_match"]:
        for vtype in ["obj", "evt"]:
            for var in variables[vtype]:
                hists[s][obj+"_"+var] = ROOT.TH1F(s+"_"+obj+"_"+var, s, variables[vtype][var]["nbins"], variables[vtype][var]["xmin"], variables[vtype][var]["xmax"])

hists2d = {}
for s in files:
    hists2d[s] = {}
    hists2d[s]["trk_eta_v_trk_pt"] = ROOT.TH2F(f"trk_eta_v_trk_pt_{s}", f"trk_eta_v_trk_pt_{s}", 30,-3,3,30,0,3000)
    hists2d[s]["trk_eta_v_trk_phi"] = ROOT.TH2F(f"trk_eta_v_trk_phi_{s}", f"trk_eta_v_trk_phi_{s}", 30,-3,3,30,-3,3)
    hists2d[s]["trk_eta_v_trk_n"] = ROOT.TH2F(f"trk_eta_v_trk_n_{s}", f"trk_eta_v_trk_n_{s}", 30,-3,3,20,0,20)
    hists2d[s]["trk_eta_v_mcp_eta"] = ROOT.TH2F(f"trk_eta_v_mcp_eta_{s}", f"trk_eta_v_mcp_eta_{s}", 30,-3,3,30,-3,3)
    hists2d[s]["trk_eta_v_mcp_pt"] = ROOT.TH2F(f"trk_eta_v_mcp_pt_{s}", f"trk_eta_v_mcp_pt_{s}", 30,-3,3,30,0,3000)
    hists2d[s]["trk_eta_v_mcp_phi"] = ROOT.TH2F(f"trk_eta_v_mcp_phi_{s}", f"trk_eta_v_mcp_phi_{s}", 30,-3,3,30,-3,3)
    hists2d[s]["trk_pt_v_mcp_pt"] = ROOT.TH2F(f"trk_pt_v_mcp_pt_{s}", f"trk_pt_v_mcp_pt_{s}", 30,0,3000,30,0,3000)
    hists2d[s]["pfo_pt_v_mcp_pt"] = ROOT.TH2F(f"pfo_pt_v_mcp_pt_{s}", f"pfo_pt_v_mcp_pt_{s}", 30,0,3000,30,0,3000)
    hists2d[s]["mcp_E_v_mcp_p"] = ROOT.TH2F(f"mcp_E_v_mcp_p_{s}", f"mcp_E_v_mcp_p_{s}", 30,0,1000,30,0,1000)



# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1 and abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp() < 0.2:
        return True
    return False

# Create a reader object to use for the rest of the time
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks_Refitted"])

# Loop over the different samples
for s in files:
    print("Working on sample", s)
    i = 0

    # Loop over the files in a sample
    for f in files[s]:
        if max_events > 0 and i >= max_events: break

        #reader = Reader(f)
        #reader.collections = ["MCParticle", "PandoraPFOs"]#, "SiTracks_Refitted"]
        #print(type(reader))
        #print(dir(reader))

        reader.open(f)

        # Loop over events in each file
        #for event in reader.get('events'):
        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i%100 == 0: print("\tProcessing event:", i)

            # Make counters and key objects to track
            n_trk = 0
            n_matched_trk = 0
            n_pfo = 0
            n_mcp = 0
            n_pfo_el = 0
            n_mcp_el = 0
            n_matched_el = 0
            my_mcp_el = None

            # Get the collections we care about
            #mcps = event.get("MCParticle")
            #pfos = event.get("PandoraPFOs")
            #trks = event.get("SiTracks_Refitted")
            #lnks = event.get("MCParticle_SiTracks_Refitted")
            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")
            trks = event.getCollection("SiTracks_Refitted")

            ######## Loop over MCPs

            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                mcp_tlv = getTLV(mcp)
                if abs(mcp_tlv.Eta())>2: continue
                fillObjHists(hists[s], "mcp", mcp_tlv)
                momentum = math.sqrt(mcp.getMomentum()[0]**2+mcp.getMomentum()[1]**2+mcp.getMomentum()[2]**2)
                hists2d[s]["mcp_E_v_mcp_p"].Fill(mcp.getEnergy(), momentum)
                n_mcp += 1


                # Look at electrons only
                if abs(mcp.getPDG()) == 11:
                    fillObjHists(hists[s], "mcp_el", mcp_tlv)
                    my_mcp_el = mcp_tlv
                    n_mcp_el += 1

            # Only consider events that had at least one electron in our fiducial region
            if n_mcp < 1: continue

            hists[s]["mcp_n"].Fill(n_mcp)
            hists[s]["mcp_el_n"].Fill(n_mcp_el)

            ######## Loop over PFOs

            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                fillObjHists(hists[s], "pfo", pfo_tlv)
                n_pfo += 1

                hists2d[s]["pfo_pt_v_mcp_pt"].Fill(pfo_tlv.Perp(), my_mcp_el.Perp())

                # Look at electrons only
                if abs(pfo.getType()) == 11:
                    fillObjHists(hists[s], "pfo_el", pfo_tlv)
                    n_pfo_el += 1

                    # Look at electrons matched to the gun electron
                    if isMatched(pfo_tlv, my_mcp_el):
                        fillObjHists(hists[s], "mcp_el_match", my_mcp_el)
                        n_matched_el += 1

            hists[s]["pfo_n"].Fill(n_pfo)
            hists[s]["pfo_el_n"].Fill(n_pfo_el)
            hists[s]["mcp_el_match_n"].Fill(n_matched_el)

            ######## Loop over tracks
            for trk in trks:
                trk_tlv = getTrackTLV(trk, m=0.0005)
                fillObjHists(hists[s], "trk", trk_tlv)

                hists2d[s]["trk_eta_v_trk_pt"].Fill(trk_tlv.Eta(), trk_tlv.Perp())
                hists2d[s]["trk_eta_v_trk_phi"].Fill(trk_tlv.Eta(), trk_tlv.Phi())
                hists2d[s]["trk_eta_v_trk_n"].Fill(trk_tlv.Eta(), len(trks))
                hists2d[s]["trk_eta_v_mcp_eta"].Fill(trk_tlv.Eta(), my_mcp_el.Eta())
                hists2d[s]["trk_eta_v_mcp_pt"].Fill(trk_tlv.Eta(), my_mcp_el.Perp())
                hists2d[s]["trk_eta_v_mcp_phi"].Fill(trk_tlv.Eta(), my_mcp_el.Phi())
                hists2d[s]["trk_pt_v_mcp_pt"].Fill(trk_tlv.Perp(), my_mcp_el.Perp())
                n_trk += 1

                if isMatched(trk_tlv, my_mcp_el):
                    fillObjHists(hists[s], "trk_el_match", my_mcp_el)
                    n_matched_trk += 1

            hists[s]["trk_n"].Fill(n_trk)
            hists[s]["trk_el_match_n"].Fill(n_matched_trk)

            # Iterate counter
            i += 1

        reader.close()

print("ntracks entries", hists[s]["trk_n"].GetEntries())
print("ntracks integral from 1", hists[s]["trk_n"].GetEntries())
print("trackpt entries", hists[s]["trk_pt"].GetEntries())

# Draw all the 1D histograms you filled
for i, h in enumerate(hists[s]):

    # Collect hists that go on a single plot
    hists_to_plot = {}
    for j, s in enumerate(hists):
        hists_to_plot[s] = hists[s][h]
    var_name = h.split("_")[-1]
    try:
        xlabel = variables["obj"][var_name]["label"]
    except:
        xlabel = variables["evt"][var_name]["label"]

    # Call plotting function
    plotHistograms(hists_to_plot, "plots/electrons/"+h+".png", xlabel, "Entries")
    plotHistograms(hists_to_plot, "plots/electrons/"+h+".root", xlabel, "Entries")
    #plotHistograms(hists_to_plot, "plots/electrons_no_el_req/"+h+".png", xlabel, "Entries")
    #plotHistograms(hists_to_plot, "plots/electrons_no_el_req/"+h+".root", xlabel, "Entries")

for s in hists2d:
    for h in hists2d[s]:
        c = ROOT.TCanvas("can", "can")
        hists2d[s][h].Draw("colz")
        hists2d[s][h].GetXaxis().SetTitle(h.split("_v_")[0])
        hists2d[s][h].GetYaxis().SetTitle(h.split("_v_")[1])
        c.SaveAs(f"plots/electrons/{hists2d[s][h].GetName()}.png")
