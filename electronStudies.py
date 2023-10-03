# Note: This script was to test out opening k4hep files with ROOT. An actual example can be found here:
# https://gitlab.desy.de/ftx-sft-key4hep/tutorials/-/blob/main/edm4hep_analysis/edm4hep_python.ipynb

# Before running this script be sure to run the commands needed to access the software:
# apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest
# apptainer run --no-home -B /collab/project/snowmass21/data/muonc:/data -B /home/$USER k4toroid.sif
# source /setup.sh

import glob
import ROOT
from ROOT import edm4hep
from podio.root_io import Reader
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Open the edm4hep files with ROOT
samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/k4reco/electronGun*")
files = {}
for s in samples:
    sname = s.split("/")[-1]
    files[sname] = glob.glob(f"{s}/*.root")

# Set up histograms
hists = {}
for s in files:
    hists[s] = {}
    for obj in ["pfo", "pfo_el", "mcp", "mcp_el", "mcp_el_match"]:
        for vtype in ["obj", "evt"]:
            for var in variables[vtype]:
                hists[s][obj+"_"+var] = ROOT.TH1F(s+"_"+obj+"_"+var, s, variables[vtype][var]["nbins"], variables[vtype][var]["xmin"], variables[vtype][var]["xmax"])

# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1 and abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp() < 0.2:
        return True
    return False

# Loop over the different samples
for s in files:
    print("Working on sample", s)
    i = 0

    # Loop over the files in a sample
    for f in files[s]:
        if max_events > 0 and i >= max_events: break

        reader = Reader(f)
        reader.collections = ["MCParticle", "PandoraPFOs"]#, "SiTracks_Refitted"]
        print(type(reader))
        print(dir(reader))


        # Loop over events in each file
        for event in reader.get('events'):
            if max_events > 0 and i >= max_events: break
            if i%100 == 0: print("\tProcessing event:", i)

            # Make counters and key objects to track
            n_pfo = 0
            n_mcp = 0
            n_pfo_el = 0
            n_mcp_el = 0
            n_matched_el = 0
            my_mcp_el = None

            # Get the collections we care about
            mcps = event.get("MCParticle")
            pfos = event.get("PandoraPFOs")
            trks = event.get("SiTracks_Refitted")
            #lnks = event.get("MCParticle_SiTracks_Refitted")

            ######## Loop over MCPs

            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                mcp_tlv = getTLV(mcp)
                fillObjHists(hists[s], "mcp", mcp_tlv)
                n_mcp += 1

                # Look at electrons only
                if abs(mcp.getPDG()) == 11:
                    fillObjHists(hists[s], "mcp_el", mcp_tlv)
                    my_mcp_el = mcp_tlv
                    n_mcp_el += 1

            hists[s]["mcp_n"].Fill(n_mcp)
            hists[s]["mcp_el_n"].Fill(n_mcp_el)

            ######## Loop over PFOs

            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                fillObjHists(hists[s], "pfo", pfo_tlv)
                n_pfo += 1

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

            #for trk in trks:

            # Iterate counter
            i += 1

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

