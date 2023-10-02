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
max_events = 10

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
        for var in variables["obj"]:
            hists[s][obj+"_"+var] = ROOT.TH1F(s+"_"+obj+"_"+var, s, variables["obj"][var]["nbins"], variables["obj"][var]["xmin"], variables["obj"][var]["xmax"])

# Loop over the different samples
for s in files:
    print("Working on sample", s)
    i = 0

    # Loop over the files in a sample
    for f in files[s]:
        if max_events > 0 and i >= max_events: break
        reader = Reader(f)

        # Loop over events in each file
        for event in reader.get('events'):
            if max_events > 0 and i >= max_events: break
            if i%100 == 0: print("\tProcessing event:", i)

            # Get the collections we care about
            mcps = event.get("MCParticle")
            pfos = event.get("PandoraPFOs")
            trks = event.get('SiTracks_Refitted')

            # Loop over MCPs
            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                mcp_tlv = getTLV(mcp)
                fillObjHists(hists[s], "mcp", mcp_tlv)

                # Look at electrons only
                if abs(mcp.getPDG()) == 11:
                    fillObjHists(hists[s], "mcp_el", mcp_tlv)

            # Loop over PFOs
            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                fillObjHists(hists[s], "pfo", pfo_tlv)

                # Look at electrons only
                if abs(pfo.getType()) == 11:
                    fillObjHists(hists[s], "pfo_el", pfo_tlv)

            # Iterate counter
            i += 1

# Draw all the 1D histograms you filled
for i, h in enumerate(hists[s]):

    c = ROOT.TCanvas("c%i"%i, "c%i"%i)
    for j, s in enumerate(hists):
        if j==0:
            hists[s][h].Draw()
            hists[s][h].GetXaxis().SetTitle(h)
            hists[s][h].GetYaxis().SetTitle("Entries")
        else:
            hists[s][h].Draw("same")
        hists[s][h].SetLineColor(colors[j])
        hists[s][h].SetMarkerColor(colors[j])

    c.BuildLegend()
    c.SaveAs("plots/electrons/"+h+".png")
