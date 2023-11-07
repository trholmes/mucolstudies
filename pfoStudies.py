# Note: This script was to test out opening k4hep files with ROOT. An actual example can be found here:
# https://gitlab.desy.de/ftx-sft-key4hep/tutorials/-/blob/main/edm4hep_analysis/edm4hep_python.ipynb

# Before running this script be sure to run the commands needed to access the software:
# apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest
# apptainer run --no-home -B /collab/project/snowmass21/data/muonc:/data -B /home/$USER k4toroid.sif
# source /setup.sh

import glob
import ROOT
import matplotlib.pyplot as plt

#from ROOT import edm4hep
#from podio.root_io import Reader
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Open the edm4hep files with ROOT
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/k4reco/electronGun*")
samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/electronGun*")
files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
#slices = ["1000_5000", "0_50", "50_250", "250_1000"]
for s in slices: files[f"electronGun_pT_{s}"] = []
for s in samples:
    sname = s.split("/")[-1]
    #files[sname] = glob.glob(f"{s}/*.root")
    files[sname] = glob.glob(f"{s}/*.slcio")

hists= {}
for s in files:
    hists[s] = {}
    hists[s]["pfo_types"] = ROOT.TH1F(f"{s}_pfo_types", f"{s}_pfo_types", 100, 0, 100)

# Create a reader object to use for the rest of the time
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks_Refitted"])

all_pfos = {}
for s in files:
    print("Working on sample", s)
    i = 0

    all_pfos[s] = []
    # Loop over the files in a sample
    for f in files[s]:
        if max_events > 0 and i >= max_events: break
        reader.open(f)

        # Loop over events in each file
        #for event in reader.get('events'):
        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i%100 == 0: print("\tProcessing event:", i)

            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")
            trks = event.getCollection("SiTracks_Refitted")

            for pfo in pfos:
                all_pfos[s].append(str(pfo.getType()))

            i+=1
        reader.close()

for s in files:
    plt.hist(all_pfos[s])
    plt.title(s)
    plt.savefig(f"plots/electrons/pfo_types_{s}.png")
    plt.clf()
