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

# Set up some options
max_events = 1

# Open the edm4hep files with ROOT
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/k4reco/muonGun_pT_250_1000/*.root")

reader = Reader(fnames[0])
for i, event in enumerate(reader.get('events')):
    print(i)
    tracks = event.get('SiTracks_Refitted')
    for trk in tracks:
        print(trk.getChi2())

