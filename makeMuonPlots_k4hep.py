import glob
import ROOT

# Open the edm4hep files with ROOT
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/edm4hep/*")
t = ROOT.TChain("events")
for f in fnames: t.Add(f)

# Check how many events are in the file
n_events = t.GetEntries()
print("Found %i events."%n_events)





