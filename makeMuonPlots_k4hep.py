import glob
import ROOT

# Set up some options
max_events = 1

# Open the edm4hep files with ROOT
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/edm4hep/*")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/photonGun_1000/reco_k4/*")
t = ROOT.TChain("events")
for f in fnames: t.Add(f)

# Check how many events are in the file
n_events = t.GetEntries()
print("Found %i events."%n_events)

#t.Draw("MCParticle.PDG")
#input("...")

# Can successfully Draw() but running into issues when doing any kind of loop. Not sure of the cause.
# Loop over events
#for e in t: print("!")

'''
for i, event in enumerate(t):
    print(i)
    #if max_events > 0 and i >= max_events: break
    #if i%100 == 0: print("Processing event %i."%i)
'''


