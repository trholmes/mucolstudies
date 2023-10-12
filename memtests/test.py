import glob
#import gc
max_events = 1000

'''
# pylcio simplified setup
import pyLCIO
files = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun_pT_250_1000/*.slcio")
i = 0

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks_Refitted"])
reader.setReadCollectionNames(["PandoraPFOs", "SiTracks_Refitted"])

for f in files:
    if i >= max_events: break
    reader.open(f)
    for event in reader:
        if i >= max_events: break
        if i%100 == 0: print("Processing event ", i)
        i += 1
    reader.close()
    #del(reader)
    #gc.collect()
'''

# key4hep simplified setup
#from podio.root_io import Reader
from ROOT import edm4hep
from ROOT import podio

files = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/k4reco/electronGun_pT_250_1000/*.root")
i = 0

for f in files:
    if i >= max_events: break
    #reader = Reader(f)
    reader = podio.ROOTFrameReader()
    #reader.collections = ["MCParticle", "PandoraPFOs", "SiTracks_Refitted"]
    reader.openFile(f)

    #for event in reader.get('events'):
    for event in podio.Frame(reader.readEntry('events', 0)):
        if i >= max_events: break
        if i%100 == 0: print("Processing event ", i)
        i += 1

print("Done.")
