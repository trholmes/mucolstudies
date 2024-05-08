import pyLCIO
import glob
import ctypes
import math

exec(open("./plotHelper.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

infile = "nuGun_reco_v0A.slcio"

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["ECalBarrelCollection", "ECalEndcapCollection"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])

hists = {
        "rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -1000, 1000, 100, 0, 1000)
        }

reader.open(infile)
for event in reader:

    ecalb = event.getCollection("ECalBarrelCollection")
    ecale = event.getCollection("ECalEndcapCollection")

    for simhit in ecalb:
        hists["rZ_ecal"].Fill(simhit.z(), math.sqrt(simhit.x()**2+simhit.y()**2))
    for simhit in ecale:
        hists["rZ_ecal"].Fill(simhit.z(), math.sqrt(simhit.x()**2+simhit.y()**2))

    break

hists["rZ_ecal"].Draw()
