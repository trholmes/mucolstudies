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
reader.setReadCollectionNames(["ECalBarrelCollection", "ECalEndcapCollection", "HCalBarrelCollection", "HCalEndcapCollection", "VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])

hists = {
        "rZ_cal": ROOT.TH2D("rZ_cal", "rZ_cal", 100, -5000, 5000, 100, 0, 5000),
        "rZ_hcal": ROOT.TH2D("rZ_hcal", "rZ_hcal", 100, -5000, 5000, 100, 0, 5000),
        "rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -3000, 3000, 100, 0, 2500),
        "rZ_tracker": ROOT.TH2D("rZ_tracker", "rZ_tracker", 500, -2300, 2300, 500, 0, 1700),
        }

reader.open(infile)
for event in reader:

    ecalb = event.getCollection("ECalBarrelCollection")
    ecale = event.getCollection("ECalEndcapCollection")
    hcalb = event.getCollection("HCalBarrelCollection")
    hcale = event.getCollection("HCalEndcapCollection")
    vtb = event.getCollection("VertexBarrelCollection")
    vte = event.getCollection("VertexEndcapCollection")
    itb = event.getCollection("InnerTrackerBarrelCollection")
    ite = event.getCollection("InnerTrackerEndcapCollection")
    otb = event.getCollection("OuterTrackerBarrelCollection")
    ote = event.getCollection("OuterTrackerEndcapCollection")

    for simhit in ecalb:
        pos = simhit.getPosition()
        hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
        hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for simhit in ecale:
        pos = simhit.getPosition()
        hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
        hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for simhit in hcalb:
        pos = simhit.getPosition()
        hists["rZ_hcal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
        hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for simhit in hcale:
        pos = simhit.getPosition()
        hists["rZ_hcal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
        hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for collection in [vtb, vte, itb, ite, otb, ote]:
        for simhit in collection:
            pos = simhit.getPosition()
            hists["rZ_tracker"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2))

    break

#ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)
for h in hists:
    print(h, hists[h].Integral())

    can = ROOT.TCanvas()
    can.SetFillColor(0)
    can.SetFillStyle(0)
    can.SetFillColorAlpha(ROOT.kWhite, 0)
    can.SetRightMargin(0.18)

    #p = ROOT.TPad()
    #p.SetFillColor(0)
    #p.SetFillStyle(0)
    #p.SetFrameFillColor(0)
    #p.SetFrameFillStyle(0)

    hists[h].Draw("colz")
    hists[h].GetXaxis().SetTitle("z [mm]")
    hists[h].GetYaxis().SetTitle("R [mm]")
    hists[h].GetXaxis().SetAxisColor(ROOT.kWhite)
    hists[h].GetYaxis().SetAxisColor(ROOT.kWhite)
    hists[h].GetZaxis().SetAxisColor(ROOT.kWhite)
    hists[h].GetXaxis().SetLabelColor(ROOT.kWhite)
    hists[h].GetYaxis().SetLabelColor(ROOT.kWhite)
    hists[h].GetZaxis().SetLabelColor(ROOT.kWhite)
    hists[h].GetXaxis().SetTitleColor(ROOT.kWhite)
    hists[h].GetYaxis().SetTitleColor(ROOT.kWhite)
    hists[h].GetZaxis().SetTitleColor(ROOT.kWhite)
    hists[h].SetFillColorAlpha(ROOT.kWhite, 0)

    if "cal" in h: hists[h].GetZaxis().SetTitle("E [GeV]")
    else:
        hists[h].GetZaxis().SetTitle("Number of Hits")
        hists[h].SetMinimum(100)
    can.SetRealAspectRatio(1)
    can.SetLogz(0)
    can.SaveAs(f"plots/bib/{h}_linear.png")
    can.SaveAs(f"plots/bib/{h}_linear.pdf")
    can.SaveAs(f"plots/bib/{h}_linear.eps")
    can.SetLogz()
    can.SaveAs(f"plots/bib/{h}.png")
    can.SaveAs(f"plots/bib/{h}.pdf")
    can.SaveAs(f"plots/bib/{h}.eps")
