import pyLCIO
import glob
import ctypes
import math

exec(open("./plotHelper.py").read())
ROOT.gStyle.SetPalette(ROOT.kGreyScale)
ROOT.TColor.InvertPalette()

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

#infile = "nuGun_reco_v0A.slcio"
infile = "/data/fmeloni/DataMuC_MuColl10_v0A/v2/digi/nuGun_digi_v0A.slcio"
#infile = "/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl10_v0A/nuGun/digi/nuGun_BIB10_digi.slcio "
append = "_closeup_phislice_gray2"
makeTransparent = False

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["ECalBarrelCollection", "ECalEndcapCollection", "HCalBarrelCollection", "HCalEndcapCollection", "VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection", "MUON"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])

hists = {
        "rZ_cal": ROOT.TH2D("rZ_cal", "rZ_cal", 100, -5000, 5000, 100, 0, 5000),
        #"rZ_hcal": ROOT.TH2D("rZ_hcal", "rZ_hcal", 100, -5000, 5000, 100, 0, 5000),
        "rZ_hcal": ROOT.TH2D("rZ_hcal", "rZ_hcal", 100, -2000, 2000, 100, 2100, 5000),
        "rZ_muon": ROOT.TH2D("rZ_muon", "rZ_muon", 100, -6500, 6500, 100, 0, 8000),
        "xy_muon": ROOT.TH2D("xy_muon", "xy_muon", 100, -2000, 2000, 100, -2000, 2000),
        "rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -2000, 2000, 100, 1800, 2300),
        #"rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -3000, 3000, 100, 0, 2500),
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
    muon = event.getCollection("MUON")

    for simhit in ecalb:
        pos = simhit.getPosition()
        #if math.atan2(pos[1],pos[0]) > 1.1 and math.atan2(pos[1],pos[0]) < 1.2:
        if math.atan2(pos[1],pos[0]) > 2.25 and math.atan2(pos[1],pos[0]) < 2.35:
            hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for simhit in ecale:
        pos = simhit.getPosition()
        #if math.atan2(pos[1],pos[0]) > 1.1 and math.atan2(pos[1],pos[0]) < 1.2:
        if math.atan2(pos[1],pos[0]) > 2.25 and math.atan2(pos[1],pos[0]) < 2.35:
            hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for simhit in hcalb:
        pos = simhit.getPosition()
        if math.atan2(pos[1],pos[0]) > -1.4 and math.atan2(pos[1],pos[0]) < -1:
            hists["rZ_hcal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for simhit in hcale:
        pos = simhit.getPosition()
        if math.atan2(pos[1],pos[0]) > -1.4 and math.atan2(pos[1],pos[0]) < -1:
            hists["rZ_hcal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
    for collection in [vtb, vte, itb, ite, otb, ote]:
        for simhit in collection:
            pos = simhit.getPosition()
            hists["rZ_tracker"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2))
    for hit in muon:
        pos = hit.getPosition()
        #print(pos[0], pos[1], pos[2])
        #print(pos[2], math.sqrt(pos[0]**2+pos[1]**2))
        hists["rZ_muon"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2))
        if pos[2]>0: hists["xy_muon"].Fill(pos[0], pos[1])

    break

for h in hists:
    print(h, hists[h].Integral())

    can = ROOT.TCanvas()
    can.SetRightMargin(0.18)
    hists[h].Draw("colz")
    hists[h].GetXaxis().SetTitle("z [mm]")
    hists[h].GetYaxis().SetTitle("R [mm]")
    if "xy" in h:
        hists[h].GetXaxis().SetTitle("x [mm]")
        hists[h].GetYaxis().SetTitle("y [mm]")

    if makeTransparent:
        can.SetFillColor(0)
        can.SetFillStyle(0)
        can.SetFillColorAlpha(ROOT.kWhite, 0)

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

    hists[h].SetMinimum(0)
    hists[h].SetMaximum(0.06)
    if "cal" in h: hists[h].GetZaxis().SetTitle("E [GeV]")
    else:
        hists[h].GetZaxis().SetTitle("Number of Hits")
        #hists[h].SetMinimum(100)
    #can.SetRealAspectRatio(1)
    can.SetLogz(0)
    #can.SaveAs(f"plots/bib/{h}_linear{append}.png")
    can.SaveAs(f"plots/bib/{h}_linear{append}.pdf")
    #can.SaveAs(f"plots/bib/{h}_linear{append}.eps")
    can.SetLogz()
    #can.SaveAs(f"plots/bib/{h}{append}.png")
    can.SaveAs(f"plots/bib/{h}{append}.pdf")
    #can.SaveAs(f"plots/bib/{h}{append}.eps")
