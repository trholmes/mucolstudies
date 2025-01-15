import pyLCIO
import glob
import ctypes
import math

exec(open("./plotHelper.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

#infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/neutronGun_E_50_250/*")
infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_50_250/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun_pT_50_250/*")
#infiles = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl10_v0A/photonGun_500/reco/*")
append = "_singlephoton"
makeTransparent = False

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["ECalBarrelCollection", "ECalEndcapCollection", "MCParticle", "HCalBarrelCollection", "HCalEndcapCollection", "VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])

hists = {
        #"rZ_cal": ROOT.TH2D("rZ_cal", "rZ_cal", 100, -5000, 5000, 100, 0, 5000),
        #"rZ_hcal": ROOT.TH2D("rZ_hcal", "rZ_hcal", 100, -2000, 2000, 100, 2100, 5000),
        #"phiZ_hcal": ROOT.TH2D("phiZ_hcal", "phiZ_hcal", 100, -2000, 2000, 100, -.8, -1.5),
        "rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -2000, 2000, 100, 1800, 2300),
        "phiZ_ecal": ROOT.TH2D("phiZ_ecal", "phiZ_ecal", 100, -2000, 2000, 100, 2, 2.5),
        #"phiZ_ecal": ROOT.TH2D("phiZ_ecal", "phiZ_ecal", 100, -2000, 2000, 100, -3.1415, 3.1415),
        #"rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -3000, 3000, 100, 0, 2500),
        #"rZ_tracker": ROOT.TH2D("rZ_tracker", "rZ_tracker", 500, -2300, 2300, 500, 0, 1700),
        }

for infile in infiles:
    reader.open(infile)
    for event in reader:

        ecalb = event.getCollection("ECalBarrelCollection")
        ecale = event.getCollection("ECalEndcapCollection")
        hcalb = event.getCollection("HCalBarrelCollection")
        hcale = event.getCollection("HCalEndcapCollection")
        mcpCollection = event.getCollection("MCParticle")

        E = 0
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)

            if mcp.getGeneratorStatus()==1:# and isGood(mcp_tlv):
                print("E", mcp_tlv.E())
                E = mcp_tlv.E()

        #if E < 200: continue

        for simhit in ecalb:
            pos = simhit.getPosition()
            if math.atan2(pos[1],pos[0]) > 2.25 and math.atan2(pos[1],pos[0]) < 2.35:
                #hists["rZ_hcal"].Fill(pos[2], 85+math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
                hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
                hists["phiZ_ecal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
                #hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
        for simhit in ecale:
            pos = simhit.getPosition()
            if math.atan2(pos[1],pos[0]) > 2.25 and math.atan2(pos[1],pos[0]) < 2.35:
                #hists["rZ_hcal"].Fill(pos[2], 85+math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
                hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
                hists["phiZ_ecal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
                #hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
        break
    break

for h in hists:
    print(h, hists[h].Integral())

    can = ROOT.TCanvas()
    can.SetRightMargin(0.18)
    hists[h].Draw("colz")
    hists[h].GetXaxis().SetTitle("z [mm]")
    hists[h].GetYaxis().SetTitle("R [mm]")
    hists[h].GetZaxis().SetTitle("E [GeV]")

    if makeTransparent:
        can.SetFillColor(0)
        can.SetFillStyle(0)
        can.SetFillColorAlpha(ROOT.kWhite, 0)
        can.SetRightMargin(0.18)

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
    hists[h].SetMaximum(0.15)
    #can.SetRealAspectRatio(1)
    can.SetLogz(0)
    #can.SaveAs(f"plots/bib/{h}_linear{append}.png")
    can.SaveAs(f"plots/bib/{h}_linear{append}.pdf")
    #can.SaveAs(f"plots/bib/{h}_linear{append}.eps")
    can.SetLogz()
    #can.SaveAs(f"plots/bib/{h}{append}.png")
    can.SaveAs(f"plots/bib/{h}{append}.pdf")
    #can.SaveAs(f"plots/bib/{h}{append}.eps")
