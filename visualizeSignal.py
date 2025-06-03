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
infiles = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v2/reco/neutronGun_E_50_250/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_50_250/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun_pT_50_250/*")
#infiles = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl10_v0A/photonGun_500/reco/*")
append = "_neutron"
makeTransparent = False

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["EcalBarrelCollectionSel", "EcalEndcapCollectionSel", "HcalBarrelCollectionConed", "HcalEndcapCollectionConed", "MCParticle"])
#reader.setReadCollectionNames(["ECalBarrelCollection", "ECalEndcapCollection", "MCParticle", "HCalBarrelCollection", "HCalEndcapCollection", "VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])

max_events = 2

i_event = 0
for infile in infiles:
    if i_event >= max_events: break
    reader.open(infile)
    for event in reader:

        if i_event >= max_events: break

        hists = {
        "rZ_cal": ROOT.TH2D("rZ_cal", "rZ_cal", 100, -5000, 5000, 100, 0, 5000),
        "xy_cal": ROOT.TH2D("xy_cal", "xy_cal", 100, -5000, 5000, 100, -5000, 5000),
        "phiZ_cal": ROOT.TH2D("phiZ_cal", "phiZ_cal", 100, -5000, 5000, 100, -3.2, 3.2),
        "phiEta_cal": ROOT.TH2D("phiEta_cal", "phiEta_cal", 100, -2.5, 2.5, 100, -3.2, 3.2),
        "rZ_hcal": ROOT.TH2D("rZ_hcal", "rZ_hcal", 100, -5000, 5000, 100, 2100, 5000),
        "xy_hcal": ROOT.TH2D("xy_hcal", "xy_hcal", 100, -5000, 5000, 100, -5000, 5000),
        "phiZ_hcal": ROOT.TH2D("phiZ_hcal", "phiZ_hcal", 100, -5000, 5000, 100, -3.2, 3.2),
        "rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -2000, 2000, 100, 1800, 2000),
        "xy_ecal": ROOT.TH2D("xy_ecal", "xy_ecal", 100, -2000, 2000, 100, -2000, 2000),
        "phiZ_ecal": ROOT.TH2D("phiZ_ecal", "phiZ_ecal", 100, -2000, 2000, 100, -3.2, 3.2),
        #"xy_ecal": ROOT.TH2D("xy_ecal", "xy_ecal", 100, -2000, 2000, 100, -3.1415, 3.1415),
        #"rZ_ecal": ROOT.TH2D("rZ_ecal", "rZ_ecal", 100, -3000, 3000, 100, 0, 2500),
        #"rZ_tracker": ROOT.TH2D("rZ_tracker", "rZ_tracker", 500, -2300, 2300, 500, 0, 1700),
        }

        try: ecalb = event.getCollection("EcalBarrelCollectionSel")
        except: ecalb = []
        try: ecale = event.getCollection("EcalEndcapCollectionSel")
        except: ecale = []
        try: hcalb = event.getCollection("HcalBarrelCollectionConed")
        except: hcalb = []
        try: hcale = event.getCollection("HcalEndcapCollectionConed")
        except: hcale = []
        mcpCollection = event.getCollection("MCParticle")

        print("Event", i_event)
        E = 0
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)
            if mcp_tlv.E() > 0.01:
                print(f"\tMCP: pdgid:{mcp.getPDG()} E: {mcp_tlv.E():.2f} GeV prodR: {(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2)**.5:.2f} status: {mcp.getGeneratorStatus()}")
                try:
                    parent = mcp.getParents()[0]
                    print(f"\t\tparent: {parent.getPDG()}, {parent.getEnergy():.2f} GeV")
                except:
                    children = mcp.getDaughters()
                    for child in children:
                        print(f"\t\tchild: {child.getPDG()}, {child.getEnergy():.2f} GeV")
            if mcp.getGeneratorStatus()==1:# and isGood(mcp_tlv):
                #print("E", mcp_tlv.E())
                E = mcp_tlv.E()
        #if E < 200: continue

        for simhit in ecalb:
            pos = simhit.getPosition()
            vpos = ROOT.TVector3(pos[0], pos[1], pos[2])
            hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_ecal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_ecal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_cal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_cal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["phiEta_cal"].Fill(vpos.Eta(), vpos.Phi(), simhit.getEnergy())
        for simhit in ecale:
            pos = simhit.getPosition()
            vpos = ROOT.TVector3(pos[0], pos[1], pos[2])
            hists["rZ_ecal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_ecal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_ecal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_cal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_cal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["phiEta_cal"].Fill(vpos.Eta(), vpos.Phi(), simhit.getEnergy())
        for simhit in hcalb:
            pos = simhit.getPosition()
            vpos = ROOT.TVector3(pos[0], pos[1], pos[2])
            hists["rZ_hcal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_hcal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_hcal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_cal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_cal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["phiEta_cal"].Fill(vpos.Eta(), vpos.Phi(), simhit.getEnergy())
        for simhit in hcale:
            pos = simhit.getPosition()
            vpos = ROOT.TVector3(pos[0], pos[1], pos[2])
            hists["rZ_hcal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_hcal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_hcal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["rZ_cal"].Fill(pos[2], math.sqrt(pos[0]**2+pos[1]**2), simhit.getEnergy())
            hists["xy_cal"].Fill(pos[0], pos[1], simhit.getEnergy())
            hists["phiZ_cal"].Fill(pos[2], math.atan2(pos[1],pos[0]), simhit.getEnergy())
            hists["phiEta_cal"].Fill(vpos.Eta(), vpos.Phi(), simhit.getEnergy())

        for h in hists:
            print(h, hists[h].Integral())

            can = ROOT.TCanvas()
            can.SetRightMargin(0.18)

            hists[h].Draw("colz")
            if h.startswith("rZ"):
                hists[h].GetXaxis().SetTitle("z [mm]")
                hists[h].GetYaxis().SetTitle("R [mm]")
            elif h.startswith("xy"):
                hists[h].GetXaxis().SetTitle("x [mm]")
                hists[h].GetYaxis().SetTitle("y [mm]")
            elif h.startswith("phiZ"):
                hists[h].GetXaxis().SetTitle("z [mm]")
                hists[h].GetYaxis().SetTitle("#phi")
            elif h.startswith("phiEta"):
                hists[h].GetXaxis().SetTitle("#eta")
                hists[h].GetYaxis().SetTitle("#phi")

            hists[h].GetZaxis().SetTitle("E [GeV]")

            latex = ROOT.TLatex()
            latex.SetNDC() # Use normalized device coordinates
            latex.SetTextAlign(31) # Align right (3), top (1)
            latex.SetTextSize(0.04) # Adjust size as needed
            latex.DrawLatex(0.75, 0.85, f"Neutron Energy: {E:.1f} GeV")
            latex.DrawLatex(0.75, 0.80, f"Calo Energy: {hists[h].Integral():.1f} GeV")

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

            #hists[h].SetMinimum(0)
            #hists[h].SetMaximum(0.15)
            if (h.startswith("phiEta") or h.startswith("rZ") or h.startswith("xy")) and "_cal" in h:
                can.SetRealAspectRatio(1)
            can.SetLogz(0)
            can.SaveAs(f"plots/bib/tests/{i_event}_{h}_linear{append}.pdf")
            #can.SetLogz()
            #can.SaveAs(f"plots/bib/tests/{h}{append}_{i_event}.pdf")
        i_event += 1
