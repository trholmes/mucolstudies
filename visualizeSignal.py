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
#infiles = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v4/neutronGun_E_50_250/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v4/pionGun_pT_0_50/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/dijet_mjj_10000/*")
infiles = glob.glob("/phchang/data/mumu_ZH/*reco_0.slcio")
#infiles = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v3/photonGun_E_50_250/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/photonGun_E_50_250/*")
#infiles = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun_pT_50_250/*")
#infiles = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl10_v0A/photonGun_500/reco/*")
append = "_phil"
makeTransparent = False

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
#reader.setReadCollectionNames(["EcalBarrelCollectionSel", "EcalEndcapCollectionSel", "HcalBarrelCollectionConed", "HcalEndcapCollectionConed", "MCParticle"])
#reader.setReadCollectionNames(["ECalBarrelCollection", "ECalEndcapCollection", "MCParticle", "HCalBarrelCollection", "HCalEndcapCollection", "VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])

max_events = 10
print_mcp_info = True
draw_mcp_lines = True
min_mcp_E = 0.01
is_jet = True

i_event = 0
for infile in infiles:
    if i_event >= max_events: break
    reader.open(infile)
    print(infile)
    for i, event in enumerate(reader):

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

        pfos = event.getCollection("PandoraPFOs")

        try: ecalb = event.getCollection("EcalBarrelCollectionRec")
        #try: ecalb = event.getCollection("EcalBarrelCollectionConed")
        except: ecalb = []
        try: ecale = event.getCollection("EcalEndcapCollectionRec")
        #try: ecale = event.getCollection("EcalEndcapCollectionConed")
        except: ecale = []
        try: hcalb = event.getCollection("HcalBarrelCollectionRec")
        #try: hcalb = event.getCollection("HcalBarrelCollectionConed")
        except: hcalb = []
        try: hcale = event.getCollection("HcalEndcapCollectionRec")
        #try: hcale = event.getCollection("HcalEndcapCollectionConed")
        except: hcale = []
        mcpCollection = event.getCollection("MCParticle")

        print("Event", i_event)
        print(infile)
        lines = []
        E = 0
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)
            if draw_mcp_lines and mcp_tlv.E() > min_mcp_E:
                lines.append([mcp.getVertex(), mcp.getEndpoint()])
            if print_mcp_info and mcp_tlv.E() > min_mcp_E:
                print(f"\tMCP: pdgid:{mcp.getPDG()} E: {mcp_tlv.E():.2f} GeV eta: {mcp_tlv.Eta():.2f} prodR: {(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2)**.5:.2f} status: {mcp.getGeneratorStatus()} xyz: {mcp.getVertex()[0]:.2f}, {mcp.getVertex()[1]:.2f}, {mcp.getVertex()[2]:.2f}, time: {mcp.getTime():.2f}")
                print(f"\t\tMCP end xyz: {mcp.getEndpoint()[0]:.2f}, {mcp.getEndpoint()[1]:.2f}, {mcp.getEndpoint()[2]:.2f}, E: {(mcp.getMomentumAtEndpoint()[0]**2 + mcp.getMomentumAtEndpoint()[1]**2 + mcp.getMomentumAtEndpoint()[2]**2)**0.5:.2f}")
                try:
                    parents = mcp.getParents()
                    if parents:
                        for parent in parents:
                            print(f"\t\tparent: pdgid:{parent.getPDG()} E: {parent.getEnergy():.2f} GeV prodR: {(parent.getVertex()[0]**2 + parent.getVertex()[1]**2)**.5:.2f} status: {mcp.getGeneratorStatus()} xyz: {parent.getVertex()[0]:.2f}, {parent.getVertex()[1]:.2f}, {parent.getVertex()[2]:.2f}")
                except:
                    pass
                #children = mcp.getDaughters()
                #for child in children:
                #    print(f"\t\tchild: {child.getPDG()}, {child.getEnergy():.2f} GeV")
            if (is_jet and mcp.getGeneratorStatus()==23) or (not is_jet and mcp.getGeneratorStatus()==1):# and isGood(mcp_tlv):
                #print("E", mcp_tlv.E())
                E += mcp_tlv.E()
        #if E < 200: continue


        ecal_total = 0
        hcal_total = 0
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
            ecal_total+=simhit.getEnergy()
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
            ecal_total+=simhit.getEnergy()
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
            hcal_total+=simhit.getEnergy()
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
            hcal_total+=simhit.getEnergy()

        #if (hcal_total+ecal_total)/E > 0.1: continue

        print(f"total ECAL: {ecal_total:.2f}, total HCAL: {hcal_total:.2f}")
        print(f"Calo E / True E: {(hcal_total+ecal_total)/E:.2f}")
        #for pfo in pfos:
        #    if True: #pfo.getEnergy() > 40:
        #        print(f"pfo E: {pfo.getEnergy()}, pdgid: {pfo.getType()}")

        for h in hists:
            #print(h, hists[h].Integral())

            if "_cal" not in h: continue

            can = ROOT.TCanvas()
            can.SetRightMargin(0.18)

            hists[h].Draw("colz")
            if h.startswith("rZ"):
                if draw_mcp_lines:
                    rlines = []
                    for i_line,line in enumerate(lines):
                        rlines.append(ROOT.TLine(line[0][2], (line[0][0]**2+line[0][1]**2)**.5, line[1][2], (line[1][0]**2+line[1][1]**2)**.5))
                        rlines[i_line].SetLineWidth(2)
                        rlines[i_line].SetLineColor(ROOT.kMagenta)
                        rlines[i_line].Draw()
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
            latex.DrawLatex(0.75, 0.85, f"True Energy: {E:.1f} GeV")
            latex.DrawLatex(0.75, 0.80, f"Calo Energy: {hists[h].Integral():.1f} GeV")
            latex.DrawLatex(0.75, 0.75, f"Fraction Reco: {hists[h].Integral()/E:.2f}")

            latex.DrawLatex(0.75, 0.70, f"ECAL Fraction: {safeDivide(ecal_total,(ecal_total+hcal_total)):.2f}")
            latex.DrawLatex(0.75, 0.65, f"File {infile.split('_')[-1].split('.')[0]}, event {i}")

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
