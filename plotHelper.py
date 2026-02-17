import ROOT
import math

ROOT.gStyle.SetOptStat(0)

print("Loading in plotHelper.py")

variables = {"obj": {}, "evt": {}}
variables["obj"]["pt"] =  {"nbins": 30, "xmin": 0,     "xmax": 3000,    "accessor": ".Perp()",  "label": "p_{T} [GeV]"}
variables["obj"]["eta"] = {"nbins": 20, "xmin": -3,    "xmax": 3,       "accessor": ".Eta()",   "label": "#eta"}
variables["obj"]["phi"] = {"nbins": 20, "xmin": -3.5,  "xmax": 3.5,     "accessor": ".Phi()",   "label": "#phi"}
variables["evt"]["n"] =   {"nbins": 20, "xmin": 0,     "xmax": 20,      "accessor": "",         "label": "number per event"}

colors =        [ROOT.TColor.GetColor("#FFA900"),   # Sunny Orange
                ROOT.TColor.GetColor("#3629AC"),    # Dark Blue
                ROOT.TColor.GetColor("#2FC494"),    # Seafoam
                ROOT.TColor.GetColor("#F65866"),    # Pink
                ROOT.TColor.GetColor("#0E81C4"),    # Light Blue
                ]

def getTLV(obj):
    obj_p = obj.getMomentum()
    obj_e = obj.getEnergy()
    obj_tlv = ROOT.TLorentzVector()
    #obj_tlv.SetPxPyPzE(obj_p.x, obj_p.y, obj_p.z, obj_e)
    obj_tlv.SetPxPyPzE(obj_p[0], obj_p[1], obj_p[2], obj_e)
    return obj_tlv

def safeDivide(a, b):
    try: return a/b
    except: return -1

# Get pT from track object
# Taken from here
# https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf
def getPt(trk, b_field = 5):
    return 3e-4*abs(b_field/trk.getOmega())
def getP(trk, b_field = 5):
    return getPt(trk)*math.sqrt(1+trk.getTanLambda()**2)
def getTrackTLV(trk, m = .106, b_field = 5):
    pt = getPt(trk, b_field)
    p = getP(trk, b_field)
    px = pt*math.cos(trk.getPhi())
    py = pt*math.sin(trk.getPhi())
    pz = pt*trk.getTanLambda()
    E = math.sqrt(p**2 + m**2) # Assume this is a muon
    trk_tlv = ROOT.TLorentzVector()
    trk_tlv.SetPxPyPzE(px, py, pz, E)
    return trk_tlv

def fillObjHists(hists, objtype, obj):
    for var in variables["obj"]:
        hists[objtype+"_"+var].Fill(eval("obj"+variables["obj"][var]["accessor"]))

def getMaximum(hists, isLog=False):
    max_factor = 1.1
    if isLog: max_factor = 10
    maximum = max(hists[i].GetMaximum() for i in range(len(hists)))
    return maximum*max_factor

def colorHists(hists):
    i=0
    for h in hists:
        h.SetLineColor(colors[i])
        h.SetMarkerColor(colors[i])
        i += 1
        i %= len(colors)
    return


# From a map of histograms (key as label name), plot them all on a canvas and save
def plotHistograms(h_map, save_name, xlabel="", ylabel="", interactive=False, logy=False, atltext="", legx=.66):

    can = ROOT.TCanvas("can", "can")
    h_keys = list(h_map.keys())
    h_values = list(h_map.values())

    if len(h_keys)<1:
        print("No histograms found in h_map. Drawing blank canvas.")
        return

    # Get maxes/mins of all hists
    maxy = 1.5*getMaximum(h_values)
    miny = 0
    if logy:
        maxy*=1e4
        miny = 1e-1

    # Draw histograms
    colorHists(h_values)
    h_map[h_keys[0]].GetXaxis().SetTitle(xlabel)
    h_map[h_keys[0]].GetYaxis().SetTitle(ylabel)
    h_values[0].SetMinimum(miny)
    h_values[0].SetMaximum(maxy)
    h_map[h_keys[0]].Draw("hist")
    if logy: ROOT.gPad.SetLogy(1)
    for k in h_keys[1:]:
        h_map[k].Draw("hist same")

    leg = ROOT.TLegend(legx, .64, .9, .88)
    for k in h_keys:
        leg.AddEntry(h_map[k], k, "l")
    leg.Draw()
    if interactive: raw_input("...")

    if atltext != "":
        ROOT.ATLASLabel(0.2,0.85, atltext[0])
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.04)
        for i, t in enumerate(atltext[1:]):
            text.DrawLatex(0.2,0.85-0.07*(i+1), t)

    can.SaveAs(save_name)
    return

def plotEfficiencies(eff_map, save_name, xlabel="", ylabel="", xrange=""):

    can = ROOT.TCanvas("can", "can")
    if len(eff_map)<1:
        return

    colorHists(eff_map.values())

    for i, k in enumerate(eff_map):
        if i==0:
            if not xrange=="":
                background = ROOT.TH1D("bkg", "bkg", 1, x_range[0], x_range[1])
                background.SetLineWidth(0)
                background.SetTitle(";%s;%s"%(xlabel,ylabel))
                background.SetMinimum(0)
                background.SetMaximum(1)
                background.Draw()
                eff_map[k].Draw("pe same")
            else:
                eff_map[k].Draw("ape")
            eff_map[k].SetTitle(";%s;%s"%(xlabel,ylabel))
            ROOT.gPad.Update()

            eff_map[k].GetPaintedGraph().GetXaxis().SetTitle(xlabel)
            eff_map[k].GetPaintedGraph().GetYaxis().SetTitle(ylabel)
            eff_map[k].GetPaintedGraph().SetMinimum(0)
            eff_map[k].GetPaintedGraph().SetMaximum(1)
            if not xrange=="":
                print("Updating axis range", xrange)
                #eff_map[k].GetPaintedGraph().GetXaxis().SetRange(x_range[0], x_range[1])
                #eff_map[k].GetPaintedGraph().GetXaxis().SetMin(x_range[0])
                #eff_map[k].GetPaintedGraph().GetXaxis().SetMax(x_range[1])
                eff_map[k].Draw("pe same")
            else:
                eff_map[k].Draw("ape")
            ROOT.gPad.Update()

        else: eff_map[k].Draw("pe same")

        ROOT.gPad.Update()

    leg = ROOT.TLegend(.66, .24, .9, .38)
    for k in eff_map:
        leg.AddEntry(eff_map[k], k, "p")
    leg.Draw()

    can.SaveAs(save_name)
    return
