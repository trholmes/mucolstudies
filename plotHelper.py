
print("Loading in plotHelper.py")

variables = {"obj": {}, "evt": {}}
variables["obj"]["pt"] =  {"nbins": 20, "xmin": 0,     "xmax": 2000,    "accessor": ".Perp()",  "label": "p_{T} [GeV]"}
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
    obj_tlv.SetPxPyPzE(obj_p.x, obj_p.y, obj_p.z, obj_e)
    return obj_tlv

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
def plotHistograms(h_map, save_name, xlabel="", ylabel="", interactive=False, logy=False, atltext=""):

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

    leg = ROOT.TLegend(.66, .64, .9, .88)
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

