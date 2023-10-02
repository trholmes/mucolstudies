
print("Loading in plotHelper.py")

variables = {"obj": {}, "event": {}}
variables["obj"]["pt"] =  {"nbins": 20, "xmin": 0,     "xmax": 2000,    "accessor": ".Perp()"}
variables["obj"]["eta"] = {"nbins": 20, "xmin": -3,    "xmax": 3,       "accessor": ".Eta()"}
variables["obj"]["phi"] = {"nbins": 20, "xmin": -3.5,  "xmax": 3.5,     "accessor": ".Phi()"}
variables["event"]["n"] = {"nbins": 20, "xmin": 0,     "xmax": 20,      "accessor": ""}

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
