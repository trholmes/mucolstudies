import pyLCIO
import glob
import ctypes

#exec(open("helpers.py").read())
exec(open("./plotHelper.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1
obj_type = "ph"
magnetic_field = 5.00
max_E = 1000

# Set up things for each object
settings = {
        "fnames": { "ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/photonGun*",
                    #"ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco_highrange/photonGun*",
                    "mu": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun*",
                    "el": "/data/fmeloni/DataMuC_MuColl10_v0A/reco_highrange/electronGun*"},
                    #"el": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*"},
        "labelname": {  "ph": "Photon",
                        "mu": "Muon",
                        "el": "Electron"},
        "plotdir":{ "ph": "photons",
                    "mu": "muons",
                    "el": "electrons"},
        "pdgid":  { "ph": 22,
                    "mu": 13,
                    "el": 11},
        "mass":   { "ph": 0,
                    "mu": 0.106,
                    "el": 0.000511}
}
print("Running on", settings["labelname"][obj_type])

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
samples = glob.glob(settings["fnames"][obj_type])
fnames = []
for s in samples:
    fnames += glob.glob(f"{s}/*.slcio")
print("Found %i files."%len(fnames))

# Get pT from track object
# Taken from here
# https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf
def getPt(trk):
    return 3e-4*abs(magnetic_field/trk.getOmega())
def getP(trk):
    return getPt(trk)*math.sqrt(1+trk.getTanLambda()**2)
def getTrackTLV(trk):
    pt = getPt(trk)
    p = getP(trk)
    px = pt*math.cos(trk.getPhi())
    py = pt*math.sin(trk.getPhi())
    pz = pt*trk.getTanLambda()
    E = math.sqrt(p**2 + settings["mass"][obj_type]**2)
    trk_tlv = ROOT.TLorentzVector()
    trk_tlv.SetPxPyPzE(px, py, pz, E)
    return trk_tlv

# Define good particle
def isGood(tlv):
    if abs(tlv.Eta()) < 2:
        return True
    return False

# Perform matching between two TLVs
def isMatched(tlv1, tlv2, req_pt = True):
    if tlv1.DeltaR(tlv2) > 0.1: return False
    if req_pt:
        drelpt = abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp()
        if drelpt > 0.2: return False # Simple 20% matching
        #if drelpt > 0.1*tlv2.Perp()/100: return False # Require 10% at 100, 20% at 200, ...
    return True

# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
# Set up histograms
# This is an algorithmic way of making a bunch of histograms and storing them in a dictionary
variables = {}
variables["pt"] =  {"nbins": 30, "xmin": 0, "xmax": max_E,   "title": "p_{T} [GeV]"}
variables["eta"] = {"nbins": 30, "xmin": -3, "xmax": 3,     "title": "#eta"}
variables["phi"] = {"nbins": 30, "xmin": -3.5, "xmax": 3.5, "title": "#phi"}
variables["n"] =   {"nbins": 20, "xmin": 0, "xmax": 20,     "title": "n"}
hists = {}

objects = {}
objects["trk"] = "Track"
objects["pfo"] = "Reconstructed"
objects["trk_ob"] = f"Tracks Matched to {settings['labelname'][obj_type]}"
objects["mcp_ob"] = f"True {settings['labelname'][obj_type]}"
objects["pfo_ob"] = f"Reconstructed {settings['labelname'][obj_type]}"
objects["mcp_ob_trkmatch"] = "Tracking Efficiency"
objects["mcp_ob_pfomatch"] = "Reconstruction Efficiency"
for obj in objects:
    for var in variables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, objects[obj], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# Making a separate set of binning conventions for plots showing resolutions
# these plots will all be filled with the difference between a pfo and a mcp object value
dvariables = {}
dvariables["dpt"] =     {"nbins": 100, "xmin": -500, "xmax": 500,       "title": "p_{T}^{meas.} - p_{T}^{true.} [GeV]"}
dvariables["drelpt"] =  {"nbins": 100, "xmin": -0.5, "xmax": 0.5,       "title": "(p_{T}^{meas.} - p_{T}^{true.})/p_{T}^{true.}"}
dvariables["drelE"] =  {"nbins": 100, "xmin": -0.5, "xmax": 0.5,       "title": "(E^{meas.} - E^{true.})/E^{true.}"}
dvariables["deta"] =    {"nbins": 100, "xmin": -0.001, "xmax": 0.001,   "title": "#eta^{meas.} - #eta^{true}"}
dvariables["dphi"] =    {"nbins": 100, "xmin": -0.001, "xmax": 0.001,   "title": "#phi^{meas.} - #phi^{true}"}
for obj in ["d_ob", "d_trk"]:
    for var in dvariables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, obj+"_"+var, dvariables[var]["nbins"], dvariables[var]["xmin"], dvariables[var]["xmax"])

hists2d = {}
for obj in ["trk_ob", "pfo_ob"]:
    for var in variables:
        hists2d[obj+"_v_mcp_ob_"+var] = ROOT.TH2F(obj+"_v_mcp_ob_"+var, obj+"_v_mcp_ob_"+var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])


# Finally making one 2D histogram non-algorithmically; this is what I'll use for a
# pT resolution vs. pT plot.
h_2d_pforelpt = ROOT.TH2F("h_2d_pforelpt", "h_2d_pforelpt", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelpt = ROOT.TH2F("h_2d_trkrelpt", "h_2d_trkrelpt", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelpt1p1 = ROOT.TH2F("h_2d_pforelpt1p1", "h_2d_pforelpt1p1", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelpt1p1 = ROOT.TH2F("h_2d_trkrelpt1p1", "h_2d_trkrelpt1p1", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelpt1p2 = ROOT.TH2F("h_2d_pforelpt1p2", "h_2d_pforelpt1p2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelpt1p2 = ROOT.TH2F("h_2d_trkrelpt1p2", "h_2d_trkrelpt1p2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelpt2 = ROOT.TH2F("h_2d_pforelpt2", "h_2d_pforelpt2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelpt2 = ROOT.TH2F("h_2d_trkrelpt2", "h_2d_trkrelpt2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelpt_eta = ROOT.TH2F("h_2d_pforelpt_eta", "h_2d_pforelpt_eta", 30, -3, 3, 500, -0.5, 0.5)
h_2d_trkrelpt_eta = ROOT.TH2F("h_2d_trkrelpt_eta", "h_2d_trkrelpt_eta", 30, -3, 3, 500, -0.5, 0.5)
h_2d_pforelE = ROOT.TH2F("h_2d_pforelE", "h_2d_pforelE", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelE_vmeas = ROOT.TH2F("h_2d_pforelE_vmeas", "h_2d_pforelE_vmeas", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelE = ROOT.TH2F("h_2d_trkrelE", "h_2d_trkrelE", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelE1p1 = ROOT.TH2F("h_2d_pforelE1p1", "h_2d_pforelE1p1", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelE1p1 = ROOT.TH2F("h_2d_trkrelE1p1", "h_2d_trkrelE1p1", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelE1p2 = ROOT.TH2F("h_2d_pforelE1p2", "h_2d_pforelE1p2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelE1p2 = ROOT.TH2F("h_2d_trkrelE1p2", "h_2d_trkrelE1p2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelE2 = ROOT.TH2F("h_2d_pforelE2", "h_2d_pforelE2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_trkrelE2 = ROOT.TH2F("h_2d_trkrelE2", "h_2d_trkrelE2", 30, 0, max_E, 500, -0.5, 0.5)
h_2d_pforelE_eta = ROOT.TH2F("h_2d_pforelE_eta", "h_2d_pforelE_eta", 30, -3, 3, 500, -0.5, 0.5)
h_2d_trkrelE_eta = ROOT.TH2F("h_2d_trkrelE_eta", "h_2d_trkrelE_eta", 30, -3, 3, 500, -0.5, 0.5)
h_2d_pfoSFs = ROOT.TH2F("h_2d_pfoSFs", "h_2d_pfoSFs", 60, 0, max_E, 1000, 0.5, 3)

# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks_Refitted"])
i = 0
for f in fnames:
    reader.open(f)
    if max_events > 0 and i >= max_events: break

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")
        trkCollection = event.getCollection("SiTracks_Refitted")

        # Make counter variables
        n_mcp_ob = 0
        n_pfo_ob = 0
        has_mcp_ob = False
        has_pfo_ob = False
        has_trk_ob = False
        my_pfo_ob = 0
        my_trk_ob = 0
        my_mcp_ob = 0

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)

            if abs(mcp.getPDG())==settings['pdgid'][obj_type] and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_ob = True
                n_mcp_ob += 1
                my_mcp_ob = mcp_tlv

        # Loop over the reconstructed objects and fill histograms
        # If there are multiple, it'll keep the one with the higher pT
        for pfo in pfoCollection:
            pfo_tlv = getTLV(pfo)

            if abs(pfo.getType())==settings['pdgid'][obj_type]:
                if has_mcp_ob: # and isMatched(pfo_tlv, my_mcp_ob, req_pt = True):
                    n_pfo_ob += 1
                    has_pfo_ob = True
                    if n_pfo_ob == 1:
                        my_pfo_ob = pfo_tlv
                    elif n_pfo_ob > 1 and pfo_tlv.Perp() > my_pfo_ob.Perp():
                        my_pfo_ob = pfo_tlv

        n_matched_tracks = 0
        # Loop over track collection and save a matched track
        # If there are multiple, it'll keep the one with the higher pT
        for trk in trkCollection:
            trk_tlv = getTrackTLV(trk)
            if has_mcp_ob: # and isMatched(trk_tlv, my_mcp_ob, req_pt = True):
                has_trk_ob = True
                n_matched_tracks += 1
                if n_matched_tracks == 1:
                    my_trk_ob = trk_tlv
                elif n_matched_tracks > 1 and trk_tlv.Perp() > my_trk_ob.Perp():
                    my_trk_ob = trk_tlv

        if n_matched_tracks > 1:
            print("Found multiple matched tracks:", n_matched_tracks)
            for j, trk in enumerate(trkCollection):
                trk.tlv = getTrackTLV(trk)
                print(f"Track {j}: pT {trk_tlv.Perp()}")

        if n_pfo_ob > 1: print("Found multiple PFOs:", n_pfo_ob)
        hists["trk_n"].Fill(len(trkCollection))
        hists["trk_ob_n"].Fill(n_matched_tracks)

        # Only make plots for events with isGood mcps
        if has_mcp_ob:
            #print("filling mcp_ob_pt with", my_mcp_ob.Perp())
            hists["mcp_ob_pt"].Fill(my_mcp_ob.Perp())
            hists["mcp_ob_eta"].Fill(my_mcp_ob.Eta())
            hists["mcp_ob_phi"].Fill(my_mcp_ob.Phi())
            if has_pfo_ob:
                hists["pfo_ob_pt"].Fill(my_pfo_ob.Perp())
                hists["pfo_ob_eta"].Fill(my_pfo_ob.Eta())
                hists["pfo_ob_phi"].Fill(my_pfo_ob.Phi())
                hists["mcp_ob_pfomatch_pt"].Fill(my_mcp_ob.Perp())
                hists["mcp_ob_pfomatch_eta"].Fill(my_mcp_ob.Eta())
                hists["mcp_ob_pfomatch_phi"].Fill(my_mcp_ob.Phi())
                hists2d["pfo_ob_v_mcp_ob_pt"].Fill(my_mcp_ob.Perp(), my_pfo_ob.Perp())
                hists2d["pfo_ob_v_mcp_ob_eta"].Fill(my_mcp_ob.Eta(), my_pfo_ob.Eta())
                hists2d["pfo_ob_v_mcp_ob_phi"].Fill(my_mcp_ob.Phi(), my_pfo_ob.Phi())
                hists["d_ob_dpt"].Fill((my_pfo_ob.Perp()-my_mcp_ob.Perp()))
                hists["d_ob_drelpt"].Fill((my_pfo_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                hists["d_ob_drelE"].Fill((my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                hists["d_ob_dphi"].Fill((my_pfo_ob.Phi()-my_mcp_ob.Phi()))
                hists["d_ob_deta"].Fill((my_pfo_ob.Eta()-my_mcp_ob.Eta()))
                h_2d_pforelpt.Fill(my_mcp_ob.Perp(), (my_pfo_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                h_2d_pforelE.Fill(my_mcp_ob.E(), (my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                h_2d_pforelE_vmeas.Fill(my_pfo_ob.E(), (my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                if abs(my_mcp_ob.Eta()) < 1.1:
                    h_2d_pforelpt1p1.Fill(my_mcp_ob.Perp(), (my_pfo_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                    h_2d_pforelE1p1.Fill(my_mcp_ob.E(), (my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                elif abs(my_mcp_ob.Eta()) < 1.2:
                    h_2d_pforelpt1p2.Fill(my_mcp_ob.Perp(), (my_pfo_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                    h_2d_pforelE1p2.Fill(my_mcp_ob.E(), (my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                else:
                    h_2d_pforelpt2.Fill(my_mcp_ob.Perp(), (my_pfo_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                    h_2d_pforelE2.Fill(my_mcp_ob.E(), (my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                h_2d_pforelpt_eta.Fill(my_mcp_ob.Eta(), (my_pfo_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                h_2d_pforelE_eta.Fill(my_mcp_ob.Eta(), (my_pfo_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                h_2d_pfoSFs.Fill(my_pfo_ob.E(), my_mcp_ob.E()/my_pfo_ob.E())

            if has_trk_ob:
                hists["trk_ob_pt"].Fill(my_trk_ob.Perp())
                hists["trk_ob_eta"].Fill(my_trk_ob.Eta())
                hists["trk_ob_phi"].Fill(my_trk_ob.Phi())
                hists["mcp_ob_trkmatch_pt"].Fill(my_mcp_ob.Perp())
                hists["mcp_ob_trkmatch_eta"].Fill(my_mcp_ob.Eta())
                hists["mcp_ob_trkmatch_phi"].Fill(my_mcp_ob.Phi())
                hists2d["trk_ob_v_mcp_ob_pt"].Fill(my_mcp_ob.Perp(), my_trk_ob.Perp())
                hists2d["trk_ob_v_mcp_ob_eta"].Fill(my_mcp_ob.Eta(), my_trk_ob.Eta())
                hists2d["trk_ob_v_mcp_ob_phi"].Fill(my_mcp_ob.Phi(), my_trk_ob.Phi())
                hists["d_trk_dpt"].Fill((my_trk_ob.Perp()-my_mcp_ob.Perp()))
                hists["d_trk_drelpt"].Fill((my_trk_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                hists["d_trk_drelE"].Fill((my_trk_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                hists["d_trk_dphi"].Fill((my_trk_ob.Phi()-my_mcp_ob.Phi()))
                hists["d_trk_deta"].Fill((my_trk_ob.Eta()-my_mcp_ob.Eta()))
                h_2d_trkrelpt.Fill(my_mcp_ob.Perp(), (my_trk_ob.Perp() - my_mcp_ob.Perp())/my_mcp_ob.Perp())
                h_2d_trkrelE.Fill(my_mcp_ob.E(), (my_trk_ob.E() - my_mcp_ob.E())/my_mcp_ob.E())
                if abs(my_mcp_ob.Eta()) < 1.1:
                    h_2d_trkrelpt1p1.Fill(my_mcp_ob.Perp(), (my_trk_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                    h_2d_trkrelE1p1.Fill(my_mcp_ob.E(), (my_trk_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                elif abs(my_mcp_ob.Eta()) < 1.2:
                    h_2d_trkrelpt1p2.Fill(my_mcp_ob.Perp(), (my_trk_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                    h_2d_trkrelE1p2.Fill(my_mcp_ob.E(), (my_trk_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                else:
                    h_2d_trkrelpt2.Fill(my_mcp_ob.Perp(), (my_trk_ob.Perp()-my_mcp_ob.Perp())/my_mcp_ob.Perp())
                    h_2d_trkrelE2.Fill(my_mcp_ob.E(), (my_trk_ob.E()-my_mcp_ob.E())/my_mcp_ob.E())
                h_2d_trkrelpt_eta.Fill(my_mcp_ob.Eta(), (my_trk_ob.Perp() - my_mcp_ob.Perp())/my_mcp_ob.Perp())
                h_2d_trkrelE_eta.Fill(my_mcp_ob.Eta(), (my_trk_ob.E() - my_mcp_ob.E())/my_mcp_ob.E())

            if has_pfo_ob and not has_trk_ob and not obj_type=="ph":
                print("Found event with a PFO but no track")

        i+=1
    reader.close()


# Fill the histograms comparing properties


# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################

print("Overall track efficiency:", hists["trk_ob_pt"].Integral()/hists["mcp_ob_pt"].Integral())
print("Overall pfo efficiency:", hists["pfo_ob_pt"].Integral()/hists["mcp_ob_pt"].Integral())

# Draw all the 1D histograms you filled
for var in dvariables:
    for obj in ["d_trk", "d_ob"]:
        c = ROOT.TCanvas("c%s%s"%(var,obj), "c")
        h = hists[obj+"_"+var]
        h.Draw()
        h.GetXaxis().SetTitle(dvariables[var]["title"])
        h.GetYaxis().SetTitle("Entries")
        f = ROOT.TF1("f%s%s"%(obj,var), "gaus")
        f.SetLineColor(colors[0])
        h.Fit("f%s%s"%(obj,var))
        c.SetLogy()
        latex = ROOT.TLatex()
        p = f.GetParameters()
        latex.DrawLatexNDC(.64, .85, "Mean: %f"%p[1])
        latex.DrawLatexNDC(.64, .78, "Sigma: %f"%p[2])
        if obj == "d_trk": latex.DrawLatexNDC(.64, .71, "Track Resolution")
        if obj == "d_ob": latex.DrawLatexNDC(.64, .71, "PFO Resolution")
        c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/fit_{obj}_{var}.png")

# Draw basic distributions
for var in ["pt", "eta", "phi"]:
    h_to_plot = {}
    for obj in ["trk_ob", "pfo_ob", "mcp_ob"]:
        h_to_plot[obj] = hists[obj+"_"+var]
    plotHistograms(h_to_plot, f"plots/{settings['plotdir'][obj_type]}/comp_{var}.png", variables[var]["title"], "Count")

# Make efficiency plots
for var in ["pt", "eta", "phi"]:
    efficiency_map = {}
    for obj in ["mcp_ob_trkmatch", "mcp_ob_pfomatch"]:
        efficiency_map[objects[obj]] = ROOT.TEfficiency(hists[obj+"_"+var], hists["mcp_ob_"+var])
        efficiency_map[objects[obj]].SetName("eff_"+obj+"_"+var)
    plotEfficiencies(efficiency_map, f"plots/{settings['plotdir'][obj_type]}/comp_eff_{var}.png", variables[var]["title"], "Efficiency")

# Make 2D plots comparing true v reco quantities
for hist in hists2d:
    c = ROOT.TCanvas("c_%s"%hist, "c")
    hists2d[hist].Draw("colz")
    var = hist.split("_")[-1]
    obj = hist.split("_")[0]
    hists2d[hist].GetXaxis().SetTitle("True "+settings['labelname'][obj_type]+" "+variables[var]["title"])
    hists2d[hist].GetYaxis().SetTitle(objects[obj]+" "+variables[var]["title"])
    c.SetRightMargin(0.18)
    c.SetLogz()
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{hist}.png")

# Make 2D plot and a TProfile to understand pT resolution v pT
for h in [h_2d_pforelpt, h_2d_trkrelpt, h_2d_trkrelpt1p1, h_2d_trkrelpt1p2, h_2d_trkrelpt2, h_2d_pforelpt1p1, h_2d_pforelpt1p2, h_2d_pforelpt2]:
    c = ROOT.TCanvas("can", "can")
    c.SetRightMargin(0.18)
    h.Draw("colz")
    h.GetXaxis().SetTitle(settings['labelname'][obj_type]+" p_{T} [GeV]")
    h.GetYaxis().SetTitle(dvariables["drelpt"]["title"])
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.root")

    c = ROOT.TCanvas("crelpt2dprof", "crelpt2dprof")
    h_prof = h.ProfileX("_pfx", 1, -1, "s")
    h_prof.Draw()
    h_prof.SetMinimum(-0.5)
    h_prof.SetMaximum(0.5)
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.root")

for h in [h_2d_pforelpt_eta, h_2d_trkrelpt_eta]:
    c = ROOT.TCanvas("can", "can")
    c.SetRightMargin(0.18)
    h.Draw("colz")
    h.GetXaxis().SetTitle(settings['labelname'][obj_type]+" #eta [GeV]")
    h.GetYaxis().SetTitle(dvariables["drelpt"]["title"])
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.root")

    c = ROOT.TCanvas("crelpt2dprof", "crelpt2dprof")
    h_prof = h.ProfileX("_pfx", 1, -1, "s")
    h_prof.Draw()
    h_prof.SetMinimum(-0.5)
    h_prof.SetMaximum(0.5)
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.root")

for h in [h_2d_pfoSFs, h_2d_pforelE, h_2d_pforelE_vmeas, h_2d_trkrelE, h_2d_pforelE1p1, h_2d_pforelE1p2, h_2d_pforelE2, h_2d_trkrelE1p1, h_2d_trkrelE1p1, h_2d_trkrelE2]:
    c = ROOT.TCanvas("can", "can")
    c.SetRightMargin(0.18)
    h.Draw("colz")
    h.GetXaxis().SetTitle(settings['labelname'][obj_type]+" E [GeV]")
    h.GetYaxis().SetTitle(dvariables["drelE"]["title"])
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.root")

    c = ROOT.TCanvas("crelE2dprof", "crelE2dprof")
    h_prof = h.ProfileX("_pfx", 1, -1, "s")
    h_prof.Draw()
    h_prof.SetMinimum(-0.5)
    h_prof.SetMaximum(0.5)
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.root")

for h in [h_2d_pforelE_eta, h_2d_trkrelE_eta]:
    c = ROOT.TCanvas("can", "can")
    c.SetRightMargin(0.18)
    h.Draw("colz")
    h.GetXaxis().SetTitle(settings['labelname'][obj_type]+" #eta [GeV]")
    h.GetYaxis().SetTitle(dvariables["drelE"]["title"])
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}.root")

    c = ROOT.TCanvas("crelE2dprof", "crelE2dprof")
    h_prof = h.ProfileX("_pfx", 1, -1, "s")
    h_prof.Draw()
    h_prof.SetMinimum(-0.5)
    h_prof.SetMaximum(0.5)
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.png")
    c.SaveAs(f"plots/{settings['plotdir'][obj_type]}/{h.GetTitle()}_prof.root")
