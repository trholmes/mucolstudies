import pyLCIO
import glob
import ctypes
import math

exec(open("./plotHelper.py").read())

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

'''
magnetic_field = 5.00
calibration_mip = 0.0001575
#calibration_mip_to_reco = 0.00641222630095
calibration_mip_to_reco = 0.0066150 # Value used in v3
sampling_scaling = calibration_mip_to_reco/calibration_mip

hcal_calibration_mip = 0.0004725
hcal_calibration_mip_to_reco = 0.024625
hcal_sampling_scaling = hcal_calibration_mip_to_reco/hcal_calibration_mip
'''
append = "v5"
file_loc = "/data/fmeloni/DataMuC_MAIA_v0/v5/reco/dijet_mjj*"
#file_loc = "/phchang/data/mumu_ZH*"

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
samples = glob.glob(file_loc)
fnames = []
for s in samples:
    print(s)
    #if "5000" in s or "1000" in s:
    #    print("\tskipping.")
    #    continue
    fnames += glob.glob(f"{s}/*.slcio")
print("Found %i files."%len(fnames))

# Define good particle
def isGood(tlv):
    if abs(tlv.Eta()) < 2:
        return True
    return False

# Perform matching between two TLVs
def isMatched(tlv1, tlv2, req_pt = False):
    if tlv1.DeltaR(tlv2) > 1: return False
    if req_pt:
        drelpt = abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp()
        if drelpt > 0.1*tlv2.Perp()/100: return False # Require 10% at 100, 20% at 200, ...
    return True

def getClusterEta(cluster):
    theta = cluster.getITheta()
    return -1*math.ln(math.tan(theta/2))

# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
# Set up histograms
# This is an algorithmic way of making a bunch of histograms and storing them in a dictionary
variables = {}
variables["pt"] =  {"nbins": 60, "xmin": 0, "xmax": 6000,   "title": "p_{T} [GeV]"}
variables["E"] =   {"nbins": 60, "xmin": 0, "xmax": 6000,   "title": "E [GeV]"}
variables["eta"] = {"nbins": 30, "xmin": -3, "xmax": 3,     "title": "#eta"}
variables["phi"] = {"nbins": 30, "xmin": -3.5, "xmax": 3.5, "title": "#phi"}
variables["dR"] = {"nbins": 50, "xmin": 0, "xmax": 1, "title": "#Delta R"}
#variables["n"] =   {"nbins": 20, "xmin": 0, "xmax": 20,     "title": "n"}
hists = {}

objects = {}
objects["mcp"] = "True Parton"
#objects["sim"] = "Sim Calorimeter"
#objects["dig"] = "Digi Calorimeter"
#objects["rec"] = "Selected Calorimeter"
#objects["clu"] = "Matched Cluster"
#objects["pfo"] = f"Reconstructed Jet}"
objects["jet"] = "Reconstructed Jet"

for obj in objects:
    for var in variables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, objects[obj], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

ranges = ["_0to1p1", "_1p1to1p2", "_1p2to2", ""]

# Initialize all the 2D histograms: the each of the above variables at each level vs the mcp value
hists2d = {}
for obj in objects:
    for var in variables:
        if obj == "mcp": continue
        for r in ranges:
            hists2d[obj+"_v_mcp_"+var+r] = ROOT.TH2F(obj+"_v_mcp_"+var+r, obj+"_v_mcp_"+var+r, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "JetOut"])
#reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])
i = 0
for f in fnames:
    reader.open(f)
    if max_events > 0 and i >= max_events: break

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        jetCollection = event.getCollection("JetOut")

        # Make counter variables
        n_mcp_ob = 0

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)
            my_jet_ob = None
            n_matched_jet_ob = 0

            # Look at initiating quarks/gluons
            if mcp.getGeneratorStatus()==23 and isGood(mcp_tlv):
                n_mcp_ob += 1

                # Loop over the reconstructed jets and fill histograms
                # If there are multiple, it'll keep the one with the higher pT
                for jet in jetCollection:
                    jet_tlv = getTLV(jet)
                    if isMatched(jet_tlv, mcp_tlv):
                        n_matched_jet_ob += 1
                        if n_matched_jet_ob == 1:
                            my_jet_ob = jet_tlv
                        elif jet_tlv.E() > my_jet_ob.E():
                            my_jet_ob = jet_tlv

               # if n_matched_jet_ob > 1: continue #FIXME

                hists["mcp_E"].Fill(mcp_tlv.E())
                hists["mcp_pt"].Fill(mcp_tlv.Perp())
                hists["mcp_eta"].Fill(mcp_tlv.Eta())
                hists["mcp_phi"].Fill(mcp_tlv.Phi())

                if n_matched_jet_ob == 0: continue
                hists["jet_E"].Fill(my_jet_ob.E())
                hists["jet_pt"].Fill(my_jet_ob.Perp())
                hists["jet_eta"].Fill(my_jet_ob.Eta())
                hists["jet_phi"].Fill(my_jet_ob.Phi())
                hists["jet_dR"].Fill(mcp_tlv.DeltaR(my_jet_ob))
                hists2d["jet_v_mcp_E"].Fill(mcp_tlv.E(), my_jet_ob.E())
                hists2d["jet_v_mcp_pt"].Fill(mcp_tlv.Perp(), my_jet_ob.Perp())
                hists2d["jet_v_mcp_eta"].Fill(mcp_tlv.Eta(), my_jet_ob.Eta())
                hists2d["jet_v_mcp_phi"].Fill(mcp_tlv.Phi(), my_jet_ob.Phi())

                # Print out 2D distributions per eta range
                for r in ranges:
                    if r == "": selection_string = "True"
                    else:
                        r1 = r.replace("p", ".").strip("_")
                        low_eta = r1.split("to")[0]
                        high_eta = r1.split("to")[1]
                        selection_string = f"mcp_tlv.Eta()>={low_eta} and mcp_tlv.Eta()<{high_eta}"
                    if eval(selection_string):
                        hists2d["jet_v_mcp_E"+r].Fill(mcp_tlv.E(), my_jet_ob.E())
                        hists2d["jet_v_mcp_pt"+r].Fill(mcp_tlv.Perp(), my_jet_ob.Perp())
                        hists2d["jet_v_mcp_eta"+r].Fill(mcp_tlv.Eta(), my_jet_ob.Eta())
                        hists2d["jet_v_mcp_phi"+r].Fill(mcp_tlv.Phi(), my_jet_ob.Phi())

        i+=1
    reader.close()


# Fill the histograms comparing properties
# Draw basic distributions
for var in variables:
    h_to_plot = {}
    for obj in objects:
        h_to_plot[obj] = hists[obj+"_"+var]
    plotHistograms(h_to_plot, f"plots/jets/comp_{var}_{append}.png", variables[var]["title"], "Count")

'''
# Make efficiency plots
for var in ["E, pt, eta, phi"]:
    efficiency_map = {}
    for obj in ["jet"]:
        efficiency_map[objects[obj]] = ROOT.TEfficiency(hists[obj+"_"+var], hists["mcp_"+var])
        efficiency_map[objects[obj]].SetName("eff_"+obj+"_"+var)
    plotEfficiencies(efficiency_map, f"plots/jets/comp_eff_{var}.png", variables[var]["title"], "Efficiency")
'''

# Make 2D plots comparing true v reco quantities
for hist in hists2d:
    c = ROOT.TCanvas("c_%s"%hist, "c")
    hists2d[hist].Draw("colz")
    var = hist.split("_")[3]
    obj = hist.split("_")[0]
    hists2d[hist].GetXaxis().SetTitle("True Parton "+variables[var]["title"])
    hists2d[hist].GetYaxis().SetTitle(objects[obj]+" "+variables[var]["title"])
    c.SetRightMargin(0.18)
    c.SetLogz()
    c.SaveAs(f"plots/jets/{hist}_{append}.png")

