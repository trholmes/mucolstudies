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
obj_type = "ph"
append = "240415"

# Set up things for each object
settings = {
        "fnames": {
                    "ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/photonGun_E_250_1000*",
                    "mu": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun*",
                    "el": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*"},
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

# Define good particle
def isGood(tlv):
    if abs(tlv.Eta()) < 2.436:
        return True
    return False

# Perform matching between two TLVs
def isMatched(tlv1, tlv2, req_pt = True):
    if tlv1.DeltaR(tlv2) > 0.1: return False
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
#variables["z"] = {"nbins":500, "xmin": -10000, "xmax": 10000, "title": "z [mm]"}
variables["scale"] = {"nbins":500, "xmin": 0, "xmax": 100, "title": "true-over-sim ratio"}
#variables["pt"] =  {"nbins": 30, "xmin": 0, "xmax": 3000,   "title": "p_{T} [GeV]"}
#variables["E"] =   {"nbins": 50, "xmin": 0, "xmax": 1000,   "title": "E [GeV]"}
variables["eta"] = {"nbins": 48, "xmin": -2.4, "xmax": 2.4,     "title": "#eta"}
variables["theta"] = {"nbins": 63, "xmin": 0, "xmax": 3.15,     "title": "#theta"}
#variables["phi"] = {"nbins": 30, "xmin": -3.5, "xmax": 3.5, "title": "#phi"}
#variables["n"] =   {"nbins": 20, "xmin": 0, "xmax": 20,     "title": "n"}
hists = {}

# Initialize all the 2D histograms: the each of the above variables at each level vs the mcp value
hists2d = {}
for var in variables:
    if var=="scale": continue
    hists2d[var] = ROOT.TH2F(var, var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"], variables["scale"]["nbins"], variables["scale"]["xmin"], variables["scale"]["xmax"])

# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "ECalBarrelCollection", "ECalEndcapCollection"])
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
        simCollection_b = event.getCollection("ECalBarrelCollection")
        simCollection_e = event.getCollection("ECalEndcapCollection")

        # Make counter variables
        n_mcp_ob = 0
        has_mcp_ob = False
        my_mcp_ob = None

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)

            if abs(mcp.getPDG())==settings['pdgid'][obj_type] and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_ob = True
                n_mcp_ob += 1
                my_mcp_ob = mcp_tlv

        # Figure out total sim energy
        values_b = [x.getEnergy() for x in simCollection_b]
        values_e = [x.getEnergy() for x in simCollection_e]
        sim_E = sum(values_b) + sum(values_e)

        # Only make plots for events with isGood mcps
        if has_mcp_ob:
            #print("filling mcp_ob_pt with", my_mcp_ob.Perp())
            hists2d["eta"].Fill(my_mcp_ob.Eta(), my_mcp_ob.E()/sim_E)
            hists2d["theta"].Fill(my_mcp_ob.Theta(), my_mcp_ob.E()/sim_E)

        i+=1
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################

# Make 2D plots comparing true v reco quantities
for hist in hists2d:
    c = ROOT.TCanvas("c_%s"%hist, "c")
    hists2d[hist].Draw("colz")
    var = hist
    hists2d[hist].GetXaxis().SetTitle(variables[var]["title"])
    hists2d[hist].GetYaxis().SetTitle(variables["scale"]["title"])
    c.SetRightMargin(0.18)
    c.SetLogz()
    c.SaveAs(f"calib/{hist}_{append}.png")
    c.SaveAs(f"calib/{hist}_{append}.root")

    c = ROOT.TCanvas("cprof", "cprof")
    h_prof = hists2d[hist].ProfileX("_pfx", 1, -1, "s")
    h_prof.Draw()
    #h_prof.SetMinimum(-0.5)
    #h_prof.SetMaximum(0.5)
    c.SaveAs(f"calib/{hist}_prof_{append}.png")
    c.SaveAs(f"calib/{hist}_prof_{append}.root")

