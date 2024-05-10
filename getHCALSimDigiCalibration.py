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
obj_type = "ne"
append = "240417"
do_ECAL = False
do_HCAL = True
ECAL_calib_file = "calib/theta_prof_240415.root"

# Set up things for each object
settings = {
        "fnames": {
                    "ne": "/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/neutronGun_E_250_1000*",
                    "pi": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/pionGun_pT_250_1000*",
                    "ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/photonGun_E_250_1000*",
                    "mu": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun*",
                    "el": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*"},
        "labelname": {  "ne": "Neutron",
                        "pi": "Pion",
                        "ph": "Photon",
                        "mu": "Muon",
                        "el": "Electron"},
        "plotdir":{ "ne": "neutrons",
                    "pi": "pions",
                    "ph": "photons",
                    "mu": "muons",
                    "el": "electrons"},
        "pdgid":  { "ne": [2112],
                    "pi": [211, 111],
                    "ph": [22],
                    "mu": [13],
                    "el": [11]},
        "mass":   { "ne": 0.940,
                    "pi": 0.135,
                    "ph": 0,
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

# Get ECAL calibration histogram
h_calib = ROOT.TFile.Open(ECAL_calib_file).Get("cprof").GetPrimitive("theta_pfx")

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
reader.setReadCollectionNames(["MCParticle", "ECalBarrelCollection", "ECalEndcapCollection", "HCalBarrelCollection", "HCalEndcapCollection"])

i = 0
for f in fnames:
    reader.open(f)
    if max_events > 0 and i >= max_events: break

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        ecal_simCollection_b = event.getCollection("ECalBarrelCollection")
        ecal_simCollection_e = event.getCollection("ECalEndcapCollection")
        hcal_simCollection_b = event.getCollection("HCalBarrelCollection")
        hcal_simCollection_e = event.getCollection("HCalEndcapCollection")

        # Make counter variables
        n_mcp_ob = 0
        has_mcp_ob = False
        my_mcp_ob = None

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)

            if abs(mcp.getPDG()) in settings['pdgid'][obj_type] and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_ob = True
                n_mcp_ob += 1
                my_mcp_ob = mcp_tlv

        # Figure out total sim energy
        e_values_b = [x.getEnergy() for x in ecal_simCollection_b]
        e_values_e = [x.getEnergy() for x in ecal_simCollection_e]
        e_sim_E = sum(e_values_b) + sum(e_values_e)

        h_values_b = [x.getEnergy() for x in hcal_simCollection_b]
        h_values_e = [x.getEnergy() for x in hcal_simCollection_e]
        h_sim_E = sum(h_values_b) + sum(h_values_e)

        # Only make plots for events with isGood mcps
        if has_mcp_ob and h_sim_E > 0:
            # Calibrate the ecal energy and subtract off
            ecal_calib_value = h_calib.GetBinContent(h_calib.FindBin(my_mcp_ob.Theta()))
            estimated_hcal_energy = my_mcp_ob.E() - e_sim_E*ecal_calib_value
            #print(ecal_calib_value, my_mcp_ob.E(), estimated_hcal_energy)

            #print("filling mcp_ob_pt with", my_mcp_ob.Perp())
            hists2d["eta"].Fill(my_mcp_ob.Eta(), estimated_hcal_energy/h_sim_E)
            hists2d["theta"].Fill(my_mcp_ob.Theta(), estimated_hcal_energy/h_sim_E)

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
    c.SaveAs(f"calib/hcal_{hist}_{append}.png")
    c.SaveAs(f"calib/hcal_{hist}_{append}.root")

    c = ROOT.TCanvas("cprof", "cprof")
    h_prof = hists2d[hist].ProfileX("_pfx", 1, -1, "s")
    h_prof.Draw()
    #h_prof.SetMinimum(-0.5)
    #h_prof.SetMaximum(0.5)
    c.SaveAs(f"calib/hcal_{hist}_prof_{append}.png")
    c.SaveAs(f"calib/hcal_{hist}_prof_{append}.root")

