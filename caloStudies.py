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
magnetic_field = 5.00
calibration_mip = 0.0001575
calibration_mip_to_reco = 0.00641222630095
sampling_scaling = calibration_mip_to_reco/calibration_mip
append = "rose3"

# Set up things for each object
settings = {
        "fnames": {
                    "ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/photonGun*",
                    #"ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco_highrange/photonGun*",
                    #"ph": "/data/fmeloni/DataMuC_MuColl10_v0A/reco/photonGun*",
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
    if abs(tlv.Eta()) < 2:
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
#variables["pt"] =  {"nbins": 30, "xmin": 0, "xmax": 3000,   "title": "p_{T} [GeV]"}
variables["E"] =   {"nbins": 50, "xmin": 0, "xmax": 1000,   "title": "E [GeV]"}
#variables["eta"] = {"nbins": 30, "xmin": -3, "xmax": 3,     "title": "#eta"}
#variables["phi"] = {"nbins": 30, "xmin": -3.5, "xmax": 3.5, "title": "#phi"}
#variables["n"] =   {"nbins": 20, "xmin": 0, "xmax": 20,     "title": "n"}
hists = {}

objects = {}
objects["mcp"] = f"True {settings['labelname'][obj_type]}"
objects["sim"] = "Sim Calorimeter"
objects["dig"] = "Digi Calorimeter"
objects["rec"] = "Reco Calorimeter"
objects["clu"] = "Matched Cluster"
objects["pfo"] = f"Reconstructed {settings['labelname'][obj_type]}"

for obj in objects:
    for var in variables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, objects[obj], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

ranges = ["_0to1p1", "_1p1to1p2", "_1p2to2"]

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
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])
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
        simCollection_b = event.getCollection("ECalBarrelCollection")
        simCollection_e = event.getCollection("ECalEndcapCollection")
        try: digCollection_b = event.getCollection("EcalBarrelCollectionDigi")
        except: digCollection_b = None
        try: digCollection_e = event.getCollection("EcalEndcapCollectionDigi")
        except: digCollection_e = None
        try: recCollection_b = event.getCollection("EcalBarrelCollectionRec")
        except: recCollection_b = None
        try: recCollection_e = event.getCollection("EcalEndcapCollectionRec")
        except: recCollection_e = None
        cluCollection = event.getCollection("PandoraClusters")

        # Make counter variables
        n_mcp_ob = 0
        n_pfo_ob = 0
        n_clu_ob = 0
        has_mcp_ob = False
        has_pfo_ob = False
        has_clu_ob = False
        my_pfo_ob = None
        my_clu_ob = None
        my_mcp_ob = None

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)

            if abs(mcp.getPDG())==settings['pdgid'][obj_type] and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_ob = True
                n_mcp_ob += 1
                my_mcp_ob = mcp_tlv

        # Loop over the reconstructed clusters and fill histograms
        # If there are multiple, it'll keep the one with the higher pT
        for clu in cluCollection:
            #clu_tlv = ROOT.TLorentzVector()
            #clu_tlv.SetE(clu.getEnergy())
            #clu_tlv.SetTheta(clu.getITheta())
            #clu_tlv.SetPhi(clu.getIPhi())
            if has_mcp_ob: # and isMatched(clu_tlv, my_mcp_ob, req_pt = False):
                n_clu_ob += 1
                has_clu_ob = True
                if n_clu_ob == 1:
                    my_clu_ob = clu
                elif n_clu_ob > 1 and clu.getEnergy() > my_clu_ob.getEnergy():
                    my_clu_ob = clu

        # Loop over the reconstructed objects and fill histograms
        # If there are multiple, it'll keep the one with the higher pT
        for pfo in pfoCollection:
            pfo_tlv = getTLV(pfo)

            if abs(pfo.getType())==settings['pdgid'][obj_type]:
                if has_mcp_ob: # and isMatched(pfo_tlv, my_mcp_ob, req_pt = False):
                    n_pfo_ob += 1
                    has_pfo_ob = True
                    if n_pfo_ob == 1:
                        my_pfo_ob = pfo_tlv
                    elif n_pfo_ob > 1 and pfo_tlv.E() > my_pfo_ob.E():
                        my_pfo_ob = pfo_tlv

        # Loop over sim calo hits and sum
        sim_E = 0
        if simCollection_b:
            #print("\nn barrel sim hits:", len(simCollection_b))
            for sim in simCollection_b: sim_E += sim.getEnergy()*sampling_scaling
        if simCollection_e:
            for sim in simCollection_e: sim_E += sim.getEnergy()*sampling_scaling

        # Loop over digi hits and sum
        dig_E = 0
        if digCollection_b:
            #print("n barrel digi hits:", len(digCollection_b))
            #print(type(digCollection_b.data()))
            #tmp_E = sum(digCollection_b.data())
            for dig in digCollection_b: dig_E += dig.getEnergy()*calibration_mip_to_reco
        if digCollection_e:
            for dig in digCollection_e: dig_E += dig.getEnergy()*calibration_mip_to_reco

        # Loop over reco hits and sum
        rec_E = 0
        if recCollection_b:
            #print("n barrel rec hits:", len(recCollection_b))
            for rec in recCollection_b: rec_E += rec.getEnergy()
        if recCollection_e:
            for rec in recCollection_e: rec_E += rec.getEnergy()

        # Print out for detailed studies
        '''
        for sim, dig, rec in zip(simCollection_b, digCollection_b, recCollection_b):
            print("sim:", sim.getEnergy(), "dig:", dig.getEnergy(), "ratio:", dig.getEnergy()/sim.getEnergy())
            print("sim:", sim.getEnergy(), "dig*calibration_mip", dig.getEnergy()*calibration_mip)
            print("rec:", rec.getEnergy(), "dig*calibration_mip_to_reco", dig.getEnergy()*calibration_mip_to_reco)
        '''

        #print("mcp", my_mcp_ob.E(), "clus", my_clu_ob.getEnergy(), "sim:", sim_E, "dig:", dig_E)
        #print("total energies - mcp:", my_mcp_ob.E(), "simx35:", sim_E, "digi*calib:", dig_E*calibration_mip, "reco:", rec_E)


        if n_pfo_ob > 1: print("Found multiple PFOs:", n_pfo_ob)

        # Only make plots for events with isGood mcps
        if has_mcp_ob:
            #print("filling mcp_ob_pt with", my_mcp_ob.Perp())
            hists["mcp_E"].Fill(my_mcp_ob.E())
            if has_pfo_ob:
                hists["pfo_E"].Fill(my_pfo_ob.E())
                #hists2d["pfo_v_mcp_E"].Fill(my_mcp_ob.E(), my_pfo_ob.E())
            if has_clu_ob:
                hists["clu_E"].Fill(my_clu_ob.getEnergy())
                #hists2d["clu_v_mcp_E"].Fill(my_mcp_ob.E(), my_clu_ob.getEnergy())
            hists["sim_E"].Fill(sim_E)
            hists["dig_E"].Fill(dig_E)
            hists["rec_E"].Fill(rec_E)
            #hists2d["sim_v_mcp_E"].Fill(my_mcp_ob.E(), sim_E)
            #hists2d["dig_v_mcp_E"].Fill(my_mcp_ob.E(), dig_E)
            #hists2d["rec_v_mcp_E"].Fill(my_mcp_ob.E(), rec_E)

            # Print out 2D distributions per eta range
            for r in ranges:
                r1 = r.replace("p", ".").strip("_")
                low_eta = r1.split("to")[0]
                high_eta = r1.split("to")[1]
                selection_string = f"my_mcp_ob.Eta()>={low_eta} and my_mcp_ob.Eta()<{high_eta}"
                if eval(selection_string):
                    if has_pfo_ob: hists2d["pfo_v_mcp_E"+r].Fill(my_mcp_ob.E(), my_pfo_ob.E())
                    if has_clu_ob: hists2d["clu_v_mcp_E"+r].Fill(my_mcp_ob.E(), my_clu_ob.getEnergy())
                    hists2d["sim_v_mcp_E"+r].Fill(my_mcp_ob.E(), sim_E)
                    hists2d["dig_v_mcp_E"+r].Fill(my_mcp_ob.E(), dig_E)
                    hists2d["rec_v_mcp_E"+r].Fill(my_mcp_ob.E(), rec_E)

        i+=1
    reader.close()


# Fill the histograms comparing properties


# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################

# Draw basic distributions
for var in variables:
    h_to_plot = {}
    for obj in objects:
        h_to_plot[obj] = hists[obj+"_"+var]
    plotHistograms(h_to_plot, f"plots/calo/comp_{var}_{append}.png", variables[var]["title"], "Count")

'''
# Make efficiency plots
for var in ["pt", "eta", "phi"]:
    efficiency_map = {}
    for obj in ["mcp_ob_trkmatch", "mcp_ob_pfomatch"]:
        efficiency_map[objects[obj]] = ROOT.TEfficiency(hists[obj+"_"+var], hists["mcp_ob_"+var])
        efficiency_map[objects[obj]].SetName("eff_"+obj+"_"+var)
    plotEfficiencies(efficiency_map, f"plots/calo/comp_eff_{var}.png", variables[var]["title"], "Efficiency")
'''

# Make 2D plots comparing true v reco quantities
for hist in hists2d:
    c = ROOT.TCanvas("c_%s"%hist, "c")
    hists2d[hist].Draw("colz")
    var = hist.split("_")[-2]
    obj = hist.split("_")[0]
    hists2d[hist].GetXaxis().SetTitle("True "+settings['labelname'][obj_type]+" "+variables[var]["title"])
    hists2d[hist].GetYaxis().SetTitle(objects[obj]+" "+variables[var]["title"])
    c.SetRightMargin(0.18)
    c.SetLogz()
    c.SaveAs(f"plots/calo/{hist}_{append}.png")

