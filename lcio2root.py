import ROOT
import pyLCIO
import argparse
import numpy as np

# ------------------------------------------------------------------------------
def parseArgs():

    parser = argparse.ArgumentParser(
        add_help=True,
        description=''
    )

    parser.add_argument("-i", "--input_files", action="store", help="Input file (LCIO)", required=True, nargs="+")
    parser.add_argument("-o", "--output_file", action="store", help="Output file (ROOT)", required=True)
    parser.add_argument("-n", "--n_events",     action="store", help="Number of events (max)", default=-1)

    args = parser.parse_args()

    return args

# -------------------------------------------------------------------------------------------------
def get_kinematics(vector_xyz):
    x, y, z = vector_xyz[0], vector_xyz[1], vector_xyz[2]
    transverse = np.hypot(x, y)
    magnitude  = np.sqrt(x**2 + y**2 + z**2)

    # Avoid division by zero
    theta = np.arccos(z / magnitude) if magnitude != 0 else 0.0
    eta = -np.log(np.tan(theta / 2.0)) if theta != 0 else 0.0
    phi = np.arctan2(y, x)

    return transverse, eta, theta, phi

# -------------------------------------------------------------------------------------------------
def update_particles( event, branches_temp, collection_name, prefix=""):

    if collection_name not in event.getCollectionNames():
        return branches_temp

    collection = event.getCollection(collection_name)

    for i in range(collection.getNumberOfElements()):
        p = collection.getElementAt(i)

        branches_temp[prefix+"pdg"].push_back(p.getPDG())

        mom = p.getMomentum()
        pt, eta, theta, phi = get_kinematics(mom)
        branches_temp[prefix+"pt"].push_back(pt)
        branches_temp[prefix+"eta"].push_back(eta)
        branches_temp[prefix+"theta"].push_back(theta)
        branches_temp[prefix+"phi"].push_back(phi)

        branches_temp[prefix+"energy"].push_back(p.getEnergy())

        vert = p.getVertex()
        branches_temp[prefix+"vx"].push_back(vert[0])
        branches_temp[prefix+"vy"].push_back(vert[1])
        branches_temp[prefix+"vz"].push_back(vert[2])
        if pt > 10: print( p.getPDG(), p.getEnergy() )
    return branches_temp

# -------------------------------------------------------------------------------------------------
def update_hits( event, branches_temp, collection_name, prefix=""):

    if collection_name not in event.getCollectionNames():
        return branches_temp

    collection = event.getCollection(collection_name)

    for i in range(collection.getNumberOfElements()):
        hit = collection.getElementAt(i)

        pos = hit.getPosition()
        branches_temp[prefix+"x"].push_back(pos[0])
        branches_temp[prefix+"y"].push_back(pos[1])
        branches_temp[prefix+"z"].push_back(pos[2])

        branches_temp[prefix+"energy"].push_back(hit.getEnergy())
        branches_temp[prefix+"time"].push_back(hit.getTime())

        depth_temp, eta_temp, theta_temp, phi_temp = get_kinematics(pos)
        branches_temp[prefix+"eta"].push_back(eta_temp)
        branches_temp[prefix+"theta"].push_back(theta_temp)
        branches_temp[prefix+"phi"].push_back(phi_temp)
        branches_temp[prefix+"depth"].push_back(depth_temp)

    return branches_temp

# -------------------------------------------------------------------------------------------------
def main():

    args = parseArgs()

    # ----- Declare Input ----- #

    input_files_lcio  = args.input_files

    print("Running over", len(input_files_lcio), "input files...")
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

    # ----- Declare Output ----- #

    output_file_root = args.output_file

    print("Creating output ROOT file:", output_file_root)
    f = ROOT.TFile(output_file_root, "RECREATE")
    tree = ROOT.TTree("Events", "Events")

    vec_float = ROOT.std.vector("float")

    branches = {}

    branches["mcp_pdg"]         = vec_float()
    branches["mcp_pt"]          = vec_float()
    branches["mcp_eta"]         = vec_float()
    branches["mcp_theta"]       = vec_float()
    branches["mcp_phi"]         = vec_float()
    branches["mcp_energy"]      = vec_float()
    branches["mcp_vx"]          = vec_float()
    branches["mcp_vy"]          = vec_float()
    branches["mcp_vz"]          = vec_float()

    branches["ecal_hit_x"]      = vec_float()
    branches["ecal_hit_y"]      = vec_float()
    branches["ecal_hit_z"]      = vec_float()
    branches["ecal_hit_energy"] = vec_float()
    branches["ecal_hit_time"]   = vec_float()
    branches["ecal_hit_theta"]  = vec_float()
    branches["ecal_hit_eta"]    = vec_float()
    branches["ecal_hit_phi"]    = vec_float()
    branches["ecal_hit_depth"]  = vec_float()

    branches["hcal_hit_x"]      = vec_float()
    branches["hcal_hit_y"]      = vec_float()
    branches["hcal_hit_z"]      = vec_float()
    branches["hcal_hit_energy"] = vec_float()
    branches["hcal_hit_time"]   = vec_float()
    branches["hcal_hit_theta"]  = vec_float()
    branches["hcal_hit_eta"]    = vec_float()
    branches["hcal_hit_phi"]    = vec_float()
    branches["hcal_hit_depth"]  = vec_float()


    for branch in branches:
        tree.Branch(branch, branches[branch])

    # ----- Loop ----- #

    do_break   = False
    max_events = int(args.n_events)
    ievt       = -1

    for input_file_lcio in input_files_lcio:
        print( input_files_lcio )

        if do_break: break

        reader.open( input_file_lcio )

        for ievt_in_file, event in enumerate(reader):

            #event = reader.readNextEvent()
            if event is None: break

            ievt += 1
            if ievt == max_events:
                do_break = True
                break

            if ievt % 100 == 0: print( "Processing event", ievt)

            # Reset vector branches
            for b in branches:
                branches[b].clear()

            branches = update_particles( event, branches, "MCParticle", "mcp_")
            branches = update_hits( event, branches, "EcalBarrelCollectionSel", "ecal_hit_")
            branches = update_hits( event, branches, "EcalEndcapCollectionSel", "ecal_hit_")
            branches = update_hits( event, branches, "HcalBarrelCollectionsConed", "hcal_hit_")
            branches = update_hits( event, branches, "HcalEndcapCollectionsConed", "hcal_hit_")

            tree.Fill()

    reader.close()

    print("Ran over", ievt, "events")
    f.Write()
    print("File written to", output_file_root)


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main()
