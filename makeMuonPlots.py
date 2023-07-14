import pyLCIO
import ROOT
import glob



# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/reco/*.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/gen_muonGun/recoBIB/*.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/muonGun_1000/recoBIB/*.slcio")
print("Found %i files."%len(fnames))



# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
# Set up histograms
# This is an algorithmic way of making a bunch of histograms and storing them in a dictionary
variables = {}
variables["pt"] =  {"nbins": 20, "xmin": 0, "xmax": 2000}
variables["eta"] = {"nbins": 20, "xmin": -3, "xmax": 3}
variables["phi"] = {"nbins": 20, "xmin": -3.5, "xmax": 3.5}
variables["n"] = {"nbins": 20, "xmin": 0, "xmax": 20}
hists = {}
for obj in ["pfo", "pfo_mu", "mcp", "mcp_mu", "mcp_mu_match"]:
    for var in variables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, obj+"_"+var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# Making a separate set of binning conventions for plots showing resolutions
# these plots will all be filled with the difference between a pfo and a mcp object value
dvariables = {}
dvariables["dpt"] =     {"nbins": 100, "xmin": -500, "xmax": 500}
dvariables["drelpt"] =  {"nbins": 100, "xmin": -0.5, "xmax": 0.5}
dvariables["dphi"] =    {"nbins": 100, "xmin": -0.001, "xmax": 0.001}
dvariables["deta"] =    {"nbins": 100, "xmin": -0.001, "xmax": 0.001}
for obj in ["d_mu"]:
    for var in dvariables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, obj+"_"+var, dvariables[var]["nbins"], dvariables[var]["xmin"], dvariables[var]["xmax"])

# Finally making one 2D histogram non-algorithmically; this is what I'll use for a
# pT resolution vs. pT plot.
h_2d_relpt = ROOT.TH2F("h_2d_relpt", "h_2d_relpt", 20, 0, 1000, 500, -0.5, 0.5)




# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
i = 0
for f in fnames:
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")

        # Make counter variables
        n_mcp_mu = 0
        n_pfo_mu = 0
        has_pfo_mu = False
        my_pfo_mu = 0

        # Loop over the reconstructed objects and fill histograms
        for pfo in pfoCollection:
            pfo_p = pfo.getMomentum()
            pfo_tlv = ROOT.TLorentzVector()
            pfo_tlv.SetPxPyPzE(pfo_p[0], pfo_p[1], pfo_p[2], pfo.getEnergy())
            hists["pfo_pt"].Fill(pfo_tlv.Perp())
            hists["pfo_eta"].Fill(pfo_tlv.Eta())
            hists["pfo_phi"].Fill(pfo_tlv.Phi())

            if abs(pfo.getType())==13:
                hists["pfo_mu_pt"].Fill(pfo_tlv.Perp())
                hists["pfo_mu_eta"].Fill(pfo_tlv.Eta())
                hists["pfo_mu_phi"].Fill(pfo_tlv.Phi())
                n_pfo_mu += 1
                has_pfo_mu = True
                my_pfo_mu = pfo_tlv     # Storing this to use for matching in the next loop

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            hists["mcp_pt"].Fill(mcp_tlv.Perp())
            hists["mcp_eta"].Fill(mcp_tlv.Eta())
            hists["mcp_phi"].Fill(mcp_tlv.Phi())

            if abs(mcp.getPDG())==13 and mcp.getGeneratorStatus()==1:
                hists["mcp_mu_pt"].Fill(mcp_tlv.Perp())
                hists["mcp_mu_eta"].Fill(mcp_tlv.Eta())
                hists["mcp_mu_phi"].Fill(mcp_tlv.Phi())
                n_mcp_mu += 1

                # For events in which a PFO mu was reconstructed, fill histograms that will
                # be used for efficiency. Both numerator and denominator must be filled with truth values!
                # Also fill resolution histograms
                if has_pfo_mu:
                    hists["mcp_mu_match_pt"].Fill(mcp_tlv.Perp())
                    hists["mcp_mu_match_eta"].Fill(mcp_tlv.Eta())
                    hists["mcp_mu_match_phi"].Fill(mcp_tlv.Phi())

                    hists["d_mu_dpt"].Fill(my_pfo_mu.Perp() - mcp_tlv.Perp())
                    hists["d_mu_drelpt"].Fill((my_pfo_mu.Perp() - mcp_tlv.Perp())/mcp_tlv.Perp())
                    hists["d_mu_deta"].Fill(my_pfo_mu.Eta() - mcp_tlv.Eta())
                    hists["d_mu_dphi"].Fill(my_pfo_mu.Phi() - mcp_tlv.Phi())
                    h_2d_relpt.Fill(mcp_tlv.Perp(), (my_pfo_mu.Perp() - mcp_tlv.Perp())/mcp_tlv.Perp())

        # This is here to check that we never reconstruct multiple muons
        # If we did, we'd have to match the correct muon to the MCP object to do eff/res plots
        # But since we don't, we can skip that step
        if n_pfo_mu > 1: print(n_pfo_mu)
        hists["mcp_n"].Fill(len(mcpCollection))
        hists["pfo_n"].Fill(len(pfoCollection))
        hists["mcp_mu_n"].Fill(n_mcp_mu)
        hists["pfo_mu_n"].Fill(n_pfo_mu)
        hists["mcp_mu_match_n"].Fill(n_pfo_mu)

        i+=1



# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
print("\t%i MCPs"%hists["mcp_pt"].GetEntries())
print("\t%i mu MCPs"%hists["mcp_mu_pt"].GetEntries())
print("\t%i PFOs"%hists["pfo_pt"].GetEntries())
print("\t%i mu PFOs"%hists["pfo_mu_pt"].GetEntries())

# Draw all the 1D histograms you filled
for i, h in enumerate(hists):
    c = ROOT.TCanvas("c%i"%i, "c%i"%i)
    hists[h].Draw()
    hists[h].GetXaxis().SetTitle(h)
    hists[h].GetYaxis().SetTitle("Entries")

    # For resolution plots, fit them and get the mean and sigma
    if h.startswith("d_mu"):
        f = ROOT.TF1("f%i"%i, "gaus")
        f.SetLineColor(ROOT.kRed)
        hists[h].Fit("f%i"%i)
        c.SetLogy()
        latex = ROOT.TLatex()
        p = f.GetParameters()
        latex.DrawLatexNDC(.64, .85, "Mean: %f"%p[1])
        latex.DrawLatexNDC(.64, .78, "Sigma: %f"%p[2])
    c.SaveAs("plots/%s.png"%h)

# Make efficiency plots
# In these files, there are at most 1 PFO mu, so matching isn't needed
for v in variables:
    if v=="n": continue
    c = ROOT.TCanvas("c%s"%v, "c%s"%v)
    eff = ROOT.TEfficiency(hists["mcp_mu_match_"+v], hists["mcp_mu_"+v])
    eff.Draw("ape")
    ROOT.gPad.Update()
    eff.SetLineWidth(2)
    eff.GetPaintedGraph().SetMinimum(0)
    eff.GetPaintedGraph().SetMaximum(1)
    eff.SetTitle(";%s;Efficiency"%v)
    c.SaveAs("plots/eff_%s.png"%v)

# Make 2D plot and a TProfile to understand pT resolution v pT
c = ROOT.TCanvas("crelpt2d", "crelpt2d")
h_2d_relpt.Draw("colz")
h_2d_relpt.GetXaxis().SetTitle("pt")
h_2d_relpt.GetYaxis().SetTitle("drelpt")
c.SaveAs("plots/d_mu_relpt_2d.png")

c = ROOT.TCanvas("crelpt2dprof", "crelpt2dprof")
h_prof = h_2d_relpt.ProfileX("_pfx", 1, -1, "s")
h_prof.GetXaxis().SetTitle("pt")
h_prof.GetYaxis().SetTitle("drelpt")
h_prof.Draw()
c.SaveAs("plots/d_mu_relpt_prof.png")
