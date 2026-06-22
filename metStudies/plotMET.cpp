// plotMET.cpp
//
// Fast MET plotter for slcio files (designed for full-BIB samples).
//
// Two things make this fast compared to the python scripts in this repo:
//   1. The event loop is compiled C++, so there's no per-PFO interpreter overhead.
//   2. We call setReadCollectionNames() so the reader ONLY unpacks the PFO
//      collection. BIB files are dominated by calo/tracker hit collections,
//      and skipping their decompression is typically a 10-50x speedup.
//
// Build:  make            (see Makefile; run inside the mucoll singularity env)
// Run:    ./plotMET [-n maxEvents] [-o out.root] [-c collectionName] file1.slcio [file2.slcio ...]
// e.g.:   ./plotMET -o met_nuGun.root /data/fmeloni/DataMuC_MuColl10_v0A/muonGun_1000/recoBIB/*.slcio
//
// Output: a ROOT file with the histograms, plus quick-look PNGs in plots/.

#include <lcio.h>
#include <IOIMPL/LCFactory.h>
#include <IO/LCReader.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <Exceptions.h>

#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int main(int argc, char** argv) {

    // ############## ARGUMENT PARSING #############################
    int max_events = -1;
    std::string out_name = "met_hists.root";
    std::string pfo_col_name = "PandoraPFOs";
    std::vector<std::string> fnames;

    for (int i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "-n") == 0 && i+1 < argc) max_events = atoi(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc) out_name = argv[++i];
        else if (strcmp(argv[i], "-c") == 0 && i+1 < argc) pfo_col_name = argv[++i];
        else fnames.push_back(argv[i]);
    }
    if (fnames.empty()) {
        printf("Usage: %s [-n maxEvents] [-o out.root] [-c collectionName] file1.slcio [file2.slcio ...]\n", argv[0]);
        return 1;
    }
    printf("Found %zu files.\n", fnames.size());

    // ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
    gROOT->SetBatch(true);
    TH1F h_met      ("met",       ";MET [GeV];Events",                100,     0, 1000);
    TH1F h_met_low  ("met_low",   ";MET [GeV];Events",                100,     0,  100);
    TH1F h_met_phi  ("met_phi",   ";MET #phi;Events",                  32, -3.2,  3.2);
    TH1F h_met_x    ("met_x",     ";MET_{x} [GeV];Events",            100,  -500,  500);
    TH1F h_met_y    ("met_y",     ";MET_{y} [GeV];Events",            100,  -500,  500);
    TH1F h_sumet    ("sumet",     ";#Sigma E_{T} [GeV];Events",       100,     0, 5000);
    TH1F h_npfo     ("npfo",      ";Number of PFOs;Events",           100,     0, 2000);
    // Full 3D missing momentum (includes pz), and total visible energy.
    // At a muon collider sqrt(s) is known, so E_miss = sqrt(s) - sum E can be
    // derived from sume for any sample.
    TH1F h_pmiss    ("pmiss",     ";p_{miss} [GeV];Events",           100,     0, 1000);
    TH1F h_pmiss_low("pmiss_low", ";p_{miss} [GeV];Events",           100,     0,  100);
    TH1F h_pmiss_z  ("pmiss_z",   ";p_{miss,z} [GeV];Events",         100,  -500,  500);
    TH1F h_sume     ("sume",      ";#Sigma E [GeV];Events",           100,     0, 5000);
    for (TH1F* h : {&h_met, &h_met_low, &h_met_phi, &h_met_x, &h_met_y, &h_sumet, &h_npfo,
                    &h_pmiss, &h_pmiss_low, &h_pmiss_z, &h_sume})
        h->Sumw2();

    // ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
    IO::LCReader* reader = IOIMPL::LCFactory::getInstance()->createLCReader();

    // The crucial line: only unpack the PFO collection, skip the (huge) hit collections
    std::vector<std::string> read_cols = { pfo_col_name };
    reader->setReadCollectionNames(read_cols);

    reader->open(fnames);

    auto t_start = std::chrono::steady_clock::now();
    int i_evt = 0;
    EVENT::LCEvent* event = nullptr;
    while ((event = reader->readNextEvent()) != nullptr) {
        if (max_events > 0 && i_evt >= max_events) break;
        if (i_evt % 100 == 0) printf("Processing event %i.\n", i_evt);
        i_evt++;

        EVENT::LCCollection* pfos = nullptr;
        try { pfos = event->getCollection(pfo_col_name); }
        catch (EVENT::DataNotAvailableException&) { continue; }

        double mex = 0, mey = 0, mez = 0, sumet = 0, sume = 0;
        int n = pfos->getNumberOfElements();
        for (int j = 0; j < n; j++) {
            auto* pfo = static_cast<EVENT::ReconstructedParticle*>(pfos->getElementAt(j));
            const double* p = pfo->getMomentum();
            mex   -= p[0];
            mey   -= p[1];
            mez   -= p[2];
            sumet += std::sqrt(p[0]*p[0] + p[1]*p[1]);
            sume  += pfo->getEnergy();
        }

        double met   = std::sqrt(mex*mex + mey*mey);
        double pmiss = std::sqrt(mex*mex + mey*mey + mez*mez);
        h_met.Fill(met);
        h_met_low.Fill(met);
        h_met_phi.Fill(std::atan2(mey, mex));
        h_met_x.Fill(mex);
        h_met_y.Fill(mey);
        h_sumet.Fill(sumet);
        h_npfo.Fill(n);
        h_pmiss.Fill(pmiss);
        h_pmiss_low.Fill(pmiss);
        h_pmiss_z.Fill(mez);
        h_sume.Fill(sume);
    }
    reader->close();

    // ############## SUMMARY, SAVE, AND QUICK-LOOK PLOTS #############################
    auto t_end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(t_end - t_start).count();
    printf("\nSummary statistics:\n");
    printf("Ran over %i events in %.1f s (%.1f events/s).\n", i_evt, seconds, i_evt/seconds);
    printf("Mean MET: %.2f GeV\n", h_met.GetMean());

    TFile fout(out_name.c_str(), "RECREATE");
    gSystem->mkdir("plots", true);
    for (TH1F* h : {&h_met, &h_met_low, &h_met_phi, &h_met_x, &h_met_y, &h_sumet, &h_npfo,
                    &h_pmiss, &h_pmiss_low, &h_pmiss_z, &h_sume}) {
        h->Write();
        TCanvas c("c", "c");
        h->SetLineWidth(2);
        h->Draw("hist");
        c.SetLogy();
        c.SaveAs(Form("plots/%s.png", h->GetName()));
    }
    fout.Close();
    printf("Wrote histograms to %s and quick-look plots to plots/.\n", out_name.c_str());

    return 0;
}
