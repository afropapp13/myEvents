#include "TChain.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include <iostream>
#include <fstream>

void lengthTruthRecoComparison() {
 
    //
    // setup
    //

    TChain *t = new TChain("analysistree/anatree");
    t->Add("/uboone/data/users/alister1/10mslifetime_v06_26_01_02.root");

    // TTree *t = (TTree*)_file0->Get("analysistree/anatree");

    // open output text file
    std::ofstream ofs;
    ofs.open ("listBadLength.txt", std::ofstream::out | std::ofstream::app);

    // Open output root file
    TFile f("lengthHistograms.root", "RECREATE");

    int nentries = 10000;
    const int maxG4Entries = 50000;
    const int maxTracks = 200;
    // aux info for tracking down badly reconstructed protons
    
    Int_t run;
    Int_t subrun;
    Int_t event;
    Float_t trkstartx[maxTracks];
    Float_t trkstartz[maxTracks];

    // g4 params
    Int_t no_primaries;
    Int_t pdg[maxG4Entries];
    Int_t status[maxG4Entries];
    Float_t pathlen[maxG4Entries];
    Int_t origin[maxG4Entries];
    Int_t TrackId[maxG4Entries];

    // reconstructed info
    Short_t ntracks;
    Float_t trklen[maxTracks];
    Int_t trkorig[maxTracks];
    Int_t trkg4id[maxTracks];

    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("subrun", &subrun);
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("trkstartx_pandoraNu", trkstartx);
    t->SetBranchAddress("trkstartz_pandoraNu", trkstartz);

    t->SetBranchAddress("geant_list_size", &no_primaries);
    t->SetBranchAddress("pdg", pdg);
    t->SetBranchAddress("status", status);
    t->SetBranchAddress("pathlen", pathlen);
    t->SetBranchAddress("origin", origin);
    t->SetBranchAddress("TrackId", TrackId);

    t->SetBranchAddress("ntracks_pandoraNu", &ntracks);
    t->SetBranchAddress("trklen_pandoraNu", trklen);
    t->SetBranchAddress("trkorig_pandoraNu", trkorig);
    t->SetBranchAddress("trkg4id_pandoraNu", trkg4id);

    //
    // Logic begins
    //

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
    TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);

    TH2D* hist = new TH2D("hist", "", 50, 0, 100, 50, 0, 100);
    TH1D* hist2 = new TH1D("hist2", "", 100, -50, 50);

    int recoID, G4ID;
    for (Long64_t i = 0; i < nentries; i++) {
    
        if (i%1000 == 0) std::cout << i << "/" << nentries << std::endl;
   
        t->GetEntry(i);

        
        // reconstructed info
        for (int recoTracks = 0; recoTracks < ntracks; recoTracks++){

            recoID = trkg4id[recoTracks];

            for (int G4Tracks = 0; G4Tracks < no_primaries; G4Tracks++ ){
            
                G4ID = TrackId[G4Tracks];

                if (recoID == G4ID){
                
                     if (pdg[G4Tracks] == 2212) {
                        
                         // std::cout << "true pathlen: " << pathlen[G4Tracks] << " reco pathlen: " << trklen[recoTracks] << std::endl;
                         hist->Fill(pathlen[G4Tracks], trklen[recoTracks]);
                         hist2->Fill(pathlen[G4Tracks] - trklen[recoTracks]);
                    
                     }
                }
            
            }

            
        }
        
    }

ofs.close();

gStyle->SetOptStat(0);

c1->cd();

hist->GetXaxis()->SetTitle("True pathlen");
hist->GetYaxis()->SetTitle("Reco tracklen");
hist->Draw("colz");

c2->cd();
hist2->GetXaxis()->SetTitle("truth - reco");
hist2->Draw();

c1->SaveAs("lengthTruthReco.eps", "eps");
c2->SaveAs("lengthTruthMinusReco.eps", "eps");
f.Write();

}

# ifndef __CINT__
int main() {

    lengthTruthRecoComparison();
    return 0;

}
# endif
