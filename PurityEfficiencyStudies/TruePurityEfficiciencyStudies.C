#define TruePurityEfficiciencyStudies_cxx
#include "TruePurityEfficiciencyStudies.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TProfile.h>

#include <iostream>
#include <iomanip>
#include <vector>

#include "ubana/myClasses/Tools.h"
#include "ubana/myClasses/STV_Tools.h"

using namespace std;

TString TrueToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}


void TruePurityEfficiciencyStudies::Loop() {

	// -----------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		Cuts = Cuts + VectorCuts[i];

		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

		// -----------------------------------------------------------------------------------------------------------------------------

	if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
	TH1D::SetDefaultSumw2();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	TString Extension = "";

	// Output Files

	TString FileName = "OutputFiles/TruthPurityEfficiencyStudies_Overlay9_Run1"+Extension+Cuts+".root";
	TFile* OutputFile = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

	// --------------------------------------------------------------------------------------------------------------------------------------------

	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	Tools tools;

	// -------------------------------------------------------------------------------------------------------------------------------------------

	int TrueCC1pCounter = 0;
	int TrueCCQElikeCounter = 0;
	
	double SumWeights = 0.;	
	
	// --------------------------------------------------------------------------------------------------------------------------------

	// POT Counting

	double POTCount = -99.;

	if (string(fWhichSample).find("Overlay") != std::string::npos) {

		TString PathToPOTFile = "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
		TFile* POTFile = TFile::Open(PathToPOTFile,"readonly");
		TH1D* POTCountHist = (TH1D*)(POTFile->Get("POTCountHist"));
		POTCount = POTCountHist->GetBinContent(1);
		POTFile->Close();

	}	
	
	// ------------------------------------------------------------------------------------------------------------------

	// POT Scaling

	double POTScalingFactor = 1.;

	double tor860_wcut = 1;
	double E1DCNT_wcut = 1.;
	double EXT = 1.;

	if (string(fWhichSample).find("Run1") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run1;
		E1DCNT_wcut = E1DCNT_wcut_Run1;
		EXT = EXT_Run1;

	}
	
	if (string(fWhichSample).find("Run2") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run2;
		E1DCNT_wcut = E1DCNT_wcut_Run2;
		EXT = EXT_Run2;

	}	
	
	if (string(fWhichSample).find("Run3") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run3;
		E1DCNT_wcut = E1DCNT_wcut_Run3;
		EXT = EXT_Run3;

	}	
	
	if (string(fWhichSample).find("Run4") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run4;
		E1DCNT_wcut = E1DCNT_wcut_Run4;
		EXT = EXT_Run4;

	}		

	if (string(fWhichSample).find("Run5") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run5;
		E1DCNT_wcut = E1DCNT_wcut_Run5;
		EXT = EXT_Run5;

	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;	

		// ------------------------------------------------------------------------------------------------------------------------------

//		double weight = 1.;
//		double T2Kweight = 1.;
		if (Weight <= 0 || Weight > 30) { continue; }
		if (T2KWeight <= 0 || T2KWeight > 30) { continue; }		
		// Weight from v3.0.4 to v.3.0.6 * weight from application of T2K tune
		POTScalingFactor = tor860_wcut / POTCount;
		double weight = POTScalingFactor * Weight * T2KWeight * ROOTinoWeight;

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Analysis over the simb::MCParticles

		std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
		std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();

		int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0;
		bool TrueCC1pEvent = false;
		bool TrueCCQElikeEvent = false;

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & 0 pi+/- (Ppi > 70 MeV / c)

		if (CC1p == 1) {
		
			// --------------------------------------------------------------------------------------------------------------------		
			
			// Containment of the vertex defined as the start point of the muon in the TPC
			// Soft fiducial volume for true vertex

			if (Muon_MCParticle_StartContainment->at(0) == 0) { continue; }					 				

			// --------------------------------------------------------------------------------------------------------------------		

			// True muon

			double TrueMuonCosTheta = Muon_MCParticle_CosTheta->at(0);

			// True CC1p event

			TrueCC1pEvent = true;
			TrueCC1pCounter++;

			TrueMuonCosThetaPlot->Fill(TrueMuonCosTheta,weight);

		} // End of signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & pi+/- (Ppi > 70 MeV / c)

	} // End of the loop over the events

	// --------------------------------------------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	std::cout << std::endl << TrueCC1pCounter << " CC1p events selected"<< std::endl << std::endl;

	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();

} // End of the program
