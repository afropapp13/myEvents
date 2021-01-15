#define PurityEfficiencyStudies_cxx
#include "PurityEfficiencyStudies.h"
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
#include <string>
#include <sstream>

using namespace std;

#include "ubana/myClasses/Tools.h"

TString ToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

void PurityEfficiencyStudies::Loop() {

	// ---------------------------------------------------------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	int OldNEventsPassingSelectionCuts = 0;
	int CC1pEventsPassingSelectionCuts = 0;
		
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
		double weight = 1.;

		// ---------------------------------------------------------------------------------------------------------------------------

		TString Extension = "";

		TString FileName = "OutputFiles/PurityEfficiencyStudies_"+fWhichSample+Extension+Cuts+".root";

		TFile* file = new TFile(FileName,"recreate");
		std::cout << std::endl << "Creating a new file: " << FileName << std::endl << std::endl << std::endl;

		// --------------------------------------------------------------------------------------------------------------------------------

		int NBinsNuScore = 20;
		double MinNuScore = 0., MaxNuScore = 1.;
		double NuScoreStep = (MaxNuScore - MinNuScore) / double(NBinsNuScore);

		int NBinsLLP = 20;
		double MinLLP = -5., MaxLLP = 5.;
		double LLPStep = (MaxLLP - MinLLP) / double(NBinsLLP);

		int NBinsLLMu = 20;
		double MinLLMu = -5., MaxLLMu = 5.;
		double LLMuStep = (MaxLLMu - MinLLMu) / double(NBinsLLMu);

		int NBinsLength = 20;
		double MinLength = -150., MaxLength = 500.;
		double LengthStep = (MaxLength - MinLength) / double(NBinsLength);

		int NBinsPMissMinus = 20;
		double MinPMissMinus = 0., MaxPMissMinus = 1.5;
		double PMissMinusStep = (MaxPMissMinus - MinPMissMinus) / double(NBinsPMissMinus);

		int NBinskMiss = 20;
		double MinkMiss = 0., MaxkMiss = 1.05;
		double kMissStep = (MaxkMiss - MinkMiss) / double(NBinskMiss);

		int NBinsDeltaTheta = 19;
		double MinDeltaTheta = 5., MaxDeltaTheta = 90.;
		double DeltaThetaStep = (MaxDeltaTheta - MinDeltaTheta) / double(NBinsDeltaTheta);

		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots

		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		TH1D* OldCutsRecoMuonCosThetaPlot = new TH1D("OldCutsRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		// --------------------------------------------------------------------------------------------------------------------------------

		// 2D NuScore vs LLp Study 

		TH1D* RecoMuonCosThetaPlotStudy[NBinsNuScore][NBinsLLP];
		TH1D* CC1pRecoMuonCosThetaPlotStudy[NBinsNuScore][NBinsLLP];

		for (int WhichXBin = 0; WhichXBin < NBinsNuScore; WhichXBin++) {

			for (int WhichLegBin = 0; WhichLegBin < NBinsLLP; WhichLegBin++) {

				TString PlotThres = "RecoMuonCosThetaPlot_NuScoreVsLLPThres_"+TString(std::to_string(WhichXBin))+"_"+TString(std::to_string(WhichLegBin));

				RecoMuonCosThetaPlotStudy[WhichXBin][WhichLegBin] = new TH1D(PlotThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotThres = "CC1pRecoMuonCosThetaPlot_NuScoreVsLLPThres_"+TString(std::to_string(WhichXBin))+"_"+TString(std::to_string(WhichLegBin));

				CC1pRecoMuonCosThetaPlotStudy[WhichXBin][WhichLegBin] = new TH1D(CC1pPlotThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

			}

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// DeltaTheta Only Study

		TH1D* RecoMuonCosThetaPlotDeltaThetaOnlyStudy[NBinsLLP];
		TH1D* CC1pRecoMuonCosThetaPlotDeltaThetaOnlyStudy[NBinsLLP];

		for (int WhichBin = 0; WhichBin < NBinsDeltaTheta; WhichBin++) {

				TString PlotOnlyThres = "RecoMuonCosThetaPlot_DeltaThetaThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotDeltaThetaOnlyStudy[WhichBin] = new TH1D(PlotOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotOnlyThres = "CC1pRecoMuonCosThetaPlot_DeltaThetaThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotDeltaThetaOnlyStudy[WhichBin] = new TH1D(CC1pPlotOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// LLP Only Study

		TH1D* RecoMuonCosThetaPlotLLPOnlyStudy[NBinsLLP];
		TH1D* CC1pRecoMuonCosThetaPlotLLPOnlyStudy[NBinsLLP];

		for (int WhichBin = 0; WhichBin < NBinsLLP; WhichBin++) {

				TString PlotLLPOnlyThres = "RecoMuonCosThetaPlot_LLPThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLLPOnlyStudy[WhichBin] = new TH1D(PlotLLPOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLLPOnlyThres = "CC1pRecoMuonCosThetaPlot_LLPThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLLPOnlyStudy[WhichBin] = new TH1D(CC1pPlotLLPOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// LLMu Only Study

		TH1D* RecoMuonCosThetaPlotLLMuOnlyStudy[NBinsLLMu];
		TH1D* CC1pRecoMuonCosThetaPlotLLMuOnlyStudy[NBinsLLMu];

		for (int WhichBin = 0; WhichBin < NBinsLLMu; WhichBin++) {

				TString PlotLLMuOnlyThres = "RecoMuonCosThetaPlot_LLMuThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLLMuOnlyStudy[WhichBin] = new TH1D(PlotLLMuOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLLMuOnlyThres = "CC1pRecoMuonCosThetaPlot_LLMuThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLLMuOnlyStudy[WhichBin] = new TH1D(CC1pPlotLLMuOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// NuScore Only Study

		TH1D* RecoMuonCosThetaPlotNuScoreOnlyStudy[NBinsNuScore];
		TH1D* CC1pRecoMuonCosThetaPlotNuScoreOnlyStudy[NBinsNuScore];

		for (int WhichBin = 0; WhichBin < NBinsNuScore; WhichBin++) {

				TString PlotNuScoreOnlyThres = "RecoMuonCosThetaPlot_NuScoreThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotNuScoreOnlyStudy[WhichBin] = new TH1D(PlotNuScoreOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotNuScoreOnlyThres = "CC1pRecoMuonCosThetaPlot_NuScoreThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotNuScoreOnlyStudy[WhichBin] = new TH1D(CC1pPlotNuScoreOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// Length Only Study

		TH1D* RecoMuonCosThetaPlotLengthOnlyStudy[NBinsLength];
		TH1D* CC1pRecoMuonCosThetaPlotLengthOnlyStudy[NBinsLength];

		for (int WhichBin = 0; WhichBin < NBinsLength; WhichBin++) {

				TString PlotLengthOnlyThres = "RecoMuonCosThetaPlot_LengthThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLengthOnlyStudy[WhichBin] = new TH1D(PlotLengthOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLengthOnlyThres = "CC1pRecoMuonCosThetaPlot_LengthThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLengthOnlyStudy[WhichBin] = new TH1D(CC1pPlotLengthOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// PMissMinus Only Study

		TH1D* RecoMuonCosThetaPlotPMissMinusOnlyStudy[NBinsPMissMinus];
		TH1D* CC1pRecoMuonCosThetaPlotPMissMinusOnlyStudy[NBinsPMissMinus];

		for (int WhichBin = 0; WhichBin < NBinsPMissMinus; WhichBin++) {

				TString PlotPMissMinusOnlyThres = "RecoMuonCosThetaPlot_PMissMinusThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotPMissMinusOnlyStudy[WhichBin] = new TH1D(PlotPMissMinusOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotPMissMinusOnlyThres = "CC1pRecoMuonCosThetaPlot_PMissMinusThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotPMissMinusOnlyStudy[WhichBin] = new TH1D(CC1pPlotPMissMinusOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		// kMiss Only Study

		TH1D* RecoMuonCosThetaPlotkMissOnlyStudy[NBinskMiss];
		TH1D* CC1pRecoMuonCosThetaPlotkMissOnlyStudy[NBinskMiss];

		for (int WhichBin = 0; WhichBin < NBinskMiss; WhichBin++) {

				TString PlotkMissOnlyThres = "RecoMuonCosThetaPlot_kMissThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotkMissOnlyStudy[WhichBin] = new TH1D(PlotkMissOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotkMissOnlyThres = "CC1pRecoMuonCosThetaPlot_kMissThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotkMissOnlyStudy[WhichBin] = new TH1D(CC1pPlotkMissOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

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

		// ------------------------------------------------------------------------------------------------------------------------------

		Tools tools;

		// ------------------------------------------------------------------------------------------------------------------

		// POT Scaling

		double POTScale = 1.;

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
		
		if (string(fWhichSample).find("ExtBNB9") != std::string::npos) { weight = E1DCNT_wcut / EXT; POTScale = weight; }

		if (string(fWhichSample).find("Overlay") != std::string::npos) { weight = tor860_wcut / POTCount; POTScale = weight; }			


		// --------------------------------------------------------------------------------------------------------------------------------

		// Loop over the events

		for (Long64_t jentry=0; jentry<nentries;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

			// ---------------------------------------------------------------------------------------------------------------------

			if (jentry%10000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

			// ------------------------------------------------------------------------------------------------------------------------

			// Demand for exactly one muon/proton candidate

			if (CandidateMu_P_MCS->size() != 1) { continue; }
			if (CandidateP_P_Range->size() != 1) { continue; }

			// ------------------------------------------------------------------------------------------------------------------------

			// Containment Demand
			// Fully contained proton candidates
			// Either fully contained or semi-contained muon candidates

			if (CandidateMu_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_EndContainment->at(0) == 0) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 
			
				if (Weight < 0 || Weight > 10) { continue; }
				if (T2KWeight < 0 || T2KWeight > 10) { continue; }				
				weight = ( tor860_wcut / POTCount) * Weight * T2KWeight * ROOTinoWeight; 
				
			}

			// -----------------------------------------------------------------------------------------------------------------------

			if ( fabs(weight) != weight) { continue; } // Securing against infinities

			// -----------------------------------------------------------------------------------------------------------------------

			// hard coded limit because something still looks weird in the dirt sample at low nu-score
			// COH interaction & infinite weight
			
//			if (NuScore < 0.04) { continue; }

			// -----------------------------------------------------------------------------------------------------------------------------

			// Muon & proton kinetic variables

			double reco_Pmu_mcs = CandidateMu_P_Range->at(0);
			if (CandidateMu_EndContainment->at(0) == 0) { reco_Pmu_mcs = CandidateMu_P_MCS->at(0); }
			double reco_Pmu_cos_theta = CandidateMu_CosTheta->at(0);
			double reco_Pmu_phi = CandidateMu_Phi->at(0) * TMath::Pi() / 180.;
			double reco_Emu = TMath::Sqrt( reco_Pmu_mcs*reco_Pmu_mcs + MuonMass_GeV*MuonMass_GeV );
			
			TVector3 TVector3CandidateMuon(-1,-1,-1);
			TVector3CandidateMuon.SetMag(reco_Pmu_mcs);
			TVector3CandidateMuon.SetTheta(TMath::ACos(reco_Pmu_cos_theta));
			TVector3CandidateMuon.SetPhi(reco_Pmu_phi);						

			double reco_Pp = CandidateP_P_Range->at(0);
			double reco_Pp_cos_theta = CandidateP_CosTheta->at(0);
			double reco_Pp_phi = CandidateP_Phi->at(0) * TMath::Pi() / 180.;
			double reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );
			double reco_Tp = reco_Ep - ProtonMass_GeV;
			
			TVector3 TVector3CandidateProton(-1,-1,-1);
			TVector3CandidateProton.SetMag(reco_Pp);
			TVector3CandidateProton.SetTheta(TMath::ACos(reco_Pp_cos_theta));
			TVector3CandidateProton.SetPhi(reco_Pp_phi);	

			// --------------------------------------------------------------------------------------------------------------------------

			// Calorimetry

			double reco_Pmu_chi2 = CandidateMu_Chi2_YPlane->at(0);
			double reco_Pp_chi2 = CandidateP_Chi2_YPlane->at(0);

			double reco_Pmu_ThreePlanechi2 = CandidateMu_ThreePlaneChi2->at(0);
			double reco_Pp_ThreePlanechi2 = CandidateP_ThreePlaneChi2->at(0);

			double reco_Pmu_ThreePlaneLogLikelihood = CandidateMu_ThreePlaneLogLikelihood->at(0);
			double reco_Pp_ThreePlaneLogLikelihood = CandidateP_ThreePlaneLogLikelihood->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			double l_muCandidate = CandidateMu_Length->at(0);
			double l_pCandidate = CandidateP_Length->at(0);
			double LengthDiff = l_muCandidate - l_pCandidate;

			// -----------------------------------------------------------------------------------------------------------------------------

			// Miss quantities

			double kMiss = Reco_kMiss->at(0);

			double PMissMinus = Reco_PMissMinus->at(0);

			double MissMomentum = Reco_PMiss->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			// STV quantities
			
			double TransMissMomentum = Reco_Pt->at(0);

			double DeltaAlphaT = Reco_DeltaAlphaT->at(0);

			double DeltaPhiT = Reco_DeltaPhiT->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction

			double ECal = Reco_ECal->at(0);

			// QE Energy Reconstruction

			double EQE = Reco_EQE->at(0);

			// Reconstructed Q2

			double reco_Q2 = Reco_Q2->at(0);

			// -------------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Make sure that the same events fill the same plots

			if (reco_Pmu_mcs < ArrayNBinsMuonMomentum[0]) { continue; }
//			if (reco_Pmu_cos_theta < ArrayNBinsMuonCosTheta[0]) { continue; }
//			if (reco_Pmu_phi < ArrayNBinsMuonPhi[0]) { continue; }

			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }
//			if (reco_Pp_cos_theta < ArrayNBinsProtonCosTheta[0]) { continue; }
//			if (reco_Pp_phi < ArrayNBinsProtonPhi[0]) { continue; }

//			if (TransMissMomentum < ArrayNBinsDeltaPT[0]) { continue; }
//			if (DeltaAlphaT < ArrayNBinsDeltaAlphaT[0]) { continue; }
//			if (DeltaPhiT < ArrayNBinsDeltaPhiT[0]) { continue; }

//			if (ECal < ArrayNBinsECal[0]) { continue; }
//			if (EQE < ArrayNBinsEQE[0]) { continue; }
//			if (reco_Q2 < ArrayNBinsQ2[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

//			if (reco_Pmu_mcs > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
//			if (reco_Pmu_cos_theta > ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]) { continue; }
//			if (reco_Pmu_phi > ArrayNBinsMuonPhi[NBinsMuonPhi]) { continue; }

//			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }
//			if (reco_Pp_cos_theta > ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]) { continue; }
//			if (reco_Pp_phi > ArrayNBinsProtonPhi[NBinsProtonPhi]) { continue; }

//			if (TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { continue; }
//			if (DeltaAlphaT > ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]) { continue; }
//			if (DeltaPhiT > ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]) { continue; }

//			if (ECal > ArrayNBinsECal[NBinsECal]) { continue; }
//			if (EQE > ArrayNBinsEQE[NBinsEQE]) { continue; }
//			if (reco_Q2 > ArrayNBinsQ2[NBinsQ2]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

			// Default Plots

			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);

			if (NuScore > 0.6 && reco_Pp_ThreePlaneLogLikelihood > -1) {

				OldNEventsPassingSelectionCuts++;

				OldCutsRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);

			}

			// -------------------------------------------------------------------------------------------------------------------------

			// Studies for selection cuts

//			if (reco_Pp_ThreePlaneLogLikelihood < 0.) { continue; }
//			if (NuScore < 0.4) { continue; }
//			if (LengthDiff < 0) { continue; }
//			if (kMiss > 0.6) { continue; }
//			if (PMissMinus < 0.6) { continue; }

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

			// Relative angles

			double DeltaThetaProtonMuon_Deg = Reco_DeltaTheta->at(0);
			double DeltaPhiProtonMuon_Deg = Reco_DeltaPhi->at(0);

			// ------------------------------------------------------------------------------------------------------------------------

			// Optical info

			double NPE = BeamFlashes_TotalPE->at(0);

			TVector3 BeamFlash(0,BeamFlashes_YCenter->at(0),BeamFlashes_ZCenter->at(0));
			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));
			double dYZ = (BeamFlash - RecoVertex).Mag();

			double distance = CandidateMuP_Distance->at(0);	

			// -------------------------------------------------------------------------------------------------------------------------

			NEventsPassingSelectionCuts++;

			// --------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;

			double true_kMiss = -1;
			double true_PMissMinus = -1;
			double true_PMiss = -1;			
			double true_TransMissMomentum = -1;
			double true_DeltaAlphaT = -1;
			double true_DeltaPhiT = -1;
			double true_ECal = -1;
			double true_EQE = -1;
			double true_Q2 = -1;			
			
			if (
				string(fWhichSample).find("Overlay9") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

				true_kMiss = True_kMiss->at(0);
				true_PMissMinus = True_PMissMinus->at(0);
				true_PMiss = True_PMiss->at(0);				
				true_TransMissMomentum = True_Pt->at(0);
				true_DeltaAlphaT = True_DeltaAlphaT->at(0);
				true_DeltaPhiT = True_DeltaPhiT->at(0);
				true_ECal = True_ECal->at(0);
				true_EQE = True_EQE->at(0);
				true_Q2 = True_Q2->at(0);				
				
			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D LLP Only Study

			for (int WhichBin = 0; WhichBin < NBinsLLP; WhichBin++) {

				double LocalThres = MinLLP + LLPStep * WhichBin;

				if (reco_Pp_ThreePlaneLogLikelihood > LocalThres) {
					
					RecoMuonCosThetaPlotLLPOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D LLMu Only Study

			for (int WhichBin = 0; WhichBin < NBinsLLMu; WhichBin++) {

				double LocalThres = MinLLMu + LLMuStep * WhichBin;

				if (reco_Pmu_ThreePlaneLogLikelihood > LocalThres) {
					
					RecoMuonCosThetaPlotLLMuOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D NuScore Only Study

			for (int WhichBin = 0; WhichBin < NBinsNuScore; WhichBin++) {

				double LocalThres = MinNuScore + NuScoreStep * WhichBin;

				if (NuScore > LocalThres) {
					
					RecoMuonCosThetaPlotNuScoreOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D Length Only Study

			for (int WhichBin = 0; WhichBin < NBinsLength; WhichBin++) {

				double LocalThres = MinLength + LengthStep * WhichBin;

				if (LengthDiff > LocalThres) {
					
					RecoMuonCosThetaPlotLengthOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D PMissMinus Only Study

			for (int WhichBin = 0; WhichBin < NBinsPMissMinus; WhichBin++) {

				double LocalThres = MinPMissMinus + PMissMinusStep * WhichBin;

				if (PMissMinus > LocalThres) {
					
					RecoMuonCosThetaPlotPMissMinusOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D kMiss Only Study

			for (int WhichBin = 0; WhichBin < NBinskMiss; WhichBin++) {

				double LocalThres = MinkMiss + kMissStep * WhichBin;

				if (kMiss > LocalThres) {
					
					RecoMuonCosThetaPlotkMissOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D DeltaTheta Only Study

			for (int WhichBin = 0; WhichBin < NBinsDeltaTheta; WhichBin++) {

				double LocalThres = MinDeltaTheta + DeltaThetaStep * WhichBin;

				if ( TMath::Abs(DeltaThetaProtonMuon_Deg - 90.) < LocalThres) {
					
					RecoMuonCosThetaPlotDeltaThetaOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// 2D LL & nu score study

			for (int WhichNuScore = 0; WhichNuScore < NBinsNuScore; WhichNuScore++) {

				double LocalNuScoreThres = MinNuScore + NuScoreStep * WhichNuScore;

				for (int WhichLL = 0; WhichLL < NBinsLLP; WhichLL++) {

					double LocalLLThres = MinLLP + LLPStep * WhichLL;

					if (NuScore > LocalNuScoreThres && reco_Pp_ThreePlaneLogLikelihood > LocalLLThres) {
					
						RecoMuonCosThetaPlotStudy[WhichNuScore][WhichLL]->Fill(reco_Pmu_cos_theta,weight);

					}

				}

			}
		

			// ---------------------------------------------------------------------------------------------------------------------------

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

				// CC1p Signal

				if (CC1p == 1 && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg 
				    && True_CandidateMu_StartContainment->at(0) == 1) {

					// ---------------------------------------------------------------------------------------------------------------------------

					// Default plots
				
					CC1pEventsPassingSelectionCuts++;

					RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);

					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D LLP only study

					for (int WhichBin = 0; WhichBin < NBinsLLP; WhichBin++) {

						double LocalThres = MinLLP + LLPStep * WhichBin;

						if (reco_Pp_ThreePlaneLogLikelihood > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotLLPOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D LLMu only study

					for (int WhichBin = 0; WhichBin < NBinsLLMu; WhichBin++) {

						double LocalThres = MinLLMu + LLMuStep * WhichBin;

						if (reco_Pmu_ThreePlaneLogLikelihood > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotLLMuOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D NuScore only study

					for (int WhichBin = 0; WhichBin < NBinsNuScore; WhichBin++) {

						double LocalThres = MinNuScore + NuScoreStep * WhichBin;

						if (NuScore > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotNuScoreOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D Length only study

					for (int WhichBin = 0; WhichBin < NBinsLength; WhichBin++) {

						double LocalThres = MinLength + LengthStep * WhichBin;

						if (LengthDiff > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotLengthOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D PMissMinus only study

					for (int WhichBin = 0; WhichBin < NBinsPMissMinus; WhichBin++) {

						double LocalThres = MinPMissMinus + PMissMinusStep * WhichBin;

						if (PMissMinus > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotPMissMinusOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D kMiss only study

					for (int WhichBin = 0; WhichBin < NBinskMiss; WhichBin++) {

						double LocalThres = MinkMiss + kMissStep * WhichBin;

						if (kMiss > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotkMissOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ----------------------------------------------------------------------------------------------------------------------

					// 1D DeltaTheta Only Study

					for (int WhichBin = 0; WhichBin < NBinsDeltaTheta; WhichBin++) {

						double LocalThres = MinDeltaTheta + DeltaThetaStep * WhichBin;

						if ( TMath::Abs(DeltaThetaProtonMuon_Deg - 90.) < LocalThres) {
							
							CC1pRecoMuonCosThetaPlotDeltaThetaOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

					// 2D LL & nu score study

					for (int WhichNuScore = 0; WhichNuScore < NBinsNuScore; WhichNuScore++) {

						double LocalNuScoreThres = MinNuScore + NuScoreStep * WhichNuScore;

						for (int WhichLL = 0; WhichLL < NBinsLLP; WhichLL++) {

							double LocalLLThres = MinLLP + LLPStep * WhichLL;

							if (NuScore > LocalNuScoreThres && reco_Pp_ThreePlaneLogLikelihood > LocalLLThres) {
							
								CC1pRecoMuonCosThetaPlotStudy[WhichNuScore][WhichLL]->Fill(reco_Pmu_cos_theta,weight);

							}

						}

					}

					// ---------------------------------------------------------------------------------------------------------------------------

				}

				// --------------------------------------------------------------------------------------------------------------------------------

			} // End of the Overlay case and the breakdown into CC1p

			// -------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the events

		std::cout << std::endl << "Created a new file: " << FileName << std::endl << std::endl << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

		// -------------------------------------------------------------------------------------------------------------------------
		
		double nentriesError = sqrt(nentries);		
		
		std::cout << std::endl << "Number of " << fWhichSample << " initial entries = " << nentries << " +/- " << nentriesError << std::endl;		
		
		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed events passing the selection criteria

		double NEventsPassingSelectionCutsError = sqrt(NEventsPassingSelectionCuts);
		double OldNEventsPassingSelectionCutsError = sqrt(OldNEventsPassingSelectionCuts);

		std::cout << std::endl << "Number of events passing our \"no cut\" selection criteria = " << NEventsPassingSelectionCuts << " +/- " 
		<< NEventsPassingSelectionCutsError << std::endl;

		std::cout << std::endl << "Number of events passing our \"old cut\" selection criteria = " << OldNEventsPassingSelectionCuts << " +/- " 
		<< OldNEventsPassingSelectionCutsError << std::endl;

		cout << endl;

		file->cd();
		file->Write();
		file->Close();

		// -------------------------------------------------------------------------------------------------------------------------	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
