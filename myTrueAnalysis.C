#define myTrueAnalysis_cxx
#include "myTrueAnalysis.h"
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

#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Tools.h"

using namespace std;

void myTrueAnalysis::Loop() {

	if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
	TH1D::SetDefaultSumw2();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Output Files

	TString FileName = "./OutputFiles/"+UBCodeVersion+"/TruthSTVAnalysis_"+fWhichSample+"_"+UBCodeVersion+".root";
	TFile* OutputFile = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

	ofstream TxtFile;
	TxtFile.open ("./TxtFiles/"+UBCodeVersion+"/TruthSTVAnalysis_"+fWhichSample+"_"+UBCodeVersion+".txt");

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// 1D True Transverse Variables

	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueEQEPlot = new TH1D("TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	Tools tools;

	// -------------------------------------------------------------------------------------------------------------------------------------------

	int TrueCC1pCounter = 0;
	int TrueCCQElikeCounter = 0;
	
	// --------------------------------------------------------------------------------------------------------------------------------

	// POT Counting

	double POTCount = -99.;

	if (string(fWhichSample).find("Overlay") != std::string::npos) {

		// Locally
		//TString PathToPOTFile = "mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
		// gpvm's
		TString PathToPOTFile = "mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
		TFile* POTFile = TFile::Open(PathToPOTFile,"readonly");
		TH1D* POTCountHist = (TH1D*)(POTFile->Get("POTCountHist"));
		POTCount = POTCountHist->GetBinContent(1);
		POTFile->Close();

	}				

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//	for (Long64_t jentry=0; jentry<50000;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;
		
		// ------------------------------------------------------------------------------------------------------------------

		// POT Scaling


		double tor860_wcut = 1;
		double E1DCNT_wcut = 1.;
		double EXT = 1.;

		if (string(fWhichSample).find("Run1") != std::string::npos) {

			tor860_wcut = tor860_wcut_Run1;
			E1DCNT_wcut = E1DCNT_wcut_Run1;
			EXT = EXT_Run1;

		}			

		// ------------------------------------------------------------------------------------------------------------------------------

//		double weight = 1.;
//		double T2Kweight = 1.;
		double weight = ( tor860_wcut / POTCount) * Weight * T2KWeight; // Weight from v3.0.4 to v.3.0.6 * weight from application of T2K tune
		if (weight <= 0 || weight > 10) { continue; }

		// -------------------------------------------------------------------------------------------------------------------------------------

		// EventWeight weights

// Fix it !!!!

		// Genie
//		if (string(fWhichSample).find("Genie_All") != std::string::npos) { weight = weight * EventWeightValues[fWhichSample][fUniverse];}

		// Flux
//		if (string(fWhichSample).find("FluxUnisim") != std::string::npos || string(fWhichSample).find("Primary") != std::string::npos) 
//			{ weight = weight * EventWeightValues[fWhichSample][fUniverse];}

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Analysis over the simb::MCParticles

		std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
		std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();

		int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0;
		bool TrueCC1pEvent = false;
		bool TrueCCQElikeEvent = false;

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Loop over the MCParticles and determine the populations

		for (int WhichMCParticle = 0; WhichMCParticle < NumberMCParticles; WhichMCParticle++) {
		
			double MCParticleMomentum = MCParticle_Mom->at(WhichMCParticle);
			int MCParticlePdg = MCParticle_Pdg->at(WhichMCParticle);

			if (MCParticlePdg == MuonPdg && MCParticleMomentum > ArrayNBinsMuonMomentum[0]) { 
				TrueMuonCounter++; 
				VectorTrueMuonIndex.push_back(WhichMCParticle);
			}

			if (MCParticlePdg == ProtonPdg && MCParticleMomentum > ArrayNBinsProtonMomentum[0]) { 
				TrueProtonCounter++; 
				VectorTrueProtonIndex.push_back(WhichMCParticle);
			}

			if ( fabs(MCParticlePdg) == AbsChargedPionPdg && MCParticleMomentum > ChargedPionMomentumThres) 
				{ TrueChargedPionCounter++; }

		} // end of the loop over the simb::MCParticles

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & 0 pi+/- (Ppi > 70 MeV / c)

		if (/*TrueMuonCounter == 1 && TrueProtonCounter == 1 && TrueChargedPionCounter == 0 &&*/ CC1p == 1) {

			// True muon

			double TrueMuonCosTheta = MCParticle_CosTheta->at(VectorTrueMuonIndex.at(0));
			double TrueMuonTheta = TMath::ACos(TrueMuonCosTheta);
			double TrueMuonPhi_Deg = MCParticle_Phi->at(VectorTrueMuonIndex.at(0));
			double TrueMuonPhi = TrueMuonPhi_Deg * TMath::Pi() / 180.;
			double TrueMuonMomentum_GeV = MCParticle_Mom->at(VectorTrueMuonIndex.at(0)); // GeV
			double TrueMuonMomentum_MeV = 1000. * TrueMuonMomentum_GeV; // MeV
			double TrueMuon_KE_MeV = tools.PToKE(MuonPdg,TrueMuonMomentum_MeV); // MeV
			double TrueMuon_KE_GeV = TrueMuon_KE_MeV / 1000.; // GeV
			double TrueMuon_E_GeV = TrueMuon_KE_GeV + MuonMass_GeV; // GeV
			int TrueMuonStartContainment = MCParticle_StartContainment->at(VectorTrueMuonIndex.at(0));
			int TrueMuonEndContainment = MCParticle_EndContainment->at(VectorTrueMuonIndex.at(0));

			TVector3 TVector3TrueMuon;
			TVector3TrueMuon.SetMagThetaPhi(TrueMuonMomentum_GeV,TrueMuonTheta,TrueMuonPhi);
			TLorentzVector TrueMuon4V(TVector3TrueMuon,TrueMuon_E_GeV);
			TVector3 TVector3TrueMuonTransverse(TVector3TrueMuon.X(),TVector3TrueMuon.Y(),0);
			double TVector3TrueMuonTransverseMag = TVector3TrueMuonTransverse.Mag();

			// True proton

			double TrueProtonCosTheta = MCParticle_CosTheta->at(VectorTrueProtonIndex.at(0));
			double TrueProtonTheta = TMath::ACos(TrueProtonCosTheta);
			double TrueProtonPhi_Deg = MCParticle_Phi->at(VectorTrueProtonIndex.at(0));
			double TrueProtonPhi = TrueProtonPhi_Deg * TMath::Pi() / 180.;
			double TrueProtonMomentum_GeV = MCParticle_Mom->at(VectorTrueProtonIndex.at(0)); // GeV
			double TrueProtonMomentum_MeV = 1000. * TrueProtonMomentum_GeV; // MeV
			double TrueProton_KE_MeV = tools.PToKE(ProtonPdg,TrueProtonMomentum_MeV); // MeV
			double TrueProton_KE_GeV = TrueProton_KE_MeV / 1000.; // GeV
			double TrueProton_E_GeV = TrueProton_KE_GeV + ProtonMass_GeV; // GeV
			int TrueProtonStartContainment = MCParticle_StartContainment->at(VectorTrueProtonIndex.at(0));
			int TrueProtonEndContainment = MCParticle_EndContainment->at(VectorTrueProtonIndex.at(0));

			TVector3 TVector3TrueProton;
			TVector3TrueProton.SetMagThetaPhi(TrueProtonMomentum_GeV,TrueProtonTheta,TrueProtonPhi);
			TLorentzVector TrueProton4V(TVector3TrueProton,TrueProton_E_GeV);
			TVector3 TVector3TrueProtonTransverse(TVector3TrueProton.X(),TVector3TrueProton.Y(),0);
			double TVector3TrueProtonTransverseMag = TVector3TrueProtonTransverse.Mag();

			double TrueDeltaPhiProtonMuon = TVector3TrueMuon.DeltaPhi(TVector3TrueProton);
			double TrueDeltaThetaProtonMuon = TVector3TrueMuon.Angle(TVector3TrueProton);
			double TrueDeltaPhiProtonMuon_Deg = TrueDeltaPhiProtonMuon * 180. / TMath::Pi();
			double TrueDeltaThetaProtonMuon_Deg = TrueDeltaThetaProtonMuon * 180. / TMath::Pi();
			if (TrueDeltaThetaProtonMuon_Deg < 0.) { TrueDeltaPhiProtonMuon_Deg += 180.; }
			if (TrueDeltaThetaProtonMuon_Deg > 180.) { TrueDeltaThetaProtonMuon_Deg -= 180.; }
			if (TrueDeltaPhiProtonMuon_Deg < 0.) { TrueDeltaPhiProtonMuon_Deg += 360.; }
			if (TrueDeltaPhiProtonMuon_Deg > 360.) { TrueDeltaPhiProtonMuon_Deg -= 360.; }

			// Reconstructed calorimetric energy using true level info / MCParticles

			double TrueRecoECal = TrueProton_KE_GeV + TrueMuon_E_GeV + BE; // GeV

			// Reconstructed QE energy using true level info / MCParticles

			double TrueRecoEQE = ( ProtonMass_GeV * BE + ProtonMass_GeV * TrueMuon_E_GeV ) / ( ProtonMass_GeV - TrueMuon_E_GeV + TrueMuonMomentum_GeV * TrueMuonCosTheta);
		
			// Reconstructed neutrino using true level info / MCParticles

			TVector3 RecoTrueNeutrinoZaxis(0,0,TrueRecoECal);
			TLorentzVector RecoTrueNeutrinoZaxis4V(0,0,TrueRecoECal,TrueRecoECal);

			TVector3 RecoTrueNeutrino = TVector3TrueMuon + TVector3TrueProton;
			double InitialRecoTrueNeutrinoMag = RecoTrueNeutrino.Mag();

			// Definition of Transverse Variables

			TVector3 TVector3TrueTransMissMomentum(RecoTrueNeutrino.X(),RecoTrueNeutrino.Y(),0);
			double TrueTransMissMomentum = TVector3TrueTransMissMomentum.Mag();

			double TrueDeltaAlphaT = 
				TMath::ACos( (- TVector3TrueMuonTransverse*TVector3TrueTransMissMomentum) / ( TVector3TrueMuonTransverseMag*TrueTransMissMomentum ) ) * 180./TMath::Pi();
			double TrueDeltaPhiT = 
			TMath::ACos( (- TVector3TrueMuonTransverse*TVector3TrueProtonTransverse) / ( TVector3TrueMuonTransverseMag*TVector3TrueProtonTransverseMag ) ) * 180./TMath::Pi();

			// Reconstructed q vector using true level info / MCParticles

			TLorentzVector Trueq4V = RecoTrueNeutrinoZaxis4V - TrueMuon4V;
			double RecoTrueQ2 = -Trueq4V.Mag2();

			// Demand that the true muon / proton start points and the true proton end point are contained
			// Angle selection cuts: collinearity to reject broken tracks
			// Momentum threshold
			// Same events fill all the plots

			if (
			    TrueMuonStartContainment == true 
			    && TrueProtonStartContainment == true && //TrueProtonEndContainment == true
			    /*&&*/ /*TrueDeltaThetaProtonMuon_Deg  < DeltaThetaCut
			    &&*/ TrueMuonMomentum_GeV > ArrayNBinsMuonMomentum[0]
			    && TrueProtonMomentum_GeV > ArrayNBinsProtonMomentum[0]
			) {

				// -----------------------------------------------------------------------------------------------------------------

				// STV analysis

				if (
				    TrueTransMissMomentum > ArrayNBinsDeltaPT[0] 
				    && TrueTransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				    && TrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] 
				    && TrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				    && TrueDeltaPhiT > ArrayNBinsDeltaPhiT[0] 
				    && TrueDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
				    
				    && TrueMuonMomentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
				    && TrueProtonMomentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				    
				    && TrueMuonCosTheta > ArrayNBinsMuonCosTheta[0]
				    && TrueMuonCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
				    && TrueProtonCosTheta > ArrayNBinsProtonCosTheta[0]
				    && TrueProtonCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
				    
				    && TrueRecoECal > ArrayNBinsECal[0]
				    && TrueRecoEQE > ArrayNBinsEQE[0]
				    && RecoTrueQ2 > ArrayNBinsQ2[0]
				    && TrueRecoECal < ArrayNBinsECal[NBinsECal]
				    && TrueRecoEQE < ArrayNBinsEQE[NBinsEQE]
				    && RecoTrueQ2 < ArrayNBinsQ2[NBinsQ2]
				) {

					// True CC1p event

					TrueCC1pEvent = true;
					TrueCC1pCounter++;

					// STV

					TrueDeltaPTPlot->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);

					// Reconstructed energy & Q2

					TrueECalPlot->Fill(TrueRecoECal,weight);
					TrueEQEPlot->Fill(TrueRecoEQE,weight);
					TrueQ2Plot->Fill(RecoTrueQ2,weight);

					// Kinematic variables

					TrueMuonMomentumPlot->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonPhiPlot->Fill(TrueMuonPhi_Deg,weight);
					TrueMuonCosThetaPlot->Fill(TrueMuonCosTheta,weight);

					TrueProtonMomentumPlot->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonPhiPlot->Fill(TrueProtonPhi_Deg,weight);
					TrueProtonCosThetaPlot->Fill(TrueProtonCosTheta,weight);
				}

				// -----------------------------------------------------------------------------------------------------------------

			} // End of the angle selection cuts && the demand that we fill the plots with the same events

		} // End of signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & pi+/- (Ppi > 70 MeV / c)

	} // End of the loop over the events

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// STV-CC1p analysis summary

	std::cout << std::endl;

	if ( string(fWhichSample).find("Overlay9") != std::string::npos ) {

		std::cout << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;
		TxtFile << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

	}

	// --------------------------------------------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();

} // End of the program
