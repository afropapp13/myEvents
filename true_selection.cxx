#define true_selection_cxx
#include "true_selection.h"
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

#include "../myClasses/Tools.h"

//--------------------------------------------------//

TString TrueToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//--------------------------------------------------//

void true_selection::Loop() {

	//--------------------------------------------------//	

	if (fChain == 0) return; 
	Long64_t nentries = fChain->GetEntriesFast(); 
	Long64_t nbytes = 0, nb = 0;

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	//--------------------------------------------------//

	Tools tools;	

	//--------------------------------------------------//

	TString Extension = "";

	// For overlays only for genie, flux and reinteraction uncertainties

	if (fUniverseIndex != -1) {

		Extension = "_"+fEventWeightLabel+"_"+TrueToStringInt(fUniverseIndex); 

	}

	// Output Files

	TString FileName = PathToFiles+fTune+"TruthSTVAnalysis_"+fWhichSample+Extension+"_"+UBCodeVersion+".root";	
	TFile* OutputFile = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

	//--------------------------------------------------//

	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueMuonCosThetaSingleBinPlot[NInte];
	TH1D* TrueThetaZPlot[NInte];

	//--------------------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//--------------------------------------------------//

		// 1D analysis

		TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TrueMuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,0.,1.);
		TrueThetaZPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);

		//--------------------------------------------------//

	} // End of the loop over the interaction processes	

	//--------------------------------------------------//

	int TrueCC1pCounter = 0;
	int TrueCCQElikeCounter = 0;
	
	double SumWeights = 0.;

	//--------------------------------------------------//

	// Loop over the events

	cout << nentries << " events included in the file" << endl;	

	//--------------------------------------------------//	

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		//--------------------------------------------------//

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;	

		//--------------------------------------------------//

		// For detector variations, the eventweight weights are -1., set them back to 1.
		if (Weight == -1.) { Weight = 1.; }
		if (T2KWeight == -1.) { T2KWeight = 1.; }

		if (Weight <= 0 || Weight > 30) { continue; }
		if (T2KWeight <= 0 || T2KWeight > 30) { continue; }		
		// Weight from v3.0.4 to v.3.0.6 * weight from application of T2K tune
		double weight = POTWeight * Weight * T2KWeight * ROOTinoWeight;

		// Fake data studies: removing the T2K tune weight
		if (fTune == "GENIEv2") { weight = POTWeight; }
		if (fTune == "NoTune") { weight = POTWeight * Weight * ROOTinoWeight; }
		// Double the MEC weight  (mode = 10)
		if (fTune == "TwiceMEC" && Muon_MCParticle_Mode->at(0) == 10) { weight = 2 * POTWeight * Weight * T2KWeight * ROOTinoWeight; }
		if (fTune == "TwiceMEC" && Muon_MCParticle_Mode->at(0) != 10) { weight = POTWeight * Weight * T2KWeight * ROOTinoWeight; }		

		//--------------------------------------------------//

		// Genie, flux & reinteraction weights for systematics

		if ( 
			   fUniverseIndex != -1 && (
			   fWhichSample == "Overlay9_Run1" 
			|| fWhichSample == "Overlay9_Run2" 
			|| fWhichSample == "Overlay9_Run3" 
			|| fWhichSample == "Overlay9_Run4a"
			|| fWhichSample == "Overlay9_Run4b" 
			|| fWhichSample == "Overlay9_Run4c" 
			|| fWhichSample == "Overlay9_Run4d" 
			|| fWhichSample == "Overlay9_Run5" 
			|| fWhichSample == "Overlay9_Combined" 
			|| fWhichSample == "OverlayDirt9_Run1" 
			|| fWhichSample == "OverlayDirt9_Run2" 
			|| fWhichSample == "OverlayDirt9_Run3" 
			|| fWhichSample == "OverlayDirt9_Run4a" 
			|| fWhichSample == "OverlayDirt9_Run4b" 
			|| fWhichSample == "OverlayDirt9_Run4c" 
			|| fWhichSample == "OverlayDirt9_Run4d" 
			|| fWhichSample == "OverlayDirt9_Run5" 
			|| fWhichSample == "OverlayDirt9_Combined"

			) 
		) {

			// Genie weights

			if (fEventWeightLabel == "All_UBGenie") { weight = weight*All_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "AxFFCCQEshape_UBGenie") { weight = weight*AxFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "DecayAngMEC_UBGenie") { weight = weight*DecayAngMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "NormCCCOH_UBGenie") { weight = weight*NormCCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "NormNCCOH_UBGenie") { weight = weight*NormNCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "RPA_CCQE_UBGenie") { weight = weight*RPA_CCQE_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "ThetaDelta2NRad_UBGenie") { weight = weight*ThetaDelta2NRad_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "Theta_Delta2Npi_UBGenie") { weight = weight*Theta_Delta2Npi_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "VecFFCCQEshape_UBGenie") { weight = weight*VecFFCCQEshape_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "XSecShape_CCMEC_UBGenie") { weight = weight*XSecShape_CCMEC_UBGenie->at(fUniverseIndex)/ T2KWeight; }

			// Flux weights
			if (fEventWeightLabel == "fluxes") { weight = weight*fluxes->at(fUniverseIndex); }

			// Reinteraction weights
			if (fEventWeightLabel == "reinteractions") { weight = weight*reinteractions->at(fUniverseIndex); }		

			// MC_Stat weights // bootstrapping
			if (fEventWeightLabel == "MC_Stat") { 

				int concat = tools.ConcatRunSubRunEvent(Run,SubRun,Event,fUniverseIndex);
				weight = weight*tools.PoissonRandomNumber(concat); 
				
			}			

		}

		//----------------------------------------//

		SumWeights += weight / POTWeight;					

		//--------------------------------------------------//

		// Analysis over the simb::MCParticles

		std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
		std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();

		int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0;
		bool TrueCC1pEvent = false;
		bool TrueCCQElikeEvent = false;

		//--------------------------------------------------//

		// Signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & 0 pi+/- (Ppi > 70 MeV / c)

		if (CC1p == 1 && NumberPi0 == 0) {
		
			//--------------------------------------------------//		
			
			// Containment of the vertex has already been demanded in PreTruthSelection
			// Soft fiducial volume for true vertex

			//if (Muon_MCParticle_StartContainment->at(0) == 0) { continue; }

			//--------------------------------------------------//	

			// True muon

			double TrueMuonCosTheta = Muon_MCParticle_CosTheta->at(0);
			double TrueMuonTheta = TMath::ACos(TrueMuonCosTheta);
			double TrueMuonPhi_Deg = Muon_MCParticle_Phi->at(0);
			double TrueMuonPhi = TrueMuonPhi_Deg * TMath::Pi() / 180.;
			double TrueMuonMomentum_GeV = Muon_MCParticle_Mom->at(0); // GeV
			double TrueMuon_E_GeV = TMath::Sqrt( TMath::Power(TrueMuonMomentum_GeV,2.) + TMath::Power(MuonMass_GeV,2.) ); // GeV

			// True proton

			double TrueProtonCosTheta = Proton_MCParticle_CosTheta->at(0);
			double TrueProtonTheta = TMath::ACos(TrueProtonCosTheta);
			double TrueProtonPhi_Deg = Proton_MCParticle_Phi->at(0);
			double TrueProtonPhi = TrueProtonPhi_Deg * TMath::Pi() / 180.;
			double TrueProtonMomentum_GeV = Proton_MCParticle_Mom->at(0); // GeV
			double TrueProton_E_GeV = TMath::Sqrt( TMath::Power(TrueProtonMomentum_GeV,2.) + TMath::Power(ProtonMass_GeV,2.) ); // GeV

			double TrueDeltaThetaProtonMuon_Deg = True_DeltaTheta->at(0);
			double TrueThetaZ = True_ThetaZ->at(0);

                       	// Underflow / overflow
                        if (TrueThetaZ < ArrayNBinsThetaZ[0]) { TrueThetaZ = (ArrayNBinsThetaZ[0] + ArrayNBinsThetaZ[1])/2.; }
                        if (TrueThetaZ > ArrayNBinsThetaZ[NBinsThetaZ]) { TrueThetaZ = (ArrayNBinsThetaZ[NBinsThetaZ] + ArrayNBinsThetaZ[NBinsThetaZ-1])/2.; }
        
			// Reconstructed calorimetric energy using true level info / MCParticles

			double TrueRecoECal = True_ECal->at(0); // GeV

			// Transverse Variables

			double TrueTransMissMomentum = True_Pt->at(0);
			double TrueDeltaAlphaT = True_DeltaAlphaT->at(0);	
			double TrueDeltaAlpha3Dq = True_DeltaAlpha3Dq->at(0);
			double TrueDeltaPhiT = True_DeltaPhiT->at(0);
			double TrueDeltaPhi3D = True_DeltaPhi3D->at(0);			

			double TruePn = True_Pn->at(0);

			//--------------------------------------------------//	

			// True Vertex

			TVector3 TrueVertex(True_Vx,True_Vy,True_Vz);
				
			//--------------------------------------------------//		

			// Demand that the true muon / proton start points and the true proton end point are contained
			// Angle selection cuts: collinearity to reject broken tracks
			// Momentum threshold
			// Same events fill all the plots

			if (
			    /*TrueMuonStartContainment == true 
			    && TrueProtonStartContainment == true 
			    &&*/ TrueMuonMomentum_GeV > ArrayNBinsMuonMomentum[0]
			    && TrueProtonMomentum_GeV > ArrayNBinsProtonMomentum[0]
			) {

				//--------------------------------------------------//

				// STV analysis

				if (
				    TrueMuonMomentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
				    && TrueProtonMomentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				) {

					//----------------------------------------//

					int genie_mode = -1;

					if (Muon_MCParticle_Mode->at(0) == 0) { genie_mode = 1; }
					else if (Muon_MCParticle_Mode->at(0) == 10) { genie_mode = 2; }
					else if (Muon_MCParticle_Mode->at(0) == 1) { genie_mode = 3; }
					else if (Muon_MCParticle_Mode->at(0) == 2) { genie_mode = 4; }
					else { genie_mode = 5; }																				

					//----------------------------------------//

					// True CC1p event

					TrueCC1pEvent = true;
					TrueCC1pCounter++;

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(TrueMuonMomentum_GeV,2.) );
					double true_Nu = True_Ev - true_MuonEnergy;

					// Playground for CC1p true momenta (longitudinal & perpendicular) ratios

					TVector3 TrueCandidateMuon(1,1,1);
					TrueCandidateMuon.SetMag(TrueMuonMomentum_GeV);
					TrueCandidateMuon.SetPhi(TrueMuonPhi);
					TrueCandidateMuon.SetTheta(TMath::ACos(TrueMuonCosTheta));

					TVector3 TrueCandidateProton(1,1,1);
					TrueCandidateProton.SetMag(TrueProtonMomentum_GeV);
					TrueCandidateProton.SetPhi(TrueProtonPhi);
					TrueCandidateProton.SetTheta(TMath::ACos(TrueProtonCosTheta));	

					//----------------------------------------//	

					// 1D analysis		

					TrueMuonCosThetaPlot[0]->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot[0]->Fill(0.5,weight);
					TrueThetaZPlot[0]->Fill(TrueThetaZ,weight);
					
					TrueMuonCosThetaPlot[genie_mode]->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot[genie_mode]->Fill(0.5,weight);
					TrueThetaZPlot[genie_mode]->Fill(TrueThetaZ,weight);
					
					//----------------------------------------//									

				} // End of the event selection

				// -----------------------------------------------------------------------------------------------------------------

			} // End of the angle selection cuts && the demand that we fill the plots with the same events

		} // End of signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 300 MeV / c) & pi+/- (Ppi > 70 MeV / c)

		// -------------------------------------------------------------------------------------------------------------------------

	} // End of the loop over the events

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// STV-CC1p analysis summary

	std::cout << std::endl;

	// -------------------------------------------------------------------------------------------------------------------------

	if ( string(fWhichSample).find("Overlay9") != std::string::npos ) {

		//std::cout << std::endl << "Samdef events = " << SamdefEvents << std::endl;
		std::cout << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

	}

	//----------------------------------------//

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();

	fFile->Close();

	//----------------------------------------//

} // End of the program
