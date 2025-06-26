#define mcc9_10_true_selection_cxx
#include "mcc9_10_true_selection.h"
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
#include <TRandom.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include "../../../generators/Tools.h"

//--------------------------------------------------//

TString TrueToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//--------------------------------------------------//

void mcc9_10_true_selection::Loop() {

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

	//--------------------------------------------------//
	
	// file for fake data study

	TFile* f_fds = new TFile("reinseghal/rs_spline.root","readonly");
	TGraph* g_fds = (TGraph*)(f_fds->Get("h_spline"));

	//--------------------------------------------------//	

	// Output Files

	TString FileName = event_selection_file_path+fTune+"Truthncpi0_"+fWhichSample+Extension+".root";	
	TFile* OutputFile = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

	//--------------------------------------------------//

	TH1D* TruePi0CosThetaPlot[NInte];
	TH1D* TrueSingleBinPlot[NInte];
	TH1D* TruePi0MomentumPlot[NInte];

	//--------------------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		TruePi0CosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TrueSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSingleBinPlot","",1,0.,1.);
		TruePi0MomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
	
	} // End of the loop over the interaction processes	

	//--------------------------------------------------//

	// Loop over the events

	cout << nentries << " events included in the file" << endl;	

	//--------------------------------------------------//	

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		//--------------------------------------------------//

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;	

		//--------------------------------------------------//

		// For detector variations runs 1-3, the eventweight weights are -1., set them back to 1.
		if (Weight == -1.) { Weight = 1.; }
		if (T2KWeight == -1.) { T2KWeight = 1.; }
			
		// For detector variations runs 4-5, the event weights are NOT -1
		// However setting the weights to 1 for consistency
		/*if (
			string(fWhichSample).find("CV") != std::string::npos || 
			string(fWhichSample).find("CVextra") != std::string::npos || 
			string(fWhichSample).find("LYDown") != std::string::npos || 
			string(fWhichSample).find("LYRayleigh") != std::string::npos || 
			string(fWhichSample).find("LYAttenuation") != std::string::npos || 
			string(fWhichSample).find("SCE") != std::string::npos || 
			string(fWhichSample).find("Recombination2") != std::string::npos || 
			string(fWhichSample).find("X") != std::string::npos || 
			string(fWhichSample).find("YZ") != std::string::npos || 
			string(fWhichSample).find("ThetaXZ") != std::string::npos || 
			string(fWhichSample).find("ThetaYZ") != std::string::npos
		) {

			Weight = 1.;
			T2KWeight = 1.;
				
		}*/

		// Set some limits to make sure that the weights are not negative or unreasonable / infinity
		if (Weight <= 0 || Weight > 30) { continue; }
		if (T2KWeight <= 0 || T2KWeight > 30) { continue; }		
		// Weight from v3.0.4 to v.3.0.6 * weight from application of T2K tune
		double weight = POTWeight * Weight * T2KWeight * ROOTinoWeight;

		// Fake data studies: reweight to Rein Sehgal (RS)
		if (fTune == "RS" && pi0_MCParticle_Mode->at(0) == 3) { 
				
			double rw = g_fds->Eval(True_Ev);
			weight = rw * weight; 
				
		}		

		//--------------------------------------------------//

		// Genie, flux & reinteraction weights for multisim systematics
		// !!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!
		// divide all the weights by 1000
		// PeLEE choice to stor integers instead of doubles	
		// Not applicable to MCStat weights		

		if ( 
			   fUniverseIndex != -1 && (
			   fWhichSample == "mcc9_10_Overlay9_Run1" 
			|| fWhichSample == "mcc9_10_Overlay9_Run2" 
			|| fWhichSample == "mcc9_10_Overlay9_Run3" 
			|| fWhichSample == "mcc9_10_Overlay9_Run4a"
			|| fWhichSample == "mcc9_10_Overlay9_Run4b" 
			|| fWhichSample == "mcc9_10_Overlay9_Run4b_unified" 
			|| fWhichSample == "mcc9_10_Overlay9_Run4b_standalone" 
			|| fWhichSample == "mcc9_10_Overlay9_Run4c" 
			|| fWhichSample == "mcc9_10_Overlay9_Run4d" 
			|| fWhichSample == "mcc9_10_Overlay9_Run5" 
			|| fWhichSample == "mcc9_10_Overlay9_Combined" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run1" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run2" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run3" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run4a" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run4b" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run4b_unified" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run4b_standalone"			
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run4c" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run4d" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Run5" 
			|| fWhichSample == "mcc9_10_OverlayDirt9_Combined"

			) 
		) {

			// Genie weights

			if (fEventWeightLabel == "All_UBGenie") { 

				if ( int(All_UBGenie->size()) > fUniverseIndex ) {

					weight = weight*All_UBGenie->at(fUniverseIndex) / T2KWeight / 1000.; 

				}

			}

			if (fEventWeightLabel == "AxFFCCQEshape_UBGenie") { 

				if ( int(AxFFCCQEshape_UBGenie->size()) > fUniverseIndex ) {
				
					weight = weight*AxFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; 

				}

			}

			if (fEventWeightLabel == "DecayAngMEC_UBGenie") { 

				if ( int(DecayAngMEC_UBGenie->size()) > fUniverseIndex ) {
				
					weight = weight*DecayAngMEC_UBGenie->at(fUniverseIndex) / T2KWeight;
 
				}

			}

			if (fEventWeightLabel == "NormCCCOH_UBGenie") { 

				if ( int(NormCCCOH_UBGenie->size()) > fUniverseIndex ) {

					weight = weight*NormCCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; 

				}

			}

			if (fEventWeightLabel == "NormNCCOH_UBGenie") { 

				if ( int(NormNCCOH_UBGenie->size()) > fUniverseIndex ) {

					weight = weight*NormNCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; 
				
				}

			}

			if (fEventWeightLabel == "RPA_CCQE_UBGenie") { 
	
				if ( int(RPA_CCQE_UBGenie->size()) > fUniverseIndex ) {

					weight = weight*RPA_CCQE_UBGenie->at(fUniverseIndex)/ T2KWeight; 

				}

			}

			if (fEventWeightLabel == "ThetaDelta2NRad_UBGenie") { 

				if ( int(ThetaDelta2NRad_UBGenie->size()) > fUniverseIndex ) {
				
					weight = weight*ThetaDelta2NRad_UBGenie->at(fUniverseIndex)/ T2KWeight; 

				}

			}

			if (fEventWeightLabel == "Theta_Delta2Npi_UBGenie") { 

				if ( int(Theta_Delta2Npi_UBGenie->size()) > fUniverseIndex ) {
				
					weight = weight*Theta_Delta2Npi_UBGenie->at(fUniverseIndex)/ T2KWeight; 

				}

			}

			if (fEventWeightLabel == "VecFFCCQEshape_UBGenie") { 
	
				if ( int(VecFFCCQEshape_UBGenie->size()) > fUniverseIndex ) {

					weight = weight*VecFFCCQEshape_UBGenie->at(fUniverseIndex)/ T2KWeight; 

				}

			}

			if (fEventWeightLabel == "XSecShape_CCMEC_UBGenie") { 

				if ( int(XSecShape_CCMEC_UBGenie->size()) > fUniverseIndex ) {

					weight = weight*XSecShape_CCMEC_UBGenie->at(fUniverseIndex)/ T2KWeight; 

				}

			}

			// Flux weights
			if (fEventWeightLabel == "fluxes") {

				if ( int(fluxes->size()) > fUniverseIndex ) {

					 weight = weight*fluxes->at(fUniverseIndex) / 1000.; 

				}

			}

			// Reinteraction weights
			if (fEventWeightLabel == "reinteractions") { 

				if ( int(reinteractions->size()) > fUniverseIndex ) {
				
					weight = weight*reinteractions->at(fUniverseIndex) / 1000.;

				} 

			}		

			// MC_Stat weights // bootstrapping
			if (fEventWeightLabel == "MC_Stat") { 

				int concat = tools.ConcatRunSubRunEvent(Run,SubRun,Event,fUniverseIndex);
				weight = weight*tools.PoissonRandomNumber(concat); 
				
			}						

		} // end of the multisim weights			

		//--------------------------------------------------//

		// Analysis over the simb::MCParticles

		std::vector<int> VectorTruePi0Index; VectorTruePi0Index.clear();

		//--------------------------------------------------//

		// Signal definition: 1 pi0 of any momentum with costheta > 0.5
		// No neutrons above 20 MeV KE
		// No protons above 20 MeV KE
		// No other heavier mesons or baryons

		if (signal) {
		
			//--------------------------------------------------//	

			// True Vertex

			TVector3 TrueVertex(True_Vx,True_Vy,True_Vz);
			if ( !tools.inFVVector(TrueVertex) ) { continue; }
	
			//----------------------------------------//

			int genie_mode = -1;

			if (pi0_MCParticle_Mode->at(0) == 0) { genie_mode = 1; } // qe
			else if (pi0_MCParticle_Mode->at(0) == 10) { genie_mode = 2; } // mec
			else if (pi0_MCParticle_Mode->at(0) == 1) { genie_mode = 3; } // res
			else if (pi0_MCParticle_Mode->at(0) == 2) { genie_mode = 4; } // dis
			else if (pi0_MCParticle_Mode->at(0) == 3) { genie_mode = 5; } // coh
			else { genie_mode = 6; } // other												
	
			//----------------------------------------//

			TVector3 TrueCandidatePi0(1,1,1);
			TrueCandidatePi0.SetMag(pi0_MCParticle_Mom->at(0));
			TrueCandidatePi0.SetPhi(pi0_MCParticle_Phi->at(0)*TMath::Pi()/180.);
			TrueCandidatePi0.SetTheta(TMath::ACos(pi0_MCParticle_CosTheta->at(0)));
	
			//----------------------------------------//	

			// all processes

			TrueSingleBinPlot[0]->Fill(0.5,weight);
			TruePi0CosThetaPlot[0]->Fill(pi0_MCParticle_CosTheta->at(0),weight);
			TruePi0MomentumPlot[0]->Fill(pi0_MCParticle_Mom->at(0),weight);

			// specific processes

			TrueSingleBinPlot[genie_mode]->Fill(0.5,weight);
			TruePi0CosThetaPlot[genie_mode]->Fill(pi0_MCParticle_CosTheta->at(0),weight);
			TruePi0MomentumPlot[genie_mode]->Fill(pi0_MCParticle_Mom->at(0),weight);
					
			//----------------------------------------//									


		} // End of the event selection

		//----------------------------------------//

	} // End of the loop over the events

	//----------------------------------------//

	std::cout << std::endl;

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();

	fFile->Close();
	f_fds->Close();

	//----------------------------------------//

} // End of the program
