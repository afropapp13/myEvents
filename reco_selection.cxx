#define reco_selection_cxx
#include "reco_selection.h"
#include <TH2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TF1.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include "../myClasses/Tools.h"
#include "../myClasses/STV_Tools.h"

using namespace std;

//----------------------------------------//

TString ToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

void reco_selection::Loop() {

	//----------------------------------------//

	// Function for proton momentum recalibration in %

	TF1* fPP = new TF1("fPP","29.354172*x-14.674918",0.3,0.5);	

	//----------------------------------------//

	int NEventsPassingSelectionCuts = 0;
	int CC1pEventsPassingSelectionCuts = 0;
	int NonCC1pEventsPassingSelectionCuts = 0;
	
	//----------------------------------------//

	Tools tools;			

	//----------------------------------------//
		
	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	VectorCuts.push_back("");
	VectorCuts.push_back("_PID_NuScore");
	VectorCuts.push_back("_CRT");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		Cuts = Cuts + VectorCuts[i];

		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

		//----------------------------------------//

		if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
		TH1D::SetDefaultSumw2();
		TH2D::SetDefaultSumw2();
		double weight = 1.;

		//----------------------------------------//

		TString Extension = "";

		// For overlays only for genie, flux and reinteraction uncertainties

		if (fUniverseIndex != -1) {

			Extension = "_"+fEventWeightLabel+"_"+ToStringInt(fUniverseIndex); 

		}

		TString FileName = PathToFiles+Cuts+"/"+fTune+"STVStudies_"+fWhichSample+Extension+Cuts+".root";
		TFile* file = new TFile(FileName,"recreate");
		std::cout << std::endl << "Creating a new file: " << FileName << std::endl << std::endl << std::endl;

		//----------------------------------------//

		// Muon CosTheta

		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);	
		TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoMuonCosThetaPlot2D = new TH2D("CC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		
		//----------------------------------------//

		// Muon CosTheta Single Bin

		TH1D* RecoMuonCosThetaSingleBinPlot = new TH1D("RecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CC1pRecoMuonCosThetaSingleBinPlot = new TH1D("CC1pRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CC1pTrueMuonCosThetaSingleBinPlot = new TH1D("CC1pTrueMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH2D* CC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("CC1pRecoMuonCosThetaSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH2D* POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH1D* NonCC1pRecoMuonCosThetaSingleBinPlot = new TH1D("NonCC1pRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCQERecoMuonCosThetaSingleBinPlot = new TH1D("CCQERecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCMECRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCRESRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCDISRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISRecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		//----------------------------------------//

		// ThetaZ

		TH1D* RecoThetaZPlot = new TH1D("RecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CC1pRecoThetaZPlot = new TH1D("CC1pRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);	
		TH1D* CC1pTrueThetaZPlot = new TH1D("CC1pTrueThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH2D* CC1pRecoThetaZPlot2D = new TH2D("CC1pRecoThetaZPlot2D",LabelXAxisThetaZ2D,NBinsThetaZ,
			ArrayNBinsThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH2D* POTScaledCC1pRecoThetaZPlot2D = new TH2D("POTScaledCC1pRecoThetaZPlot2D",LabelXAxisThetaZ2D,NBinsThetaZ,
			ArrayNBinsThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* NonCC1pRecoThetaZPlot = new TH1D("NonCC1pRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCQERecoThetaZPlot = new TH1D("CCQERecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCMECRecoThetaZPlot = new TH1D("CCMECRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCRESRecoThetaZPlot = new TH1D("CCRESRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCDISRecoThetaZPlot = new TH1D("CCDISRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
	
		//----------------------------------------//

		TH1D* CC1pThetaZDiff_ECalSlicesPlot[TwoDNBinsECal];	
		TH1D* CC1pThetaZReso_ECalSlicesPlot[TwoDNBinsECal];	

		// Loop over the ECal slices
		
		for (int iecal = 0; iecal < TwoDNBinsECal; iecal++) {

			CC1pThetaZDiff_ECalSlicesPlot[iecal] = new TH1D("CC1pThetaZDiff_ECalSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,-30,30);
			CC1pThetaZReso_ECalSlicesPlot[iecal] = new TH1D("CC1pThetaZReso_ECalSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,-100,100);
		
		} // End of the loop over ECal slices

		//----------------------------------------//

		TH1D* CC1pThetaZDiff_MuonMomentumSlicesPlot[TwoDNBinsMuonMomentum];	
		TH1D* CC1pThetaZReso_MuonMomentumSlicesPlot[TwoDNBinsMuonMomentum];	

		// Loop over the MuonMomentum slices
		
		for (int iecal = 0; iecal < TwoDNBinsMuonMomentum; iecal++) {

			CC1pThetaZDiff_MuonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZDiff_MuonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,-30,30);
			CC1pThetaZReso_MuonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZReso_MuonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,-100,100);
		
		} // End of the loop over MuonMomentum slices
	
		//----------------------------------------//

		TH1D* CC1pThetaZDiff_ProtonMomentumSlicesPlot[TwoDNBinsProtonMomentum];	
		TH1D* CC1pThetaZReso_ProtonMomentumSlicesPlot[TwoDNBinsProtonMomentum];	

		// Loop over the ProtonMomentum slices
		
		for (int iecal = 0; iecal < TwoDNBinsProtonMomentum; iecal++) {

			CC1pThetaZDiff_ProtonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZDiff_ProtonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,-30,30);
			CC1pThetaZReso_ProtonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZReso_ProtonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,-100,100);
		
		} // End of the loop over ProtonMomentum slices
	
		//----------------------------------------//
		//----------------------------------------//

		// Loop over the events

		cout << nentries << " events included in the file" << endl;

		for (Long64_t jentry=0; jentry<nentries;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

			// ---------------------------------------------------------------------------------------------------------------------

			if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;

			// ------------------------------------------------------------------------------------------------------------------------

			// Demand for exactly one muon/proton candidate

			if (CandidateMu_P_MCS->size() != 1) { continue; }
			if (CandidateP_P_Range->size() != 1) { continue; }

			// ------------------------------------------------------------------------------------------------------------------------

			// Containment demand
			// Fully contained muon / proton candidates
			// Just to be sure, that has already been done at a preselection level 

			if (CandidateMu_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_StartContainment->at(0) == 0) { continue; }
			if (CandidateMu_EndContainment->at(0) == 0) { continue; }
			if (CandidateP_EndContainment->at(0) == 0) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------

			weight = POTWeight;

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

				// For detector variations, the eventweight weights are -1., set them back to 1.
				if (Weight == -1.) { Weight = 1.; }
				if (T2KWeight == -1.) { T2KWeight = 1.; }
			
				if (Weight <= 0 || Weight > 30) { continue; } // bug fix weight 
				if (T2KWeight <= 0 || T2KWeight > 30) { continue; }	// T2K tune weight	
				// For the detector variations, Weight (bug fix) = 1		
				weight = POTWeight * Weight * T2KWeight * ROOTinoWeight; 

				// Fake data studies: removing the T2K tune weight
				if (fTune == "GENIEv2") { weight = POTWeight; }
				if (fTune == "NoTune") { weight = POTWeight * Weight * ROOTinoWeight; }
				// Double the MEC weight (mode = 10)
				if (fTune == "TwiceMEC" && MCParticle_Mode == 10) { weight = 2 * POTWeight * Weight * T2KWeight * ROOTinoWeight; }
				if (fTune == "TwiceMEC" && MCParticle_Mode != 10) { weight = POTWeight * Weight * T2KWeight * ROOTinoWeight; }									
				
			}
			
			// --------------------------------------------------------------------------------------------------------------------------------

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

				// Watch out: The EventWeight weights already include the weight for the tune

				// Genie weights
				
				if (fUniverseIndex < (int)(All_UBGenie->size())) {

					if (fEventWeightLabel == "All_UBGenie") { weight = weight*All_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "AxFFCCQEshape_UBGenie") { weight = weight*AxFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "DecayAngMEC_UBGenie") { weight = weight*DecayAngMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "NormCCCOH_UBGenie") { weight = weight*NormCCCOH_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "NormNCCOH_UBGenie") { weight = weight*NormNCCOH_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "RPA_CCQE_UBGenie") { weight = weight*RPA_CCQE_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "ThetaDelta2NRad_UBGenie") { weight = weight*ThetaDelta2NRad_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "Theta_Delta2Npi_UBGenie") { weight = weight*Theta_Delta2Npi_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "VecFFCCQEshape_UBGenie") { weight = weight*VecFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
					if (fEventWeightLabel == "XSecShape_CCMEC_UBGenie") { weight = weight*XSecShape_CCMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }

				}
                                else {

                                	cout <<  ", run = " << Run << " subrun = " << SubRun << " event = " << Event << endl;

                                }

				// Flux weights
				if (fEventWeightLabel == "fluxes") { 

					if ( fUniverseIndex < (int)(fluxes->size()) ) {

						weight = weight*fluxes->at(fUniverseIndex); 

					}
					else {

						cout << ", run = " << Run << " subrun = " << SubRun << " event = " << Event << endl;

					}

				}

				// Reinteraction weights
				if (fEventWeightLabel == "reinteractions") { 
				
					if ( fUniverseIndex < (int)(reinteractions->size()) ) {

						weight = weight*reinteractions->at(fUniverseIndex); 

					}

				}

				// MC_Stat weights // bootstrapping
				if (fEventWeightLabel == "MC_Stat") { 

					int concat = tools.ConcatRunSubRunEvent(Run,SubRun,Event,fUniverseIndex);
					weight = weight*tools.PoissonRandomNumber(concat); 
				
				}				

			}	

			// -----------------------------------------------------------------------------------------------------------------------

			if ( fabs(weight) != weight) { continue; } // Securing against infinities

			// -----------------------------------------------------------------------------------------------------------------------------

			// Contained Reconstructed Vertex

			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));

			if ( !tools.inFVVector(RecoVertex) ) { continue; }

			// -----------------------------------------------------------------------------------------------------------------------------

			// Muon info

			double reco_Pmu = CandidateMu_P_Range->at(0);

			// Quality cut suggested by Anne Schukraft // May 7 2021
			// For contained tracks, we have two methods for the muon momentum reconstruction
			// MCS & range, we can require that the two are withing 25%
			// to reject misreconstructed tracks 
	
			if (CandidateMu_EndContainment->at(0) == 1) { 

				double Reso =  TMath::Abs(CandidateMu_P_MCS->at(0) - CandidateMu_P_Range->at(0) ) / CandidateMu_P_Range->at(0) ; 
				if (Reso > MuRangeMCSAgreeValue) { continue; }

			}

			// -------------------------------------------------------------------------------------------------------------------------

			// Ensure that we don't have flipped tracks 
			// by demanding that muon start - vertex distance < muon end - vertex distance,
			// that proton start - vertex distance < proton end - vertex distance

			if (CandidateMuStartVertexDistance->at(0) > CandidateMuEndVertexDistance->at(0)) { continue; }
			if (CandidatePStartVertexDistance->at(0) > CandidatePEndVertexDistance->at(0)) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------			

			TVector3 CandidateMuonStart(CandidateMu_StartX->at(0),CandidateMu_StartY->at(0),CandidateMu_StartZ->at(0));
			TVector3 CandidateMuonEnd(CandidateMu_EndX->at(0),CandidateMu_EndY->at(0),CandidateMu_EndZ->at(0));

			TVector3 CandidateProtonStart(CandidateP_StartX->at(0),CandidateP_StartY->at(0),CandidateP_StartZ->at(0));
			TVector3 CandidateProtonEnd(CandidateP_EndX->at(0),CandidateP_EndY->at(0),CandidateP_EndZ->at(0));

			double DistanceStartPoints = StartToStartDistance->at(0);
			double DistanceEndPoints = EndToEndDistance->at(0); 

			if (DistanceStartPoints > DistanceEndPoints) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------

			// Muon info

			double reco_Pmu_cos_theta = CandidateMu_CosTheta->at(0);
			double reco_Pmu_phi = CandidateMu_Phi->at(0); // deg
			double reco_Emu = TMath::Sqrt( reco_Pmu*reco_Pmu + MuonMass_GeV*MuonMass_GeV );
			
			TVector3 TVector3CandidateMuon(-1,-1,-1);
			TVector3CandidateMuon.SetMag(reco_Pmu);
			TVector3CandidateMuon.SetTheta(TMath::ACos(reco_Pmu_cos_theta));
			TVector3CandidateMuon.SetPhi(reco_Pmu_phi * TMath::Pi() / 180.);	

			// --------------------------------------------------------------------------------------------------------------------------

			// Proton info

			double reco_Pp = CandidateP_P_Range->at(0);
			double reco_Pp_cos_theta = CandidateP_CosTheta->at(0);
			double reco_Pp_phi = CandidateP_Phi->at(0); // deg 
			double reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );

			TVector3 TVector3CandidateProton(-1,-1,-1);
			TVector3CandidateProton.SetMag(reco_Pp);
			TVector3CandidateProton.SetTheta(TMath::ACos(reco_Pp_cos_theta));
			TVector3CandidateProton.SetPhi(reco_Pp_phi * TMath::Pi() / 180.);	

			// --------------------------------------------------------------------------------------------------------------------------

			// Calorimetry

			double reco_p_LLR_Score = CandidateP_LLR_PID->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction
	
			double ECal = Reco_ECal->at(0);
			
			double ThetaZ = Reco_ThetaZ->at(0);

			// -------------------------------------------------------------------------------------------------------------------------
			
			STV_Tools reco_stv_tool(TVector3CandidateMuon,TVector3CandidateProton,reco_Emu,reco_Ep);
	
			// Redefinition if P_p < 0.5 where the biases have been observed

			if (reco_Pp < 0.5) { 

				reco_Pp = ( 1.-0.01*fPP->Eval(reco_Pp) ) * reco_Pp ;
				TVector3CandidateProton.SetMag(reco_Pp);
				reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );

				ECal = reco_stv_tool.ReturnECal();
				ThetaZ = reco_stv_tool.ReturnThetaZ();

			}

			// Underflow / overflow
			if (ThetaZ < ArrayNBinsThetaZ[0]) { ThetaZ = (ArrayNBinsThetaZ[0] + ArrayNBinsThetaZ[1])/2.; }
			if (ThetaZ > ArrayNBinsThetaZ[NBinsThetaZ]) { ThetaZ = (ArrayNBinsThetaZ[NBinsThetaZ] + ArrayNBinsThetaZ[NBinsThetaZ-1])/2.; }
	
			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

				if (VectorCuts[i] == "_PID_NuScore" && !(reco_p_LLR_Score < ProtonLLRPIDScore) ) 
					{ PassedSelection = false; }

				if (VectorCuts[i] == "_PID_NuScore_CRT" && !(crtveto == 0) ) 
					{ PassedSelection = false; }


			}

			if (PassedSelection == false) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Make sure that the same events fill the same plots

			if (reco_Pmu < ArrayNBinsMuonMomentum[0]) { continue; }
			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

			if (reco_Pmu > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }

			NEventsPassingSelectionCuts++;

			// --------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;

			double true_ECal = -1;
			double true_ThetaZ = -1;

			// Only for MC to obtain true vales			
			
			if (
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

				true_ECal = True_ECal->at(0);
				true_ThetaZ = True_ThetaZ->at(0);
	
                        	// Underflow / overflow
                                if (true_ThetaZ < ArrayNBinsThetaZ[0]) { true_ThetaZ = (ArrayNBinsThetaZ[0] + ArrayNBinsThetaZ[1])/2.; }
                                if (true_ThetaZ > ArrayNBinsThetaZ[NBinsThetaZ]) { true_ThetaZ = (ArrayNBinsThetaZ[NBinsThetaZ] + ArrayNBinsThetaZ[NBinsThetaZ-1])/2.; }
                        	
			} // End of if statement: Only for MC to obtain true vales

			//----------------------------------------//

			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			RecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
			RecoThetaZPlot->Fill(ThetaZ,weight);

			//------------------------------//

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

				bool isSignal = false;

				int inte_mode = -1;
				if (genie_mode == 0) { inte_mode = 1; } // QE 
				else if (genie_mode == 10) { inte_mode = 2; } // MEC 
				else if (genie_mode == 1) { inte_mode = 3; } // RES 
				else if (genie_mode == 2) { inte_mode = 4; } // DIS 
				else { inte_mode = 5; } // COH / Other 

				// CC1p Signal

				if (    CC1p == 1 && NumberPi0 == 0 && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg 
				     && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg 
				     && True_CandidateMu_StartContainment->at(0) == 1
 
				     && True_CandidateMu_P->at(0) > ArrayNBinsMuonMomentum[0] 
				     && True_CandidateP_P->at(0) > ArrayNBinsProtonMomentum[0]

				     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
				     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]

				) {

					// --------------------------------------------------------------------------------------------------
				
					CC1pEventsPassingSelectionCuts++;
					isSignal = true;

					// ---------------------------------------------------------------------------------------------------------------------------
					// ---------------------------------------------------------------------------------------------------------------------------
					// 1D Plots using True level info for selected CC1p events  

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(True_CandidateMu_P->at(0),2.) );
					double true_Nu = True_Ev - true_MuonEnergy;

					CC1pTrueMuonCosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
					CC1pTrueMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CC1pTrueThetaZPlot->Fill(true_ThetaZ,weight);

					//----------------------------------------//

					// 1D Reco Plots for the selected CC1p events 

					CC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CC1pRecoThetaZPlot->Fill(ThetaZ,weight);
					
					//------------------------------//

					CC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoMuonCosThetaSingleBinPlot2D->Fill(0.5,0.5);
					CC1pRecoThetaZPlot2D->Fill(true_ThetaZ,ThetaZ);

					POTScaledCC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
					POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D->Fill(0.5,0.5,weight);
					POTScaledCC1pRecoThetaZPlot2D->Fill(true_ThetaZ,ThetaZ,weight);
					
					//-----------------------------------------------------------------------------

					// Atmospherics
				
					double diff = ThetaZ - true_ThetaZ; // deg
					double reso = diff / true_ThetaZ * 100.; // deg

					// ECal slices

					int ECalTwoDIndex = tools.ReturnIndex(ECal, TwoDArrayNBinsECal);

					CC1pThetaZDiff_ECalSlicesPlot[ECalTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_ECalSlicesPlot[ECalTwoDIndex]->Fill(reso,weight);

					// MuonMomentum slices

					int MuonMomentumTwoDIndex = tools.ReturnIndex(reco_Pmu, TwoDArrayNBinsMuonMomentum);

					CC1pThetaZDiff_MuonMomentumSlicesPlot[MuonMomentumTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_MuonMomentumSlicesPlot[MuonMomentumTwoDIndex]->Fill(reso,weight);

					// ProtonMomentum slices

					int ProtonMomentumTwoDIndex = tools.ReturnIndex(reco_Pp, TwoDArrayNBinsProtonMomentum);

					CC1pThetaZDiff_ProtonMomentumSlicesPlot[ProtonMomentumTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_ProtonMomentumSlicesPlot[ProtonMomentumTwoDIndex]->Fill(reso,weight);


				} // End of the CC1p signal

				// -------------------------------------------------------------------------------------------------------------

				// Non-CC1p beam related background or EXT BNB

				else {

					NonCC1pEventsPassingSelectionCuts++;

					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					NonCC1pRecoThetaZPlot->Fill(ThetaZ,weight);

					//------------------------------//

				} // End of the Non-CC1p beam related background

				// -------------------------------------------------------------------------------------------------------------------------
				// ------------------------------------------------------------------------------------------------------------------------

				// CCQE

				if (genie_mode == 0) {

					CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCQERecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCQERecoThetaZPlot->Fill(ThetaZ,weight);

				} // End of CCQE selection

				// ------------------------------------------------------------------------------------------------------------------------

				// CCMEC

				if (genie_mode == 10) {

					CCMECRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCMECRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCMECRecoThetaZPlot->Fill(ThetaZ,weight);
			
				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCRES

				if (genie_mode == 1) {

					CCRESRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCRESRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCRESRecoThetaZPlot->Fill(ThetaZ,weight);
	
				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCDIS

				if (genie_mode == 2) {

					CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCDISRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCDISRecoThetaZPlot->Fill(ThetaZ,weight);

				}

				// --------------------------------------------------------------------------------------------------------------------------

				// Overlay particle breakdown using the Backtracker for PID studies

				if ( !(CandidateMu_MCParticle_Pdg->size() > 0 && CandidateP_MCParticle_Pdg->size() > 0) ) { 

					cout << "muon candidate true pdg = " << CandidateMu_MCParticle_Pdg->at(0) << endl;
					cout << "proton candidate true pdg = " << CandidateP_MCParticle_Pdg->at(0) << endl;

				}

			} // End of the Overlay case and the breakdown into CC1p/NonCC1p & QE,MEC,RES,DIS

			// -------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the events

		std::cout << std::endl << "Created file: " << FileName << std::endl << std::endl;

		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

		//----------------------------------------//

		cout << fWhichSample << "  " << NEventsPassingSelectionCuts << " events passing selection criteria" << endl << endl;

		//----------------------------------------//	

		file->cd();
		file->Write();
		file->Close();
		fFile->Close();

		//----------------------------------------//	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
