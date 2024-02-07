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

	// Txt file to keep track of the event reduction at each stage

	TString TxtName = "/uboone/data/users/apapadop/myEvents/my3DTxtFiles/"+UBCodeVersion+"/TxtmyTrueEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);

	// Txt file to keep track of the run/subrun/event of the candidate events

	TString RunTxtName = "/uboone/data/users/apapadop/myEvents/my3DTxtFiles/"+UBCodeVersion+"/TxtmyTrueRunSubRunEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
	ofstream myRunTxtFile;
	myRunTxtFile.open(RunTxtName);
	myRunTxtFile << std::fixed << std::setprecision(2);
	myRunTxtFile << fWhichSample;

	//--------------------------------------------------//

	// 1D analysis

	TH1D* TrueDeltaPTPlot[NInte];
	TH1D* TrueDeltaAlphaTPlot[NInte];
	TH1D* TrueDeltaAlpha3DqPlot[NInte];
	TH1D* TrueDeltaPnPlot[NInte];
	TH1D* TrueMuonMomentumPlot[NInte];
	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueMuonCosThetaSingleBinPlot[NInte];
	TH1D* TrueProtonMomentumPlot[NInte];
	TH1D* TrueProtonCosThetaPlot[NInte];
	TH1D* TrueECalPlot[NInte];

	//--------------------------------------------------//	

	// 3D analysis (uncorrelated)	

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
	TH1D* TrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[NInte][TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
	TH1D* TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[NInte][TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
	TH1D* TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[NInte][TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];	

	//--------------------------------------------------//	

	// 3D analysis in 1D grid		
	
	TH1D* SerialTrueECal_InDeltaPTDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPnDeltaAlpha3DqPlot[NInte];
	TH1D* SerialTrueECal_InMuonCosThetaMuonMomentumPlot[NInte];
	TH1D* SerialTrueECal_InProtonCosThetaProtonMomentumPlot[NInte];

	//--------------------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//--------------------------------------------------//

		// 1D analysis

		TrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TrueDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);
		TrueDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TrueMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TrueMuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,0.,1.);
		TrueProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TrueECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		//--------------------------------------------------//	

		// 3D analysis (uncorrelated)

		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

			for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

				TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[inte][WhichDeltaPT][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

			}

		}

		for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

			for (int WhichDeltaAlpha3Dq = 0; WhichDeltaAlpha3Dq < TwoDNBinsDeltaAlpha3Dq; WhichDeltaAlpha3Dq++) {	

				TString ECalTwoDInDeltaPnDeltaAlpha3DqLabel = "ECal_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"_DeltaAlpha3Dq_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq+1])+"Plot";
				TrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[inte][WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);

			}

		}


		for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

			for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

				TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";
				TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[inte][WhichMuonCosTheta][WhichMuonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);

			}

		}	

		for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

			for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

				TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
				TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[inte][WhichProtonCosTheta][WhichProtonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);

			}

		}

		//--------------------------------------------------//	

		// 3D analysis in 1D grid		
		
		SerialTrueECal_InDeltaPTDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		SerialTrueECal_InDeltaPnDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		SerialTrueECal_InMuonCosThetaMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		SerialTrueECal_InProtonCosThetaProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);

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

		// For detector variations runs 1-3, the eventweight weights are -1., set them back to 1.
		if (Weight == -1.) { Weight = 1.; }
		if (T2KWeight == -1.) { T2KWeight = 1.; }
			
		// For detector variations runs 4-5, the event weights are NOT -1
		// However setting the weights to 1 for consistency
		if (
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
				
		}


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

			// Overflow bins
			// Affects

			// DeltaPT
			// DeltaPn
			// ECal

			if (TrueTransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { TrueTransMissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
			if (TruePn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { TruePn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }
			if (TrueRecoECal > ArrayNBinsECal[NBinsECal]) { TrueRecoECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }

			//--------------------------------------------------//

			// Underflow bins
			// Affects

			// ECal
			
			if (TrueRecoECal < ArrayNBinsECal[0]) { TrueRecoECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			
			//--------------------------------------------------//
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
				    TrueTransMissMomentum > ArrayNBinsDeltaPT[0] 
				    //&& TrueTransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				    && TrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] 
				    && TrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				    
				    && TrueMuonMomentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
				    && TrueProtonMomentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				    
				    && TrueMuonCosTheta > ArrayNBinsMuonCosTheta[0]
				    && TrueMuonCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
				    && TrueProtonCosTheta > ArrayNBinsProtonCosTheta[0]
				    && TrueProtonCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
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

					// Indices for 2D & 3D analysis

					int DeltaPTTwoDIndex = tools.ReturnIndex(TrueTransMissMomentum, TwoDArrayNBinsDeltaPT);
					int DeltaPnTwoDIndex = tools.ReturnIndex(TruePn, TwoDArrayNBinsDeltaPn);					
					int DeltaAlphaTTwoDIndex = tools.ReturnIndex(TrueDeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaAlpha3DqTwoDIndex = tools.ReturnIndex(TrueDeltaAlpha3Dq, TwoDArrayNBinsDeltaAlpha3Dq);
					int MuonCosThetaTwoDIndex = tools.ReturnIndex(TrueMuonCosTheta, TwoDArrayNBinsMuonCosTheta);
					int ProtonCosThetaTwoDIndex = tools.ReturnIndex(TrueProtonCosTheta, TwoDArrayNBinsProtonCosTheta);
					int MuonMomentumTwoDIndex = tools.ReturnIndex(TrueMuonMomentum_GeV, TwoDArrayNBinsMuonMomentum);
					int ProtonMomentumTwoDIndex = tools.ReturnIndex(TrueProtonMomentum_GeV, TwoDArrayNBinsProtonMomentum);



					int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,TrueRecoECal);
					int SerialECalInDeltaPnDeltaAlpha3DqIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices,DeltaPnTwoDIndex,DeltaAlpha3DqTwoDIndex,TrueRecoECal);
					int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,TrueRecoECal);
					int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,TrueRecoECal);																							

					//----------------------------------------//	

					// 1D analysis		

					TrueDeltaPTPlot[0]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot[0]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlpha3DqPlot[0]->Fill(TrueDeltaAlpha3Dq,weight);
					TrueDeltaPnPlot[0]->Fill(TruePn,weight);
					TrueECalPlot[0]->Fill(TrueRecoECal,weight);
					TrueMuonMomentumPlot[0]->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonCosThetaPlot[0]->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot[0]->Fill(0.5,weight);
					TrueProtonMomentumPlot[0]->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonCosThetaPlot[0]->Fill(TrueProtonCosTheta,weight);

					TrueDeltaPTPlot[genie_mode]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot[genie_mode]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlpha3DqPlot[genie_mode]->Fill(TrueDeltaAlpha3Dq,weight);
					TrueDeltaPnPlot[genie_mode]->Fill(TruePn,weight);
					TrueECalPlot[genie_mode]->Fill(TrueRecoECal,weight);
					TrueMuonMomentumPlot[genie_mode]->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonCosThetaPlot[genie_mode]->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot[genie_mode]->Fill(0.5,weight);
					TrueProtonMomentumPlot[genie_mode]->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonCosThetaPlot[genie_mode]->Fill(TrueProtonCosTheta,weight);					
					//----------------------------------------//										

					// 3D analysis (uncorrelated)

					TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[0][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[0][DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[0][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[0][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);

					TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[genie_mode][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[genie_mode][DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[genie_mode][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);					

					//----------------------------------------//						

					// 3D analysis treated in 1D grid

					SerialTrueECal_InDeltaPTDeltaAlphaTPlot[0]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialTrueECal_InDeltaPnDeltaAlpha3DqPlot[0]->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialTrueECal_InMuonCosThetaMuonMomentumPlot[0]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialTrueECal_InProtonCosThetaProtonMomentumPlot[0]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);

					SerialTrueECal_InDeltaPTDeltaAlphaTPlot[genie_mode]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialTrueECal_InDeltaPnDeltaAlpha3DqPlot[genie_mode]->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialTrueECal_InMuonCosThetaMuonMomentumPlot[genie_mode]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialTrueECal_InProtonCosThetaProtonMomentumPlot[genie_mode]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);					

					//----------------------------------------//									
				} // End of the event selection

				// -----------------------------------------------------------------------------------------------------------------

			} // End of the angle selection cuts && the demand that we fill the plots with the same events

		} // End of signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 300 MeV / c) & pi+/- (Ppi > 70 MeV / c)

		// -------------------------------------------------------------------------------------------------------------------------

		// Storing the run/subrun/event of the candidate events 

		myRunTxtFile << "Run = " << Run << ", SubRun = " << SubRun << ", Event = " << Event ;

	} // End of the loop over the events

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// STV-CC1p analysis summary

	std::cout << std::endl;

	// -------------------------------------------------------------------------------------------------------------------------

	if ( string(fWhichSample).find("Overlay9") != std::string::npos ) {

		//std::cout << std::endl << "Samdef events = " << SamdefEvents << std::endl;
		std::cout << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

		//myTxtFile << "Samdef events = " << SamdefEvents << std::endl;
		myTxtFile << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

	}

	//----------------------------------------//

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();
	myTxtFile.close();

	fFile->Close();

	//----------------------------------------//

} // End of the program
