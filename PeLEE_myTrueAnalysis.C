#define PeLEE_myTrueAnalysis_cxx
#include "PeLEE_myTrueAnalysis.h"
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

//--------------------------------------------------//

TString TrueToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//--------------------------------------------------//

void PeLEE_myTrueAnalysis::Loop() {

	//--------------------------------------------------//	

	if (fChain == 0) return; 
	Long64_t nentries = fChain->GetEntriesFast(); 
	Long64_t nbytes = 0, nb = 0;

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

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

	TString TxtName = "/uboone/data/users/apapadop/myEvents/myTxtFiles/"+UBCodeVersion+"/TxtmyTrueEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);

	// Txt file to keep track of the run/subrun/event of the candidate events

	TString RunTxtName = "/uboone/data/users/apapadop/myEvents/myTxtFiles/"+UBCodeVersion+"/TxtmyTrueRunSubRunEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
	ofstream myRunTxtFile;
	myRunTxtFile.open(RunTxtName);
	myRunTxtFile << std::fixed << std::setprecision(2);
	myRunTxtFile << fWhichSample;

	//--------------------------------------------------//

	// 1D True Variables

	TH1D* TruePi0Plot = new TH1D("TruePi0Plot",";# #pi^{0}",4,-0.5,3.5);
	TH1D* TrueNeutronPlot = new TH1D("TrueNeutronPlot",";# neutrons",6,-0.5,5.5);
	
	TH1D* TruekMissPlot = new TH1D("TruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
	TH1D* TruePMissMinusPlot = new TH1D("TruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
	TH1D* TruePMissPlot = new TH1D("TruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);	

	TH1D* TrueDeltaPLPlot = new TH1D("TrueDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
	TH1D* TrueDeltaPnPlot = new TH1D("TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
	TH1D* TrueDeltaPtxPlot = new TH1D("TrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
	TH1D* TrueDeltaPtyPlot = new TH1D("TrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
	TH1D* TrueAPlot = new TH1D("TrueAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	TH1D* TrueMuonCosThetaSingleBinPlot = new TH1D("TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

	TH1D* TrueContainedEndMuonCosThetaPlot = new TH1D("TrueContainedEndMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueEQEPlot = new TH1D("TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	TH1D* TrueCCQEMuonCosThetaPlot = new TH1D("TrueCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
	TH1D* TrueCCQEMuonMomentumPlot = new TH1D("TrueCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
	TH1D* TrueCCQEMuonPhiPlot = new TH1D("TrueCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

	TH1D* TrueCCQEProtonCosThetaPlot = new TH1D("TrueCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
	TH1D* TrueCCQEProtonMomentumPlot = new TH1D("TrueCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
	TH1D* TrueCCQEProtonPhiPlot = new TH1D("TrueCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

	TH1D* TrueCCQEECalPlot = new TH1D("TrueCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
	TH1D* TrueCCQEQ2Plot = new TH1D("TrueCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);	

	TH1D* TrueVertexXPlot = new TH1D("TrueVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
	TH1D* TrueVertexYPlot = new TH1D("TrueVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
	TH1D* TrueVertexZPlot = new TH1D("TrueVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

	TH1D* TrueEvPlot = new TH1D("TrueEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
	TH1D* TrueNuPlot = new TH1D("TrueNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);

	//--------------------------------------------------//

	Tools tools;	

	//--------------------------------------------------//

	// Now let's get serious with the 2D analysis

	// Ecal in DeltaPT and DeltaAlphaT bins		

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];	

	//--------------------------------------------------//

	// DeltaAlphaT in DeltaPT slices

	TH1D* TrueDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//

	// ECal in DeltaPT slices

	TH1D* TrueECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//

	// DeltaPhiT in DeltaPT slices

	TH1D* TrueDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//

	// DeltaPn in DeltaPT slices

	TH1D* TrueDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//	

	for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

		//------------------------------//

		// DeltaAlphaT in DeltaPT slices

		TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);

		//------------------------------//

		// ECal in DeltaPT slices

		TString ECalTwoDInDeltaPTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);		

		//------------------------------//

		// DeltaPhiT in DeltaPT slices

		TString DeltaPhiTTwoDInDeltaPTLabel = "DeltaPhiT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	

		//------------------------------//

		// DeltaPn in DeltaPT slices

		TString DeltaPnTwoDInDeltaPTLabel = "DeltaPn_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	

		//------------------------------//

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

			TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
			TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

		}

	}

	//--------------------------------------------------//

	// Ecal in MuonCosTheta and MuonMomentum bins		

	TH1D* TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];	

	// DeltaAlphaT in muon cos theta slices

	TH1D* TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];

	// DeltaPT in muon cos theta slices

	TH1D* TrueDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];	

	// Pmu in muon cos theta slices

	TH1D* TrueMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];

	// proton cos theta in muon cos theta slices

	TH1D* TrueProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];	

	for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

		TString DeltaAlphaTTwoDInMuonCosThetaLabel = "DeltaAlphaT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);

		TString DeltaPTTwoDInMuonCosThetaLabel = "DeltaPT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);		

		TString MuonMomentumTwoDInMuonCosThetaLabel = "MuonMomentum_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);

		TString ProtonCosThetaTwoDInMuonCosThetaLabel = "ProtonCosTheta_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);		

		for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

			TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";
			TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("True"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);

		}

	}	

	//--------------------------------------------------//

	// Ecal in ProtonCosTheta and ProtonMomentum bins		

	TH1D* TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];	

	// Pp in proton cos theta slices

	TH1D* TrueProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];

	// DeltaAlphaT in proton cos theta slices

	TH1D* TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];	

	// DeltaPT in proton cos theta slices

	TH1D* TrueDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];		

	for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

		TString ProtonMomentumTwoDInProtonCosThetaLabel = "ProtonMomentum_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);	

		TString DeltaAlphaTTwoDInProtonCosThetaLabel = "DeltaAlphaT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

		TString DeltaPTTwoDInProtonCosThetaLabel = "DeltaPT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

		for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

			TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
			TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("True"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);

		}

	}

	//--------------------------------------------------//

	// DeltaPty in DeltaPtx slices

	TH1D* TrueDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];

	for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

		TString DeltaPtyTwoDInDeltaPtxLabel = "DeltaPty_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			
		TrueDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("True"+DeltaPtyTwoDInDeltaPtxLabel,LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);	

	}	

	//--------------------------------------------------//

	// DeltaPtx in DeltaPty slices

	TH1D* TrueDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];

	for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

		TString DeltaPtxTwoDInDeltaPtyLabel = "DeltaPtx_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";			
		TrueDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("True"+DeltaPtxTwoDInDeltaPtyLabel,LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);	

	}		

	//--------------------------------------------------//

	// DeltaPT in DeltaAlphaT slices

	TH1D* TrueDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];

	// ECal in DeltaAlphaT slices

	TH1D* TrueECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];	

	for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

		TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
		TrueDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("True"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

		TString ECalTwoDInDeltaAlphaTLabel = "ECal_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
		TrueECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("True"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

	}

	//--------------------------------------------------//

	// 2D now let's try to do it right with correlations

	// DeltaPT in cos theta mu bins
	TH1D* SerialTrueDeltaPT_InMuonCosThetaPlot = new TH1D("TrueSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);

	// DeltaPT in cos theta p bins
	TH1D* SerialTrueDeltaPT_InProtonCosThetaPlot = new TH1D("TrueSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);	

	// DeltaPT in DeltaAlphaT bins
	TH1D* SerialTrueDeltaPT_InDeltaAlphaTPlot = new TH1D("TrueSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);

	// Pmu in cos theta mu bins
	TH1D* SerialTrueMuonMomentum_InMuonCosThetaPlot = new TH1D("TrueSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);

	// Pp in cos theta p bins
	TH1D* SerialTrueProtonMomentum_InProtonCosThetaPlot = new TH1D("TrueSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);	

	// DeltaAlphaT in cos theta mu bins
	TH1D* SerialTrueDeltaAlphaT_InMuonCosThetaPlot = new TH1D("TrueSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);

	// DeltaAlphaT in cos theta p bins
	TH1D* SerialTrueDeltaAlphaT_InProtonCosThetaPlot = new TH1D("TrueSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);

	// DeltaAlphaT in DeltaPT bins
	TH1D* SerialTrueDeltaAlphaT_InDeltaPTPlot = new TH1D("TrueSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);		

	// DeltaPhiT in DeltaPT bins
	TH1D* SerialTrueDeltaPhiT_InDeltaPTPlot = new TH1D("TrueSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);

	// DeltaPn in DeltaPT bins
	TH1D* SerialTrueDeltaPn_InDeltaPTPlot = new TH1D("TrueSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);

	// ECal in DeltaAlphaT bins
	TH1D* SerialTrueECal_InDeltaAlphaTPlot = new TH1D("TrueSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);

	// ECal in DeltaPT bins
	TH1D* SerialTrueECal_InDeltaPTPlot = new TH1D("TrueSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);	

	// cos theta p in cos theta mu bins
	TH1D* SerialTrueProtonCosTheta_InMuonCosThetaPlot = new TH1D("TrueSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);

	// DeltaPty in DeltaPtx bins
	TH1D* SerialTrueDeltaPty_InDeltaPtxPlot = new TH1D("TrueSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);

	// DeltaPtx in DeltaPty bins
	TH1D* SerialTrueDeltaPtx_InDeltaPtyPlot = new TH1D("TrueSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);

	// ECal in DeltaPT & DeltaAlphaT bins
	TH1D* SerialTrueECal_InDeltaPTDeltaAlphaTPlot = new TH1D("TrueSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);

	// ECal in MuonCosTheta & MuonMomentum bins
	TH1D* SerialTrueECal_InMuonCosThetaMuonMomentumPlot = new TH1D("TrueSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);

	// ECal in ProtonCosTheta & ProtonMomentum bins
	TH1D* SerialTrueECal_InProtonCosThetaProtonMomentumPlot = new TH1D("TrueSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);

	//--------------------------------------------------//	
	
	// 2D Analysis

	// For now and until box opening
	int NBins2DAnalysis = 4;
		
	TH2D* TrueCosThetaMuPmuPlot = new TH2D("TrueCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
	TH2D* TrueCosThetaPPpPlot = new TH2D("TrueCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);

	//--------------------------------------------------//

	TH1D* TrueMuonTrueMomentumLongitudinalRatio = new TH1D("TrueMuonTrueMomentumLongitudinalRatio",";P^{true}_{#mu,||}/P^{true}_{#mu}",25,0.,1.);		
	TH1D* TrueProtonTrueMomentumLongitudinalRatio = new TH1D("TrueProtonTrueMomentumLongitudinalRatio",";P^{true}_{p,||}/P^{true}_{p}",25,0.,1.);

	TH1D* TrueMuonTrueMomentumTransverseRatio = new TH1D("TrueMuonTrueMomentumTransverseRatio",";P^{true}_{#mu,T}/P^{true}_{#mu}",25,0.,1.);		
	TH1D* TrueProtonTrueMomentumTransverseRatio = new TH1D("TrueProtonTrueMomentumTransverseRatio",";P^{true}_{p,T}/P^{true}_{p}",25,0.,1.);	

	//--------------------------------------------------//

	int TrueCC1pCounter = 0;
	int TrueCCQElikeCounter = 0;
	
	double SumWeights = 0.;
	
	//--------------------------------------------------//

	TH1D* POTScalePlot = new TH1D("POTScalePlot","",1,0,1);
	TH1D* NEventsPlot = new TH1D("NEventsPlot","",1,0,1);
	TH1D* NSelectedPlot = new TH1D("NSelectedPlot","",1,0,1);

	//--------------------------------------------------//

	// Loop over the events

	cout << nentries << " events included in the file" << endl;	

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		//--------------------------------------------------//

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;	

		//--------------------------------------------------//

//		double weight = 1.;
//		double T2Kweight = 1.;

		// For detector variations, the eventweight weights are -1., set them back to 1.
		if (Weight == -1.) { Weight = 1.; }
		if (T2KWeight == -1.) { T2KWeight = 1.; }

		if (Weight <= 0 || Weight > 30) { continue; }
		if (T2KWeight <= 0 || T2KWeight > 30) { continue; }		
		// Weight from v3.0.4 to v.3.0.6 * weight from application of T2K tune
		double weight = POTWeight * Weight * T2KWeight * ROOTinoWeight;

		// Fake data studies: removing the T2K tune weight
		if (fTune == "NoTune") { weight = POTWeight * Weight * ROOTinoWeight; }
		// Double the MEC weight  (mode = 10)
		if (fTune == "TwiceMEC" && Muon_MCParticle_Mode->at(0) == 10) { weight = 2 * POTWeight * Weight * T2KWeight * ROOTinoWeight; }
		if (fTune == "TwiceMEC" && Muon_MCParticle_Mode->at(0) != 10) { weight = POTWeight * Weight * T2KWeight * ROOTinoWeight; }		

		// -------------------------------------------------------------------------------------------------------------------------------------

		// Genie, flux & reinteraction weights for systematics

		if ( 
			fUniverseIndex != -1 && (fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || 
			fWhichSample == "Overlay9_Run3" || fWhichSample == "Overlay9_Run4a"
		|| fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5" || fWhichSample == "Overlay9_Combined" 
		|| fWhichSample == "OverlayDirt9_Run1" || fWhichSample == "OverlayDirt9_Run2" || 
			fWhichSample == "OverlayDirt9_Run3" || fWhichSample == "OverlayDirt9_Run4a"
		|| fWhichSample == "OverlayDirt9_Run4" || fWhichSample == "OverlayDirt9_Run5" || fWhichSample == "OverlayDirt9_Combined") ) {

			// Genie weights

			if (fEventWeightLabel == "All_UBGenie") { weight = weight*All_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "AxFFCCQEshape_UBGenie") { weight = weight*AxFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "DecayAngMEC_UBGenie") { weight = weight*DecayAngMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "NormCCCOH_UBGenie") { weight = weight*NormCCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "NormNCCOH_UBGenie") { weight = weight*NormNCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; }
//			if (fEventWeightLabel == "RPA_CCQE_Reduced_UBGenie") { weight = weight*RPA_CCQE_Reduced_UBGenie->at(fUniverseIndex)/ T2KWeight; }
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

				int concat = tools.ConcatRunSubRunEvent(Run,SubRun,Event);
				weight = weight*tools.PoissonRandomNumber(concat); 
				
			}			

		}

		// ---------------------------------------------------------------------------------------------------------------------------------

		SumWeights += weight / POTWeight;					

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Analysis over the simb::MCParticles

		std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
		std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();

		int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0;
		bool TrueCC1pEvent = false;
		bool TrueCCQElikeEvent = false;

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & 0 pi+/- (Ppi > 70 MeV / c)

		if (CC1p == 1 && NumberPi0 == 0) {
		
			// --------------------------------------------------------------------------------------------------------------------		
			
			// Containment of the vertex has already been demanded in PreTruthSelection
			// Soft fiducial volume for true vertex

			//if (Muon_MCParticle_StartContainment->at(0) == 0) { continue; }

			// --------------------------------------------------------------------------------------------------------------------		

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

			double TrueDeltaPhiProtonMuon_Deg = True_DeltaPhi->at(0);
			double TrueDeltaThetaProtonMuon_Deg = True_DeltaTheta->at(0);

			// Reconstructed calorimetric energy using true level info / MCParticles

			double TrueRecoECal = True_ECal->at(0); // GeV

			// Reconstructed QE energy using true level info / MCParticles

			double TrueRecoEQE = True_EQE->at(0);

			// Transverse Variables

			double TruekMiss = True_kMiss->at(0);
			double TruePMissMinus = True_PMissMinus->at(0);
			double TrueMissMomentum = True_PMiss->at(0);

			double TrueTransMissMomentum = True_Pt->at(0);
			double TrueDeltaAlphaT = True_DeltaAlphaT->at(0);
			double TrueDeltaPhiT = True_DeltaPhiT->at(0);

			double TruePL = True_PL->at(0);
			double TruePn = True_Pn->at(0);
			double TruePtx = True_Ptx->at(0);
			double TruePty = True_Pty->at(0);
			double TrueA = True_A->at(0);

			// -----------------------------------------------------------------------------------------------------------------	

			// True Vertex

			TVector3 TrueVertex(True_Vx,True_Vy,True_Vz);
				
			// Reconstructed q vector using true level info / MCParticles

			double RecoTrueQ2 = True_Q2->at(0);

			// -----------------------------------------------------------------------------------------------------------------	

			// Overflow bins
			// Affects

			// DeltaPT
			// DeltaPtx
			// DeltaPty
			// DeltaPL
			// DeltaPn
			// Q2
			// ECal
			// EQE
			// alpha
			// kMiss
			// PMiss
			// PMissMinus


			if (TrueTransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { TrueTransMissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
			if (TruePtx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { TruePtx = 0.5 * (ArrayNBinsDeltaPtx[NBinsDeltaPtx] + ArrayNBinsDeltaPtx[NBinsDeltaPtx-1]); }
			if (TruePty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { TruePty = 0.5 * (ArrayNBinsDeltaPty[NBinsDeltaPty] + ArrayNBinsDeltaPty[NBinsDeltaPty-1]); }
			if (TruePL > ArrayNBinsDeltaPL[NBinsDeltaPL]) { TruePL = 0.5 * (ArrayNBinsDeltaPL[NBinsDeltaPL] + ArrayNBinsDeltaPL[NBinsDeltaPL-1]); }						
			if (TruePn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { TruePn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

			if (TrueRecoECal > ArrayNBinsECal[NBinsECal]) { TrueRecoECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
			if (TrueRecoEQE > ArrayNBinsEQE[NBinsEQE]) { TrueRecoEQE = 0.5 * (ArrayNBinsEQE[NBinsEQE] + ArrayNBinsEQE[NBinsEQE-1]); }
			if (RecoTrueQ2 > ArrayNBinsQ2[NBinsQ2]) { RecoTrueQ2 = 0.5 * (ArrayNBinsQ2[NBinsQ2] + ArrayNBinsQ2[NBinsQ2-1]); }

			if (TrueA > ArrayNBinsA[NBinsA]) { TrueA = 0.5 * (ArrayNBinsA[NBinsA] + ArrayNBinsA[NBinsA-1]); }	
			if (TruekMiss > ArrayNBinskMiss[NBinskMiss]) { TruekMiss = 0.5 * (ArrayNBinskMiss[NBinskMiss] + ArrayNBinskMiss[NBinskMiss-1]); }														
			if (TrueMissMomentum > ArrayNBinsPMiss[NBinsPMiss]) { TrueMissMomentum = 0.5 * (ArrayNBinsPMiss[NBinsPMiss] + ArrayNBinsPMiss[NBinsPMiss-1]); }
			if (TruePMissMinus > ArrayNBinsPMissMinus[NBinsPMissMinus]) { TruePMissMinus = 0.5 * (ArrayNBinsPMissMinus[NBinsPMissMinus] + ArrayNBinsPMissMinus[NBinsPMissMinus-1]); }

			// ---------------------------------------------------------------------------------------------------------------------------

			// Underflow bins
			// Affects

			// ECal
			// EQE
			// DeltaPtx
			// DeltaPty
			// DeltaPL
			// alpha
			// PMissMinus
			
			if (TrueRecoECal < ArrayNBinsECal[0]) { TrueRecoECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			
			if (TrueRecoEQE < ArrayNBinsEQE[0]) { TrueRecoEQE = 0.5 * (ArrayNBinsEQE[0] + ArrayNBinsEQE[1]); }			
			if (TruePtx < ArrayNBinsDeltaPtx[0]) { TruePtx = 0.5 * (ArrayNBinsDeltaPtx[0] + ArrayNBinsDeltaPtx[1]); }
			if (TruePty < ArrayNBinsDeltaPty[0]) { TruePty = 0.5 * (ArrayNBinsDeltaPty[0] + ArrayNBinsDeltaPty[1]); }
			if (TruePL < ArrayNBinsDeltaPL[0]) { TruePL = 0.5 * (ArrayNBinsDeltaPL[0] + ArrayNBinsDeltaPL[1]); }						
			if (TrueA < ArrayNBinsA[0]) { TrueA = 0.5 * (ArrayNBinsA[0] + ArrayNBinsA[1]); }
			if (TruePMissMinus < ArrayNBinsPMissMinus[0]) { TruePMissMinus = 0.5 * (ArrayNBinsPMissMinus[0] + ArrayNBinsPMissMinus[1]); }			

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------									

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

				// -----------------------------------------------------------------------------------------------------------------

				// STV analysis

				if (
				    TrueTransMissMomentum > ArrayNBinsDeltaPT[0] 
				    //&& TrueTransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
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
				    
//				    && TrueRecoECal > ArrayNBinsECal[0]
//				    && TrueRecoEQE > ArrayNBinsEQE[0]
//				    && RecoTrueQ2 > ArrayNBinsQ2[0]
//				    && TrueRecoECal < ArrayNBinsECal[NBinsECal]
//				    && TrueRecoEQE < ArrayNBinsEQE[NBinsEQE]
//				    && RecoTrueQ2 < ArrayNBinsQ2[NBinsQ2]
				) {

					//----------------------------------------//

					// Indices for 2D analysis

					int DeltaPTTwoDIndex = tools.ReturnIndex(TrueTransMissMomentum, TwoDArrayNBinsDeltaPT);
					int DeltaAlphaTTwoDIndex = tools.ReturnIndex(TrueDeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int MuonCosThetaTwoDIndex = tools.ReturnIndex(TrueMuonCosTheta, TwoDArrayNBinsMuonCosTheta);
					int ProtonCosThetaTwoDIndex = tools.ReturnIndex(TrueProtonCosTheta, TwoDArrayNBinsProtonCosTheta);
					int DeltaPtxTwoDIndex = tools.ReturnIndex(TruePtx, TwoDArrayNBinsDeltaPtx);	
					int DeltaPtyTwoDIndex = tools.ReturnIndex(TruePty, TwoDArrayNBinsDeltaPty);	
					int MuonMomentumTwoDIndex = tools.ReturnIndex(TrueMuonMomentum_GeV, TwoDArrayNBinsMuonMomentum);
					int ProtonMomentumTwoDIndex = tools.ReturnIndex(TrueProtonMomentum_GeV, TwoDArrayNBinsProtonMomentum);

					// 2D analysis treated in 1D grid

					int SerialMuonMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueMuonMomentum_GeV);
					int SerialProtonMomentumInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueProtonMomentum_GeV);				
					int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueDeltaAlphaT);
					int SerialDeltaAlphaTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueDeltaAlphaT);
					int SerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,DeltaPTTwoDIndex,TrueDeltaAlphaT);
					int SerialDeltaPhiTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices,DeltaPTTwoDIndex,TrueDeltaPhiT);
					int SerialDeltaPnInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaPTSlices,DeltaPTTwoDIndex,TruePn);	
					int SerialECalInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPTSlices,DeltaPTTwoDIndex,TrueRecoECal);															
					int SerialECalInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TrueRecoECal);
					int SerialProtonCosThetaInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueProtonCosTheta);
					int SerialDeltaPtyInDeltaPtxIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices,DeltaPtxTwoDIndex,TruePty);				
					int SerialDeltaPtxInDeltaPtyIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtxInDeltaPtySlices,DeltaPtyTwoDIndex,TruePtx);
					int SerialDeltaPTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueTransMissMomentum);
					int SerialDeltaPTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueTransMissMomentum);
					int SerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TrueTransMissMomentum);					

					//----------------------------------------//

					// 3D analysis	

					TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);	

					// 2D analysis

					TrueDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);					
					TrueDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
					TrueDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TruePn,weight);
					TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);					
					TrueDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TrueTransMissMomentum,weight);					
					TrueMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TrueMuonMomentum_GeV,weight);
					TrueProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TrueProtonCosTheta,weight);					
					TrueProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TrueProtonMomentum_GeV,weight);
					TrueDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(TruePty,weight);
					TrueDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(TruePtx,weight);	
					TrueDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TrueTransMissMomentum,weight);
					TrueECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueRecoECal,weight);															

					// 2D analysis in 1D grid

					SerialTrueMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialTrueProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialTrueDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialTrueDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialTrueECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialTrueProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialTrueDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialTrueDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialTrueDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialTrueDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);
					SerialTrueDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialTrueECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);	

					// 3D analysis treated in 1D grid

					int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,TrueRecoECal);
					int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,TrueRecoECal);
					int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,TrueRecoECal);																																																									

					SerialTrueECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialTrueECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialTrueECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);

					//----------------------------------------//									

					// No weight to be applied in the multiplicity plots

					TruePi0Plot->Fill(NumberPi0);
					TrueNeutronPlot->Fill(NumberNeutrons);

					// True CC1p event

					TrueCC1pEvent = true;
					TrueCC1pCounter++;

					// STV

					TruekMissPlot->Fill(TruekMiss,weight);
					TruePMissMinusPlot->Fill(TruePMissMinus,weight);
					TruePMissPlot->Fill(TrueMissMomentum,weight);

					TrueDeltaPTPlot->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);

					TrueDeltaPLPlot->Fill(TruePL,weight);
					TrueDeltaPnPlot->Fill(TruePn,weight);
					TrueDeltaPtxPlot->Fill(TruePtx,weight);
					TrueDeltaPtyPlot->Fill(TruePty,weight);
					TrueAPlot->Fill(TrueA,weight);

					// Reconstructed energy & Q2

					TrueECalPlot->Fill(TrueRecoECal,weight);
					TrueEQEPlot->Fill(TrueRecoEQE,weight);
					TrueQ2Plot->Fill(RecoTrueQ2,weight);

					// Kinematic variables

					TrueMuonMomentumPlot->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonPhiPlot->Fill(TrueMuonPhi_Deg,weight);
					TrueMuonCosThetaPlot->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot->Fill(TrueMuonCosTheta,weight);

					// Given that we keep only contained reco tracks, we want to also calculate the efficiency with respect to that part of the MC
					if (Muon_MCParticle_EndContainment->at(0) == 1) { TrueContainedEndMuonCosThetaPlot->Fill(TrueMuonCosTheta,weight); }

					TrueProtonMomentumPlot->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonPhiPlot->Fill(TrueProtonPhi_Deg,weight);
					TrueProtonCosThetaPlot->Fill(TrueProtonCosTheta,weight);

					// True Vertex

					TrueVertexXPlot->Fill(TrueVertex.X(),weight);
					TrueVertexYPlot->Fill(TrueVertex.Y(),weight);
					TrueVertexZPlot->Fill(TrueVertex.Z(),weight);

					// True Energy

					TrueEvPlot->Fill(True_Ev,weight);

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(TrueMuonMomentum_GeV,2.) );
					double true_Nu = True_Ev - true_MuonEnergy;
					TrueNuPlot->Fill(true_Nu,weight);
					
					// 2D Analysis
		
					TrueCosThetaMuPmuPlot->Fill(TrueMuonCosTheta,TrueMuonMomentum_GeV,weight);
					TrueCosThetaPPpPlot->Fill(TrueProtonCosTheta,TrueProtonMomentum_GeV,weight);

					// Playground for CC1p true momenta (longitudinal & perpendicular) ratios

					TVector3 TrueCandidateMuon(1,1,1);
					TrueCandidateMuon.SetMag(TrueMuonMomentum_GeV);
					TrueCandidateMuon.SetPhi(TrueMuonPhi);
					TrueCandidateMuon.SetTheta(TMath::ACos(TrueMuonCosTheta));

					TVector3 TrueCandidateProton(1,1,1);
					TrueCandidateProton.SetMag(TrueProtonMomentum_GeV);
					TrueCandidateProton.SetPhi(TrueProtonPhi);
					TrueCandidateProton.SetTheta(TMath::ACos(TrueProtonCosTheta));	

					TrueMuonTrueMomentumLongitudinalRatio->Fill(TrueCandidateMuon.Z() / TrueCandidateMuon.Mag(),weight);
					TrueMuonTrueMomentumTransverseRatio->Fill(TrueCandidateMuon.Pt() / TrueCandidateMuon.Mag(),weight);

					TrueProtonTrueMomentumLongitudinalRatio->Fill(TrueCandidateProton.Z() / TrueCandidateProton.Mag(),weight);
					TrueProtonTrueMomentumTransverseRatio->Fill(TrueCandidateProton.Pt() / TrueCandidateProton.Mag(),weight);				
					
					// CCQElike

					if (
						TrueTransMissMomentum < 0.35
						&& TMath::Abs(TrueDeltaThetaProtonMuon_Deg - 90.) < 55.
						&& TMath::Abs(TrueDeltaPhiProtonMuon_Deg - 180.) < 35.						

						&& TrueMuonMomentum_GeV > ArrayNBinsMuonMomentum[0]
						&& TrueProtonMomentum_GeV > ArrayNBinsProtonMomentum[0]						
						&& TrueMuonMomentum_GeV < ArrayNBinsMuonMomentum[CCQENBinsMuonMomentum]
						&& TrueProtonMomentum_GeV < ArrayNBinsProtonMomentum[CCQENBinsProtonMomentum]
						
						&& TrueMuonCosTheta > CCQEArrayNBinsMuonCosTheta[0]
						&& TrueMuonCosTheta < CCQEArrayNBinsMuonCosTheta[CCQENBinsMuonCosTheta]
						&& TrueProtonCosTheta > CCQEArrayNBinsProtonCosTheta[0]
						&& TrueProtonCosTheta < CCQEArrayNBinsProtonCosTheta[CCQENBinsProtonCosTheta]
						
					) {

						TrueCCQEECalPlot->Fill(TrueRecoECal,weight);
						TrueCCQEQ2Plot->Fill(RecoTrueQ2,weight);

						TrueCCQEMuonMomentumPlot->Fill(TrueMuonMomentum_GeV,weight);
						TrueCCQEMuonPhiPlot->Fill(TrueMuonPhi_Deg,weight);
						TrueCCQEMuonCosThetaPlot->Fill(TrueMuonCosTheta,weight);

						TrueCCQEProtonMomentumPlot->Fill(TrueProtonMomentum_GeV,weight);
						TrueCCQEProtonPhiPlot->Fill(TrueProtonPhi_Deg,weight);
						TrueCCQEProtonCosThetaPlot->Fill(TrueProtonCosTheta,weight);

					}					


				}

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
	
	POTScalePlot->SetBinContent(1,POTWeight);
	POTScalePlot->SetBinError(1,0);

	NEventsPlot->SetBinContent(1,nentries);
	NEventsPlot->SetBinError(1,TMath::Sqrt(nentries));

	NSelectedPlot->SetBinContent(1,TrueCC1pCounter);
	NSelectedPlot->SetBinError(1,sqrt(TrueCC1pCounter));	

	//----------------------------------------//

//	double ScaleDueToWeights = double(fChain->GetEntries()) / double(SumWeights);

	//----------------------------------------//

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();
	myTxtFile.close();

	fFile->Close();

	//----------------------------------------//

} // End of the program
