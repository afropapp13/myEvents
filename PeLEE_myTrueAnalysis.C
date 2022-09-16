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

	// Reference plots

	TH1D* TruePi0Plot = new TH1D("TruePi0Plot",";# #pi^{0}",4,-0.5,3.5);
	TH1D* TrueNeutronPlot = new TH1D("TrueNeutronPlot",";# neutrons",6,-0.5,5.5);
	TH1D* TrueMuonTrueMomentumLongitudinalRatio = new TH1D("TrueMuonTrueMomentumLongitudinalRatio",";P^{true}_{#mu,||}/P^{true}_{#mu}",25,0.,1.);		
	TH1D* TrueProtonTrueMomentumLongitudinalRatio = new TH1D("TrueProtonTrueMomentumLongitudinalRatio",";P^{true}_{p,||}/P^{true}_{p}",25,0.,1.);
	TH1D* TrueMuonTrueMomentumTransverseRatio = new TH1D("TrueMuonTrueMomentumTransverseRatio",";P^{true}_{#mu,T}/P^{true}_{#mu}",25,0.,1.);		
	TH1D* TrueProtonTrueMomentumTransverseRatio = new TH1D("TrueProtonTrueMomentumTransverseRatio",";P^{true}_{p,T}/P^{true}_{p}",25,0.,1.);	
	TH1D* TrueEvPlot = new TH1D("TrueEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
	TH1D* TrueNuPlot = new TH1D("TrueNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
	TH1D* POTScalePlot = new TH1D("POTScalePlot","",1,0,1);
	TH1D* NEventsPlot = new TH1D("NEventsPlot","",1,0,1);
	TH1D* NSelectedPlot = new TH1D("NSelectedPlot","",1,0,1);	

	//--------------------------------------------------//

	// 1D analysis

	TH1D* TrueVertexXPlot[NInte];
	TH1D* TrueVertexYPlot[NInte];
	TH1D* TrueVertexZPlot[NInte];	

	TH1D* TrueDeltaPTPlot[NInte];
	TH1D* TrueDeltaAlphaTPlot[NInte];
	TH1D* TrueDeltaAlpha3DqPlot[NInte];
	TH1D* TrueDeltaAlpha3DMuPlot[NInte];		
	TH1D* TrueDeltaPhiTPlot[NInte];
	TH1D* TrueDeltaPLPlot[NInte];
	TH1D* TrueDeltaPnPlot[NInte];
	TH1D* TrueDeltaPtxPlot[NInte];
	TH1D* TrueDeltaPtyPlot[NInte];
	TH1D* TrueAPlot[NInte];
	TH1D* TrueMuonMomentumPlot[NInte];
	TH1D* TrueMuonPhiPlot[NInte];
	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueMuonCosThetaSingleBinPlot[NInte];
	TH1D* TrueProtonMomentumPlot[NInte];
	TH1D* TrueProtonPhiPlot[NInte];
	TH1D* TrueProtonCosThetaPlot[NInte];
	TH1D* TrueECalPlot[NInte];
	TH1D* TrueEQEPlot[NInte];
	TH1D* TrueQ2Plot[NInte];
	TH1D* TrueCCQEMuonMomentumPlot[NInte];
	TH1D* TrueCCQEMuonPhiPlot[NInte];
	TH1D* TrueCCQEMuonCosThetaPlot[NInte];	
	TH1D* TrueCCQEProtonMomentumPlot[NInte];
	TH1D* TrueCCQEProtonPhiPlot[NInte];
	TH1D* TrueCCQEProtonCosThetaPlot[NInte];
	TH1D* TrueCCQEECalPlot[NInte];
	TH1D* TrueCCQEQ2Plot[NInte];	
	TH1D* TruekMissPlot[NInte];
	TH1D* TruePMissMinusPlot[NInte];
	TH1D* TruePMissPlot[NInte];

	//--------------------------------------------------//

	// 2D analysis (uncorrelated)

	TH1D* TrueDeltaPT_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsMuonCosTheta];
	TH1D* TrueDeltaPT_InProtonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];	
	TH1D* TrueDeltaPT_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
	TH1D* TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[NInte][TwoDNBinsDeltaAlpha3Dq];
	TH1D* TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[NInte][TwoDNBinsDeltaAlpha3DMu];		
	TH1D* TrueProtonMomentum_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];	
	TH1D* TrueMuonMomentum_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsMuonCosTheta];
	TH1D* TrueProtonMomentum_InProtonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];			
	TH1D* TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsMuonCosTheta];
	TH1D* TrueDeltaAlphaT_InMuonMomentumTwoDPlot[NInte][TwoDNBinsMuonMomentum];	
	TH1D* TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];
	TH1D* TrueDeltaAlphaT_InProtonMomentumTwoDPlot[NInte][TwoDNBinsProtonMomentum];			
	TH1D* TrueDeltaAlphaT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
	TH1D* TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
	TH1D* TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];					
	TH1D* TrueProtonCosTheta_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];	
	TH1D* TrueDeltaPhiT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];	
	TH1D* TrueDeltaPn_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
	TH1D* TrueDeltaPn_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
	TH1D* TrueDeltaPty_InDeltaPtxTwoDPlot[NInte][TwoDNBinsDeltaPtx];
	TH1D* TrueDeltaPtx_InDeltaPtyTwoDPlot[NInte][TwoDNBinsDeltaPty];	
	TH1D* TrueECal_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];	
	TH1D* TrueECal_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
	TH1D* TrueECal_InDeltaPtxTwoDPlot[NInte][TwoDNBinsDeltaPtx];
	TH1D* TrueECal_InDeltaPtyTwoDPlot[NInte][TwoDNBinsDeltaPty];		

	//--------------------------------------------------//	

	// 3D analysis (uncorrelated)	

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
	TH1D* TrueECal_InDeltaPtxDeltaPtyTwoDPlot[NInte][TwoDNBinsDeltaPtx][TwoDNBinsDeltaPty];	
	TH1D* TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[NInte][TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
	TH1D* TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[NInte][TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];	

	//--------------------------------------------------//

	// 2D analysis in 1D grid	

	TH1D* SerialTrueDeltaPT_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaPT_InProtonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaPT_InDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaAlpha3DqPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[NInte];		
	TH1D* SerialTrueProtonMomentum_InDeltaAlphaTPlot[NInte];	
	TH1D* SerialTrueMuonMomentum_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueProtonMomentum_InProtonCosThetaPlot[NInte];	
	TH1D* SerialTrueDeltaAlphaT_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaAlphaT_InMuonMomentumPlot[NInte];	
	TH1D* SerialTrueDeltaAlphaT_InProtonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaAlphaT_InProtonMomentumPlot[NInte];	
	TH1D* SerialTrueDeltaAlphaT_InDeltaPTPlot[NInte];
	TH1D* SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[NInte];
	TH1D* SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[NInte];			
	TH1D* SerialTrueDeltaPhiT_InDeltaPTPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaPTPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaAlphaTPlot[NInte];	
	TH1D* SerialTrueProtonCosTheta_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaPty_InDeltaPtxPlot[NInte];
	TH1D* SerialTrueDeltaPtx_InDeltaPtyPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPTPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPtxPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPtyPlot[NInte];		
	TH1D* SerialTrueECal_InDeltaAlphaTPlot[NInte];

	//--------------------------------------------------//	

	// 3D analysis in 1D grid		
	
	TH1D* SerialTrueECal_InDeltaPTDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPtxDeltaPtyPlot[NInte];	
	TH1D* SerialTrueECal_InMuonCosThetaMuonMomentumPlot[NInte];
	TH1D* SerialTrueECal_InProtonCosThetaProtonMomentumPlot[NInte];

	//--------------------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//--------------------------------------------------//

		// 1D analysis

		TrueVertexXPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TrueVertexYPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TrueVertexZPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

		TrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TrueDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);
		TrueDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlpha3DMuPlot",LabelXAxisDeltaAlpha3DMu,NBinsDeltaAlpha3DMu,ArrayNBinsDeltaAlpha3DMu);				
		TrueDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TrueDeltaPLPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TrueDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TrueDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TrueDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TrueAPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);
		TrueMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TrueMuonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TrueMuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,0.,1.);
		TrueProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TrueProtonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
		TrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TrueECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TrueEQEPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TrueQ2Plot[inte] = new TH1D(InteractionLabels[inte]+"TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		TrueCCQEMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TrueCCQEMuonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);
		TrueCCQEMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);	
		TrueCCQEProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TrueCCQEProtonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);
		TrueCCQEProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TrueCCQEECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TrueCCQEQ2Plot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);	
		TruekMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TruePMissMinusPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TruePMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		//--------------------------------------------------//	

		// 2D & 3D analysis (uncorrelated)

		for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

			TString DeltaAlpha3DqTwoDInDeltaPnLabel = "DeltaAlpha3Dq_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";			
			TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlpha3DqTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3Dq,TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][0]);

			TString DeltaAlpha3DMuTwoDInDeltaPnLabel = "DeltaAlpha3DMu_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";			
			TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlpha3DMuTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3DMu,TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][0]);			

		}

		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

			TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);			

			TString ECalTwoDInDeltaPTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueECal_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);		


			TString DeltaPhiTTwoDInDeltaPTLabel = "DeltaPhiT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueDeltaPhiT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	

			TString DeltaPnTwoDInDeltaPTLabel = "DeltaPn_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueDeltaPn_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	

			for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

				TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[inte][WhichDeltaPT][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

			}

		}

		for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {

			TString DeltaAlphaTTwoDInMuonMomentumLabel = "DeltaAlphaT_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";			
			TrueDeltaAlphaT_InMuonMomentumTwoDPlot[inte][WhichMuonMomentum] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInMuonMomentumLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonMomentumSlices[WhichMuonMomentum].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonMomentumSlices[WhichMuonMomentum][0]);

		}

		for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

			TString DeltaAlphaTTwoDInMuonCosThetaLabel = "DeltaAlphaT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			TString DeltaPTTwoDInMuonCosThetaLabel = "DeltaPT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueDeltaPT_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);		

			TString MuonMomentumTwoDInMuonCosThetaLabel = "MuonMomentum_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueMuonMomentum_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			TString ProtonCosThetaTwoDInMuonCosThetaLabel = "ProtonCosTheta_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueProtonCosTheta_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

				TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";
				TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[inte][WhichMuonCosTheta][WhichMuonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);

			}

		}	

		for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {

			TString DeltaAlphaTTwoDInProtonMomentumLabel = "DeltaAlphaT_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";			
			TrueDeltaAlphaT_InProtonMomentumTwoDPlot[inte][WhichProtonMomentum] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInProtonMomentumLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonMomentumSlices[WhichProtonMomentum].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonMomentumSlices[WhichProtonMomentum][0]);

		}			

		for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

			TString ProtonMomentumTwoDInProtonCosThetaLabel = "ProtonMomentum_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
			TrueProtonMomentum_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			TString DeltaAlphaTTwoDInProtonCosThetaLabel = "DeltaAlphaT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
			TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			TString DeltaPTTwoDInProtonCosThetaLabel = "DeltaPT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
			TrueDeltaPT_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

				TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
				TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[inte][WhichProtonCosTheta][WhichProtonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);

			}

		}

		for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

			TString DeltaPtyTwoDInDeltaPtxLabel = "DeltaPty_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			
			TrueDeltaPty_InDeltaPtxTwoDPlot[inte][WhichDeltaPtx] = new TH1D(InteractionLabels[inte]+"True"+DeltaPtyTwoDInDeltaPtxLabel,LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);	

			TString ECalTwoDInDeltaPtxLabel = "ECal_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			
			TrueECal_InDeltaPtxTwoDPlot[inte][WhichDeltaPtx] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPtxLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsECalInDeltaPtxSlices[WhichDeltaPtx][0]);

			for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {	

				TString ECalTwoDInDeltaPtxDeltaPtyLabel = "ECal_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";
				TrueECal_InDeltaPtxDeltaPtyTwoDPlot[inte][WhichDeltaPtx][WhichDeltaPty] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPtxDeltaPtyLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[WhichDeltaPtx][WhichDeltaPty].size()-1,&TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[WhichDeltaPtx][WhichDeltaPty][0]);

			}

		}	

		for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

			TString DeltaPtxTwoDInDeltaPtyLabel = "DeltaPtx_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";			
			TrueDeltaPtx_InDeltaPtyTwoDPlot[inte][WhichDeltaPty] = new TH1D(InteractionLabels[inte]+"True"+DeltaPtxTwoDInDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);	

			TString ECalTwoDInDeltaPtyLabel = "ECal_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";			
			TrueECal_InDeltaPtyTwoDPlot[inte][WhichDeltaPty] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPtyLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsECalInDeltaPtySlices[WhichDeltaPty][0]);

		}		

		for (int WhichDeltaAlpha3Dq = 0; WhichDeltaAlpha3Dq < TwoDNBinsDeltaAlpha3Dq; WhichDeltaAlpha3Dq++) {

			TString DeltaPnTwoDInDeltaAlpha3DqLabel = "DeltaPn_DeltaAlpha3Dq_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq+1])+"Plot";			
			TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaAlpha3DqLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][0]);			

		}

		for (int WhichDeltaAlpha3DMu = 0; WhichDeltaAlpha3DMu < TwoDNBinsDeltaAlpha3DMu; WhichDeltaAlpha3DMu++) {

			TString DeltaPnTwoDInDeltaAlpha3DMuLabel = "DeltaPn_DeltaAlpha3DMu_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3DMu[WhichDeltaAlpha3DMu])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3DMu[WhichDeltaAlpha3DMu+1])+"Plot";			
			TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaAlpha3DMuLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][0]);			

		}		

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

			TString ProtonMomentumTwoDInDeltaAlphaTLabel = "ProtonMomentum_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueProtonMomentum_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ProtonMomentumTwoDInDeltaAlphaTLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsProtonMomentumInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

			TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

			TString ECalTwoDInDeltaAlphaTLabel = "ECal_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueECal_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

			TString DeltaPnTwoDInDeltaAlphaTLabel = "DeltaPn_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueDeltaPn_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

		}

		//--------------------------------------------------//	

		// 2D analysis in 1D grid	

		SerialTrueDeltaPT_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		SerialTrueDeltaPT_InProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);	
		SerialTrueDeltaPT_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		SerialTrueDeltaPn_InDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPn_DeltaAlpha3DqPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices)[0]);
		SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPn_DeltaAlpha3DMuPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices)[0]);				
		SerialTrueProtonMomentum_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialProtonMomentum_DeltaAlphaTPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInDeltaAlphaTSlices)[0]);		
		SerialTrueDeltaPn_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);		
		SerialTrueMuonMomentum_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		SerialTrueProtonMomentum_InProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);	
		SerialTrueDeltaAlphaT_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		SerialTrueDeltaAlphaT_InMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_MuonMomentumPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonMomentumSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonMomentumSlices)[0]);
		SerialTrueDeltaAlphaT_InProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		SerialTrueDeltaAlphaT_InProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_ProtonMomentumPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonMomentumSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonMomentumSlices)[0]);		
		SerialTrueDeltaAlphaT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlpha3Dq_DeltaPnPlot",LabelXAxisDeltaAlpha3Dq,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices)[0]);
		SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlpha3DMu_DeltaPnPlot",LabelXAxisDeltaAlpha3DMu,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices)[0]);						
		SerialTrueDeltaPhiT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		SerialTrueDeltaPn_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		SerialTrueProtonCosTheta_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		SerialTrueDeltaPty_InDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		SerialTrueDeltaPtx_InDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		SerialTrueECal_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		SerialTrueECal_InDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPtxPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPtxSlices)[0]);
		SerialTrueECal_InDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPtyPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPtySlices)[0]);				
		SerialTrueECal_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);

		//--------------------------------------------------//	

		// 3D analysis in 1D grid		
		
		SerialTrueECal_InDeltaPTDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		SerialTrueECal_InDeltaPtxDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPtxDeltaPtyPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices)[0]);		
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
			fUniverseIndex != -1 && (fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || 
			fWhichSample == "Overlay9_Run3" || fWhichSample == "Overlay9_Run4a" ||
			fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5" || fWhichSample == "Overlay9_Combined" ||
			fWhichSample == "OverlayDirt9_Run1" || fWhichSample == "OverlayDirt9_Run2" || 
			fWhichSample == "OverlayDirt9_Run3" || fWhichSample == "OverlayDirt9_Run4a" ||
			fWhichSample == "OverlayDirt9_Run4" || fWhichSample == "OverlayDirt9_Run5" || fWhichSample == "OverlayDirt9_Combined") 
		) {

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

				int concat = tools.ConcatRunSubRunEvent(Run,SubRun,Event,fUniverseIndex);
				weight = weight*tools.PoissonRandomNumber(concat); 
				
			}			

		}

		//--------------------------------------------------//

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
			double TrueDeltaAlpha3Dq = True_DeltaAlpha3Dq->at(0);
			double TrueDeltaAlpha3DMu = True_DeltaAlpha3DMu->at(0);						
			double TrueDeltaPhiT = True_DeltaPhiT->at(0);

			double TruePL = True_PL->at(0);
			double TruePn = True_Pn->at(0);
			double TruePtx = True_Ptx->at(0);
			double TruePty = True_Pty->at(0);
			double TrueA = True_A->at(0);

			//--------------------------------------------------//	

			// True Vertex

			TVector3 TrueVertex(True_Vx,True_Vy,True_Vz);
				
			// Reconstructed q vector using true level info / MCParticles

			double RecoTrueQ2 = True_Q2->at(0);

			//--------------------------------------------------//	

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

			//--------------------------------------------------//

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

					int genie_mode = -1;

					if (Muon_MCParticle_Mode->at(0) == 0) { genie_mode = 1; }
					else if (Muon_MCParticle_Mode->at(0) == 10) { genie_mode = 2; }
					else if (Muon_MCParticle_Mode->at(0) == 1) { genie_mode = 3; }
					else if (Muon_MCParticle_Mode->at(0) == 2) { genie_mode = 4; }
					else { genie_mode = 5; }																				

					//----------------------------------------//

					// No weight to be applied in the multiplicity plots

					TruePi0Plot->Fill(NumberPi0);
					TrueNeutronPlot->Fill(NumberNeutrons);

					// True CC1p event

					TrueCC1pEvent = true;
					TrueCC1pCounter++;

					// True Energy

					TrueEvPlot->Fill(True_Ev,weight);

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(TrueMuonMomentum_GeV,2.) );
					double true_Nu = True_Ev - true_MuonEnergy;
					TrueNuPlot->Fill(true_Nu,weight);

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

					//----------------------------------------//

					// Indices for 2D & 3D analysis

					int DeltaPTTwoDIndex = tools.ReturnIndex(TrueTransMissMomentum, TwoDArrayNBinsDeltaPT);
					int DeltaPnTwoDIndex = tools.ReturnIndex(TruePn, TwoDArrayNBinsDeltaPn);					
					int DeltaAlphaTTwoDIndex = tools.ReturnIndex(TrueDeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int DeltaAlpha3DqTwoDIndex = tools.ReturnIndex(TrueDeltaAlpha3Dq, TwoDArrayNBinsDeltaAlpha3Dq);
					int DeltaAlpha3DMuTwoDIndex = tools.ReturnIndex(TrueDeltaAlpha3DMu, TwoDArrayNBinsDeltaAlpha3DMu);									
					int MuonCosThetaTwoDIndex = tools.ReturnIndex(TrueMuonCosTheta, TwoDArrayNBinsMuonCosTheta);
					int ProtonCosThetaTwoDIndex = tools.ReturnIndex(TrueProtonCosTheta, TwoDArrayNBinsProtonCosTheta);
					int DeltaPtxTwoDIndex = tools.ReturnIndex(TruePtx, TwoDArrayNBinsDeltaPtx);	
					int DeltaPtyTwoDIndex = tools.ReturnIndex(TruePty, TwoDArrayNBinsDeltaPty);	
					int MuonMomentumTwoDIndex = tools.ReturnIndex(TrueMuonMomentum_GeV, TwoDArrayNBinsMuonMomentum);
					int ProtonMomentumTwoDIndex = tools.ReturnIndex(TrueProtonMomentum_GeV, TwoDArrayNBinsProtonMomentum);

					int SerialMuonMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueMuonMomentum_GeV);
					int SerialProtonMomentumInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueProtonMomentum_GeV);				
					int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueDeltaAlphaT);
					int SerialDeltaAlphaTInMuonMomentumIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInMuonMomentumSlices,MuonMomentumTwoDIndex,TrueDeltaAlphaT);					
					int SerialDeltaAlphaTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueDeltaAlphaT);
					int SerialDeltaAlphaTInProtonMomentumIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInProtonMomentumSlices,ProtonMomentumTwoDIndex,TrueDeltaAlphaT);					
					int SerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,DeltaPTTwoDIndex,TrueDeltaAlphaT);
					int SerialDeltaAlpha3DqInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices,DeltaPnTwoDIndex,TrueDeltaAlpha3Dq);
					int SerialDeltaAlpha3DMuInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices,DeltaPnTwoDIndex,TrueDeltaAlpha3DMu);										
					int SerialDeltaPhiTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices,DeltaPTTwoDIndex,TrueDeltaPhiT);
					int SerialDeltaPnInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaPTSlices,DeltaPTTwoDIndex,TruePn);	
					int SerialECalInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPTSlices,DeltaPTTwoDIndex,TrueRecoECal);
					int SerialECalInDeltaPtxIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPtxSlices,DeltaPtxTwoDIndex,TrueRecoECal);
					int SerialECalInDeltaPtyIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPtySlices,DeltaPtyTwoDIndex,TrueRecoECal);																									
					int SerialECalInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TrueRecoECal);
					int SerialProtonCosThetaInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueProtonCosTheta);
					int SerialDeltaPtyInDeltaPtxIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices,DeltaPtxTwoDIndex,TruePty);				
					int SerialDeltaPtxInDeltaPtyIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtxInDeltaPtySlices,DeltaPtyTwoDIndex,TruePtx);
					int SerialDeltaPTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueTransMissMomentum);
					int SerialDeltaPTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueTransMissMomentum);
					int SerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TrueTransMissMomentum);
					int SerialDeltaPnInDeltaAlpha3DqIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices,DeltaAlpha3DqTwoDIndex,TruePn);
					int SerialDeltaPnInDeltaAlpha3DMuIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices,DeltaAlpha3DMuTwoDIndex,TruePn);										
					int SerialProtonMomentumInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonMomentumInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TrueProtonMomentum_GeV);					
					int SerialDeltaPnInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TruePn);										

					int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,TrueRecoECal);
					int SerialECalInDeltaPtxDeltaPtyIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices,DeltaPtxTwoDIndex,DeltaPtyTwoDIndex,TrueRecoECal);					
					int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,TrueRecoECal);
					int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,TrueRecoECal);																																																									

					//----------------------------------------//	

					// 1D analysis		

					// True Vertex

					TrueVertexXPlot[0]->Fill(TrueVertex.X(),weight);
					TrueVertexYPlot[0]->Fill(TrueVertex.Y(),weight);
					TrueVertexZPlot[0]->Fill(TrueVertex.Z(),weight);											

					TruekMissPlot[0]->Fill(TruekMiss,weight);
					TruePMissMinusPlot[0]->Fill(TruePMissMinus,weight);
					TruePMissPlot[0]->Fill(TrueMissMomentum,weight);
					TrueDeltaPTPlot[0]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot[0]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlpha3DqPlot[0]->Fill(TrueDeltaAlpha3Dq,weight);
					TrueDeltaAlpha3DMuPlot[0]->Fill(TrueDeltaAlpha3DMu,weight);										
					TrueDeltaPhiTPlot[0]->Fill(TrueDeltaPhiT,weight);
					TrueDeltaPLPlot[0]->Fill(TruePL,weight);
					TrueDeltaPnPlot[0]->Fill(TruePn,weight);
					TrueDeltaPtxPlot[0]->Fill(TruePtx,weight);
					TrueDeltaPtyPlot[0]->Fill(TruePty,weight);
					TrueAPlot[0]->Fill(TrueA,weight);
					TrueECalPlot[0]->Fill(TrueRecoECal,weight);
					TrueEQEPlot[0]->Fill(TrueRecoEQE,weight);
					TrueQ2Plot[0]->Fill(RecoTrueQ2,weight);
					TrueMuonMomentumPlot[0]->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonPhiPlot[0]->Fill(TrueMuonPhi_Deg,weight);
					TrueMuonCosThetaPlot[0]->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot[0]->Fill(0.5,weight);
					TrueProtonMomentumPlot[0]->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonPhiPlot[0]->Fill(TrueProtonPhi_Deg,weight);
					TrueProtonCosThetaPlot[0]->Fill(TrueProtonCosTheta,weight);

					TrueVertexXPlot[genie_mode]->Fill(TrueVertex.X(),weight);
					TrueVertexYPlot[genie_mode]->Fill(TrueVertex.Y(),weight);
					TrueVertexZPlot[genie_mode]->Fill(TrueVertex.Z(),weight);

					TruekMissPlot[genie_mode]->Fill(TruekMiss,weight);
					TruePMissMinusPlot[genie_mode]->Fill(TruePMissMinus,weight);
					TruePMissPlot[genie_mode]->Fill(TrueMissMomentum,weight);
					TrueDeltaPTPlot[genie_mode]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot[genie_mode]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlpha3DqPlot[genie_mode]->Fill(TrueDeltaAlpha3Dq,weight);
					TrueDeltaAlpha3DMuPlot[genie_mode]->Fill(TrueDeltaAlpha3DMu,weight);										
					TrueDeltaPhiTPlot[genie_mode]->Fill(TrueDeltaPhiT,weight);
					TrueDeltaPLPlot[genie_mode]->Fill(TruePL,weight);
					TrueDeltaPnPlot[genie_mode]->Fill(TruePn,weight);
					TrueDeltaPtxPlot[genie_mode]->Fill(TruePtx,weight);
					TrueDeltaPtyPlot[genie_mode]->Fill(TruePty,weight);
					TrueAPlot[genie_mode]->Fill(TrueA,weight);
					TrueECalPlot[genie_mode]->Fill(TrueRecoECal,weight);
					TrueEQEPlot[genie_mode]->Fill(TrueRecoEQE,weight);
					TrueQ2Plot[genie_mode]->Fill(RecoTrueQ2,weight);
					TrueMuonMomentumPlot[genie_mode]->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonPhiPlot[genie_mode]->Fill(TrueMuonPhi_Deg,weight);
					TrueMuonCosThetaPlot[genie_mode]->Fill(TrueMuonCosTheta,weight);
					TrueMuonCosThetaSingleBinPlot[genie_mode]->Fill(0.5,weight);
					TrueProtonMomentumPlot[genie_mode]->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonPhiPlot[genie_mode]->Fill(TrueProtonPhi_Deg,weight);
					TrueProtonCosThetaPlot[genie_mode]->Fill(TrueProtonCosTheta,weight);					

					//----------------------------------------//

					// 2D analysis (uncorrelated)

					TrueDeltaAlphaT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(TrueDeltaAlpha3Dq,weight);
					TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(TrueDeltaAlpha3DMu,weight);															
					TrueDeltaPhiT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
					TrueDeltaPn_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TruePn,weight);
					TrueDeltaPn_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(TruePn,weight);					
					TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlphaT_InMuonMomentumTwoDPlot[0][MuonMomentumTwoDIndex]->Fill(TrueDeltaAlphaT,weight);					
					TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[0][ProtonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlphaT_InProtonMomentumTwoDPlot[0][ProtonMomentumTwoDIndex]->Fill(TrueDeltaAlphaT,weight);										
					TrueDeltaPT_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaPT_InProtonCosThetaTwoDPlot[0][ProtonCosThetaTwoDIndex]->Fill(TrueTransMissMomentum,weight);					
					TrueMuonMomentum_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(TrueMuonMomentum_GeV,weight);
					TrueProtonCosTheta_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(TrueProtonCosTheta,weight);					
					TrueProtonMomentum_InProtonCosThetaTwoDPlot[0][ProtonCosThetaTwoDIndex]->Fill(TrueProtonMomentum_GeV,weight);
					TrueDeltaPty_InDeltaPtxTwoDPlot[0][DeltaPtxTwoDIndex]->Fill(TruePty,weight);
					TrueDeltaPtx_InDeltaPtyTwoDPlot[0][DeltaPtyTwoDIndex]->Fill(TruePtx,weight);	
					TrueDeltaPT_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[0][DeltaAlpha3DqTwoDIndex]->Fill(TruePn,weight);
					TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[0][DeltaAlpha3DMuTwoDIndex]->Fill(TruePn,weight);										
					TrueProtonMomentum_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(TrueProtonMomentum_GeV,weight);					
					TrueECal_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPtxTwoDPlot[0][DeltaPtxTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPtyTwoDPlot[0][DeltaPtyTwoDIndex]->Fill(TrueRecoECal,weight);										

					TrueDeltaAlphaT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(TrueDeltaAlpha3Dq,weight);
					TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(TrueDeltaAlpha3DMu,weight);															
					TrueDeltaPhiT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
					TrueDeltaPn_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TruePn,weight);
					TrueDeltaPn_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(TruePn,weight);					
					TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlphaT_InMuonMomentumTwoDPlot[genie_mode][MuonMomentumTwoDIndex]->Fill(TrueDeltaAlphaT,weight);					
					TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaAlphaT_InProtonMomentumTwoDPlot[genie_mode][ProtonMomentumTwoDIndex]->Fill(TrueDeltaAlphaT,weight);										
					TrueDeltaPT_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaPT_InProtonCosThetaTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex]->Fill(TrueTransMissMomentum,weight);					
					TrueMuonMomentum_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(TrueMuonMomentum_GeV,weight);
					TrueProtonCosTheta_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(TrueProtonCosTheta,weight);					
					TrueProtonMomentum_InProtonCosThetaTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex]->Fill(TrueProtonMomentum_GeV,weight);
					TrueDeltaPty_InDeltaPtxTwoDPlot[genie_mode][DeltaPtxTwoDIndex]->Fill(TruePty,weight);
					TrueDeltaPtx_InDeltaPtyTwoDPlot[genie_mode][DeltaPtyTwoDIndex]->Fill(TruePtx,weight);	
					TrueDeltaPT_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(TrueTransMissMomentum,weight);
					TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[genie_mode][DeltaAlpha3DqTwoDIndex]->Fill(TruePn,weight);
					TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[genie_mode][DeltaAlpha3DMuTwoDIndex]->Fill(TruePn,weight);										
					TrueProtonMomentum_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(TrueProtonMomentum_GeV,weight);					
					TrueECal_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TrueRecoECal,weight);	
					TrueECal_InDeltaPtxTwoDPlot[genie_mode][DeltaPtxTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPtyTwoDPlot[genie_mode][DeltaPtyTwoDIndex]->Fill(TrueRecoECal,weight);									

					//----------------------------------------//										

					// 3D analysis (uncorrelated)

					TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[0][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPtxDeltaPtyTwoDPlot[0][DeltaPtxTwoDIndex][DeltaPtyTwoDIndex]->Fill(TrueRecoECal,weight);					
					TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[0][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[0][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);

					TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[genie_mode][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InDeltaPtxDeltaPtyTwoDPlot[genie_mode][DeltaPtxTwoDIndex][DeltaPtyTwoDIndex]->Fill(TrueRecoECal,weight);					
					TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[genie_mode][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);
					TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(TrueRecoECal,weight);					

					//----------------------------------------//																					

					// 2D analysis in 1D grid

					SerialTrueMuonMomentum_InMuonCosThetaPlot[0]->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialTrueProtonMomentum_InProtonCosThetaPlot[0]->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InMuonCosThetaPlot[0]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InMuonMomentumPlot[0]->Fill(SerialDeltaAlphaTInMuonMomentumIndex,weight);					
					SerialTrueDeltaAlphaT_InProtonCosThetaPlot[0]->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InProtonMomentumPlot[0]->Fill(SerialDeltaAlphaTInProtonMomentumIndex,weight);					
					SerialTrueDeltaAlphaT_InDeltaPTPlot[0]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
					SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);										
					SerialTrueDeltaPhiT_InDeltaPTPlot[0]->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialTrueDeltaPn_InDeltaPTPlot[0]->Fill(SerialDeltaPnInDeltaPTIndex,weight);
					SerialTrueDeltaPn_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);						
					SerialTrueECal_InDeltaPTPlot[0]->Fill(SerialECalInDeltaPTIndex,weight);
					SerialTrueECal_InDeltaPtxPlot[0]->Fill(SerialECalInDeltaPtxIndex,weight);
					SerialTrueECal_InDeltaPtyPlot[0]->Fill(SerialECalInDeltaPtyIndex,weight);																																						
					SerialTrueProtonCosTheta_InMuonCosThetaPlot[0]->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialTrueDeltaPty_InDeltaPtxPlot[0]->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialTrueDeltaPtx_InDeltaPtyPlot[0]->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialTrueDeltaPT_InMuonCosThetaPlot[0]->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialTrueDeltaPT_InProtonCosThetaPlot[0]->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);
					SerialTrueDeltaPT_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialTrueDeltaPn_InDeltaAlpha3DqPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
					SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);										
					SerialTrueProtonMomentum_InDeltaAlphaTPlot[0]->Fill(SerialProtonMomentumInDeltaAlphaTIndex,weight);					
					SerialTrueECal_InDeltaAlphaTPlot[0]->Fill(SerialECalInDeltaAlphaTIndex,weight);

					SerialTrueMuonMomentum_InMuonCosThetaPlot[genie_mode]->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialTrueProtonMomentum_InProtonCosThetaPlot[genie_mode]->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InMuonCosThetaPlot[genie_mode]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InMuonMomentumPlot[genie_mode]->Fill(SerialDeltaAlphaTInMuonMomentumIndex,weight);					
					SerialTrueDeltaAlphaT_InProtonCosThetaPlot[genie_mode]->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialTrueDeltaAlphaT_InProtonMomentumPlot[genie_mode]->Fill(SerialDeltaAlphaTInProtonMomentumIndex,weight);					
					SerialTrueDeltaAlphaT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
					SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);										
					SerialTrueDeltaPhiT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialTrueDeltaPn_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaPnInDeltaPTIndex,weight);
					SerialTrueDeltaPn_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);						
					SerialTrueECal_InDeltaPTPlot[genie_mode]->Fill(SerialECalInDeltaPTIndex,weight);
					SerialTrueECal_InDeltaPtxPlot[genie_mode]->Fill(SerialECalInDeltaPtxIndex,weight);
					SerialTrueECal_InDeltaPtyPlot[genie_mode]->Fill(SerialECalInDeltaPtyIndex,weight);																												
					SerialTrueProtonCosTheta_InMuonCosThetaPlot[genie_mode]->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialTrueDeltaPty_InDeltaPtxPlot[genie_mode]->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialTrueDeltaPtx_InDeltaPtyPlot[genie_mode]->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialTrueDeltaPT_InMuonCosThetaPlot[genie_mode]->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialTrueDeltaPT_InProtonCosThetaPlot[genie_mode]->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);
					SerialTrueDeltaPT_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialTrueDeltaPn_InDeltaAlpha3DqPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
					SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);										
					SerialTrueProtonMomentum_InDeltaAlphaTPlot[genie_mode]->Fill(SerialProtonMomentumInDeltaAlphaTIndex,weight);					
					SerialTrueECal_InDeltaAlphaTPlot[genie_mode]->Fill(SerialECalInDeltaAlphaTIndex,weight);					

					//----------------------------------------//						

					// 3D analysis treated in 1D grid

					SerialTrueECal_InDeltaPTDeltaAlphaTPlot[0]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialTrueECal_InDeltaPtxDeltaPtyPlot[0]->Fill(SerialECalInDeltaPtxDeltaPtyIndex,weight);					
					SerialTrueECal_InMuonCosThetaMuonMomentumPlot[0]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialTrueECal_InProtonCosThetaProtonMomentumPlot[0]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);

					SerialTrueECal_InDeltaPTDeltaAlphaTPlot[genie_mode]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialTrueECal_InDeltaPtxDeltaPtyPlot[genie_mode]->Fill(SerialECalInDeltaPtxDeltaPtyIndex,weight);					
					SerialTrueECal_InMuonCosThetaMuonMomentumPlot[genie_mode]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialTrueECal_InProtonCosThetaProtonMomentumPlot[genie_mode]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);					

					//----------------------------------------//									
					
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

						TrueCCQEECalPlot[0]->Fill(TrueRecoECal,weight);
						TrueCCQEQ2Plot[0]->Fill(RecoTrueQ2,weight);
						TrueCCQEMuonMomentumPlot[0]->Fill(TrueMuonMomentum_GeV,weight);
						TrueCCQEMuonPhiPlot[0]->Fill(TrueMuonPhi_Deg,weight);
						TrueCCQEMuonCosThetaPlot[0]->Fill(TrueMuonCosTheta,weight);
						TrueCCQEProtonMomentumPlot[0]->Fill(TrueProtonMomentum_GeV,weight);
						TrueCCQEProtonPhiPlot[0]->Fill(TrueProtonPhi_Deg,weight);
						TrueCCQEProtonCosThetaPlot[0]->Fill(TrueProtonCosTheta,weight);

						TrueCCQEECalPlot[genie_mode]->Fill(TrueRecoECal,weight);
						TrueCCQEQ2Plot[genie_mode]->Fill(RecoTrueQ2,weight);
						TrueCCQEMuonMomentumPlot[genie_mode]->Fill(TrueMuonMomentum_GeV,weight);
						TrueCCQEMuonPhiPlot[genie_mode]->Fill(TrueMuonPhi_Deg,weight);
						TrueCCQEMuonCosThetaPlot[genie_mode]->Fill(TrueMuonCosTheta,weight);
						TrueCCQEProtonMomentumPlot[genie_mode]->Fill(TrueProtonMomentum_GeV,weight);
						TrueCCQEProtonPhiPlot[genie_mode]->Fill(TrueProtonPhi_Deg,weight);
						TrueCCQEProtonCosThetaPlot[genie_mode]->Fill(TrueProtonCosTheta,weight);						

					}					

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