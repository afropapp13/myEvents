#define PeLEE_myRecoAnalysis_cxx
#include "PeLEE_myRecoAnalysis.h"
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
#include <TProfile2D.h>
#include <TF1.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include "ubana/myClasses/Tools.h"
#include "ubana/myClasses/STV_Tools.h"

//----------------------------------------//

TString ToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

void PeLEE_myRecoAnalysis::Loop() {

	//----------------------------------------//

	// Function for proton momentum recalibration in %

	TF1* fPP = new TF1("fPP","29.354172*x-14.674918",0.3,0.5);	

	//----------------------------------------//

	int TotalCounter = 0;
	int ContainmentCounter = 0;
	int NoFlippedTrackCounter = 0;
	int ContainedVertexCounter = 0;
	int MuonQualityCutCounter = 0;
	int PIDCounter = 0;
	int NuScoreCounter = 0;
	int KinematicsCounter = 0;

	//----------------------------------------//

	int NEventsPassingSelectionCuts = 0;
	int CC1pEventsPassingSelectionCuts = 0;
	int NonCC1pEventsPassingSelectionCuts = 0;
	int ManualNonCC1pEventsPassingSelectionCuts = 0;
	int CC1p1piEventsPassingSelectionCuts = 0;
	int CC2pEventsPassingSelectionCuts = 0;
	int CC2p1piEventsPassingSelectionCuts = 0;
	int CCNpXpiEventsPassingSelectionCuts = 0;
	int CC3pEventsPassingSelectionCuts = 0;
	int CC3p1piEventsPassingSelectionCuts = 0;
	int CC3p2piEventsPassingSelectionCuts = 0;
	int CC3p3piEventsPassingSelectionCuts = 0;
	int CC4p0piEventsPassingSelectionCuts = 0;
	int CC0pXpiEventsPassingSelectionCuts = 0;
	int CC5pXpiEventsPassingSelectionCuts = 0;
	int CC6pXpiEventsPassingSelectionCuts = 0;
	int CC7pXpiEventsPassingSelectionCuts = 0;
	int CC8pXpiEventsPassingSelectionCuts = 0;
	int MisIndetifiedMuonAsPion = 0;
	int BrokenProtonTrack = 0;
	int MisIndetifiedMuonAsProton = 0;
	int MisIndetifiedMuonAsAntiMuon = 0;
	int MisIndetifiedProtonAsPion = 0;
	int MisIndetifiedProtonAsDeuterium = 0;
	int DeuteriumPionPairs = 0;
	int MisIndetifiedMuonAsDeuterium = 0;
	int MisIndetifiedMuonAsHelium = 0;
	int MisIndetifiedProtonAsHelium = 0;
	int MisIndetifiedProtonAsElectron = 0;
	int MisIndetifiedMuPToElectronElectron = 0;
	int MisIndetifiedMuPToMuMu = 0;
	int BrokenMuonTrack = 0;
	int MisIndetifiedMuPToPiPi = 0;
	int CandidateMuon_MCParticle_OutFV = 0;
	int CandidateProton_MCParticle_OutFV = 0;
	int MultipleVertices = 0;
	int InTimeCosmics = 0;
	int NeutronPairCounter = 0;
	int NCEvents = 0;
	int DoubleMuonEvents = 0;
	int FlippedMuPEvents = 0;
	int ArArEvents = 0;
	int ProtonNeutronEvents = 0;
	int ProtonKaonEvents = 0;
	int pi0Included = 0;
	int OutCommonRange = 0;
	int OtherMCBkg = 0;

	//----------------------------------------//
		
	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID_NuScore");

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

		// Txt file to keep track of the run/subrun/event of the candidate events

		TString RunTxtName = "/uboone/data/users/apapadop/myEvents/myTxtFiles/"+UBCodeVersion+"/TxtmyRunSubRunEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
		ofstream myRunTxtFile;
		myRunTxtFile.open(RunTxtName);
		myRunTxtFile << std::fixed << std::setprecision(2);
		myRunTxtFile << fWhichSample << endl << endl;

		//----------------------------------------//

		// 1D Reco Level Plots

		TH1D* RecoPi0Plot = new TH1D("RecoPi0Plot",";# #pi^{0}",4,-0.5,3.5);
		TH1D* RecoNeutronPlot = new TH1D("RecoNeutronPlot",";# neutrons",6,-0.5,5.5);

		TH1D* RecoNuPlot = new TH1D("RecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* RecoEvPlot = new TH1D("RecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* RecoNuScorePlot = new TH1D("RecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* RecoFlashScorePlot = new TH1D("RecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* RecoDistancePlot = new TH1D("RecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* RecoLengthDifferencePlot = new TH1D("RecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* RecoVertexXPlot = new TH1D("RecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* RecoVertexYPlot = new TH1D("RecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* RecoVertexZPlot = new TH1D("RecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* RecoCandidateMuEndYPlot = new TH1D("RecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* RecoMuonLLRPIDPlot = new TH1D("RecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* RecoProtonLLRPIDPlot = new TH1D("RecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* RecoMuonLengthPlot = new TH1D("RecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* RecodMuonTracksScorePlot = new TH1D("RecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* RecodMuonVertexDistancePlot = new TH1D("RecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* RecoContainedMuonLengthPlot = new TH1D("RecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* RecoUncontainedMuonLengthPlot = new TH1D("RecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* RecoMuonCosThetaSingleBinPlot = new TH1D("RecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* RecoMuonPhiPlot = new TH1D("RecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* RecoCCQEMuonMomentumPlot = new TH1D("RecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* RecoCCQEMuonCosThetaPlot = new TH1D("RecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* RecoCCQEMuonPhiPlot = new TH1D("RecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* RecoContainedMuonMomentumPlot = new TH1D("RecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* RecoUncontainedMuonMomentumPlot = new TH1D("RecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);

		TH1D* RecoProtonLengthPlot = new TH1D("RecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);
//		TH1D* RecodProtonTracksScorePlot = new TH1D("RecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* RecodProtonVertexDistancePlot = new TH1D("RecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);

		TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* RecoProtonPhiPlot = new TH1D("RecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* RecoCCQEProtonMomentumPlot = new TH1D("RecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* RecoCCQEProtonCosThetaPlot = new TH1D("RecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* RecoCCQEProtonPhiPlot = new TH1D("RecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* RecoDeltaThetaPlot = new TH1D("RecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* RecoDeltaForwardThetaPlot = new TH1D("RecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* RecoDeltaBackwardThetaPlot = new TH1D("RecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* RecoDeltaPhiPlot = new TH1D("RecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* RecokMissPlot = new TH1D("RecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* RecoPMissMinusPlot = new TH1D("RecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* RecoPMissPlot = new TH1D("RecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* RecoDeltaPLPlot = new TH1D("RecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* RecoDeltaPnPlot = new TH1D("RecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* RecoDeltaPtxPlot = new TH1D("RecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* RecoDeltaPtyPlot = new TH1D("RecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* RecoAPlot = new TH1D("RecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* RecoDeltaPTPlot = new TH1D("RecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* RecoDeltaAlphaTPlot = new TH1D("RecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* RecoDeltaPhiTPlot = new TH1D("RecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* RecoEQEPlot = new TH1D("RecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* RecoQ2Plot = new TH1D("RecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* RecoCCQEECalPlot = new TH1D("RecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* RecoCCQEQ2Plot = new TH1D("RecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
		// 2D Analysis

		// For now and until MicroBooNE box opening
		int NBins2DAnalysis = 4;
		
		TH2D* RecoCosThetaMuPmuPlot = new TH2D("RecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* RecoCosThetaPPpPlot = new TH2D("RecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);					

		// 2D Plot for Default Chi2 vs 3-Plane Chi2

		TH2D* RecoChi2vsThreePlaneChi2TPlot = new TH2D("RecoChi2vsThreePlaneChi2TPlot",RecoLabelXAxisChi2+RecoLabelXAxisThreePlaneChi2,
			NBinsChi2,MinChi2,MaxChi2,NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);

		// Momentum vs Length

		TH2D* RecoMuonMomentumVsLengthPlot = new TH2D("RecoMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);
		TH2D* RecoContainedMuonMomentumVsLengthPlot = new TH2D("RecoContainedMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);
		TH2D* RecoUncontainedMuonMomentumVsLengthPlot = new TH2D("RecoUncontainedMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D True Level Plots for Signal CC1p reconstructed candidates

		TH1D* CC1pTrueNuPlot = new TH1D("CC1pTrueNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CC1pTrueEvPlot = new TH1D("CC1pTrueEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CC1pTrueVertexXPlot = new TH1D("CC1pTrueVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CC1pTrueVertexYPlot = new TH1D("CC1pTrueVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CC1pTrueVertexZPlot = new TH1D("CC1pTrueVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

		TH1D* CC1pTrueMuonMomentumPlot = new TH1D("CC1pTrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pTrueMuonCosThetaSingleBinPlot = new TH1D("CC1pTrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* CC1pTrueMuonPhiPlot = new TH1D("CC1pTrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CC1pTrueCCQEMuonMomentumPlot = new TH1D("CC1pTrueCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CC1pTrueCCQEMuonCosThetaPlot = new TH1D("CC1pTrueCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CC1pTrueCCQEMuonPhiPlot = new TH1D("CC1pTrueCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* CC1pTrueProtonMomentumPlot = new TH1D("CC1pTrueProtonMomentumPlot",LabelXAxisProtonMomentum, NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CC1pTrueProtonCosThetaPlot = new TH1D("CC1pTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta, NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CC1pTrueProtonPhiPlot = new TH1D("CC1pTrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CC1pTrueCCQEProtonMomentumPlot = new TH1D("CC1pTrueCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CC1pTrueCCQEProtonCosThetaPlot = new TH1D("CC1pTrueCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CC1pTrueCCQEProtonPhiPlot = new TH1D("CC1pTrueCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CC1pTrueDeltaPTPlot = new TH1D("CC1pTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pTrueDeltaAlphaTPlot = new TH1D("CC1pTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pTrueDeltaPhiTPlot = new TH1D("CC1pTrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CC1pTrueDeltaPLPlot = new TH1D("CC1pTrueDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* CC1pTrueDeltaPnPlot = new TH1D("CC1pTrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CC1pTrueDeltaPtxPlot = new TH1D("CC1pTrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* CC1pTrueDeltaPtyPlot = new TH1D("CC1pTrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* CC1pTrueAPlot = new TH1D("CC1pTrueAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);
		TH1D* CC1pTruekMissPlot = new TH1D("CC1pTruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CC1pTruePMissPlot = new TH1D("CC1pTruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);
		TH1D* CC1pTruePMissMinusPlot = new TH1D("CC1pTruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);

		TH1D* CC1pTrueECalPlot = new TH1D("CC1pTrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CC1pTrueEQEPlot = new TH1D("CC1pTrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CC1pTrueQ2Plot = new TH1D("CC1pTrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CC1pTrueCCQEECalPlot = new TH1D("CC1pTrueCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CC1pTrueCCQEQ2Plot = new TH1D("CC1pTrueCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for Signal CC1p

		TH1D* CC1pRecoNuPlot = new TH1D("CC1pRecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CC1pRecoEvPlot = new TH1D("CC1pRecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CC1pRecoNuScorePlot = new TH1D("CC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CC1pRecoFlashScorePlot = new TH1D("CC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* CC1pRecoDistancePlot = new TH1D("CC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* CC1pRecoLengthDifferencePlot = new TH1D("CC1pRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* CC1pRecoVertexXPlot = new TH1D("CC1pRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CC1pRecoVertexYPlot = new TH1D("CC1pRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CC1pRecoVertexZPlot = new TH1D("CC1pRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* CC1pRecoCandidateMuEndYPlot = new TH1D("CC1pRecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* CC1pRecoMuonLLRPIDPlot = new TH1D("CC1pRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* CC1pRecoProtonLLRPIDPlot = new TH1D("CC1pRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CC1pRecoMuonLengthPlot = new TH1D("CC1pRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CC1pRecodMuonTracksScorePlot = new TH1D("CC1pRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CC1pRecodMuonVertexDistancePlot = new TH1D("CC1pRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);	

		TH1D* CC1pRecoContainedMuonLengthPlot = new TH1D("CC1pRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CC1pRecoUncontainedMuonLengthPlot = new TH1D("CC1pRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoMuonCosThetaSingleBinPlot = new TH1D("CC1pRecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* CC1pRecoMuonPhiPlot = new TH1D("CC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CC1pRecoCCQEMuonMomentumPlot = new TH1D("CC1pRecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CC1pRecoCCQEMuonCosThetaPlot = new TH1D("CC1pRecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoCCQEMuonPhiPlot = new TH1D("CC1pRecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* CC1pRecoContainedMuonMomentumPlot = new TH1D("CC1pRecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pRecoUncontainedMuonMomentumPlot = new TH1D("CC1pRecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);

		TH1D* CC1pRecoProtonLengthPlot = new TH1D("CC1pRecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);
//		TH1D* CC1pRecodProtonTracksScorePlot = new TH1D("CC1pRecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CC1pRecodProtonVertexDistancePlot = new TH1D("CC1pRecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);


		TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CC1pRecoProtonPhiPlot = new TH1D("CC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CC1pRecoCCQEProtonMomentumPlot = new TH1D("CC1pRecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,
			CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CC1pRecoCCQEProtonCosThetaPlot = new TH1D("CC1pRecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CC1pRecoCCQEProtonPhiPlot = new TH1D("CC1pRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CC1pRecoDeltaThetaPlot = new TH1D("CC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CC1pRecoDeltaForwardThetaPlot = new TH1D("CC1pRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CC1pRecoDeltaBackwardThetaPlot = new TH1D("CC1pRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CC1pRecoDeltaPhiPlot = new TH1D("CC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CC1pRecokMissPlot = new TH1D("CC1pRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CC1pRecoPMissMinusPlot = new TH1D("CC1pRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CC1pRecoPMissPlot = new TH1D("CC1pRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CC1pRecoDeltaPLPlot = new TH1D("CC1pRecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* CC1pRecoDeltaPnPlot = new TH1D("CC1pRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CC1pRecoDeltaPtxPlot = new TH1D("CC1pRecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* CC1pRecoDeltaPtyPlot = new TH1D("CC1pRecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* CC1pRecoAPlot = new TH1D("CC1pRecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CC1pRecoEQEPlot = new TH1D("CC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CC1pRecoQ2Plot = new TH1D("CC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CC1pRecoCCQEECalPlot = new TH1D("CC1pRecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CC1pRecoCCQEQ2Plot = new TH1D("CC1pRecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);

		TProfile2D* CC1pTProfProtonMomTheta = new TProfile2D("CC1pTProfProtonMomTheta",LabelXAxisProtonMomentum+";#theta_{p}; P_{p} Fractional Bias [%]",NBinsProtonMomentum,ArrayNBinsProtonMomentum,12,0.,180.);
		TProfile2D* CC1pTProfProtonMomPhi = new TProfile2D("CC1pTProfProtonMomPhi",LabelXAxisProtonMomentum+";#phi_{p}; P_{p} Fractional Bias [%]",NBinsProtonMomentum,ArrayNBinsProtonMomentum,12,-180.,180.);

		// 2D Analysis
		
		TH2D* CC1pRecoCosThetaMuPmuPlot = new TH2D("CC1pRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CC1pRecoCosThetaPPpPlot = new TH2D("CC1pRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// 2D Reco Level Plots for Signal CC1p, unweighted

		TH2D* CC1pRecoMuonMomentumPlot2D = new TH2D("CC1pRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
			ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH2D* CC1pRecoProtonMomentumPlot2D = new TH2D("CC1pRecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,NBinsProtonMomentum,
			ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH2D* CC1pRecoCCQEMuonMomentumPlot2D = new TH2D("CC1pRecoCCQEMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,CCQENBinsMuonMomentum,
			CCQEArrayNBinsMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH2D* CC1pRecoCCQEProtonMomentumPlot2D = new TH2D("CC1pRecoCCQEProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,CCQENBinsProtonMomentum,
			CCQEArrayNBinsProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);

		TH2D* CC1pRecoMuonCosThetaPlot2D = new TH2D("CC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("CC1pRecoMuonCosThetaSingleBinPlot2D",LabelXAxisMuonCosTheta2D,1,-1.,1.,1,-1.,1.);
		TH2D* CC1pRecoProtonCosThetaPlot2D = new TH2D("CC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
			ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH2D* CC1pRecoCCQEMuonCosThetaPlot2D = new TH2D("CC1pRecoCCQEMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,CCQENBinsMuonCosTheta,
			CCQEArrayNBinsMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoCCQEProtonCosThetaPlot2D = new TH2D("CC1pRecoCCQEProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,CCQENBinsProtonCosTheta,
			CCQEArrayNBinsProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);

		TH2D* CC1pRecoMuonPhiPlot2D = new TH2D("CC1pRecoMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,
			ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH2D* CC1pRecoProtonPhiPlot2D = new TH2D("CC1pRecoProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,
			ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH2D* CC1pRecoCCQEMuonPhiPlot2D = new TH2D("CC1pRecoCCQEMuonPhiPlot2D",LabelXAxisMuonPhi2D,CCQENBinsMuonPhi,
			CCQEArrayNBinsMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);
		TH2D* CC1pRecoCCQEProtonPhiPlot2D = new TH2D("CC1pRecoCCQEProtonPhiPlot2D",LabelXAxisProtonPhi2D,CCQENBinsProtonPhi,
			CCQEArrayNBinsProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH2D* CC1pRecoDeltaPTPlot2D = new TH2D("CC1pRecoDeltaPTPlot2D",LabelXAxisDeltaPT2D,NBinsDeltaPT,
			ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH2D* CC1pRecoDeltaAlphaTPlot2D = new TH2D("CC1pRecoDeltaAlphaTPlot2D",LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,
			ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH2D* CC1pRecoDeltaPhiTPlot2D = new TH2D("CC1pRecoDeltaPhiTPlot2D",LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,
			ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH2D* CC1pRecoECalPlot2D = new TH2D("CC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
		TH2D* CC1pRecoEQEPlot2D = new TH2D("CC1pRecoEQEPlot2D",LabelXAxisEQE2D,NBinsEQE,ArrayNBinsEQE,NBinsEQE,ArrayNBinsEQE);
		TH2D* CC1pRecoQ2Plot2D = new TH2D("CC1pRecoQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

		TH2D* CC1pRecoCCQEECalPlot2D = new TH2D("CC1pRecoCCQEECalPlot2D",LabelXAxisECal2D,CCQENBinsECal,CCQEArrayNBinsECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH2D* CC1pRecoCCQEQ2Plot2D = new TH2D("CC1pRecoCCQEQ2Plot2D",LabelXAxisQ22D,CCQENBinsQ2,CCQEArrayNBinsQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);

		TH2D* CC1pRecoDeltaPLPlot2D = new TH2D("CC1pRecoDeltaPLPlot2D",LabelXAxisDeltaPL2D,NBinsDeltaPL,ArrayNBinsDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH2D* CC1pRecoDeltaPnPlot2D = new TH2D("CC1pRecoDeltaPnPlot2D",LabelXAxisDeltaPn2D,NBinsDeltaPn,ArrayNBinsDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH2D* CC1pRecoDeltaPtxPlot2D = new TH2D("CC1pRecoDeltaPtxPlot2D",LabelXAxisDeltaPtx2D,NBinsDeltaPtx,ArrayNBinsDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH2D* CC1pRecoDeltaPtyPlot2D = new TH2D("CC1pRecoDeltaPtyPlot2D",LabelXAxisDeltaPty2D,NBinsDeltaPty,ArrayNBinsDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH2D* CC1pRecoAPlot2D = new TH2D("CC1pRecoAPlot2D",LabelXAxisA2D,NBinsA,ArrayNBinsA,NBinsA,ArrayNBinsA);

		TH2D* CC1pRecokMissPlot2D = new TH2D("CC1pRecokMissPlot2D",LabelXAxiskMiss2D,NBinskMiss,ArrayNBinskMiss,NBinskMiss,ArrayNBinskMiss);
		TH2D* CC1pRecoPMissPlot2D = new TH2D("CC1pRecoPMissPlot2D",LabelXAxisPMiss2D,NBinsPMiss,ArrayNBinsPMiss,NBinsPMiss,ArrayNBinsPMiss);
		TH2D* CC1pRecoPMissMinusPlot2D = new TH2D("CC1pRecoPMissMinusPlot2D",LabelXAxisPMissMinus2D,NBinsPMissMinus,ArrayNBinsPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);

		// -------------------------------------------------------------------------------------------------------------------------------------

		// 2D Reco Level Plots for Signal CC1p, POT Scaled for response matrices

		TH2D* POTScaledCC1pRecoMuonMomentumPlot2D = new TH2D("POTScaledCC1pRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
			ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH2D* POTScaledCC1pRecoProtonMomentumPlot2D = new TH2D("POTScaledCC1pRecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,NBinsProtonMomentum,
			ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH2D* POTScaledCC1pRecoCCQEMuonMomentumPlot2D = new TH2D("POTScaledCC1pRecoCCQEMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,CCQENBinsMuonMomentum,
			CCQEArrayNBinsMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH2D* POTScaledCC1pRecoCCQEProtonMomentumPlot2D = new TH2D("POTScaledCC1pRecoCCQEProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,CCQENBinsProtonMomentum,
			CCQEArrayNBinsProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);

		TH2D* POTScaledCC1pRecoMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D",LabelXAxisMuonCosTheta2D,1,-1.,1.,1,-1.,1.);
		TH2D* POTScaledCC1pRecoProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
			ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH2D* POTScaledCC1pRecoCCQEMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoCCQEMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,CCQENBinsMuonCosTheta,
			CCQEArrayNBinsMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoCCQEProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoCCQEProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,CCQENBinsProtonCosTheta,
			CCQEArrayNBinsProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);

		TH2D* POTScaledCC1pRecoMuonPhiPlot2D = new TH2D("POTScaledCC1pRecoMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,
			ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH2D* POTScaledCC1pRecoProtonPhiPlot2D = new TH2D("POTScaledCC1pRecoProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,
			ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH2D* POTScaledCC1pRecoCCQEMuonPhiPlot2D = new TH2D("POTScaledCC1pRecoCCQEMuonPhiPlot2D",LabelXAxisMuonPhi2D,CCQENBinsMuonPhi,
			CCQEArrayNBinsMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);
		TH2D* POTScaledCC1pRecoCCQEProtonPhiPlot2D = new TH2D("POTScaledCC1pRecoCCQEProtonPhiPlot2D",LabelXAxisProtonPhi2D,CCQENBinsProtonPhi,
			CCQEArrayNBinsProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH2D* POTScaledCC1pRecoDeltaPTPlot2D = new TH2D("POTScaledCC1pRecoDeltaPTPlot2D",LabelXAxisDeltaPT2D,NBinsDeltaPT,
			ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH2D* POTScaledCC1pRecoDeltaAlphaTPlot2D = new TH2D("POTScaledCC1pRecoDeltaAlphaTPlot2D",LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,
			ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH2D* POTScaledCC1pRecoDeltaPhiTPlot2D = new TH2D("POTScaledCC1pRecoDeltaPhiTPlot2D",LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,
			ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH2D* POTScaledCC1pRecoECalPlot2D = new TH2D("POTScaledCC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
		TH2D* POTScaledCC1pRecoEQEPlot2D = new TH2D("POTScaledCC1pRecoEQEPlot2D",LabelXAxisEQE2D,NBinsEQE,ArrayNBinsEQE,NBinsEQE,ArrayNBinsEQE);
		TH2D* POTScaledCC1pRecoQ2Plot2D = new TH2D("POTScaledCC1pRecoQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

		TH2D* POTScaledCC1pRecoCCQEECalPlot2D = new TH2D("POTScaledCC1pRecoCCQEECalPlot2D",LabelXAxisECal2D,CCQENBinsECal,CCQEArrayNBinsECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH2D* POTScaledCC1pRecoCCQEQ2Plot2D = new TH2D("POTScaledCC1pRecoCCQEQ2Plot2D",LabelXAxisQ22D,CCQENBinsQ2,CCQEArrayNBinsQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);

		TH2D* POTScaledCC1pRecoDeltaPLPlot2D = new TH2D("POTScaledCC1pRecoDeltaPLPlot2D",LabelXAxisDeltaPL2D,NBinsDeltaPL,ArrayNBinsDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH2D* POTScaledCC1pRecoDeltaPnPlot2D = new TH2D("POTScaledCC1pRecoDeltaPnPlot2D",LabelXAxisDeltaPn2D,NBinsDeltaPn,ArrayNBinsDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH2D* POTScaledCC1pRecoDeltaPtxPlot2D = new TH2D("POTScaledCC1pRecoDeltaPtxPlot2D",LabelXAxisDeltaPtx2D,NBinsDeltaPtx,ArrayNBinsDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH2D* POTScaledCC1pRecoDeltaPtyPlot2D = new TH2D("POTScaledCC1pRecoDeltaPtyPlot2D",LabelXAxisDeltaPty2D,NBinsDeltaPty,ArrayNBinsDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH2D* POTScaledCC1pRecoAPlot2D = new TH2D("POTScaledCC1pRecoAPlot2D",LabelXAxisA2D,NBinsA,ArrayNBinsA,NBinsA,ArrayNBinsA);

		TH2D* POTScaledCC1pRecokMissPlot2D = new TH2D("POTScaledCC1pRecokMissPlot2D",LabelXAxiskMiss2D,NBinskMiss,ArrayNBinskMiss,NBinskMiss,ArrayNBinskMiss);
		TH2D* POTScaledCC1pRecoPMissPlot2D = new TH2D("POTScaledCC1pRecoPMissPlot2D",LabelXAxisPMiss2D,NBinsPMiss,ArrayNBinsPMiss,NBinsPMiss,ArrayNBinsPMiss);
		TH2D* POTScaledCC1pRecoPMissMinusPlot2D = new TH2D("POTScaledCC1pRecoPMissMinusPlot2D",LabelXAxisPMissMinus2D,NBinsPMissMinus,ArrayNBinsPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);

		// -------------------------------------------------------------------------------------------------------------------------------------

		TH2D* CC1pRecoMuonMomentumVsLengthPlot = new TH2D("CC1pRecoMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);
		TH2D* CC1pRecoContainedMuonMomentumVsLengthPlot = new TH2D("CC1pRecoContainedMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);
		TH2D* CC1pRecoUncontainedMuonMomentumVsLengthPlot = new TH2D("CC1pRecoUncontainedMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);

		// -------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for non-CC1p

		TH1D* NonCC1pRecoNuPlot = new TH1D("NonCC1pRecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* NonCC1pRecoEvPlot = new TH1D("NonCC1pRecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* NonCC1pRecoNuScorePlot = new TH1D("NonCC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* NonCC1pRecoFlashScorePlot = new TH1D("NonCC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,
			NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* NonCC1pRecoDistancePlot = new TH1D("NonCC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* NonCC1pRecoLengthDifferencePlot = new TH1D("NonCC1pRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* NonCC1pRecoVertexXPlot = new TH1D("NonCC1pRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* NonCC1pRecoVertexYPlot = new TH1D("NonCC1pRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* NonCC1pRecoVertexZPlot = new TH1D("NonCC1pRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* NonCC1pRecoCandidateMuEndYPlot = new TH1D("NonCC1pRecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* NonCC1pRecoMuonLLRPIDPlot = new TH1D("NonCC1pRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* NonCC1pRecoProtonLLRPIDPlot = new TH1D("NonCC1pRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* NonCC1pRecoMuonLengthPlot = new TH1D("NonCC1pRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* NonCC1pRecodMuonTracksScorePlot = new TH1D("NonCC1pRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* NonCC1pRecodMuonVertexDistancePlot = new TH1D("NonCC1pRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* NonCC1pRecoContainedMuonLengthPlot = new TH1D("NonCC1pRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* NonCC1pRecoUncontainedMuonLengthPlot = new TH1D("NonCC1pRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* NonCC1pRecoPionMomentumPlot = new TH1D("NonCC1pRecoPionMomentumPlot",";P_{#pi} [GeV/c]",20,0.1,0.6); // GeV/c
		TH1D* NonCC1pRecoMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoMuonCosThetaSingleBinPlot = new TH1D("NonCC1pRecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* NonCC1pRecoMuonPhiPlot = new TH1D("NonCC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* NonCC1pRecoCCQEMuonMomentumPlot = new TH1D("NonCC1pRecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum); // GeV/c
		TH1D* NonCC1pRecoCCQEMuonCosThetaPlot = new TH1D("NonCC1pRecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoCCQEMuonPhiPlot = new TH1D("NonCC1pRecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* NonCC1pRecoContainedMuonMomentumPlot = new TH1D("NonCC1pRecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* NonCC1pRecoUncontainedMuonMomentumPlot = new TH1D("NonCC1pRecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);	

		TH1D* NonCC1pRecoProtonLengthPlot = new TH1D("NonCC1pRecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);
//		TH1D* NonCC1pRecodProtonTracksScorePlot = new TH1D("NonCC1pRecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* NonCC1pRecodProtonVertexDistancePlot = new TH1D("NonCC1pRecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);

		TH1D* NonCC1pRecoProtonMomentumPlot = new TH1D("NonCC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c
		TH1D* NonCC1pRecoProtonCosThetaPlot = new TH1D("NonCC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* NonCC1pRecoProtonPhiPlot = new TH1D("NonCC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* NonCC1pRecoCCQEProtonMomentumPlot = new TH1D("NonCC1pRecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,
			CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum); // GeV/c
		TH1D* NonCC1pRecoCCQEProtonCosThetaPlot = new TH1D("NonCC1pRecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* NonCC1pRecoCCQEProtonPhiPlot = new TH1D("NonCC1pRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* NonCC1pRecoDeltaThetaPlot = new TH1D("NonCC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* NonCC1pRecoDeltaForwardThetaPlot = new TH1D("NonCC1pRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* NonCC1pRecoDeltaBackwardThetaPlot = new TH1D("NonCC1pRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* NonCC1pRecoDeltaPhiPlot = new TH1D("NonCC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* NonCC1pRecokMissPlot = new TH1D("NonCC1pRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* NonCC1pRecoPMissMinusPlot = new TH1D("NonCC1pRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* NonCC1pRecoPMissPlot = new TH1D("NonCC1pRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* NonCC1pRecoDeltaPLPlot = new TH1D("NonCC1pRecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* NonCC1pRecoDeltaPnPlot = new TH1D("NonCC1pRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* NonCC1pRecoDeltaPtxPlot = new TH1D("NonCC1pRecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* NonCC1pRecoDeltaPtyPlot = new TH1D("NonCC1pRecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* NonCC1pRecoAPlot = new TH1D("NonCC1pRecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* NonCC1pRecoDeltaPTPlot = new TH1D("NonCC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* NonCC1pRecoDeltaAlphaTPlot = new TH1D("NonCC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* NonCC1pRecoDeltaPhiTPlot = new TH1D("NonCC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* NonCC1pRecoECalPlot = new TH1D("NonCC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* NonCC1pRecoEQEPlot = new TH1D("NonCC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* NonCC1pRecoQ2Plot = new TH1D("NonCC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* NonCC1pRecoCCQEECalPlot = new TH1D("NonCC1pRecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* NonCC1pRecoCCQEQ2Plot = new TH1D("NonCC1pRecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* NonCC1pRecoCosThetaMuPmuPlot = new TH2D("NonCC1pRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* NonCC1pRecoCosThetaPPpPlot = new TH2D("NonCC1pRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBinsProtonMomentum,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);

		TH2D* NonCC1pRecoMuonMomentumVsLengthPlot = new TH2D("NonCC1pRecoMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);
		TH2D* NonCC1pRecoContainedMuonMomentumVsLengthPlot = new TH2D("NonCC1pRecoContainedMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);
		TH2D* NonCC1pRecoUncontainedMuonMomentumVsLengthPlot = new TH2D("NonCC1pRecoUncontainedMuonMomentumVsLengthPlot",";l_{#mu} [cm];P_{#mu} [GeV/cm]",100,0,100,50,0,0.5);		

		// ---------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCQE

		TH1D* CCQERecoNuPlot = new TH1D("CCQERecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CCQERecoEvPlot = new TH1D("CCQERecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CCQERecoNuScorePlot = new TH1D("CCQERecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCQERecoFlashScorePlot = new TH1D("CCQERecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* CCQERecoDistancePlot = new TH1D("CCQERecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* CCQERecoLengthDifferencePlot = new TH1D("CCQERecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* CCQERecoVertexXPlot = new TH1D("CCQERecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCQERecoVertexYPlot = new TH1D("CCQERecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCQERecoVertexZPlot = new TH1D("CCQERecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* CCQERecoCandidateMuEndYPlot = new TH1D("CCQERecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* CCQERecoMuonLLRPIDPlot = new TH1D("CCQERecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* CCQERecoProtonLLRPIDPlot = new TH1D("CCQERecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCQERecoMuonLengthPlot = new TH1D("CCQERecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCQERecodMuonTracksScorePlot = new TH1D("CCQERecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* CCQERecodMuonVertexDistancePlot = new TH1D("CCQERecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCQERecoContainedMuonLengthPlot = new TH1D("CCQERecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCQERecoUncontainedMuonLengthPlot = new TH1D("CCQERecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCQERecoMuonMomentumPlot = new TH1D("CCQERecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCQERecoMuonCosThetaSingleBinPlot = new TH1D("CCQERecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* CCQERecoMuonPhiPlot = new TH1D("CCQERecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CCQERecoCCQEMuonMomentumPlot = new TH1D("CCQERecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CCQERecoCCQEMuonCosThetaPlot = new TH1D("CCQERecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CCQERecoCCQEMuonPhiPlot = new TH1D("CCQERecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* CCQERecoContainedMuonMomentumPlot = new TH1D("CCQERecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCQERecoUncontainedMuonMomentumPlot = new TH1D("CCQERecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);

		TH1D* CCQERecoProtonLengthPlot = new TH1D("CCQERecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);
//		TH1D* CCQERecodProtonTracksScorePlot = new TH1D("CCQERecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* CCQERecodProtonVertexDistancePlot = new TH1D("CCQERecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);

		TH1D* CCQERecoProtonMomentumPlot = new TH1D("CCQERecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCQERecoProtonCosThetaPlot = new TH1D("CCQERecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCQERecoProtonPhiPlot = new TH1D("CCQERecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCQERecoCCQEProtonMomentumPlot = new TH1D("CCQERecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,
			CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CCQERecoCCQEProtonCosThetaPlot = new TH1D("CCQERecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CCQERecoCCQEProtonPhiPlot = new TH1D("CCQERecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CCQERecoDeltaThetaPlot = new TH1D("CCQERecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCQERecoDeltaForwardThetaPlot = new TH1D("CCQERecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCQERecoDeltaBackwardThetaPlot = new TH1D("CCQERecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCQERecoDeltaPhiPlot = new TH1D("CCQERecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCQERecokMissPlot = new TH1D("CCQERecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCQERecoPMissMinusPlot = new TH1D("CCQERecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCQERecoPMissPlot = new TH1D("CCQERecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCQERecoDeltaPLPlot = new TH1D("CCQERecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* CCQERecoDeltaPnPlot = new TH1D("CCQERecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CCQERecoDeltaPtxPlot = new TH1D("CCQERecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* CCQERecoDeltaPtyPlot = new TH1D("CCQERecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* CCQERecoAPlot = new TH1D("CCQERecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* CCQERecoDeltaPTPlot = new TH1D("CCQERecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQERecoDeltaAlphaTPlot = new TH1D("CCQERecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQERecoDeltaPhiTPlot = new TH1D("CCQERecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCQERecoECalPlot = new TH1D("CCQERecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCQERecoEQEPlot = new TH1D("CCQERecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCQERecoQ2Plot = new TH1D("CCQERecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CCQERecoCCQEECalPlot = new TH1D("CCQERecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CCQERecoCCQEQ2Plot = new TH1D("CCQERecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCQERecoCosThetaMuPmuPlot = new TH2D("CCQERecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCQERecoCosThetaPPpPlot = new TH2D("CCQERecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ---------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCMEC

		TH1D* CCMECRecoNuPlot = new TH1D("CCMECRecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CCMECRecoEvPlot = new TH1D("CCMECRecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CCMECRecoNuScorePlot = new TH1D("CCMECRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCMECRecoFlashScorePlot = new TH1D("CCMECRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* CCMECRecoDistancePlot = new TH1D("CCMECRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* CCMECRecoLengthDifferencePlot = new TH1D("CCMECRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* CCMECRecoVertexXPlot = new TH1D("CCMECRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCMECRecoVertexYPlot = new TH1D("CCMECRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCMECRecoVertexZPlot = new TH1D("CCMECRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* CCMECRecoCandidateMuEndYPlot = new TH1D("CCMECRecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* CCMECRecoMuonLLRPIDPlot = new TH1D("CCMECRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* CCMECRecoProtonLLRPIDPlot = new TH1D("CCMECRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCMECRecoMuonLengthPlot = new TH1D("CCMECRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCMECRecodMuonTracksScorePlot = new TH1D("CCMECRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CCMECRecodMuonVertexDistancePlot = new TH1D("CCMECRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCMECRecoContainedMuonLengthPlot = new TH1D("CCMECRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCMECRecoUncontainedMuonLengthPlot = new TH1D("CCMECRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCMECRecoMuonMomentumPlot = new TH1D("CCMECRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECRecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* CCMECRecoMuonPhiPlot = new TH1D("CCMECRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CCMECRecoCCQEMuonMomentumPlot = new TH1D("CCMECRecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CCMECRecoCCQEMuonCosThetaPlot = new TH1D("CCMECRecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoCCQEMuonPhiPlot = new TH1D("CCMECRecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* CCMECRecoContainedMuonMomentumPlot = new TH1D("CCMECRecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCMECRecoUncontainedMuonMomentumPlot = new TH1D("CCMECRecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);

		TH1D* CCMECRecoProtonLengthPlot = new TH1D("CCMECRecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);	
//		TH1D* CCMECRecodProtonTracksScorePlot = new TH1D("CCMECRecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* CCMECRecodProtonVertexDistancePlot = new TH1D("CCMECRecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);

		TH1D* CCMECRecoProtonMomentumPlot = new TH1D("CCMECRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCMECRecoProtonCosThetaPlot = new TH1D("CCMECRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCMECRecoProtonPhiPlot = new TH1D("CCMECRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCMECRecoCCQEProtonMomentumPlot = new TH1D("CCMECRecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,
			CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CCMECRecoCCQEProtonCosThetaPlot = new TH1D("CCMECRecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CCMECRecoCCQEProtonPhiPlot = new TH1D("CCMECRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CCMECRecoDeltaThetaPlot = new TH1D("CCMECRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaForwardThetaPlot = new TH1D("CCMECRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaBackwardThetaPlot = new TH1D("CCMECRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaPhiPlot = new TH1D("CCMECRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCMECRecokMissPlot = new TH1D("CCMECRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCMECRecoPMissMinusPlot = new TH1D("CCMECRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCMECRecoPMissPlot = new TH1D("CCMECRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCMECRecoDeltaPLPlot = new TH1D("CCMECRecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* CCMECRecoDeltaPnPlot = new TH1D("CCMECRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CCMECRecoDeltaPtxPlot = new TH1D("CCMECRecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* CCMECRecoDeltaPtyPlot = new TH1D("CCMECRecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* CCMECRecoAPlot = new TH1D("CCMECRecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* CCMECRecoDeltaPTPlot = new TH1D("CCMECRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECRecoDeltaAlphaTPlot = new TH1D("CCMECRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECRecoDeltaPhiTPlot = new TH1D("CCMECRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCMECRecoECalPlot = new TH1D("CCMECRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECRecoEQEPlot = new TH1D("CCMECRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCMECRecoQ2Plot = new TH1D("CCMECRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CCMECRecoCCQEECalPlot = new TH1D("CCMECRecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CCMECRecoCCQEQ2Plot = new TH1D("CCMECRecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCMECRecoCosThetaMuPmuPlot = new TH2D("CCMECRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCMECRecoCosThetaPPpPlot = new TH2D("CCMECRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCRES

		TH1D* CCRESRecoNuPlot = new TH1D("CCRESRecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CCRESRecoEvPlot = new TH1D("CCRESRecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CCRESRecoNuScorePlot = new TH1D("CCRESRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCRESRecoFlashScorePlot = new TH1D("CCRESRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* CCRESRecoDistancePlot = new TH1D("CCRESRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* CCRESRecoLengthDifferencePlot = new TH1D("CCRESRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* CCRESRecoVertexXPlot = new TH1D("CCRESRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCRESRecoVertexYPlot = new TH1D("CCRESRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCRESRecoVertexZPlot = new TH1D("CCRESRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* CCRESRecoCandidateMuEndYPlot = new TH1D("CCRESRecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* CCRESRecoMuonLLRPIDPlot = new TH1D("CCRESRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* CCRESRecoProtonLLRPIDPlot = new TH1D("CCRESRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCRESRecoMuonLengthPlot = new TH1D("CCRESRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCRESRecodMuonTracksScorePlot = new TH1D("CCRESRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CCRESRecodMuonVertexDistancePlot = new TH1D("CCRESRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCRESRecoContainedMuonLengthPlot = new TH1D("CCRESRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCRESRecoUncontainedMuonLengthPlot = new TH1D("CCRESRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCRESRecoMuonMomentumPlot = new TH1D("CCRESRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESRecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* CCRESRecoMuonPhiPlot = new TH1D("CCRESRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CCRESRecoCCQEMuonMomentumPlot = new TH1D("CCRESRecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CCRESRecoCCQEMuonCosThetaPlot = new TH1D("CCRESRecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoCCQEMuonPhiPlot = new TH1D("CCRESRecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);	

		TH1D* CCRESRecoContainedMuonMomentumPlot = new TH1D("CCRESRecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCRESRecoUncontainedMuonMomentumPlot = new TH1D("CCRESRecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);

//		TH1D* CCRESRecodProtonTracksScorePlot = new TH1D("CCRESRecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* CCRESRecoProtonLengthPlot = new TH1D("CCRESRecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);
		TH1D* CCRESRecodProtonVertexDistancePlot = new TH1D("CCRESRecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);

		TH1D* CCRESRecoProtonMomentumPlot = new TH1D("CCRESRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCRESRecoProtonCosThetaPlot = new TH1D("CCRESRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCRESRecoProtonPhiPlot = new TH1D("CCRESRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCRESRecoCCQEProtonMomentumPlot = new TH1D("CCRESRecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,
			CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CCRESRecoCCQEProtonCosThetaPlot = new TH1D("CCRESRecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CCRESRecoCCQEProtonPhiPlot = new TH1D("CCRESRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CCRESRecoDeltaThetaPlot = new TH1D("CCRESRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCRESRecoDeltaForwardThetaPlot = new TH1D("CCRESRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCRESRecoDeltaBackwardThetaPlot = new TH1D("CCRESRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCRESRecoDeltaPhiPlot = new TH1D("CCRESRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCRESRecokMissPlot = new TH1D("CCRESRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCRESRecoPMissMinusPlot = new TH1D("CCRESRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCRESRecoPMissPlot = new TH1D("CCRESRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCRESRecoDeltaPLPlot = new TH1D("CCRESRecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* CCRESRecoDeltaPnPlot = new TH1D("CCRESRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CCRESRecoDeltaPtxPlot = new TH1D("CCRESRecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* CCRESRecoDeltaPtyPlot = new TH1D("CCRESRecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* CCRESRecoAPlot = new TH1D("CCRESRecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* CCRESRecoDeltaPTPlot = new TH1D("CCRESRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESRecoDeltaAlphaTPlot = new TH1D("CCRESRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESRecoDeltaPhiTPlot = new TH1D("CCRESRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCRESRecoECalPlot = new TH1D("CCRESRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESRecoEQEPlot = new TH1D("CCRESRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCRESRecoQ2Plot = new TH1D("CCRESRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CCRESRecoCCQEECalPlot = new TH1D("CCRESRecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CCRESRecoCCQEQ2Plot = new TH1D("CCRESRecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCRESRecoCosThetaMuPmuPlot = new TH2D("CCRESRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCRESRecoCosThetaPPpPlot = new TH2D("CCRESRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ---------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCDIS

		TH1D* CCDISRecoNuPlot = new TH1D("CCDISRecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CCDISRecoEvPlot = new TH1D("CCDISRecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CCDISRecoNuScorePlot = new TH1D("CCDISRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCDISRecoFlashScorePlot = new TH1D("CCDISRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* CCDISRecoDistancePlot = new TH1D("CCDISRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* CCDISRecoLengthDifferencePlot = new TH1D("CCDISRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);

		TH1D* CCDISRecoVertexXPlot = new TH1D("CCDISRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCDISRecoVertexYPlot = new TH1D("CCDISRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCDISRecoVertexZPlot = new TH1D("CCDISRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);
		TH1D* CCDISRecoCandidateMuEndYPlot = new TH1D("CCDISRecoCandidateMuEndYPlot",";Muon End Y [cm]",NBinsVertexY,MinVertexY,MaxVertexY);		

		TH1D* CCDISRecoMuonLLRPIDPlot = new TH1D("CCDISRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
		TH1D* CCDISRecoProtonLLRPIDPlot = new TH1D("CCDISRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCDISRecoMuonLengthPlot = new TH1D("CCDISRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCDISRecodMuonTracksScorePlot = new TH1D("CCDISRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CCDISRecodMuonVertexDistancePlot = new TH1D("CCDISRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCDISRecoContainedMuonLengthPlot = new TH1D("CCDISRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCDISRecoUncontainedMuonLengthPlot = new TH1D("CCDISRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCDISRecoMuonMomentumPlot = new TH1D("CCDISRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISRecoMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TH1D* CCDISRecoMuonPhiPlot = new TH1D("CCDISRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CCDISRecoCCQEMuonMomentumPlot = new TH1D("CCDISRecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CCDISRecoCCQEMuonCosThetaPlot = new TH1D("CCDISRecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoCCQEMuonPhiPlot = new TH1D("CCDISRecoCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* CCDISRecoContainedMuonMomentumPlot = new TH1D("CCDISRecoContainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCDISRecoUncontainedMuonMomentumPlot = new TH1D("CCDISRecoUncontainedMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);

		TH1D* CCDISRecoProtonLengthPlot = new TH1D("CCDISRecoProtonLengthPlot",RecoLabelXAxisProtonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);	
//		TH1D* CCDISRecodProtonTracksScorePlot = new TH1D("CCDISRecodProtonTracksScorePlot",RecoLabelXAxisProtonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* CCDISRecodProtonVertexDistancePlot = new TH1D("CCDISRecodProtonVertexDistancePlot",RecoLabelXAxisProtonVertexDistanceTrackScore,NBinsProtonVertexDistance,MinProtonVertexDistance,MaxProtonVertexDistance);

		TH1D* CCDISRecoProtonMomentumPlot = new TH1D("CCDISRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCDISRecoProtonCosThetaPlot = new TH1D("CCDISRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCDISRecoProtonPhiPlot = new TH1D("CCDISRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCDISRecoCCQEProtonMomentumPlot = new TH1D("CCDISRecoCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,
			CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CCDISRecoCCQEProtonCosThetaPlot = new TH1D("CCDISRecoCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CCDISRecoCCQEProtonPhiPlot = new TH1D("CCDISRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CCDISRecoDeltaThetaPlot = new TH1D("CCDISRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaForwardThetaPlot = new TH1D("CCDISRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaBackwardThetaPlot = new TH1D("CCDISRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaPhiPlot = new TH1D("CCDISRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCDISRecokMissPlot = new TH1D("CCDISRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCDISRecoPMissMinusPlot = new TH1D("CCDISRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCDISRecoPMissPlot = new TH1D("CCDISRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCDISRecoDeltaPLPlot = new TH1D("CCDISRecoDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TH1D* CCDISRecoDeltaPnPlot = new TH1D("CCDISRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CCDISRecoDeltaPtxPlot = new TH1D("CCDISRecoDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TH1D* CCDISRecoDeltaPtyPlot = new TH1D("CCDISRecoDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TH1D* CCDISRecoAPlot = new TH1D("CCDISRecoAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);

		TH1D* CCDISRecoDeltaPTPlot = new TH1D("CCDISRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISRecoDeltaAlphaTPlot = new TH1D("CCDISRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISRecoDeltaPhiTPlot = new TH1D("CCDISRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCDISRecoECalPlot = new TH1D("CCDISRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISRecoEQEPlot = new TH1D("CCDISRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCDISRecoQ2Plot = new TH1D("CCDISRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CCDISRecoCCQEECalPlot = new TH1D("CCDISRecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CCDISRecoCCQEQ2Plot = new TH1D("CCDISRecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCDISRecoCosThetaMuPmuPlot = new TH2D("CCDISRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCDISRecoCosThetaPPpPlot = new TH2D("CCDISRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// --------------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------------------

		// Chi2 PID Studies

		// 1D Reco Level Plots

		TH1D* RecoLLRPIDPlot = new TH1D("RecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* RecoChi2Plot = new TH1D("RecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);

		// Muon 1D Reco Level Plots

		TH1D* MuonRecoLLRPIDPlot = new TH1D("MuonRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* MuonRecoChi2Plot = new TH1D("MuonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);

		// Proton 1D Reco Level Plots

		TH1D* ProtonRecoLLRPIDPlot = new TH1D("ProtonRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* ProtonRecoChi2Plot = new TH1D("ProtonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);

		// Pion 1D Reco Level Plots

		TH1D* PionRecoLLRPIDPlot = new TH1D("PionRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* PionRecoChi2Plot = new TH1D("PionRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);

		// Cosmic 1D Reco Level Plots

		TH1D* CosmicRecoLLRPIDPlot = new TH1D("CosmicRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CosmicRecoChi2Plot = new TH1D("CosmicRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);

		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		// Optical & calorimetric info (only Reco & CC1pReco) for selection cuts

		TH2D* RecoChi2PvsChi2Mu2D = new TH2D("RecoChi2PvsChi2Mu2D",";#chi^{2}_{#mu};#chi^{2}_{p}",NBinsChi2,MinChi2,MaxChi2,NBinsChi2,MinChi2,MaxChi2);
		TH2D* CC1pRecoChi2PvsChi2Mu2D = new TH2D("CC1pRecoChi2PvsChi2Mu2D",";#chi^{2}_{#mu};#chi^{2}_{p}",NBinsChi2,MinChi2,MaxChi2,NBinsChi2,MinChi2,MaxChi2);

		TH2D* RecoLengthPvsLengthMu2D = new TH2D("RecoLengthPvsLengthMu2D",RecoLabelXAxisMuonLength+RecoLabelXAxisProtonLength,
							 NBinsMuonLength,MinMuonLength,MaxMuonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);
		TH2D* CC1pRecoLengthPvsLengthMu2D = new TH2D("CC1pRecoLengthPvsLengthMu2D",RecoLabelXAxisMuonLength+RecoLabelXAxisProtonLength,
							 NBinsMuonLength,MinMuonLength,MaxMuonLength,NBinsProtonLength,MinProtonLength,MaxProtonLength);

		// --------------------------------------------------------------------------------------------------------------------------------

		TH1D* POTScalePlot = new TH1D("POTScalePlot","",1,0,1);
		TH1D* NEventsPlot = new TH1D("NEventsPlot","",1,0,1);
		TH1D* NSelectedPlot = new TH1D("NSelectedPlot","",1,0,1);
		TH1D* NCC1pPlot = new TH1D("NCC1pPlot","",1,0,1);
		TH1D* NCC1p1piPlot = new TH1D("NCC1p1piPlot","",1,0,1);
		TH1D* NCC2pPlot = new TH1D("NCC2pPlot","",1,0,1);
		TH1D* NCC2p1piPlot = new TH1D("NCC2p1piPlot","",1,0,1);
		TH1D* NCC3pPlot = new TH1D("NCC3pPlot","",1,0,1);

		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		int NBinsCC1pSTVReso = 200; 
		double MinCC1pSTVMuMom = -1.,MaxCC1pSTVMuMom = 1.;
		double MinCC1pSTVPt = -1.,MaxCC1pSTVPt = 1.;
		double MinCC1pSTVDeltaAlphaT = -180.,MaxCC1pSTVDeltaAlphaT = 180.;
		double MinCC1pSTVDeltaPhiT = -60.,MaxCC1pSTVDeltaPhiT = 60.;

		TString LabelMuonMomentumResolution = ";P_{#mu} (Reco - True) [GeV/c]";
		TString LabelDeltaPTResolution = ";#deltaP_{T} (Reco - True) [GeV/c]";
		TString LabelDeltaAlphaTResolution = ";#delta#alpha_{T} (Reco - True) [deg]";
		TString LabelDeltaPhiTResolution = ";#delta#phi_{T} (Reco - True) [deg]";

		// --------------------------------------------------------------------------------------------------------------------------------

		TH1D* CC1pRecoMuonTrueMomentumLongitudinalRatio = new TH1D("CC1pRecoMuonTrueMomentumLongitudinalRatio",";P^{true}_{#mu,||}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pRecoProtonTrueMomentumLongitudinalRatio = new TH1D("CC1pRecoProtonTrueMomentumLongitudinalRatio",";P^{true}_{p,||}/P^{true}_{p}",25,0.,1.);

		TH1D* CC1pRecoMuonTrueMomentumTransverseRatio = new TH1D("CC1pRecoMuonTrueMomentumTransverseRatio",";P^{true}_{#mu,T}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pRecoProtonTrueMomentumTransverseRatio = new TH1D("CC1pRecoProtonTrueMomentumTransverseRatio",";P^{true}_{p,T}/P^{true}_{p}",25,0.,1.);

		TH1D* CC1pTrueMuonTrueMomentumLongitudinalRatio = new TH1D("CC1pTrueMuonTrueMomentumLongitudinalRatio",";P^{true}_{#mu,||}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pTrueProtonTrueMomentumLongitudinalRatio = new TH1D("CC1pTrueProtonTrueMomentumLongitudinalRatio",";P^{true}_{p,||}/P^{true}_{p}",25,0.,1.);

		TH1D* CC1pTrueMuonTrueMomentumTransverseRatio = new TH1D("CC1pTrueMuonTrueMomentumTransverseRatio",";P^{true}_{#mu,T}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pTrueProtonTrueMomentumTransverseRatio = new TH1D("CC1pTrueProtonTrueMomentumTransverseRatio",";P^{true}_{p,T}/P^{true}_{p}",25,0.,1.);		
		
		//----------------------------------------//

		TH1D* RecoMuonCosThetaPlotLLPOnlyStudy[NBinsLLRPID];
		TH1D* CC1pRecoMuonCosThetaPlotLLPOnlyStudy[NBinsLLRPID];

		for (int WhichBin = 0; WhichBin < NBinsLLRPID; WhichBin++) {

				TString PlotLLPOnlyThres = "RecoMuonCosThetaPlot_LLPThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLLPOnlyStudy[WhichBin] = new TH1D(PlotLLPOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLLPOnlyThres = "CC1pRecoMuonCosThetaPlot_LLPThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLLPOnlyStudy[WhichBin] = new TH1D(CC1pPlotLLPOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		//----------------------------------------//

		// NuScore Only Study

		TH1D* RecoMuonCosThetaPlotNuScoreOnlyStudy[NBinsNuScore];
		TH1D* CC1pRecoMuonCosThetaPlotNuScoreOnlyStudy[NBinsNuScore];

		for (int WhichBin = 0; WhichBin < NBinsNuScore; WhichBin++) {

				TString PlotNuScoreOnlyThres = "RecoMuonCosThetaPlot_NuScoreThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotNuScoreOnlyStudy[WhichBin] = new TH1D(PlotNuScoreOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotNuScoreOnlyThres = "CC1pRecoMuonCosThetaPlot_NuScoreThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotNuScoreOnlyStudy[WhichBin] = new TH1D(CC1pPlotNuScoreOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		//----------------------------------------//

		// Length Only Study

		TH1D* RecoMuonCosThetaPlotLengthOnlyStudy[NBinsMuonLength];
		TH1D* CC1pRecoMuonCosThetaPlotLengthOnlyStudy[NBinsMuonLength];

		for (int WhichBin = 0; WhichBin < NBinsMuonLength; WhichBin++) {

				TString PlotLengthOnlyThres = "RecoMuonCosThetaPlot_LengthThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLengthOnlyStudy[WhichBin] = new TH1D(PlotLengthOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLengthOnlyThres = "CC1pRecoMuonCosThetaPlot_LengthThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLengthOnlyStudy[WhichBin] = new TH1D(CC1pPlotLengthOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		//----------------------------------------//

		Tools tools;		

		//----------------------------------------//

		// 2D analysis
		// Ecal in DeltaPT and DeltaAlphaT bins		

		TH1D* RecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
		TH1D* CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
		TH1D* CC1pTrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];				
		TH1D* NonCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];		
		TH1D* CCQERecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
		TH1D* CCMECRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
		TH1D* CCRESRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
		TH1D* CCDISRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
		TH2D* CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];	
		TH2D* POTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];	

		//----------------------------------------//

		TH1D* RecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pTrueDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];				
		TH1D* NonCC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];		
		TH1D* CCQERecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCMECRecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCRESRecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCDISRecoDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH2D* CC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];	
		TH2D* POTScaledCC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];

		//----------------------------------------//

		TH1D* RecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pRecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pTrueECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];				
		TH1D* NonCC1pRecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];		
		TH1D* CCQERecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCMECRecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCRESRecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCDISRecoECal_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH2D* CC1pRecoECal_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];	
		TH2D* POTScaledCC1pRecoECal_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];		

		//----------------------------------------//

		TH1D* RecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pRecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pTrueDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];				
		TH1D* NonCC1pRecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];		
		TH1D* CCQERecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCMECRecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCRESRecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCDISRecoDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH2D* CC1pRecoDeltaPhiT_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];	
		TH2D* POTScaledCC1pRecoDeltaPhiT_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];		

		//----------------------------------------//

		TH1D* RecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pRecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CC1pTrueDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];				
		TH1D* NonCC1pRecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];		
		TH1D* CCQERecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCMECRecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCRESRecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH1D* CCDISRecoDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
		TH2D* CC1pRecoDeltaPn_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];	
		TH2D* POTScaledCC1pRecoDeltaPn_InDeltaPTTwoDPlot2D[TwoDNBinsDeltaPT];			

		//----------------------------------------//			

		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

			//------------------------------//

			// DeltaAlphaT in DeltaPT slices

			TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			

			RecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("Reco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pReco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pTrueDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pTrue"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);					
			NonCC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("NonCC1pReco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);
			CCQERecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCQEReco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);
			CCMECRecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCMECReco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);
			CCRESRecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCRESReco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);
			CCDISRecoDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCDISReco"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);	
			CC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("CC1pReco"+DeltaAlphaTTwoDInDeltaPTLabel+"2D",LabelXAxisDeltaAlphaT2D,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);	
			POTScaledCC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("POTScaledCC1pReco"+DeltaAlphaTTwoDInDeltaPTLabel+"2D",LabelXAxisDeltaAlphaT2D,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);																									

			//------------------------------//

			// ECal in DeltaPT slices

			TString ECalTwoDInDeltaPTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			

			RecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("Reco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pRecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pReco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pTrueECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pTrue"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);					
			NonCC1pRecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("NonCC1pReco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);
			CCQERecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCQEReco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);
			CCMECRecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCMECReco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);
			CCRESRecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCRESReco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);
			CCDISRecoECal_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCDISReco"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);	
			CC1pRecoECal_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("CC1pReco"+ECalTwoDInDeltaPTLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);	
			POTScaledCC1pRecoECal_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("POTScaledCC1pReco"+ECalTwoDInDeltaPTLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);																									

			//------------------------------//

			// DeltaPhiT in DeltaPT slices

			TString DeltaPhiTTwoDInDeltaPTLabel = "DeltaPhiT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			

			RecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("Reco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pRecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pReco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pTrueDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pTrue"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);					
			NonCC1pRecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("NonCC1pReco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);
			CCQERecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCQEReco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);
			CCMECRecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCMECReco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);
			CCRESRecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCRESReco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);
			CCDISRecoDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCDISReco"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	
			CC1pRecoDeltaPhiT_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("CC1pReco"+DeltaPhiTTwoDInDeltaPTLabel+"2D",LabelXAxisDeltaPhiT2D,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	
			POTScaledCC1pRecoDeltaPhiT_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("POTScaledCC1pReco"+DeltaPhiTTwoDInDeltaPTLabel+"2D",LabelXAxisDeltaPhiT2D,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);																									

			//------------------------------//

			// DeltaPn in DeltaPT slices

			TString DeltaPnTwoDInDeltaPTLabel = "DeltaPn_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			

			RecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("Reco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pRecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pReco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);
			CC1pTrueDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CC1pTrue"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);					
			NonCC1pRecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("NonCC1pReco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);
			CCQERecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCQEReco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);
			CCMECRecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCMECReco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);
			CCRESRecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCRESReco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);
			CCDISRecoDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("CCDISReco"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	
			CC1pRecoDeltaPn_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("CC1pReco"+DeltaPnTwoDInDeltaPTLabel+"2D",LabelXAxisDeltaPn2D,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	
			POTScaledCC1pRecoDeltaPn_InDeltaPTTwoDPlot2D[WhichDeltaPT] = new TH2D("POTScaledCC1pReco"+DeltaPnTwoDInDeltaPTLabel+"2D",LabelXAxisDeltaPn2D,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);																									

			//------------------------------//

			for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

				TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
		
				RecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("Reco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);
				CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("CC1pReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);
				CC1pTrueECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("CC1pTrue"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);					
				NonCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("NonCC1pReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);
				CCQERecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("CCQEReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);
				CCMECRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("CCMECReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);
				CCRESRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("CCRESReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);
				CCDISRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("CCDISReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);	
				CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[WhichDeltaPT][WhichDeltaAlphaT] = new TH2D("CC1pReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0],TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);	
				POTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[WhichDeltaPT][WhichDeltaAlphaT] = new TH2D("POTScaledCC1pReco"+ECalTwoDInDeltaPTDeltaAlphaTLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0],TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);																									

			} // End of the loop over the 2D DeltaAlphaT bins

		} // End of the loop over the 2D DeltaPT bins

		//----------------------------------------//

		// Ecal in MuonCosTheta and MuonMomentum bins		

		TH1D* RecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
		TH1D* CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
		TH1D* CC1pTrueECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];				
		TH1D* NonCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];		
		TH1D* CCQERecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
		TH1D* CCMECRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
		TH1D* CCRESRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
		TH1D* CCDISRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
		TH2D* CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];	
		TH2D* POTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];	

		//----------------------------------------//														

		// DeltaAlphaT in mu cos theta slices

		TH1D* RecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pTrueDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];				
		TH1D* NonCC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];		
		TH1D* CCQERecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCMECRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCRESRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCDISRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH2D* CC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];	
		TH2D* POTScaledCC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];

		//----------------------------------------//

		// DeltaPT in mu cos theta slices

		TH1D* RecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pRecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pTrueDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];				
		TH1D* NonCC1pRecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];		
		TH1D* CCQERecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCMECRecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCRESRecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCDISRecoDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH2D* CC1pRecoDeltaPT_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];	
		TH2D* POTScaledCC1pRecoDeltaPT_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];

		//----------------------------------------//

		// Pmu in mu cos theta slices

		TH1D* RecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pTrueMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];				
		TH1D* NonCC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];		
		TH1D* CCQERecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCMECRecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCRESRecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCDISRecoMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH2D* CC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];	
		TH2D* POTScaledCC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];	

		//----------------------------------------//

		// p cos theta in mu cos theta slices

		TH1D* RecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CC1pTrueProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];				
		TH1D* NonCC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];		
		TH1D* CCQERecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCMECRecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCRESRecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH1D* CCDISRecoProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
		TH2D* CC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];	
		TH2D* POTScaledCC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot2D[TwoDNBinsMuonCosTheta];						

		//----------------------------------------//			

		for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

			TString DeltaAlphaTTwoDInMuonCosThetaLabel = "DeltaAlphaT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			

			RecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("Reco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pReco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pTrueDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pTrue"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);					
			NonCC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("NonCC1pReco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCQERecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCQEReco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCMECRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCMECReco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCRESRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCRESReco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCDISRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCDISReco"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			CC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("CC1pReco"+DeltaAlphaTTwoDInMuonCosThetaLabel+"2D",LabelXAxisDeltaAlphaT2D,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			POTScaledCC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("POTScaledCC1pReco"+DeltaAlphaTTwoDInMuonCosThetaLabel+"2D",LabelXAxisDeltaAlphaT2D,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);									

			TString DeltaPTTwoDInMuonCosThetaLabel = "DeltaPT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			

			RecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("Reco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pRecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pReco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pTrueDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pTrue"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);					
			NonCC1pRecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("NonCC1pReco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCQERecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCQEReco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCMECRecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCMECReco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCRESRecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCRESReco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCDISRecoDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCDISReco"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			CC1pRecoDeltaPT_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("CC1pReco"+DeltaPTTwoDInMuonCosThetaLabel+"2D",LabelXAxisDeltaPT2D,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			POTScaledCC1pRecoDeltaPT_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("POTScaledCC1pReco"+DeltaPTTwoDInMuonCosThetaLabel+"2D",LabelXAxisDeltaPT2D,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);																									

			TString MuonMomentumTwoDInMuonCosThetaLabel = "MuonMomentum_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			

			RecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("Reco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pReco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pTrueMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pTrue"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);					
			NonCC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("NonCC1pReco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCQERecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCQEReco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCMECRecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCMECReco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCRESRecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCRESReco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCDISRecoMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCDISReco"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			CC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("CC1pReco"+MuonMomentumTwoDInMuonCosThetaLabel+"2D",LabelXAxisMuonMomentum2D,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			POTScaledCC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("POTScaledCC1pReco"+MuonMomentumTwoDInMuonCosThetaLabel+"2D",LabelXAxisMuonMomentum2D,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			TString ProtonCosThetaTwoDInMuonCosThetaLabel = "ProtonCosTheta_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			

			RecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("Reco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pReco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CC1pTrueProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CC1pTrue"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);					
			NonCC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("NonCC1pReco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCQERecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCQEReco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCMECRecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCMECReco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCRESRecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCRESReco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);
			CCDISRecoProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("CCDISReco"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			CC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("CC1pReco"+ProtonCosThetaTwoDInMuonCosThetaLabel+"2D",LabelXAxisProtonCosTheta2D,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);	
			POTScaledCC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot2D[WhichMuonCosTheta] = new TH2D("POTScaledCC1pReco"+ProtonCosThetaTwoDInMuonCosThetaLabel+"2D",LabelXAxisProtonCosTheta2D,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0],TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);																												

			//------------------------------//

			for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

				TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";

				RecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("Reco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);
				CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("CC1pReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);
				CC1pTrueECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("CC1pTrue"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);					
				NonCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("NonCC1pReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);
				CCQERecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("CCQEReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);
				CCMECRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("CCMECReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);
				CCRESRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("CCRESReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);
				CCDISRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("CCDISReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);	
				CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[WhichMuonCosTheta][WhichMuonMomentum] = new TH2D("CC1pReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0],TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);	
				POTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[WhichMuonCosTheta][WhichMuonMomentum] = new TH2D("POTScaledCC1pReco"+ECalTwoDInMuonCosThetaMuonMomentumLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0],TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);																									

			} // End of the loop over the 2D MuonMomentum bins

		}	

		//----------------------------------------//

		// Ecal in ProtonCosTheta and ProtonMomentum bins		

		TH1D* RecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];
		TH1D* CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];
		TH1D* CC1pTrueECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];				
		TH1D* NonCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];		
		TH1D* CCQERecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];
		TH1D* CCMECRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];
		TH1D* CCRESRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];
		TH1D* CCDISRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];
		TH2D* CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];	
		TH2D* POTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];		

		//----------------------------------------//														

		// DeltaAlphaT in p cos theta slices

		TH1D* RecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CC1pTrueDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];				
		TH1D* NonCC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];		
		TH1D* CCQERecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCMECRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCRESRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCDISRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH2D* CC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot2D[TwoDNBinsProtonCosTheta];	
		TH2D* POTScaledCC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot2D[TwoDNBinsProtonCosTheta];		

		//----------------------------------------//

		// DeltaPT in p cos theta slices

		TH1D* RecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CC1pRecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CC1pTrueDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];				
		TH1D* NonCC1pRecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];		
		TH1D* CCQERecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCMECRecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCRESRecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCDISRecoDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH2D* CC1pRecoDeltaPT_InProtonCosThetaTwoDPlot2D[TwoDNBinsProtonCosTheta];	
		TH2D* POTScaledCC1pRecoDeltaPT_InProtonCosThetaTwoDPlot2D[TwoDNBinsProtonCosTheta];		

		//----------------------------------------//		

		// Pp in p cos theta slices

		TH1D* RecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CC1pTrueProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];				
		TH1D* NonCC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];		
		TH1D* CCQERecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCMECRecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCRESRecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH1D* CCDISRecoProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];
		TH2D* CC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot2D[TwoDNBinsProtonCosTheta];	
		TH2D* POTScaledCC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot2D[TwoDNBinsProtonCosTheta];				

		//----------------------------------------//			

		for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {	

			//------------------------------//

			TString ProtonMomentumTwoDInProtonCosThetaLabel = "ProtonMomentum_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			

			RecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("Reco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CC1pReco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CC1pTrueProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CC1pTrue"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);					
			NonCC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("NonCC1pReco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCQERecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCQEReco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCMECRecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCMECReco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCRESRecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCRESReco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCDISRecoProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCDISReco"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);	
			CC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot2D[WhichProtonCosTheta] = new TH2D("CC1pReco"+ProtonMomentumTwoDInProtonCosThetaLabel+"2D",LabelXAxisProtonMomentum2D,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0],TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);	
			POTScaledCC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot2D[WhichProtonCosTheta] = new TH2D("POTScaledCC1pReco"+ProtonMomentumTwoDInProtonCosThetaLabel+"2D",LabelXAxisProtonMomentum2D,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0],TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);																									

			//------------------------------//

			TString DeltaAlphaTTwoDInProtonCosThetaLabel = "DeltaAlphaT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			

			RecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("Reco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CC1pReco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CC1pTrueDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CC1pTrue"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);					
			NonCC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("NonCC1pReco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCQERecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCQEReco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCMECRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCMECReco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCRESRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCRESReco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCDISRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCDISReco"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);	
			CC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot2D[WhichProtonCosTheta] = new TH2D("CC1pReco"+DeltaAlphaTTwoDInProtonCosThetaLabel+"2D",LabelXAxisDeltaAlphaT2D,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0],TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);	
			POTScaledCC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot2D[WhichProtonCosTheta] = new TH2D("POTScaledCC1pReco"+DeltaAlphaTTwoDInProtonCosThetaLabel+"2D",LabelXAxisDeltaAlphaT2D,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0],TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			//------------------------------//

			TString DeltaPTTwoDInProtonCosThetaLabel = "DeltaPT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			

			RecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("Reco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CC1pRecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CC1pReco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CC1pTrueDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CC1pTrue"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);					
			NonCC1pRecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("NonCC1pReco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCQERecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCQEReco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCMECRecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCMECReco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCRESRecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCRESReco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);
			CCDISRecoDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("CCDISReco"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);	
			CC1pRecoDeltaPT_InProtonCosThetaTwoDPlot2D[WhichProtonCosTheta] = new TH2D("CC1pReco"+DeltaPTTwoDInProtonCosThetaLabel+"2D",LabelXAxisDeltaPT2D,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0],TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);	
			POTScaledCC1pRecoDeltaPT_InProtonCosThetaTwoDPlot2D[WhichProtonCosTheta] = new TH2D("POTScaledCC1pReco"+DeltaPTTwoDInProtonCosThetaLabel+"2D",LabelXAxisDeltaPT2D,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0],TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);													

			//------------------------------//

			for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

				TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
		
				RecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("Reco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);
				CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("CC1pReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);
				CC1pTrueECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("CC1pTrue"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);					
				NonCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("NonCC1pReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);
				CCQERecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("CCQEReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);
				CCMECRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("CCMECReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);
				CCRESRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("CCRESReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);
				CCDISRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("CCDISReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);	
				CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[WhichProtonCosTheta][WhichProtonMomentum] = new TH2D("CC1pReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0],TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);	
				POTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[WhichProtonCosTheta][WhichProtonMomentum] = new TH2D("POTScaledCC1pReco"+ECalTwoDInProtonCosThetaProtonMomentumLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0],TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);																									

			} // End of the loop over the 2D ProtonMomentum bins

		}	

		//----------------------------------------//

		// DeltaPty in DeltaPtx slices

		TH1D* RecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
		TH1D* CC1pRecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
		TH1D* CC1pTrueDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];				
		TH1D* NonCC1pRecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];		
		TH1D* CCQERecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
		TH1D* CCMECRecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
		TH1D* CCRESRecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
		TH1D* CCDISRecoDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
		TH2D* CC1pRecoDeltaPty_InDeltaPtxTwoDPlot2D[TwoDNBinsDeltaPtx];	
		TH2D* POTScaledCC1pRecoDeltaPty_InDeltaPtxTwoDPlot2D[TwoDNBinsDeltaPtx];				

		//----------------------------------------//			

		for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {	

			TString DeltaPtyTwoDDeltaPtxLabel = "DeltaPty_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			

			RecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("Reco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);
			CC1pRecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("CC1pReco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);
			CC1pTrueDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("CC1pTrue"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);					
			NonCC1pRecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("NonCC1pReco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);
			CCQERecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("CCQEReco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);
			CCMECRecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("CCMECReco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);
			CCRESRecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("CCRESReco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);
			CCDISRecoDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("CCDISReco"+DeltaPtyTwoDDeltaPtxLabel,LabelXAxisDeltaPty,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);	
			CC1pRecoDeltaPty_InDeltaPtxTwoDPlot2D[WhichDeltaPtx] = new TH2D("CC1pReco"+DeltaPtyTwoDDeltaPtxLabel+"2D",LabelXAxisDeltaPty2D,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0],TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);	
			POTScaledCC1pRecoDeltaPty_InDeltaPtxTwoDPlot2D[WhichDeltaPtx] = new TH2D("POTScaledCC1pReco"+DeltaPtyTwoDDeltaPtxLabel+"2D",LabelXAxisDeltaPty2D,TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0],TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx].size()-1,&TwoDArrayNBinsDeltaPtyInDeltaPtxSlices[WhichDeltaPtx][0]);																									

		}	

		//----------------------------------------//

		// DeltaPtx in DeltaPty slices

		TH1D* RecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];
		TH1D* CC1pRecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];
		TH1D* CC1pTrueDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];				
		TH1D* NonCC1pRecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];		
		TH1D* CCQERecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];
		TH1D* CCMECRecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];
		TH1D* CCRESRecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];
		TH1D* CCDISRecoDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];
		TH2D* CC1pRecoDeltaPtx_InDeltaPtyTwoDPlot2D[TwoDNBinsDeltaPty];	
		TH2D* POTScaledCC1pRecoDeltaPtx_InDeltaPtyTwoDPlot2D[TwoDNBinsDeltaPty];				

		//----------------------------------------//			

		for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {	

			TString DeltaPtxTwoDDeltaPtyLabel = "DeltaPtx_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";			

			RecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("Reco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);
			CC1pRecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("CC1pReco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);
			CC1pTrueDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("CC1pTrue"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);					
			NonCC1pRecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("NonCC1pReco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);
			CCQERecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("CCQEReco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);
			CCMECRecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("CCMECReco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);
			CCRESRecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("CCRESReco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);
			CCDISRecoDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("CCDISReco"+DeltaPtxTwoDDeltaPtyLabel,LabelXAxisDeltaPtx,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);	
			CC1pRecoDeltaPtx_InDeltaPtyTwoDPlot2D[WhichDeltaPty] = new TH2D("CC1pReco"+DeltaPtxTwoDDeltaPtyLabel+"2D",LabelXAxisDeltaPtx2D,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0],TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);	
			POTScaledCC1pRecoDeltaPtx_InDeltaPtyTwoDPlot2D[WhichDeltaPty] = new TH2D("POTScaledCC1pReco"+DeltaPtxTwoDDeltaPtyLabel+"2D",LabelXAxisDeltaPtx2D,TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0],TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty].size()-1,&TwoDArrayNBinsDeltaPtxInDeltaPtySlices[WhichDeltaPty][0]);																									

		}		

		//----------------------------------------//

		// DeltaPT in DeltaAlphaT slices

		TH1D* RecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CC1pTrueDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];				
		TH1D* NonCC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];		
		TH1D* CCQERecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCMECRecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCRESRecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCDISRecoDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH2D* CC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaAlphaT];	
		TH2D* POTScaledCC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaAlphaT];

		// DeltaPn in DeltaAlphaT slices

		TH1D* RecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CC1pTrueDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];				
		TH1D* NonCC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];		
		TH1D* CCQERecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCMECRecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCRESRecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCDISRecoDeltaPn_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH2D* CC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaAlphaT];	
		TH2D* POTScaledCC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaAlphaT];				

		// ECal in DeltaAlphaT slices

		TH1D* RecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CC1pRecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CC1pTrueECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];				
		TH1D* NonCC1pRecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];		
		TH1D* CCQERecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCMECRecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCRESRecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH1D* CCDISRecoECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
		TH2D* CC1pRecoECal_InDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaAlphaT];	
		TH2D* POTScaledCC1pRecoECal_InDeltaAlphaTTwoDPlot2D[TwoDNBinsDeltaAlphaT];			

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

			TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			RecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("Reco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CC1pReco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CC1pTrueDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CC1pTrue"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);					
			NonCC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("NonCC1pReco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCQERecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCQEReco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCMECRecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCMECReco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCRESRecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCRESReco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCDISRecoDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCDISReco"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);	
			CC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot2D[WhichDeltaAlphaT] = new TH2D("CC1pReco"+DeltaPTTwoDInDeltaAlphaTLabel+"2D",LabelXAxisDeltaPT2D,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);	
			POTScaledCC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot2D[WhichDeltaAlphaT] = new TH2D("POTScaledCC1pReco"+DeltaPTTwoDInDeltaAlphaTLabel+"2D",LabelXAxisDeltaPT2D,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);			

			TString DeltaPnTwoDInDeltaAlphaTLabel = "DeltaPn_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			RecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("Reco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CC1pReco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CC1pTrueDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CC1pTrue"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);					
			NonCC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("NonCC1pReco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCQERecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCQEReco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCMECRecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCMECReco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCRESRecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCRESReco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCDISRecoDeltaPn_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCDISReco"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);	
			CC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot2D[WhichDeltaAlphaT] = new TH2D("CC1pReco"+DeltaPnTwoDInDeltaAlphaTLabel+"2D",LabelXAxisDeltaPn2D,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);	
			POTScaledCC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot2D[WhichDeltaAlphaT] = new TH2D("POTScaledCC1pReco"+DeltaPnTwoDInDeltaAlphaTLabel+"2D",LabelXAxisDeltaPn2D,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);			


			TString ECalTwoDInDeltaAlphaTLabel = "ECal_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			RecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("Reco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CC1pRecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CC1pReco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CC1pTrueECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CC1pTrue"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);					
			NonCC1pRecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("NonCC1pReco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCQERecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCQEReco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCMECRecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCMECReco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCRESRecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCRESReco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);
			CCDISRecoECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("CCDISReco"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);	
			CC1pRecoECal_InDeltaAlphaTTwoDPlot2D[WhichDeltaAlphaT] = new TH2D("CC1pReco"+ECalTwoDInDeltaAlphaTLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);	
			POTScaledCC1pRecoECal_InDeltaAlphaTTwoDPlot2D[WhichDeltaAlphaT] = new TH2D("POTScaledCC1pReco"+ECalTwoDInDeltaAlphaTLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

		}						

		//----------------------------------------//

		// DeltaPT in cos theta mu bins
		TH1D* SerialRecoDeltaPT_InMuonCosThetaPlot = new TH1D("RecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pRecoDeltaPT_InMuonCosThetaPlot = new TH1D("CC1pRecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPT_InMuonCosThetaPlot = new TH1D("CC1pTrueSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPT_InMuonCosThetaPlot = new TH1D("NonCC1pRecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoDeltaPT_InMuonCosThetaPlot = new TH1D("CCQERecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPT_InMuonCosThetaPlot = new TH1D("CCMECRecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPT_InMuonCosThetaPlot = new TH1D("CCRESRecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPT_InMuonCosThetaPlot = new TH1D("CCDISRecoSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPT_InMuonCosThetaPlot2D = new TH2D("CC1pRecoSerialDeltaPT_MuonCosThetaPlot2D",LabelXAxisDeltaPT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPT_InMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPT_MuonCosThetaPlot2D",LabelXAxisDeltaPT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);				

		// DeltaPT in cos theta p bins
		TH1D* SerialRecoDeltaPT_InProtonCosThetaPlot = new TH1D("RecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);	
		TH1D* SerialCC1pRecoDeltaPT_InProtonCosThetaPlot = new TH1D("CC1pRecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPT_InProtonCosThetaPlot = new TH1D("CC1pTrueSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPT_InProtonCosThetaPlot = new TH1D("NonCC1pRecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoDeltaPT_InProtonCosThetaPlot = new TH1D("CCQERecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPT_InProtonCosThetaPlot = new TH1D("CCMECRecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPT_InProtonCosThetaPlot = new TH1D("CCRESRecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPT_InProtonCosThetaPlot = new TH1D("CCDISRecoSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPT_InProtonCosThetaPlot2D = new TH2D("CC1pRecoSerialDeltaPT_ProtonCosThetaPlot2D",LabelXAxisDeltaPT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPT_InProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPT_ProtonCosThetaPlot2D",LabelXAxisDeltaPT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);

		// Pmu in cos theta mu bins
		TH1D* SerialRecoMuonMomentum_InMuonCosThetaPlot = new TH1D("RecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pRecoMuonMomentum_InMuonCosThetaPlot = new TH1D("CC1pRecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueMuonMomentum_InMuonCosThetaPlot = new TH1D("CC1pTrueSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoMuonMomentum_InMuonCosThetaPlot = new TH1D("NonCC1pRecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoMuonMomentum_InMuonCosThetaPlot = new TH1D("CCQERecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoMuonMomentum_InMuonCosThetaPlot = new TH1D("CCMECRecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoMuonMomentum_InMuonCosThetaPlot = new TH1D("CCRESRecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoMuonMomentum_InMuonCosThetaPlot = new TH1D("CCDISRecoSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoMuonMomentum_InMuonCosThetaPlot2D = new TH2D("CC1pRecoSerialMuonMomentum_MuonCosThetaPlot2D",LabelXAxisMuonMomentum2D,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoMuonMomentum_InMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialMuonMomentum_MuonCosThetaPlot2D",LabelXAxisMuonMomentum2D,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);

		// Pp in cos theta p bins
		TH1D* SerialRecoProtonMomentum_InProtonCosThetaPlot = new TH1D("RecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);	
		TH1D* SerialCC1pRecoProtonMomentum_InProtonCosThetaPlot = new TH1D("CC1pRecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueProtonMomentum_InProtonCosThetaPlot = new TH1D("CC1pTrueSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoProtonMomentum_InProtonCosThetaPlot = new TH1D("NonCC1pRecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoProtonMomentum_InProtonCosThetaPlot = new TH1D("CCQERecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoProtonMomentum_InProtonCosThetaPlot = new TH1D("CCMECRecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoProtonMomentum_InProtonCosThetaPlot = new TH1D("CCRESRecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoProtonMomentum_InProtonCosThetaPlot = new TH1D("CCDISRecoSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoProtonMomentum_InProtonCosThetaPlot2D = new TH2D("CC1pRecoSerialProtonMomentum_ProtonCosThetaPlot2D",LabelXAxisProtonMomentum2D,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoProtonMomentum_InProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialProtonMomentum_ProtonCosThetaPlot2D",LabelXAxisProtonMomentum2D,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);

		// DeltaAlphaT in cos theta mu bins
		TH1D* SerialRecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("RecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pRecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("CC1pRecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueDeltaAlphaT_InMuonCosThetaPlot = new TH1D("CC1pTrueSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("NonCC1pRecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("CCQERecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("CCMECRecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("CCRESRecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoDeltaAlphaT_InMuonCosThetaPlot = new TH1D("CCDISRecoSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoDeltaAlphaT_InMuonCosThetaPlot2D = new TH2D("CC1pRecoSerialDeltaAlphaT_MuonCosThetaPlot2D",LabelXAxisDeltaAlphaT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaAlphaT_InMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaAlphaT_MuonCosThetaPlot2D",LabelXAxisDeltaAlphaT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);

		// DeltaAlphaT in cos theta p bins
		TH1D* SerialRecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("RecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH1D* SerialCC1pRecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("CC1pRecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueDeltaAlphaT_InProtonCosThetaPlot = new TH1D("CC1pTrueSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("NonCC1pRecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("CCQERecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("CCMECRecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("CCRESRecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoDeltaAlphaT_InProtonCosThetaPlot = new TH1D("CCDISRecoSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoDeltaAlphaT_InProtonCosThetaPlot2D = new TH2D("CC1pRecoSerialDeltaAlphaT_ProtonCosThetaPlot2D",LabelXAxisDeltaAlphaT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaAlphaT_InProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaAlphaT_ProtonCosThetaPlot2D",LabelXAxisDeltaAlphaT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);

		// DeltaAlphaT in DeltaPT bins
		TH1D* SerialRecoDeltaAlphaT_InDeltaPTPlot = new TH1D("RecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);		
		TH1D* SerialCC1pRecoDeltaAlphaT_InDeltaPTPlot = new TH1D("CC1pRecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		TH1D* SerialCC1pTrueDeltaAlphaT_InDeltaPTPlot = new TH1D("CC1pTrueSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaAlphaT_InDeltaPTPlot = new TH1D("NonCC1pRecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		TH1D* SerialCCQERecoDeltaAlphaT_InDeltaPTPlot = new TH1D("CCQERecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		TH1D* SerialCCMECRecoDeltaAlphaT_InDeltaPTPlot = new TH1D("CCMECRecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaAlphaT_InDeltaPTPlot = new TH1D("CCRESRecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		TH1D* SerialCCDISRecoDeltaAlphaT_InDeltaPTPlot = new TH1D("CCDISRecoSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		TH2D* SerialCC1pRecoDeltaAlphaT_InDeltaPTPlot2D = new TH2D("CC1pRecoSerialDeltaAlphaT_DeltaPTPlot2D",LabelXAxisDeltaAlphaT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaAlphaT_InDeltaPTPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaAlphaT_DeltaPTPlot2D",LabelXAxisDeltaAlphaT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);

		// DeltaPT in DeltaAlphaT bins
		TH1D* SerialRecoDeltaPT_InDeltaAlphaTPlot = new TH1D("RecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);		
		TH1D* SerialCC1pRecoDeltaPT_InDeltaAlphaTPlot = new TH1D("CC1pRecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPT_InDeltaAlphaTPlot = new TH1D("CC1pTrueSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPT_InDeltaAlphaTPlot = new TH1D("NonCC1pRecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCQERecoDeltaPT_InDeltaAlphaTPlot = new TH1D("CCQERecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPT_InDeltaAlphaTPlot = new TH1D("CCMECRecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPT_InDeltaAlphaTPlot = new TH1D("CCRESRecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPT_InDeltaAlphaTPlot = new TH1D("CCDISRecoSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPT_InDeltaAlphaTPlot2D = new TH2D("CC1pRecoSerialDeltaPT_DeltaAlphaTPlot2D",LabelXAxisDeltaPT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPT_InDeltaAlphaTPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPT_DeltaAlphaTPlot2D",LabelXAxisDeltaPT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);

		// DeltaPn in DeltaAlphaT bins
		TH1D* SerialRecoDeltaPn_InDeltaAlphaTPlot = new TH1D("RecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);		
		TH1D* SerialCC1pRecoDeltaPn_InDeltaAlphaTPlot = new TH1D("CC1pRecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPn_InDeltaAlphaTPlot = new TH1D("CC1pTrueSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPn_InDeltaAlphaTPlot = new TH1D("NonCC1pRecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCQERecoDeltaPn_InDeltaAlphaTPlot = new TH1D("CCQERecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPn_InDeltaAlphaTPlot = new TH1D("CCMECRecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPn_InDeltaAlphaTPlot = new TH1D("CCRESRecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPn_InDeltaAlphaTPlot = new TH1D("CCDISRecoSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPn_InDeltaAlphaTPlot2D = new TH2D("CC1pRecoSerialDeltaPn_DeltaAlphaTPlot2D",LabelXAxisDeltaPn2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPn_InDeltaAlphaTPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPn_DeltaAlphaTPlot2D",LabelXAxisDeltaPn2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);				

		// DeltaPhiT in DeltaPT bins
		TH1D* SerialRecoDeltaPhiT_InDeltaPTPlot = new TH1D("RecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH1D* SerialCC1pRecoDeltaPhiT_InDeltaPTPlot = new TH1D("CC1pRecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPhiT_InDeltaPTPlot = new TH1D("CC1pTrueSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPhiT_InDeltaPTPlot = new TH1D("NonCC1pRecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH1D* SerialCCQERecoDeltaPhiT_InDeltaPTPlot = new TH1D("CCQERecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPhiT_InDeltaPTPlot = new TH1D("CCMECRecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPhiT_InDeltaPTPlot = new TH1D("CCRESRecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPhiT_InDeltaPTPlot = new TH1D("CCDISRecoSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPhiT_InDeltaPTPlot2D = new TH2D("CC1pRecoSerialDeltaPhiT_DeltaPTPlot2D",LabelXAxisDeltaPhiT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPhiT_InDeltaPTPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPhiT_DeltaPTPlot2D",LabelXAxisDeltaPhiT2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);

		// DeltaPn in DeltaPT bins
		TH1D* SerialRecoDeltaPn_InDeltaPTPlot = new TH1D("RecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH1D* SerialCC1pRecoDeltaPn_InDeltaPTPlot = new TH1D("CC1pRecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPn_InDeltaPTPlot = new TH1D("CC1pTrueSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPn_InDeltaPTPlot = new TH1D("NonCC1pRecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH1D* SerialCCQERecoDeltaPn_InDeltaPTPlot = new TH1D("CCQERecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPn_InDeltaPTPlot = new TH1D("CCMECRecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPn_InDeltaPTPlot = new TH1D("CCRESRecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPn_InDeltaPTPlot = new TH1D("CCDISRecoSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPn_InDeltaPTPlot2D = new TH2D("CC1pRecoSerialDeltaPn_DeltaPTPlot2D",LabelXAxisDeltaPn2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPn_InDeltaPTPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPn_DeltaPTPlot2D",LabelXAxisDeltaPn2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);

		// ECal in DeltaPT bins
		TH1D* SerialRecoECal_InDeltaPTPlot = new TH1D("RecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH1D* SerialCC1pRecoECal_InDeltaPTPlot = new TH1D("CC1pRecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH1D* SerialCC1pTrueECal_InDeltaPTPlot = new TH1D("CC1pTrueSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);		
		TH1D* SerialNonCC1pRecoECal_InDeltaPTPlot = new TH1D("NonCC1pRecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH1D* SerialCCQERecoECal_InDeltaPTPlot = new TH1D("CCQERecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH1D* SerialCCMECRecoECal_InDeltaPTPlot = new TH1D("CCMECRecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);								
		TH1D* SerialCCRESRecoECal_InDeltaPTPlot = new TH1D("CCRESRecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH1D* SerialCCDISRecoECal_InDeltaPTPlot = new TH1D("CCDISRecoSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH2D* SerialCC1pRecoECal_InDeltaPTPlot2D = new TH2D("CC1pRecoSerialECal_DeltaPTPlot2D",LabelXAxisECal2D,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoECal_InDeltaPTPlot2D = new TH2D("POTScaledCC1pRecoSerialECal_DeltaPTPlot2D",LabelXAxisECal2D,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);

		// ECal in DeltaAlphaT bins
		TH1D* SerialRecoECal_InDeltaAlphaTPlot = new TH1D("RecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH1D* SerialCC1pRecoECal_InDeltaAlphaTPlot = new TH1D("CC1pRecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH1D* SerialCC1pTrueECal_InDeltaAlphaTPlot = new TH1D("CC1pTrueSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);		
		TH1D* SerialNonCC1pRecoECal_InDeltaAlphaTPlot = new TH1D("NonCC1pRecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCQERecoECal_InDeltaAlphaTPlot = new TH1D("CCQERecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCMECRecoECal_InDeltaAlphaTPlot = new TH1D("CCMECRecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);								
		TH1D* SerialCCRESRecoECal_InDeltaAlphaTPlot = new TH1D("CCRESRecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH1D* SerialCCDISRecoECal_InDeltaAlphaTPlot = new TH1D("CCDISRecoSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH2D* SerialCC1pRecoECal_InDeltaAlphaTPlot2D = new TH2D("CC1pRecoSerialECal_DeltaAlphaTPlot2D",LabelXAxisECal2D,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoECal_InDeltaAlphaTPlot2D = new TH2D("POTScaledCC1pRecoSerialECal_DeltaAlphaTPlot2D",LabelXAxisECal2D,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0],tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);		

		// cos theta p in cos theta mu bins
		TH1D* SerialRecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("RecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pRecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("CC1pRecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH1D* SerialCC1pTrueProtonCosTheta_InMuonCosThetaPlot = new TH1D("CC1pTrueSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);		
		TH1D* SerialNonCC1pRecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("NonCC1pRecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH1D* SerialCCQERecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("CCQERecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH1D* SerialCCMECRecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("CCMECRecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);								
		TH1D* SerialCCRESRecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("CCRESRecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH1D* SerialCCDISRecoProtonCosTheta_InMuonCosThetaPlot = new TH1D("CCDISRecoSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH2D* SerialCC1pRecoProtonCosTheta_InMuonCosThetaPlot2D = new TH2D("CC1pRecoSerialProtonCosTheta_MuonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoProtonCosTheta_InMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoSerialProtonCosTheta_MuonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0],tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);

		// DeltaPty in DeltaPtx bins
		TH1D* SerialRecoDeltaPty_InDeltaPtxPlot = new TH1D("RecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH1D* SerialCC1pRecoDeltaPty_InDeltaPtxPlot = new TH1D("CC1pRecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH1D* SerialCC1pTrueDeltaPty_InDeltaPtxPlot = new TH1D("CC1pTrueSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPty_InDeltaPtxPlot = new TH1D("NonCC1pRecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH1D* SerialCCQERecoDeltaPty_InDeltaPtxPlot = new TH1D("CCQERecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH1D* SerialCCMECRecoDeltaPty_InDeltaPtxPlot = new TH1D("CCMECRecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPty_InDeltaPtxPlot = new TH1D("CCRESRecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH1D* SerialCCDISRecoDeltaPty_InDeltaPtxPlot = new TH1D("CCDISRecoSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH2D* SerialCC1pRecoDeltaPty_InDeltaPtxPlot2D = new TH2D("CC1pRecoSerialDeltaPty_DeltaPtxPlot2D",LabelXAxisDeltaPty2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPty_InDeltaPtxPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPty_DeltaPtxPlot2D",LabelXAxisDeltaPty2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);

		// DeltaPtx in DeltaPty bins
		TH1D* SerialRecoDeltaPtx_InDeltaPtyPlot = new TH1D("RecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH1D* SerialCC1pRecoDeltaPtx_InDeltaPtyPlot = new TH1D("CC1pRecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH1D* SerialCC1pTrueDeltaPtx_InDeltaPtyPlot = new TH1D("CC1pTrueSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);		
		TH1D* SerialNonCC1pRecoDeltaPtx_InDeltaPtyPlot = new TH1D("NonCC1pRecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH1D* SerialCCQERecoDeltaPtx_InDeltaPtyPlot = new TH1D("CCQERecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH1D* SerialCCMECRecoDeltaPtx_InDeltaPtyPlot = new TH1D("CCMECRecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);								
		TH1D* SerialCCRESRecoDeltaPtx_InDeltaPtyPlot = new TH1D("CCRESRecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH1D* SerialCCDISRecoDeltaPtx_InDeltaPtyPlot = new TH1D("CCDISRecoSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH2D* SerialCC1pRecoDeltaPtx_InDeltaPtyPlot2D = new TH2D("CC1pRecoSerialDeltaPtx_DeltaPtyPlot2D",LabelXAxisDeltaPtx2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoDeltaPtx_InDeltaPtyPlot2D = new TH2D("POTScaledCC1pRecoSerialDeltaPtx_DeltaPtyPlot2D",LabelXAxisDeltaPtx2D,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0],tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);

		// ECal in DeltaPT & DeltaAlphaT bins
		TH1D* SerialRecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("RecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH1D* SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("CC1pRecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH1D* SerialCC1pTrueECal_InDeltaPTDeltaAlphaTPlot = new TH1D("CC1pTrueSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);		
		TH1D* SerialNonCC1pRecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("NonCC1pRecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH1D* SerialCCQERecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("CCQERecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH1D* SerialCCMECRecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("CCMECRecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);								
		TH1D* SerialCCRESRecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("CCRESRecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH1D* SerialCCDISRecoECal_InDeltaPTDeltaAlphaTPlot = new TH1D("CCDISRecoSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH2D* SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot2D = new TH2D("CC1pRecoSerialECal_DeltaPTDeltaAlphaTPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTPlot2D = new TH2D("POTScaledCC1pRecoSerialECal_DeltaPTDeltaAlphaTPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);

		// ECal in MuonCosTheta & MuonMomentum bins
		TH1D* SerialRecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("RecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH1D* SerialCC1pRecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("CC1pRecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH1D* SerialCC1pTrueECal_InMuonCosThetaMuonMomentumPlot = new TH1D("CC1pTrueSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);		
		TH1D* SerialNonCC1pRecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("NonCC1pRecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH1D* SerialCCQERecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("CCQERecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH1D* SerialCCMECRecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("CCMECRecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);								
		TH1D* SerialCCRESRecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("CCRESRecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH1D* SerialCCDISRecoECal_InMuonCosThetaMuonMomentumPlot = new TH1D("CCDISRecoSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH2D* SerialCC1pRecoECal_InMuonCosThetaMuonMomentumPlot2D = new TH2D("CC1pRecoSerialECal_MuonCosThetaMuonMomentumPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumPlot2D = new TH2D("POTScaledCC1pRecoSerialECal_MuonCosThetaMuonMomentumPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);

		// ECal in ProtonCosTheta & ProtonMomentum bins
		TH1D* SerialRecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("RecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH1D* SerialCC1pRecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("CC1pRecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH1D* SerialCC1pTrueECal_InProtonCosThetaProtonMomentumPlot = new TH1D("CC1pTrueSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);		
		TH1D* SerialNonCC1pRecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("NonCC1pRecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH1D* SerialCCQERecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("CCQERecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH1D* SerialCCMECRecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("CCMECRecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);								
		TH1D* SerialCCRESRecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("CCRESRecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH1D* SerialCCDISRecoECal_InProtonCosThetaProtonMomentumPlot = new TH1D("CCDISRecoSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH2D* SerialCC1pRecoECal_InProtonCosThetaProtonMomentumPlot2D = new TH2D("CC1pRecoSerialECal_ProtonCosThetaProtonMomentumPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumPlot2D = new TH2D("POTScaledCC1pRecoSerialECal_ProtonCosThetaProtonMomentumPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);

		//----------------------------------------//
		//----------------------------------------//

		// Loop over the events

		cout << nentries << " events included in the file" << endl;

		for (Long64_t jentry=0; jentry<nentries;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

			TotalCounter++;

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

			ContainmentCounter++;

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
				if (fTune == "NoTune") { weight = POTWeight * Weight * ROOTinoWeight; }
				// Double the MEC weight (mode = 10)
				if (fTune == "TwiceMEC" && MCParticle_Mode == 10) { weight = 2 * POTWeight * Weight * T2KWeight * ROOTinoWeight; }
				if (fTune == "TwiceMEC" && MCParticle_Mode != 10) { weight = POTWeight * Weight * T2KWeight * ROOTinoWeight; }									
				
			}
			
			// --------------------------------------------------------------------------------------------------------------------------------

			// Genie, flux & reinteraction weights for systematics

			if ( 
				fUniverseIndex != -1 && 
				(fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || 
				fWhichSample == "Overlay9_Run3" || fWhichSample == "Overlay9_Run4a"
			|| fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5" || fWhichSample == "Overlay9_Combined" 
			|| fWhichSample == "OverlayDirt9_Run1" || fWhichSample == "OverlayDirt9_Run2" || 
				fWhichSample == "OverlayDirt9_Run3" || fWhichSample == "OverlayDirt9_Run4a"
			|| fWhichSample == "OverlayDirt9_Run4" || fWhichSample == "OverlayDirt9_Run5" || fWhichSample == "OverlayDirt9_Combined") 
			) {

				// Watch out: The EventWeight weights already include the weight for the tune

				// Genie weights
				if (fEventWeightLabel == "All_UBGenie") { weight = weight*All_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "AxFFCCQEshape_UBGenie") { weight = weight*AxFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "DecayAngMEC_UBGenie") { weight = weight*DecayAngMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "NormCCCOH_UBGenie") { weight = weight*NormCCCOH_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "NormNCCOH_UBGenie") { weight = weight*NormNCCOH_UBGenie->at(fUniverseIndex) / T2KWeight; }
//				if (fEventWeightLabel == "RPA_CCQE_Reduced_UBGenie") { weight = weight*RPA_CCQE_Reduced_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "RPA_CCQE_UBGenie") { weight = weight*RPA_CCQE_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "ThetaDelta2NRad_UBGenie") { weight = weight*ThetaDelta2NRad_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "Theta_Delta2Npi_UBGenie") { weight = weight*Theta_Delta2Npi_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "VecFFCCQEshape_UBGenie") { weight = weight*VecFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
				if (fEventWeightLabel == "XSecShape_CCMEC_UBGenie") { weight = weight*XSecShape_CCMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }

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

			// -----------------------------------------------------------------------------------------------------------------------

			if ( fabs(weight) != weight) { continue; } // Securing against infinities

			// -----------------------------------------------------------------------------------------------------------------------------

			// Contained Reconstructed Vertex

			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));

			if ( !tools.inFVVector(RecoVertex) ) { continue; }
			ContainedVertexCounter++;

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

			MuonQualityCutCounter++;

			// -------------------------------------------------------------------------------------------------------------------------

			// Ensure that we don't have flipped tracks 
			// by demanding that muon start - vertex distance < muon end - vertex distance,
			// that proton start - vertex distance < proton end - vertex distance

			if (CandidateMuStartVertexDistance->at(0) > CandidateMuEndVertexDistance->at(0)) { continue; }
			if (CandidatePStartVertexDistance->at(0) > CandidatePEndVertexDistance->at(0)) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------

			// Targeting the remaining muon cosmic contamination which piles up at high Y

			//if (CandidateMu_EndY->at(0) > 85) { continue; } 

			// -------------------------------------------------------------------------------------------------------------------------			

			TVector3 CandidateMuonStart(CandidateMu_StartX->at(0),CandidateMu_StartY->at(0),CandidateMu_StartZ->at(0));
			TVector3 CandidateMuonEnd(CandidateMu_EndX->at(0),CandidateMu_EndY->at(0),CandidateMu_EndZ->at(0));

			TVector3 CandidateProtonStart(CandidateP_StartX->at(0),CandidateP_StartY->at(0),CandidateP_StartZ->at(0));
			TVector3 CandidateProtonEnd(CandidateP_EndX->at(0),CandidateP_EndY->at(0),CandidateP_EndZ->at(0));

			double DistanceStartPoints = StartToStartDistance->at(0);
			double DistanceEndPoints = EndToEndDistance->at(0); 

			if (DistanceStartPoints > DistanceEndPoints) { continue; }

			NoFlippedTrackCounter++;

			// -------------------------------------------------------------------------------------------------------------------------

			double l_muCandidate = CandidateMu_Length->at(0);
			double l_pCandidate = CandidateP_Length->at(0);
			double LengthDifference = l_muCandidate - l_pCandidate;
//			double MuonTrackScore = CandidateMu_TrackScore->at(0);
//			double ProtonTrackScore = CandidateP_TrackScore->at(0);
			double distance = CandidateMuP_Distance->at(0);	

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
//			double reco_Tp = reco_Ep - ProtonMass_GeV;

			TVector3 TVector3CandidateProton(-1,-1,-1);
			TVector3CandidateProton.SetMag(reco_Pp);
			TVector3CandidateProton.SetTheta(TMath::ACos(reco_Pp_cos_theta));
			TVector3CandidateProton.SetPhi(reco_Pp_phi * TMath::Pi() / 180.);	

			// --------------------------------------------------------------------------------------------------------------------------

			// Calorimetry

			double reco_mu_LLR_Score = CandidateMu_LLR_PID->at(0);
			double reco_p_LLR_Score = CandidateP_LLR_PID->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			// STV variables
			
			double TransMissMomentum = Reco_Pt->at(0);
			double DeltaAlphaT = Reco_DeltaAlphaT->at(0);
			double DeltaPhiT = Reco_DeltaPhiT->at(0);

			// Light cone variables
			double kMiss = Reco_kMiss->at(0);
			double PMissMinus = Reco_PMissMinus->at(0);
			double MissMomentum = Reco_PMiss->at(0);
			double reco_A = Reco_A->at(0);

			// Minerva variables
			double reco_PL = Reco_PL->at(0);
			double reco_Pn = Reco_Pn->at(0);
			double reco_Ptx = Reco_Ptx->at(0);
			double reco_Pty = Reco_Pty->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction
			double ECal = Reco_ECal->at(0);

			// QE Energy Reconstruction
			double EQE = Reco_EQE->at(0);

			// Reconstructed Q2
			double reco_Q2 = Reco_Q2->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// STV redefinition if P_p < 0.5 where the biases have been observed
			if (reco_Pp < 0.5) { 

				reco_Pp = ( 1.-0.01*fPP->Eval(reco_Pp) ) * reco_Pp ;
				TVector3CandidateProton.SetMag(reco_Pp);
				reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );

				STV_Tools reco_stv_tool(TVector3CandidateMuon,TVector3CandidateProton,reco_Emu,reco_Ep);
	
				TransMissMomentum = reco_stv_tool.ReturnPt();
				DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
				DeltaPhiT = reco_stv_tool.ReturnDeltaPhiT();
				ECal = reco_stv_tool.ReturnECal();

				kMiss = reco_stv_tool.ReturnkMiss();
				PMissMinus = reco_stv_tool.ReturnPMissMinus();
				MissMomentum = reco_stv_tool.ReturnPMiss();

				reco_PL = reco_stv_tool.ReturnPL();
				reco_Pn = reco_stv_tool.ReturnPn();
				reco_Ptx = reco_stv_tool.ReturnPtx();
				reco_Pty = reco_stv_tool.ReturnPty();
				reco_A = reco_stv_tool.ReturnA();

			}

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

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

			if (TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { TransMissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
			if (reco_Ptx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { reco_Ptx = 0.5 * (ArrayNBinsDeltaPtx[NBinsDeltaPtx] + ArrayNBinsDeltaPtx[NBinsDeltaPtx-1]); }
			if (reco_Pty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { reco_Pty = 0.5 * (ArrayNBinsDeltaPty[NBinsDeltaPty] + ArrayNBinsDeltaPty[NBinsDeltaPty-1]); }
			if (reco_PL > ArrayNBinsDeltaPL[NBinsDeltaPL]) { reco_PL = 0.5 * (ArrayNBinsDeltaPL[NBinsDeltaPL] + ArrayNBinsDeltaPL[NBinsDeltaPL-1]); }						
			if (reco_Pn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { reco_Pn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

			if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
			if (EQE > ArrayNBinsEQE[NBinsEQE]) { EQE = 0.5 * (ArrayNBinsEQE[NBinsEQE] + ArrayNBinsEQE[NBinsEQE-1]); }
			if (reco_Q2 > ArrayNBinsQ2[NBinsQ2]) { reco_Q2 = 0.5 * (ArrayNBinsQ2[NBinsQ2] + ArrayNBinsQ2[NBinsQ2-1]); }

			if (reco_A > ArrayNBinsA[NBinsA]) { reco_A = 0.5 * (ArrayNBinsA[NBinsA] + ArrayNBinsA[NBinsA-1]); }	
			if (kMiss > ArrayNBinskMiss[NBinskMiss]) { kMiss = 0.5 * (ArrayNBinskMiss[NBinskMiss] + ArrayNBinskMiss[NBinskMiss-1]); }														
			if (MissMomentum > ArrayNBinsPMiss[NBinsPMiss]) { MissMomentum = 0.5 * (ArrayNBinsPMiss[NBinsPMiss] + ArrayNBinsPMiss[NBinsPMiss-1]); }
			if (PMissMinus > ArrayNBinsPMissMinus[NBinsPMissMinus]) { PMissMinus = 0.5 * (ArrayNBinsPMissMinus[NBinsPMissMinus] + ArrayNBinsPMissMinus[NBinsPMissMinus-1]); }

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

			if (ECal < ArrayNBinsECal[0]) { ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			
			if (EQE < ArrayNBinsEQE[0]) { EQE = 0.5 * (ArrayNBinsEQE[0] + ArrayNBinsEQE[1]); }			
			if (reco_Ptx < ArrayNBinsDeltaPtx[0]) { reco_Ptx = 0.5 * (ArrayNBinsDeltaPtx[0] + ArrayNBinsDeltaPtx[1]); }
			if (reco_Pty < ArrayNBinsDeltaPty[0]) { reco_Pty = 0.5 * (ArrayNBinsDeltaPty[0] + ArrayNBinsDeltaPty[1]); }
			if (reco_PL < ArrayNBinsDeltaPL[0]) { reco_PL = 0.5 * (ArrayNBinsDeltaPL[0] + ArrayNBinsDeltaPL[1]); }						
			if (reco_A < ArrayNBinsA[0]) { reco_A = 0.5 * (ArrayNBinsA[0] + ArrayNBinsA[1]); }
			if (PMissMinus < ArrayNBinsPMissMinus[0]) { PMissMinus = 0.5 * (ArrayNBinsPMissMinus[0] + ArrayNBinsPMissMinus[1]); }			

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

			// Relative angles

			double DeltaThetaProtonMuon_Deg = Reco_DeltaTheta->at(0);
			double DeltaPhiProtonMuon_Deg = Reco_DeltaPhi->at(0);		

			// ------------------------------------------------------------------------------------------------------------------------

			double MuonVertexDistance = (CandidateMuonStart-RecoVertex).Mag();	
			double ProtonVertexDistance = (CandidateProtonStart-RecoVertex).Mag();	

			// -------------------------------------------------------------------------------------------------------------------------

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

				if (VectorCuts[i] == "_PID_NuScore" && !(reco_p_LLR_Score < ProtonLLRPIDScore) ) 
					{ PassedSelection = false; }
				if (VectorCuts[i] == "_PID_NuScore" && (reco_p_LLR_Score < ProtonLLRPIDScore) ) { PIDCounter ++; }

//				if (VectorCuts[i] == "_NuScore" && !( NuScore > MinimumNuScore) )  { PassedSelection = false; }
//				if (VectorCuts[i] == "_NuScore" && ( NuScore > MinimumNuScore) && (reco_p_LLR_Score > ProtonThreePlaneChi2LogLikelihoodCut) ) { NuScoreCounter ++; }


			}

			if (PassedSelection == false) { continue; }

			NEventsPassingSelectionCuts++;

			// -------------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Make sure that the same events fill the same plots

			if (reco_Pmu < ArrayNBinsMuonMomentum[0]) { continue; }
			if (reco_Pmu_cos_theta < ArrayNBinsMuonCosTheta[0]) { continue; }
			if (reco_Pmu_phi < ArrayNBinsMuonPhi[0]) { continue; }

			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }
			if (reco_Pp_cos_theta < ArrayNBinsProtonCosTheta[0]) { continue; }
			if (reco_Pp_phi < ArrayNBinsProtonPhi[0]) { continue; }

//			if (MissMomentum < ArrayNBinsPMiss[0]) { continue; }
//			if (PMissMinus < ArrayNBinsPMissMinus[0]) { continue; }
//			if (kMiss < ArrayNBinskMiss[0]) { continue; }

			if (TransMissMomentum < ArrayNBinsDeltaPT[0]) { continue; }
			if (DeltaAlphaT < ArrayNBinsDeltaAlphaT[0]) { continue; }
			if (DeltaPhiT < ArrayNBinsDeltaPhiT[0]) { continue; }

//			if (ECal < ArrayNBinsECal[0]) { continue; }
//			if (EQE < ArrayNBinsEQE[0]) { continue; }
//			if (reco_Q2 < ArrayNBinsQ2[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

			if (reco_Pmu > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
			if (reco_Pmu_cos_theta > ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]) { continue; }
			if (reco_Pmu_phi > ArrayNBinsMuonPhi[NBinsMuonPhi]) { continue; }

			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }
			if (reco_Pp_cos_theta > ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]) { continue; }
			if (reco_Pp_phi > ArrayNBinsProtonPhi[NBinsProtonPhi]) { continue; }

//			if (MissMomentum > ArrayNBinsPMiss[NBinsPMiss]) { continue; }
//			if (PMissMinus > ArrayNBinsPMissMinus[NBinsPMissMinus]) { continue; }
//			if (kMiss > ArrayNBinskMiss[NBinskMiss]) { continue; }

//			if (TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { continue; }
			if (DeltaAlphaT > ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]) { continue; }
			if (DeltaPhiT > ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]) { continue; }

//			if (ECal > ArrayNBinsECal[NBinsECal]) { continue; }
//			if (EQE > ArrayNBinsEQE[NBinsEQE]) { continue; }
//			if (reco_Q2 > ArrayNBinsQ2[NBinsQ2]) { continue; }

			KinematicsCounter++;				

			// --------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;

			double true_kMiss = -1;
			double true_PMissMinus = -1;
			double true_PMiss = -1;	

			double true_PL = -1;
			double true_Pn = -1;
			double true_Ptx = -1;
			double true_Pty = -1;
			double true_A = -1;	
		
			double true_TransMissMomentum = -1;
			double true_DeltaAlphaT = -1;
			double true_DeltaPhiT = -1;
			double true_ECal = -1;
			double true_EQE = -1;
			double true_Q2 = -1;			
			double true_nu = -1;

			// Only for MC to obtain true vales			
			
			if (
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

				true_kMiss = True_kMiss->at(0);
				true_PMissMinus = True_PMissMinus->at(0);
				true_PMiss = True_PMiss->at(0);	

				true_PL = True_PL->at(0);
				true_Pn = True_Pn->at(0);
				true_Ptx = True_Ptx->at(0);
				true_Pty = True_Pty->at(0);
				true_A = True_A->at(0);
			
				true_TransMissMomentum = True_Pt->at(0);
				true_DeltaAlphaT = True_DeltaAlphaT->at(0);
				true_DeltaPhiT = True_DeltaPhiT->at(0);
				true_ECal = True_ECal->at(0);
				true_EQE = True_EQE->at(0);
				true_Q2 = True_Q2->at(0);			
				true_nu = True_Ev - True_CandidateMu_P->at(0);	

				// ----------------------------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------------------------

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


				if (true_TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { true_TransMissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
				if (true_Ptx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { true_Ptx = 0.5 * (ArrayNBinsDeltaPtx[NBinsDeltaPtx] + ArrayNBinsDeltaPtx[NBinsDeltaPtx-1]); }
				if (true_Pty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { true_Pty = 0.5 * (ArrayNBinsDeltaPty[NBinsDeltaPty] + ArrayNBinsDeltaPty[NBinsDeltaPty-1]); }
				if (true_PL > ArrayNBinsDeltaPL[NBinsDeltaPL]) { true_PL = 0.5 * (ArrayNBinsDeltaPL[NBinsDeltaPL] + ArrayNBinsDeltaPL[NBinsDeltaPL-1]); }						
				if (true_Pn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { true_Pn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

				if (true_ECal > ArrayNBinsECal[NBinsECal]) { true_ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
				if (true_EQE > ArrayNBinsEQE[NBinsEQE]) { true_EQE = 0.5 * (ArrayNBinsEQE[NBinsEQE] + ArrayNBinsEQE[NBinsEQE-1]); }
				if (true_Q2 > ArrayNBinsQ2[NBinsQ2]) { true_Q2 = 0.5 * (ArrayNBinsQ2[NBinsQ2] + ArrayNBinsQ2[NBinsQ2-1]); }

				if (true_A > ArrayNBinsA[NBinsA]) { true_A = 0.5 * (ArrayNBinsA[NBinsA] + ArrayNBinsA[NBinsA-1]); }	
				if (true_kMiss > ArrayNBinskMiss[NBinskMiss]) { true_kMiss = 0.5 * (ArrayNBinskMiss[NBinskMiss] + ArrayNBinskMiss[NBinskMiss-1]); }														
				if (true_PMiss > ArrayNBinsPMiss[NBinsPMiss]) { true_PMiss = 0.5 * (ArrayNBinsPMiss[NBinsPMiss] + ArrayNBinsPMiss[NBinsPMiss-1]); }
				if (true_PMissMinus > ArrayNBinsPMissMinus[NBinsPMissMinus]) { true_PMissMinus = 0.5 * (ArrayNBinsPMissMinus[NBinsPMissMinus] + ArrayNBinsPMissMinus[NBinsPMissMinus-1]); }

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
				
				if (true_ECal < ArrayNBinsECal[0]) { true_ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			
				if (true_EQE < ArrayNBinsEQE[0]) { true_EQE = 0.5 * (ArrayNBinsEQE[0] + ArrayNBinsEQE[1]); }				
				if (true_Ptx < ArrayNBinsDeltaPtx[0]) { true_Ptx = 0.5 * (ArrayNBinsDeltaPtx[0] + ArrayNBinsDeltaPtx[1]); }
				if (true_Pty < ArrayNBinsDeltaPty[0]) { true_Pty = 0.5 * (ArrayNBinsDeltaPty[0] + ArrayNBinsDeltaPty[1]); }
				if (true_PL < ArrayNBinsDeltaPL[0]) { true_PL = 0.5 * (ArrayNBinsDeltaPL[0] + ArrayNBinsDeltaPL[1]); }						
				if (true_A < ArrayNBinsA[0]) { true_A = 0.5 * (ArrayNBinsA[0] + ArrayNBinsA[1]); }
				if (true_PMissMinus < ArrayNBinsPMissMinus[0]) { true_PMissMinus = 0.5 * (ArrayNBinsPMissMinus[0] + ArrayNBinsPMissMinus[1]); }				

				// ----------------------------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------------------------						
				
			} // End of if statement: Only for MC to obtain true vales

			// ------------------------------------------------------------------------------------------------------------------------
			
			// CCQElike analysis
			// Use the DeltaTheta, DeltaPhi & Pt cuts

			bool RecoCCQElike = false;
				
			if (

				TMath::Abs(DeltaPhiProtonMuon_Deg - 180.) < 35.
				&& TMath::Abs(DeltaThetaProtonMuon_Deg - 90.) < 55.			
				&& TransMissMomentum < 0.35

				&& reco_Pmu > CCQEArrayNBinsMuonMomentum[0]
				&& reco_Pmu < CCQEArrayNBinsMuonMomentum[CCQENBinsMuonMomentum]
				&& reco_Pp > CCQEArrayNBinsProtonMomentum[0]
				&& reco_Pp < CCQEArrayNBinsProtonMomentum[CCQENBinsProtonMomentum]	

				&& reco_Pmu_cos_theta > CCQEArrayNBinsMuonCosTheta[0]
				&& reco_Pmu_cos_theta < CCQEArrayNBinsMuonCosTheta[CCQENBinsMuonCosTheta]
				&& reco_Pp_cos_theta > CCQEArrayNBinsProtonCosTheta[0]
				&& reco_Pp_cos_theta < CCQEArrayNBinsProtonCosTheta[CCQENBinsProtonCosTheta]							
				
			) { RecoCCQElike = true; }			


			bool TrueCCQElike = false;
				
			if ( string(fWhichSample).find("Overlay") != std::string::npos ) {

				if (

					TMath::Abs(True_DeltaPhi->at(0) - 180.) < 35.
					&& TMath::Abs(True_DeltaTheta->at(0) - 90.) < 55.			
					&& true_TransMissMomentum < 0.35

					&& True_CandidateMu_P->at(0) > CCQEArrayNBinsMuonMomentum[0]
					&& True_CandidateMu_P->at(0) < CCQEArrayNBinsMuonMomentum[CCQENBinsMuonMomentum]
					&& True_CandidateP_P->at(0) > CCQEArrayNBinsProtonMomentum[0]
					&& True_CandidateP_P->at(0) < CCQEArrayNBinsProtonMomentum[CCQENBinsProtonMomentum]	

					&& True_CandidateMu_CosTheta->at(0) > CCQEArrayNBinsMuonCosTheta[0]
					&& True_CandidateMu_CosTheta->at(0) < CCQEArrayNBinsMuonCosTheta[CCQENBinsMuonCosTheta]
					&& True_CandidateP_CosTheta->at(0) > CCQEArrayNBinsProtonCosTheta[0]
					&& True_CandidateP_CosTheta->at(0) < CCQEArrayNBinsProtonCosTheta[CCQENBinsProtonCosTheta]							
					
				) { TrueCCQElike = true; }

			}			

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D LLP Only Study // Purity & Efficiency Study

			for (int WhichBin = 0; WhichBin < NBinsLLRPID; WhichBin++) {

				double LocalThres = MinLLRPID + LLRPIDStep * WhichBin;

				if (reco_p_LLR_Score < LocalThres) {
					
					RecoMuonCosThetaPlotLLPOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

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

			//----------------------------------------//

			// 1D Length Only Study

			for (int WhichBin = 0; WhichBin < NBinsMuonLength; WhichBin++) {

				double LocalThres = MinMuonLength + MuonLengthStep * WhichBin;

				if (l_muCandidate > LocalThres) {
					
					RecoMuonCosThetaPlotLengthOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			//----------------------------------------//

			// Indices for 2D analysis

			int DeltaPTTwoDIndex = tools.ReturnIndex(TransMissMomentum, TwoDArrayNBinsDeltaPT);
			int DeltaAlphaTTwoDIndex = tools.ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
			int MuonCosThetaTwoDIndex = tools.ReturnIndex(reco_Pmu_cos_theta, TwoDArrayNBinsMuonCosTheta);
			int ProtonCosThetaTwoDIndex = tools.ReturnIndex(reco_Pp_cos_theta, TwoDArrayNBinsProtonCosTheta);
			int DeltaPtxTwoDIndex = tools.ReturnIndex(reco_Ptx, TwoDArrayNBinsDeltaPtx);
			int DeltaPtyTwoDIndex = tools.ReturnIndex(reco_Pty, TwoDArrayNBinsDeltaPty);	
			int MuonMomentumTwoDIndex = tools.ReturnIndex(reco_Pmu, TwoDArrayNBinsMuonMomentum);
			int ProtonMomentumTwoDIndex = tools.ReturnIndex(reco_Pp, TwoDArrayNBinsProtonMomentum);			

			// 2D analysis treated in 1D grid

			int SerialMuonMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices,MuonCosThetaTwoDIndex,reco_Pmu);
			int SerialProtonMomentumInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,reco_Pp);				
			int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,DeltaAlphaT);
			int SerialDeltaAlphaTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,DeltaAlphaT);
			int SerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,DeltaPTTwoDIndex,DeltaAlphaT);
			int SerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TransMissMomentum);
			int SerialDeltaPnInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,reco_Pn);						
			int SerialDeltaPhiTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices,DeltaPTTwoDIndex,DeltaPhiT);
			int SerialDeltaPnInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaPTSlices,DeltaPTTwoDIndex,reco_Pn);	
			int SerialECalInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPTSlices,DeltaPTTwoDIndex,ECal);
			int SerialECalInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,ECal);																		
			int SerialProtonCosThetaInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices,MuonCosThetaTwoDIndex,reco_Pp_cos_theta);
			int SerialDeltaPtyInDeltaPtxIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices,DeltaPtxTwoDIndex,reco_Pty);				
			int SerialDeltaPtxInDeltaPtyIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtxInDeltaPtySlices,DeltaPtyTwoDIndex,reco_Ptx);
			int SerialDeltaPTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TransMissMomentum);
			int SerialDeltaPTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TransMissMomentum);														

			// 3D analysis treated in 1D grid

			int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,ECal);
			int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,ECal);
			int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,ECal);								

			//----------------------------------------//

			// No weight to be applied in the multiplicity plots
			RecoPi0Plot->Fill(NumberPi0); 
			RecoNeutronPlot->Fill(NumberNeutrons); 

			RecoNuPlot->Fill(true_nu,weight);
			RecoEvPlot->Fill(True_Ev,weight);
			RecoNuScorePlot->Fill(NuScore,weight);
			RecoFlashScorePlot->Fill(FlashScore,weight);
			RecoDistancePlot->Fill(distance,weight);
			RecoLengthDifferencePlot->Fill(LengthDifference,weight);

			RecoVertexXPlot->Fill(Vertex_X->at(0),weight);
			RecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
			RecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
			RecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);			

			RecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
			RecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

			RecoMuonLengthPlot->Fill(l_muCandidate,weight);
//			RecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
			RecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

			if (CandidateMu_EndContainment->at(0) == 1) { RecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
			if (CandidateMu_EndContainment->at(0) == 0) { RecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

			RecoMuonMomentumPlot->Fill(reco_Pmu,weight);
			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			RecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
			RecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);

			RecoProtonLengthPlot->Fill(l_pCandidate,weight);
//			RecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
			RecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

			if (CandidateMu_EndContainment->at(0) == 1) { RecoContainedMuonMomentumPlot->Fill(reco_Pmu,weight); }
			if (CandidateMu_EndContainment->at(0) == 0) { RecoUncontainedMuonMomentumPlot->Fill(reco_Pmu,weight); }

			RecoProtonMomentumPlot->Fill(reco_Pp,weight);
			RecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
			RecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

			RecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			if (reco_Pmu_cos_theta > 0.9) { RecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
			else { RecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
			RecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			RecokMissPlot->Fill(kMiss,weight);
			RecoPMissMinusPlot->Fill(PMissMinus,weight);
			RecoPMissPlot->Fill(MissMomentum,weight);

			RecoDeltaPLPlot->Fill(reco_PL,weight);
			RecoDeltaPnPlot->Fill(reco_Pn,weight);
			RecoDeltaPtxPlot->Fill(reco_Ptx,weight);
			RecoDeltaPtyPlot->Fill(reco_Pty,weight);
			RecoAPlot->Fill(reco_A,weight);

			RecoDeltaPTPlot->Fill(TransMissMomentum,weight);
			RecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
			RecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

			RecoECalPlot->Fill(ECal,weight);
			RecoEQEPlot->Fill(EQE,weight);
			RecoQ2Plot->Fill(reco_Q2,weight);

			RecoMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight);
			if (CandidateMu_EndContainment->at(0) == 1) { RecoContainedMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight); }
			if (CandidateMu_EndContainment->at(0) == 0) { RecoUncontainedMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight); }
			
			// 2D Analysis
		
			RecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
			RecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);

			RecoLengthPvsLengthMu2D->Fill(l_muCandidate,l_pCandidate,weight);

			//------------------------------//

			// Ecal in DeltaPT & DeltaAlphaT bins

			RecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
			RecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
			RecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);							

			// 2D DeltaAlphaT in DeltaPT bins

			RecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);
			RecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);
			RecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);			
			RecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
			RecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);									
			RecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);
			RecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);			
			RecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);
			RecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);			
			RecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
			RecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);			
			RecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);														
			RecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);
			RecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);	
			RecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);
			RecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);

			// 2D analysis in 1D grid

			SerialRecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
			SerialRecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
			SerialRecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
			SerialRecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
			SerialRecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
			SerialRecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
			SerialRecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);			
			SerialRecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);						
			SerialRecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
			SerialRecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
			SerialRecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
			SerialRecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
			SerialRecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
			SerialRecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
			SerialRecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
			SerialRecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);

			// 3D analysis in 1D grid

			SerialRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
			SerialRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
			SerialRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);						

			// -------------------------------------------------------------------------------------------------------------------------

			// Chi2 PID Studies

			 RecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
			 RecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);

			if (RecoCCQElike) {

				RecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
				RecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
				RecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

				RecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
				RecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
				RecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

				RecoCCQEECalPlot->Fill(ECal,weight);
				RecoCCQEQ2Plot->Fill(reco_Q2,weight);								

			}

			// -------------------------------------------------------------------------------------------------------------------------

			// Playground for Reco plots

			

			// ---------------------------------------------------------------------------------------------------------------------------

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

//				TVector3 True_CandidateMuonVertex(True_CandidateMu_StartX->at(0),True_CandidateMu_StartY->at(0),True_CandidateMu_StartZ->at(0));
//				TVector3 True_CandidateProtonVertex(True_CandidateP_StartX->at(0),True_CandidateP_StartY->at(0),True_CandidateP_StartZ->at(0));

				// CC1p Signal

				if ( CC1p == 1 && NumberPi0 == 0 && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg 
					 && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg 
				     && True_CandidateMu_StartContainment->at(0) == 1
 
				     && True_CandidateMu_P->at(0) > ArrayNBinsMuonMomentum[0] 
				     && True_CandidateP_P->at(0) > ArrayNBinsProtonMomentum[0]
				     && True_CandidateMu_CosTheta->at(0) > ArrayNBinsMuonCosTheta[0] 
                     && True_CandidateP_CosTheta->at(0) > ArrayNBinsProtonCosTheta[0]
				     && True_CandidateMu_Phi->at(0) > ArrayNBinsMuonPhi[0] 
				     && True_CandidateP_Phi->at(0) > ArrayNBinsProtonPhi[0]
				     && true_TransMissMomentum > ArrayNBinsDeltaPT[0]
				     && true_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]
				     && true_DeltaPhiT > ArrayNBinsDeltaPhiT[0]

				     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
				     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				     && True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta] 
				     && True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
				     && True_CandidateMu_Phi->at(0) < ArrayNBinsMuonPhi[NBinsMuonPhi] 
				     && True_CandidateP_Phi->at(0) < ArrayNBinsProtonPhi[NBinsProtonPhi]
				     //&& true_TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				     && true_DeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]

				) {

					// --------------------------------------------------------------------------------------------------
				
					CC1pEventsPassingSelectionCuts++;
					myRunTxtFile << endl << "CC1p0pi signal event" << endl;

					// ---------------------------------------------------------------------------------------------------------------------------
					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D LLP only study

					for (int WhichBin = 0; WhichBin < NBinsLLRPID; WhichBin++) {

						double LocalThres = MinLLRPID + LLRPIDStep * WhichBin;

						if (reco_p_LLR_Score < LocalThres) {
							
							CC1pRecoMuonCosThetaPlotLLPOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

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

					for (int WhichBin = 0; WhichBin < NBinsMuonLength; WhichBin++) {

						double LocalThres = MinMuonLength + MuonLengthStep * WhichBin;

						if (l_muCandidate > LocalThres) {
							
							CC1pRecoMuonCosThetaPlotLengthOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

						}

					}

					// --------------------------------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------------------

					// 1D Plots using True level info for selected CC1p events  

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(True_CandidateMu_P->at(0),2.) );
					double true_Nu = True_Ev - true_MuonEnergy;

					CC1pTrueNuPlot->Fill(true_Nu,weight);
					CC1pTrueEvPlot->Fill(True_Ev,weight);
					CC1pTrueVertexXPlot->Fill(Vertex_X->at(0),weight);
					CC1pTrueVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CC1pTrueVertexZPlot->Fill(Vertex_Z->at(0),weight);

					CC1pTrueMuonMomentumPlot->Fill(True_CandidateMu_P->at(0),weight);
					CC1pTrueMuonCosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
					CC1pTrueMuonCosThetaSingleBinPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
					CC1pTrueMuonPhiPlot->Fill(True_CandidateMu_Phi->at(0),weight);

					CC1pTrueProtonMomentumPlot->Fill(True_CandidateP_P->at(0),weight);
					CC1pTrueProtonCosThetaPlot->Fill(True_CandidateP_CosTheta->at(0),weight);
					CC1pTrueProtonPhiPlot->Fill(True_CandidateP_Phi->at(0),weight);

					CC1pTrueDeltaPTPlot->Fill(true_TransMissMomentum,weight);
					CC1pTrueDeltaAlphaTPlot->Fill(true_DeltaAlphaT,weight);
					CC1pTrueDeltaPhiTPlot->Fill(true_DeltaPhiT,weight);

					CC1pTrueDeltaPLPlot->Fill(true_PL,weight);
					CC1pTrueDeltaPnPlot->Fill(true_Pn,weight);
					CC1pTrueDeltaPtxPlot->Fill(true_Ptx,weight);
					CC1pTrueDeltaPtyPlot->Fill(true_Pty,weight);
					CC1pTrueAPlot->Fill(true_A,weight);
					CC1pTruekMissPlot->Fill(true_kMiss,weight);
					CC1pTruePMissPlot->Fill(true_PMiss,weight);
					CC1pTruePMissMinusPlot->Fill(true_PMissMinus,weight);

					CC1pTrueECalPlot->Fill(true_ECal,weight);
					CC1pTrueEQEPlot->Fill(true_EQE,weight);
					CC1pTrueQ2Plot->Fill(true_Q2,weight);

					//------------------------------//

					int TrueDeltaPTTwoDIndex = tools.ReturnIndex(true_TransMissMomentum, TwoDArrayNBinsDeltaPT);
					int TrueDeltaAlphaTTwoDIndex = tools.ReturnIndex(true_DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int TrueMuonCosThetaTwoDIndex = tools.ReturnIndex(True_CandidateMu_CosTheta->at(0), TwoDArrayNBinsMuonCosTheta);
					int TrueProtonCosThetaTwoDIndex = tools.ReturnIndex(True_CandidateP_CosTheta->at(0), TwoDArrayNBinsProtonCosTheta);
					int TrueDeltaPtxTwoDIndex = tools.ReturnIndex(true_Ptx, TwoDArrayNBinsDeltaPtx);
					int TrueDeltaPtyTwoDIndex = tools.ReturnIndex(true_Pty, TwoDArrayNBinsDeltaPty);	
					int TrueMuonMomentumTwoDIndex = tools.ReturnIndex(True_CandidateMu_P->at(0), TwoDArrayNBinsMuonMomentum);
					int TrueProtonMomentumTwoDIndex = tools.ReturnIndex(True_CandidateP_P->at(0), TwoDArrayNBinsProtonMomentum);						

					// 2D analysis treated in 1D grid

					int TrueSerialMuonMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices,TrueMuonCosThetaTwoDIndex,True_CandidateMu_P->at(0));
					int TrueSerialProtonMomentumInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices,TrueProtonCosThetaTwoDIndex,True_CandidateP_P->at(0));				
					int TrueSerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,TrueMuonCosThetaTwoDIndex,true_DeltaAlphaT);
					int TrueSerialDeltaAlphaTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices,TrueProtonCosThetaTwoDIndex,true_DeltaAlphaT);
					int TrueSerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,TrueDeltaPTTwoDIndex,true_DeltaAlphaT);
					int TrueSerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,TrueDeltaAlphaTTwoDIndex,true_TransMissMomentum);
					int TrueSerialDeltaPnInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices,TrueDeltaAlphaTTwoDIndex,true_Pn);										
					int TrueSerialDeltaPhiTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices,TrueDeltaPTTwoDIndex,true_DeltaPhiT);
					int TrueSerialDeltaPnInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaPTSlices,TrueDeltaPTTwoDIndex,true_Pn);	
					int TrueSerialECalInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPTSlices,TrueDeltaPTTwoDIndex,true_ECal);															
					int TrueSerialECalInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaAlphaTSlices,TrueDeltaAlphaTTwoDIndex,true_ECal);
					int TrueSerialProtonCosThetaInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices,TrueMuonCosThetaTwoDIndex,True_CandidateP_CosTheta->at(0));
					int TrueSerialDeltaPtyInDeltaPtxIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices,TrueDeltaPtxTwoDIndex,true_Pty);				
					int TrueSerialDeltaPtxInDeltaPtyIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtxInDeltaPtySlices,TrueDeltaPtyTwoDIndex,true_Ptx);
					int TrueSerialDeltaPTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices,TrueMuonCosThetaTwoDIndex,true_TransMissMomentum);
					int TrueSerialDeltaPTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices,TrueProtonCosThetaTwoDIndex,true_TransMissMomentum);														

					SerialCC1pTrueDeltaPT_InMuonCosThetaPlot->Fill(TrueSerialDeltaPTInMuonCosThetaIndex,weight);
					SerialCC1pTrueDeltaPT_InProtonCosThetaPlot->Fill(TrueSerialDeltaPTInProtonCosThetaIndex,weight);
					SerialCC1pTrueDeltaPT_InDeltaAlphaTPlot->Fill(TrueSerialDeltaPTInDeltaAlphaTIndex,weight);					
					SerialCC1pTrueMuonMomentum_InMuonCosThetaPlot->Fill(TrueSerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialCC1pTrueProtonMomentum_InProtonCosThetaPlot->Fill(TrueSerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialCC1pTrueDeltaAlphaT_InMuonCosThetaPlot->Fill(TrueSerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialCC1pTrueDeltaAlphaT_InProtonCosThetaPlot->Fill(TrueSerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialCC1pTrueDeltaAlphaT_InDeltaPTPlot->Fill(TrueSerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialCC1pTrueDeltaPhiT_InDeltaPTPlot->Fill(TrueSerialDeltaPhiTInDeltaPTIndex,weight);
					SerialCC1pTrueDeltaPn_InDeltaPTPlot->Fill(TrueSerialDeltaPnInDeltaPTIndex,weight);
					SerialCC1pTrueDeltaPn_InDeltaAlphaTPlot->Fill(TrueSerialDeltaPnInDeltaAlphaTIndex,weight);																													
					SerialCC1pTrueProtonCosTheta_InMuonCosThetaPlot->Fill(TrueSerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialCC1pTrueDeltaPty_InDeltaPtxPlot->Fill(TrueSerialDeltaPtyInDeltaPtxIndex,weight);
					SerialCC1pTrueDeltaPtx_InDeltaPtyPlot->Fill(TrueSerialDeltaPtxInDeltaPtyIndex,weight);
					SerialCC1pTrueECal_InDeltaPTPlot->Fill(TrueSerialECalInDeltaPTIndex,weight);
					SerialCC1pTrueECal_InDeltaAlphaTPlot->Fill(TrueSerialECalInDeltaAlphaTIndex,weight);											

					// 3D analysis treated in 1D grid

					int TrueSerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,TrueDeltaPTTwoDIndex,TrueDeltaAlphaTTwoDIndex,true_ECal);
					int TrueSerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,TrueMuonCosThetaTwoDIndex,TrueMuonMomentumTwoDIndex,true_ECal);
					int TrueSerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,TrueProtonCosThetaTwoDIndex,TrueProtonMomentumTwoDIndex,true_ECal);						

					// 3D analysis in 1D grid

					SerialCC1pTrueECal_InDeltaPTDeltaAlphaTPlot->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCC1pTrueECal_InMuonCosThetaMuonMomentumPlot->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCC1pTrueECal_InProtonCosThetaProtonMomentumPlot->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,weight);																																						

					//----------------------------------------//

					// 2D analysis (uncorrelated)

					CC1pTrueDeltaPT_InMuonCosThetaTwoDPlot[TrueMuonCosThetaTwoDIndex]->Fill(true_TransMissMomentum,weight);		
					CC1pTrueDeltaPT_InProtonCosThetaTwoDPlot[TrueProtonCosThetaTwoDIndex]->Fill(true_TransMissMomentum,weight);
					CC1pTrueDeltaPT_InDeltaAlphaTTwoDPlot[TrueDeltaAlphaTTwoDIndex]->Fill(true_TransMissMomentum,weight);						
					CC1pTrueMuonMomentum_InMuonCosThetaTwoDPlot[TrueMuonCosThetaTwoDIndex]->Fill(True_CandidateMu_P->at(0),weight);	
					CC1pTrueProtonMomentum_InProtonCosThetaTwoDPlot[TrueProtonCosThetaTwoDIndex]->Fill(True_CandidateP_P->at(0),weight);					
					CC1pTrueDeltaAlphaT_InDeltaPTTwoDPlot[TrueDeltaPTTwoDIndex]->Fill(true_DeltaAlphaT,weight);	
					CC1pTrueDeltaAlphaT_InMuonCosThetaTwoDPlot[TrueMuonCosThetaTwoDIndex]->Fill(true_DeltaAlphaT,weight);
					CC1pTrueDeltaAlphaT_InProtonCosThetaTwoDPlot[TrueProtonCosThetaTwoDIndex]->Fill(true_DeltaAlphaT,weight);															
					CC1pTrueDeltaPhiT_InDeltaPTTwoDPlot[TrueDeltaPTTwoDIndex]->Fill(true_DeltaPhiT,weight);	
					CC1pTrueDeltaPn_InDeltaPTTwoDPlot[TrueDeltaPTTwoDIndex]->Fill(true_Pn,weight);
					CC1pTrueDeltaPn_InDeltaAlphaTTwoDPlot[TrueDeltaAlphaTTwoDIndex]->Fill(true_Pn,weight);					
					CC1pTrueProtonCosTheta_InMuonCosThetaTwoDPlot[TrueMuonCosThetaTwoDIndex]->Fill(True_CandidateP_CosTheta->at(0),weight);						
					CC1pTrueDeltaPty_InDeltaPtxTwoDPlot[TrueDeltaPtxTwoDIndex]->Fill(true_Pty,weight);
					CC1pTrueDeltaPtx_InDeltaPtyTwoDPlot[TrueDeltaPtyTwoDIndex]->Fill(true_Ptx,weight);										
					CC1pTrueECal_InDeltaAlphaTTwoDPlot[TrueDeltaAlphaTTwoDIndex]->Fill(true_ECal,weight);
					CC1pTrueECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(true_ECal,weight);					

					// 3D analysis (uncorrelated)

					CC1pTrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TrueDeltaPTTwoDIndex][TrueDeltaAlphaTTwoDIndex]->Fill(true_ECal,weight);	
					CC1pTrueECal_InMuonCosThetaMuonMomentumTwoDPlot[TrueMuonCosThetaTwoDIndex][TrueMuonMomentumTwoDIndex]->Fill(ECal,weight);
					CC1pTrueECal_InProtonCosThetaProtonMomentumTwoDPlot[TrueProtonCosThetaTwoDIndex][TrueProtonMomentumTwoDIndex]->Fill(ECal,weight);										

					//----------------------------------------//

					// 1D Reco Plots for the selected CC1p events 

					CC1pRecoNuPlot->Fill(true_nu,weight);
					CC1pRecoEvPlot->Fill(True_Ev,weight);
					CC1pRecoNuScorePlot->Fill(NuScore,weight);
					CC1pRecoFlashScorePlot->Fill(FlashScore,weight);
					CC1pRecoDistancePlot->Fill(distance,weight);
					CC1pRecoLengthDifferencePlot->Fill(LengthDifference,weight);

					CC1pRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CC1pRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CC1pRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
					CC1pRecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);					

					CC1pRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					CC1pRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CC1pRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CC1pRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CC1pRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CC1pRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CC1pRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CC1pRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CC1pRecoContainedMuonMomentumPlot->Fill(reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CC1pRecoUncontainedMuonMomentumPlot->Fill(reco_Pmu,weight); }

					CC1pRecoProtonLengthPlot->Fill(l_pCandidate,weight);
//					CC1pRecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
					CC1pRecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

					CC1pRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					if (reco_Pmu_cos_theta > 0.9) { CC1pRecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
					else { CC1pRecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
					CC1pRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CC1pRecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
					CC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);
					CC1pRecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

					CC1pRecokMissPlot->Fill(kMiss,weight);
					CC1pRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CC1pRecoPMissPlot->Fill(MissMomentum,weight);

					CC1pRecoDeltaPLPlot->Fill(reco_PL,weight);
					CC1pRecoDeltaPnPlot->Fill(reco_Pn,weight);
					CC1pRecoDeltaPtxPlot->Fill(reco_Ptx,weight);
					CC1pRecoDeltaPtyPlot->Fill(reco_Pty,weight);
					CC1pRecoAPlot->Fill(reco_A,weight);

					CC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CC1pRecoECalPlot->Fill(ECal,weight);
					CC1pRecoEQEPlot->Fill(EQE,weight);
					CC1pRecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					CC1pRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
					CC1pRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);

					CC1pRecoMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight);
					if (CandidateMu_EndContainment->at(0) == 1) { CC1pRecoContainedMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CC1pRecoUncontainedMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight); }

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 2D DeltaAlphaT in DeltaPT bins
					
					CC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);	
					CC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);
					CC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);					
					CC1pRecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
					CC1pRecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);															
					CC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);
					CC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);					
					CC1pRecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);
					CC1pRecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);					
					CC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
					CC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);					
					CC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);	
					CC1pRecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);
					CC1pRecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);	
					CC1pRecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);
					CC1pRecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);	

					// 2D analysis in 1D grid

					SerialCC1pRecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialCC1pRecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialCC1pRecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialCC1pRecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialCC1pRecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialCC1pRecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialCC1pRecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialCC1pRecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);										
					SerialCC1pRecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialCC1pRecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialCC1pRecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialCC1pRecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialCC1pRecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialCC1pRecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialCC1pRecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialCC1pRecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);		

					// 3D analysis in 1D grid

					SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCC1pRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCC1pRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);																																															

					// --------------------------------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------------------

					// 2D Plots Kinematic Variables

					CC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu);
					CC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp);

					CC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoMuonCosThetaSingleBinPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta);

					CC1pRecoMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi);
					CC1pRecoProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi);

					POTScaledCC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu,weight);
					POTScaledCC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp,weight);

					POTScaledCC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
					POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
					POTScaledCC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta,weight);

					POTScaledCC1pRecoMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi,weight);
					POTScaledCC1pRecoProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi,weight);

					CC1pTProfProtonMomTheta->Fill(reco_Pp,TMath::ACos(reco_Pp_cos_theta)*180./TMath::Pi(), (reco_Pp - True_CandidateP_P->at(0)) / True_CandidateP_P->at(0) * 100.);	
					CC1pTProfProtonMomPhi->Fill(reco_Pp,reco_Pp_phi, (reco_Pp - True_CandidateP_P->at(0)) / True_CandidateP_P->at(0) * 100.);										

					// 2D analysis in 1D grid

					SerialCC1pRecoMuonMomentum_InMuonCosThetaPlot2D->Fill(TrueSerialMuonMomentumInMuonCosThetaIndex,SerialMuonMomentumInMuonCosThetaIndex);			
					SerialCC1pRecoProtonMomentum_InProtonCosThetaPlot2D->Fill(TrueSerialProtonMomentumInProtonCosThetaIndex,SerialProtonMomentumInProtonCosThetaIndex);
					SerialCC1pRecoDeltaAlphaT_InMuonCosThetaPlot2D->Fill(TrueSerialDeltaAlphaTInMuonCosThetaIndex,SerialDeltaAlphaTInMuonCosThetaIndex);
					SerialCC1pRecoDeltaAlphaT_InProtonCosThetaPlot2D->Fill(TrueSerialDeltaAlphaTInProtonCosThetaIndex,SerialDeltaAlphaTInProtonCosThetaIndex);
					SerialCC1pRecoDeltaAlphaT_InDeltaPTPlot2D->Fill(TrueSerialDeltaAlphaTInDeltaPTIndex,SerialDeltaAlphaTInDeltaPTIndex);
					SerialCC1pRecoDeltaPT_InDeltaAlphaTPlot2D->Fill(TrueSerialDeltaPTInDeltaAlphaTIndex,SerialDeltaPTInDeltaAlphaTIndex);
					SerialCC1pRecoDeltaPn_InDeltaAlphaTPlot2D->Fill(TrueSerialDeltaPnInDeltaAlphaT,SerialDeltaPnInDeltaAlphaTIndex);						
					SerialCC1pRecoECal_InDeltaAlphaTPlot2D->Fill(TrueSerialECalInDeltaAlphaTIndex,SerialECalInDeltaAlphaTIndex);									
					SerialCC1pRecoDeltaPhiT_InDeltaPTPlot2D->Fill(TrueSerialDeltaPhiTInDeltaPTIndex,SerialDeltaPhiTInDeltaPTIndex);
					SerialCC1pRecoDeltaPn_InDeltaPTPlot2D->Fill(TrueSerialDeltaPnInDeltaPTIndex,SerialDeltaPnInDeltaPTIndex);	
					SerialCC1pRecoECal_InDeltaPTPlot2D->Fill(TrueSerialECalInDeltaPTIndex,SerialECalInDeltaPTIndex);																							
					SerialCC1pRecoProtonCosTheta_InMuonCosThetaPlot2D->Fill(TrueSerialProtonCosThetaInMuonCosThetaIndex,SerialProtonCosThetaInMuonCosThetaIndex);
					SerialCC1pRecoDeltaPty_InDeltaPtxPlot2D->Fill(TrueSerialDeltaPtyInDeltaPtxIndex,SerialDeltaPtyInDeltaPtxIndex);
					SerialCC1pRecoDeltaPtx_InDeltaPtyPlot2D->Fill(TrueSerialDeltaPtxInDeltaPtyIndex,SerialDeltaPtxInDeltaPtyIndex);	
					SerialCC1pRecoDeltaPT_InMuonCosThetaPlot2D->Fill(TrueSerialDeltaPTInMuonCosThetaIndex,SerialDeltaPTInMuonCosThetaIndex);
					SerialCC1pRecoDeltaPT_InProtonCosThetaPlot2D->Fill(TrueSerialDeltaPTInProtonCosThetaIndex,SerialDeltaPTInProtonCosThetaIndex);

					SerialPOTScaledCC1pRecoMuonMomentum_InMuonCosThetaPlot2D->Fill(TrueSerialMuonMomentumInMuonCosThetaIndex,SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialPOTScaledCC1pRecoProtonMomentum_InProtonCosThetaPlot2D->Fill(TrueSerialProtonMomentumInProtonCosThetaIndex,SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialPOTScaledCC1pRecoDeltaAlphaT_InMuonCosThetaPlot2D->Fill(TrueSerialDeltaAlphaTInMuonCosThetaIndex,SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialPOTScaledCC1pRecoDeltaAlphaT_InProtonCosThetaPlot2D->Fill(TrueSerialDeltaAlphaTInProtonCosThetaIndex,SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialPOTScaledCC1pRecoDeltaAlphaT_InDeltaPTPlot2D->Fill(TrueSerialDeltaAlphaTInDeltaPTIndex,SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialPOTScaledCC1pRecoDeltaPT_InDeltaAlphaTPlot2D->Fill(TrueSerialDeltaPTInDeltaAlphaTIndex,SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialPOTScaledCC1pRecoDeltaPn_InDeltaAlphaTPlot->Fill(TrueSerialDeltaPnInDeltaAlphaTIndex,SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialPOTScaledCC1pRecoECal_InDeltaAlphaTPlot2D->Fill(TrueSerialECalInDeltaAlphaTIndex,SerialECalInDeltaAlphaTIndex,weight);					
					SerialPOTScaledCC1pRecoDeltaPhiT_InDeltaPTPlot2D->Fill(TrueSerialDeltaPhiTInDeltaPTIndex,SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialPOTScaledCC1pRecoDeltaPn_InDeltaPTPlot2D->Fill(TrueSerialDeltaPnInDeltaPTIndex,SerialDeltaPnInDeltaPTIndex,weight);	
					SerialPOTScaledCC1pRecoECal_InDeltaPTPlot2D->Fill(TrueSerialECalInDeltaPTIndex,SerialECalInDeltaPTIndex,weight);																							
					SerialPOTScaledCC1pRecoProtonCosTheta_InMuonCosThetaPlot2D->Fill(TrueSerialProtonCosThetaInMuonCosThetaIndex,SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialPOTScaledCC1pRecoDeltaPty_InDeltaPtxPlot2D->Fill(TrueSerialDeltaPtyInDeltaPtxIndex,SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialPOTScaledCC1pRecoDeltaPtx_InDeltaPtyPlot2D->Fill(TrueSerialDeltaPtxInDeltaPtyIndex,SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialPOTScaledCC1pRecoDeltaPT_InMuonCosThetaPlot2D->Fill(TrueSerialDeltaPTInMuonCosThetaIndex,SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialPOTScaledCC1pRecoDeltaPT_InProtonCosThetaPlot2D->Fill(TrueSerialDeltaPTInProtonCosThetaIndex,SerialDeltaPTInProtonCosThetaIndex,weight);					

					// 3D analysis in 1D grid

					SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot2D->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,SerialECalInDeltaPTDeltaAlphaTIndex);
					SerialCC1pRecoECal_InMuonCosThetaMuonMomentumPlot2D->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,SerialECalInMuonCosThetaMuonMomentumIndex);
					SerialCC1pRecoECal_InProtonCosThetaProtonMomentumPlot2D->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,SerialECalInProtonCosThetaProtonMomentumIndex);

					// 3D analysis in 1D grid

					SerialPOTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTPlot2D->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialPOTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumPlot2D->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialPOTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumPlot2D->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,SerialECalInProtonCosThetaProtonMomentumIndex,weight);

					// -----------------------------------------------------------------------------------------------------

					// True Level STV

					CC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum);
					CC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT);
					CC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT);

					CC1pRecoDeltaPLPlot2D->Fill(true_PL,reco_PL);
					CC1pRecoDeltaPnPlot2D->Fill(true_Pn,reco_Pn);
					CC1pRecoDeltaPtxPlot2D->Fill(true_Ptx,reco_Ptx);
					CC1pRecoDeltaPtyPlot2D->Fill(true_Pty,reco_Pty);
					CC1pRecoAPlot2D->Fill(true_A,reco_A);
					CC1pRecokMissPlot2D->Fill(true_kMiss,kMiss);
					CC1pRecoPMissPlot2D->Fill(true_PMiss,MissMomentum);
					CC1pRecoPMissMinusPlot2D->Fill(true_PMissMinus,PMissMinus);

					POTScaledCC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum,weight);
					POTScaledCC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);
					POTScaledCC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT,weight);

					POTScaledCC1pRecoDeltaPLPlot2D->Fill(true_PL,reco_PL,weight);
					POTScaledCC1pRecoDeltaPnPlot2D->Fill(true_Pn,reco_Pn,weight);
					POTScaledCC1pRecoDeltaPtxPlot2D->Fill(true_Ptx,reco_Ptx,weight);
					POTScaledCC1pRecoDeltaPtyPlot2D->Fill(true_Pty,reco_Pty,weight);
					POTScaledCC1pRecoAPlot2D->Fill(true_A,reco_A,weight);
					POTScaledCC1pRecokMissPlot2D->Fill(true_kMiss,kMiss,weight);
					POTScaledCC1pRecoPMissPlot2D->Fill(true_PMiss,MissMomentum,weight);
					POTScaledCC1pRecoPMissMinusPlot2D->Fill(true_PMissMinus,PMissMinus,weight);

					// -----------------------------------------------------------------------------------------------------

					// True level energy reconstruction & Q2

					CC1pRecoECalPlot2D->Fill(true_ECal,ECal);
					CC1pRecoEQEPlot2D->Fill(true_EQE,EQE);
					CC1pRecoQ2Plot2D->Fill(true_Q2,reco_Q2);

					POTScaledCC1pRecoECalPlot2D->Fill(true_ECal,ECal,weight);
					POTScaledCC1pRecoEQEPlot2D->Fill(true_EQE,EQE,weight);
					POTScaledCC1pRecoQ2Plot2D->Fill(true_Q2,reco_Q2,weight);

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal,weight);

					CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(true_ECal,ECal);
										
					CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(true_ECal,ECal);	
					POTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(true_ECal,ECal,weight);														

					// 2D DeltaAlphaT in DeltaPT bins
					
					CC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_DeltaAlphaT,DeltaAlphaT);
					POTScaledCC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);	

					CC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot2D[DeltaAlphaTTwoDIndex]->Fill(true_TransMissMomentum,TransMissMomentum);
					POTScaledCC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot2D[DeltaAlphaTTwoDIndex]->Fill(true_TransMissMomentum,TransMissMomentum,weight);											

					CC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot2D[DeltaAlphaTTwoDIndex]->Fill(true_Pn,reco_Pn);
					POTSacledCC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot2D[DeltaAlphaTTwoDIndex]->Fill(true_Pn,reco_Pn,weight);					

					CC1pRecoECal_InDeltaAlphaTTwoDPlot2D[DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InDeltaAlphaTTwoDPlot2D[DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal,weight);	

					CC1pRecoECal_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_ECal,ECal);	
					POTScaledCC1pRecoECal_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_ECal,ECal,weight);													

					CC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(true_DeltaAlphaT,DeltaAlphaT);
					POTScaledCC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);

					CC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot2D[ProtonCosThetaTwoDIndex]->Fill(true_DeltaAlphaT,DeltaAlphaT);
					POTScaledCC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot2D[ProtonCosThetaTwoDIndex]->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);											

					CC1pRecoDeltaPT_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(true_TransMissMomentum,TransMissMomentum);
					POTScaledCC1pRecoDeltaPT_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(true_TransMissMomentum,TransMissMomentum,weight);																												

					CC1pRecoDeltaPT_InProtonCosThetaTwoDPlot2D[ProtonCosThetaTwoDIndex]->Fill(true_TransMissMomentum,TransMissMomentum);
					POTScaledCC1pRecoDeltaPT_InProtonCosThetaTwoDPlot2D[ProtonCosThetaTwoDIndex]->Fill(true_TransMissMomentum,TransMissMomentum,weight);					

					CC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(True_CandidateMu_P->at(0),reco_Pmu);
					POTScaledCC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(True_CandidateMu_P->at(0),reco_Pmu,weight);

					CC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta);
					POTScaledCC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot2D[MuonCosThetaTwoDIndex]->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta,weight);

					CC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot2D[ProtonCosThetaTwoDIndex]->Fill(True_CandidateP_P->at(0),reco_Pp);					
					POTScaledCC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot2D[ProtonCosThetaTwoDIndex]->Fill(True_CandidateP_P->at(0),reco_Pp,weight);

					CC1pRecoDeltaPhiT_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_DeltaPhiT,DeltaPhiT);
					POTScaledCC1pRecoDeltaPhiT_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_DeltaPhiT,DeltaPhiT,weight);

					CC1pRecoDeltaPn_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_Pn,reco_Pn);
					POTScaledCC1pRecoDeltaPn_InDeltaPTTwoDPlot2D[DeltaPTTwoDIndex]->Fill(true_Pn,reco_Pn,weight);	

					CC1pRecoDeltaPty_InDeltaPtxTwoDPlot2D[DeltaPtxTwoDIndex]->Fill(true_Pty,reco_Pty);	
					POTScaledCC1pRecoDeltaPty_InDeltaPtxTwoDPlot2D[DeltaPtxTwoDIndex]->Fill(true_Pty,reco_Pty,weight);

					CC1pRecoDeltaPtx_InDeltaPtyTwoDPlot2D[DeltaPtyTwoDIndex]->Fill(true_Ptx,reco_Ptx);
					POTScaledCC1pRecoDeltaPtx_InDeltaPtyTwoDPlot2D[DeltaPtyTwoDIndex]->Fill(true_Ptx,reco_Ptx,weight);																		

					// -----------------------------------------------------------------------------------------------------

					CC1pRecoLengthPvsLengthMu2D->Fill(l_muCandidate,l_pCandidate,weight);

					// -------------------------------------------------------------------------------------------------------------------------

					// Playground for CC1p STV resolution plots
					// Warning: they are called resolution plots but they are just reco - true differences

					double CC1pMuonMomentumReso = ( reco_Pmu - True_CandidateMu_P->at(0) );
					double CC1pDeltaPTReso = ( TransMissMomentum - true_TransMissMomentum );
					double CC1pDeltaAlphaTReso = ( DeltaAlphaT - true_DeltaAlphaT );
					double CC1pDeltaPhiTReso = ( DeltaPhiT - true_DeltaPhiT );

					// -------------------------------------------------------------------------------------------------------------------------

					// Playground for CC1p true momenta (longitudinal & perpendicular) ratios

					TVector3 TrueCandidateMuon(1,1,1);
					TrueCandidateMuon.SetMag(True_CandidateMu_P->at(0));
					TrueCandidateMuon.SetPhi(True_CandidateMu_Phi->at(0));
					TrueCandidateMuon.SetTheta(TMath::ACos(True_CandidateMu_CosTheta->at(0)));

					TVector3 TrueCandidateProton(1,1,1);
					TrueCandidateProton.SetMag(True_CandidateP_P->at(0));
					TrueCandidateProton.SetPhi(True_CandidateP_Phi->at(0));
					TrueCandidateProton.SetTheta(TMath::ACos(True_CandidateP_CosTheta->at(0)));

					CC1pRecoMuonTrueMomentumLongitudinalRatio->Fill(TVector3CandidateMuon.Z() / TVector3CandidateMuon.Mag(),weight);
					CC1pRecoMuonTrueMomentumTransverseRatio->Fill(TVector3CandidateMuon.Pt() / TVector3CandidateMuon.Mag(),weight);

					CC1pRecoProtonTrueMomentumLongitudinalRatio->Fill(TVector3CandidateProton.Z() / TVector3CandidateProton.Mag(),weight);
					CC1pRecoProtonTrueMomentumTransverseRatio->Fill(TVector3CandidateProton.Pt() / TVector3CandidateProton.Mag(),weight);

					CC1pTrueMuonTrueMomentumLongitudinalRatio->Fill(TrueCandidateMuon.Z() / TrueCandidateMuon.Mag(),weight);
					CC1pTrueMuonTrueMomentumTransverseRatio->Fill(TrueCandidateMuon.Pt() / TrueCandidateMuon.Mag(),weight);

					CC1pTrueProtonTrueMomentumLongitudinalRatio->Fill(TrueCandidateProton.Z() / TrueCandidateProton.Mag(),weight);
					CC1pTrueProtonTrueMomentumTransverseRatio->Fill(TrueCandidateProton.Pt() / TrueCandidateProton.Mag(),weight);

					// -------------------------------------------------------------------------------------------------------------------------

					if (RecoCCQElike && TrueCCQElike) {

						CC1pRecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						CC1pRecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						CC1pRecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						CC1pRecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						CC1pRecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						CC1pRecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						CC1pRecoCCQEECalPlot->Fill(ECal,weight);
						CC1pRecoCCQEQ2Plot->Fill(reco_Q2,weight);

						CC1pRecoCCQEMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu);
						CC1pRecoCCQEProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp);

						CC1pRecoCCQEMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
						CC1pRecoCCQEProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta);

						CC1pRecoCCQEMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi);
						CC1pRecoCCQEProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi);								

						CC1pRecoCCQEECalPlot2D->Fill(true_ECal,ECal);
						CC1pRecoCCQEQ2Plot2D->Fill(true_Q2,reco_Q2);

						POTScaledCC1pRecoCCQEMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu,weight);
						POTScaledCC1pRecoCCQEProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp,weight);

						POTScaledCC1pRecoCCQEMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
						POTScaledCC1pRecoCCQEProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta,weight);

						POTScaledCC1pRecoCCQEMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi,weight);
						POTScaledCC1pRecoCCQEProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi,weight);

						POTScaledCC1pRecoCCQEECalPlot2D->Fill(true_ECal,ECal,weight);
						POTScaledCC1pRecoCCQEQ2Plot2D->Fill(true_Q2,reco_Q2,weight);

						CC1pTrueCCQEMuonMomentumPlot->Fill(True_CandidateMu_P->at(0),weight);
						CC1pTrueCCQEMuonCosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
						CC1pTrueCCQEProtonMomentumPlot->Fill(True_CandidateP_P->at(0),weight);
						CC1pTrueCCQEProtonCosThetaPlot->Fill(True_CandidateP_CosTheta->at(0),weight);
						CC1pTrueCCQEECalPlot->Fill(true_ECal,weight);
						CC1pTrueCCQEQ2Plot->Fill(true_Q2,weight);

					} 
					
					if (RecoCCQElike && !TrueCCQElike) {

						NonCC1pRecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						NonCC1pRecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						NonCC1pRecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						NonCC1pRecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						NonCC1pRecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						NonCC1pRecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						NonCC1pRecoCCQEECalPlot->Fill(ECal,weight);
						NonCC1pRecoCCQEQ2Plot->Fill(reco_Q2,weight);

					}

					if ( TMath::Abs(CC1pDeltaAlphaTReso) > 150) { myRunTxtFile << endl << "Bad delta alphaT case" << endl; }

				} // End of the CC1p signal

				// -------------------------------------------------------------------------------------------------------------

				// Non-CC1p beam related background or EXT BNB

				else {

					NonCC1pEventsPassingSelectionCuts++;

					if ( 
						(TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg ) ||
						(TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg )
					) { 
							MisIndetifiedMuonAsPion++; 
							ManualNonCC1pEventsPassingSelectionCuts++; 
							myRunTxtFile << endl << "pi-p event" << endl; 
							NonCC1pRecoPionMomentumPlot->Fill(reco_Pmu,weight);
					}
					else if (CC1p1pi == 1) { CC1p1piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC1p1pi event" << endl; }
					else if (CC2p1pi == 1) { CC2p1piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC2p1pi event" << endl;}
					else if (CC2p == 1) { CC2pEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC2p event" << endl; }
					else if (CC3p == 1) { CC3pEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC3p event" << endl; }
					else if (CC3p1pi == 1) { CC3p1piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								 myRunTxtFile << endl << "CC3p1pi event" << endl; }
					else if (CC3p2pi == 1) { CC3p2piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								 myRunTxtFile << endl << "CC3p2pi event" << endl;}
					//else if (CC3p3pi == 1) { CC3p3piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; }
					else if (CC4p0pi == 1) { CC4p0piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC4p event" << endl; }
					else if (NumberPi0 > 0) { pi0Included++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "pi0 event" << endl;}
					else if ( 
						( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == MuonPdg ) ||
						( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == MuonPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg )
						) 
						{ MisIndetifiedProtonAsPion++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "mu-pi event" << endl; }
					else if ( 
						( CandidateP_MCParticle_Pdg->at(0) == DeuteriumPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == MuonPdg ) ||
						( CandidateP_MCParticle_Pdg->at(0) == MuonPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == DeuteriumPdg )
						) 
						{ MisIndetifiedProtonAsDeuterium++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "mu-D event" << endl;}
					else if ( 
						(CandidateMu_MCParticle_Pdg->at(0) == DeuteriumPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ProtonPdg )  ||
						(CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == DeuteriumPdg ) )			
						{ MisIndetifiedMuonAsDeuterium++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "p-D event" << endl;}
					else if ( 
						( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == DeuteriumPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg ) ||
						( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == DeuteriumPdg )
						) 
						{ DeuteriumPionPairs++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "pi-D event" << endl; }
					else if ( CandidateP_MCParticle_Pdg->at(0) == HeliumPdg && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) 
						{ MisIndetifiedProtonAsHelium++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "mu-He event" << endl; }
					else if (  
						(CandidateMu_MCParticle_Pdg->at(0) == HeliumPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) ||
						(CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg && CandidateP_MCParticle_Pdg->at(0) == HeliumPdg)
						) 
						{ MisIndetifiedMuonAsHelium++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "p-He event" << endl; }
					else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ElectronPdg && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) 			
						{ MisIndetifiedProtonAsElectron++;  ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "e-mu event" << endl; }
					else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ElectronPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == ElectronPdg) 
						{ MisIndetifiedMuPToElectronElectron++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "e-e event" << endl; }
					else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == MuonPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == MuonPdg) 
						{ 
 
							if (True_CandidateMu_StartZ->at(0) == True_CandidateP_StartZ->at(0)) { BrokenMuonTrack++; }
							else { MisIndetifiedMuPToMuMu++; OtherMCBkg++; }
							ManualNonCC1pEventsPassingSelectionCuts++;
							myRunTxtFile << endl << "Broken muon track" << endl;
					}
					else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg) 
						{ MisIndetifiedMuPToPiPi++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "pi-pi event" << endl; }
					else if ( CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) 
						{ 
							if (True_CandidateMu_StartZ->at(0) == True_CandidateP_StartZ->at(0)) { BrokenProtonTrack++; }
							else { MisIndetifiedMuonAsProton++; }  
							ManualNonCC1pEventsPassingSelectionCuts++; 
							myRunTxtFile << endl << "Broken proton track" << endl;
					}
					else if ( CandidateMu_MCParticle_Pdg->at(0) == -MuonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) 
						{ MisIndetifiedMuonAsAntiMuon++;  ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "mu+ - p event" << endl;}
					else if ( True_CandidateMu_StartContainment->at(0) == 0 ) 
						{ CandidateMuon_MCParticle_OutFV++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "True vertex out FV" << endl;}
					else if ( CandidateMu_MCParticle_Pdg->at(0) == -99. || CandidateP_MCParticle_Pdg->at(0) == -99.) 
						{ InTimeCosmics++;  ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "In-time cosmic event" << endl; }
					else if ( CandidateMu_MCParticle_Pdg->at(0) == NeutronPdg && CandidateP_MCParticle_Pdg->at(0) == NeutronPdg) 
						{ NeutronPairCounter++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "n-n event" << endl; }
					else if ( NC == 1) { NCEvents++;  ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "NC event" << endl; }
					else if (NumberMuons == 2) { DoubleMuonEvents++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "mu-mu event" << endl; }
					else if (NumberProtons == 5) { CC5pXpiEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								       myRunTxtFile << endl << "CC5p event" << endl; }
					else if (NumberProtons == 6) { CC6pXpiEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								       myRunTxtFile << endl << "CC6p event" << endl; }
					else if (NumberProtons == 7) { CC7pXpiEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								       myRunTxtFile << endl << "CC7p event" << endl;}
					else if (NumberProtons == 8) { CC8pXpiEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								       myRunTxtFile << endl << "CC8p event" << endl; }
					//else if (NumberChargedPions > 1) { CCNpXpiEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; }
					else if (NumberProtons == 0) { CC0pXpiEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								       myRunTxtFile << endl << "CC0p event" << endl; }
					else if ( CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg && CandidateP_MCParticle_Pdg->at(0) == MuonPdg) 
						{ FlippedMuPEvents++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "flipped mu-p event" << endl; }
					else if ( CandidateMu_MCParticle_Pdg->at(0) == ArgonPdg && CandidateP_MCParticle_Pdg->at(0) == ArgonPdg) 
						{ ArArEvents++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "Ar-Ar event" << endl; }
					else if ( 
						( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == ProtonPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == NeutronPdg ) ||
						( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == NeutronPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ProtonPdg )
						) 
						{ ProtonNeutronEvents++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "p-n event" << endl; }
					else if ( 
						( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == ProtonPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == KaonPdg ) ||
						( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == KaonPdg && TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ProtonPdg )
						) 
						{ ProtonKaonEvents++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "p-K event" << endl; }
					else if (
					     !  (True_CandidateMu_P->at(0) > ArrayNBinsMuonMomentum[0] 
					     && True_CandidateP_P->at(0) > ArrayNBinsProtonMomentum[0]
					     && True_CandidateMu_CosTheta->at(0) > ArrayNBinsMuonCosTheta[0] 
		                 && True_CandidateP_CosTheta->at(0) > ArrayNBinsProtonCosTheta[0]
					     && True_CandidateMu_Phi->at(0) > ArrayNBinsMuonPhi[0] 
					     && True_CandidateP_Phi->at(0) > ArrayNBinsProtonPhi[0]
					     && true_TransMissMomentum > ArrayNBinsDeltaPT[0]
					     && true_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]
					     && true_DeltaPhiT > ArrayNBinsDeltaPhiT[0]

					     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
					     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
					     && True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta] 
					     && True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
					     && True_CandidateMu_Phi->at(0) < ArrayNBinsMuonPhi[NBinsMuonPhi] 
					     && True_CandidateP_Phi->at(0) < ArrayNBinsProtonPhi[NBinsProtonPhi]
					     //&& true_TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
					     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
					     && true_DeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT])

					) {

						OutCommonRange++; ManualNonCC1pEventsPassingSelectionCuts++;
						myRunTxtFile << endl << "Out of range event" << endl;
					}

					else if ( True_CandidateMu_StartX->at(0) != True_CandidateP_StartX->at(0) 
					|| True_CandidateMu_StartY->at(0) != True_CandidateP_StartY->at(0)
					|| True_CandidateMu_StartZ->at(0) != True_CandidateP_StartZ->at(0) ) 
						{ MultipleVertices++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "Multiple vertices event" << endl; }

					else { 

						ManualNonCC1pEventsPassingSelectionCuts++;
						//std::cout << "new MC background topology";
						//std::cout << " CandidateMu_MCParticle_Pdg = " << CandidateMu_MCParticle_Pdg->at(0);
						//std::cout << " CandidateP_MCParticle_Pdg = " << CandidateP_MCParticle_Pdg->at(0) << endl; 
						//std::cout << " NumberProtons = " << NumberProtons << endl; 
						//std::cout << " NumberChargedPions = " << NumberChargedPions << endl; 
						//std::cout << " NumberMuons = " << NumberMuons << endl; 
						//std::cout << " NumberPi0 = " << NumberPi0 << endl; 
						//std::cout << " NumberNeutrons = " << NumberNeutrons << endl; 

					}

					NonCC1pRecoNuPlot->Fill(true_nu,weight);
					NonCC1pRecoEvPlot->Fill(True_Ev,weight);
					NonCC1pRecoNuScorePlot->Fill(NuScore,weight);
					NonCC1pRecoFlashScorePlot->Fill(FlashScore,weight);
					NonCC1pRecoDistancePlot->Fill(distance,weight);
					NonCC1pRecoLengthDifferencePlot->Fill(LengthDifference,weight);

					NonCC1pRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					NonCC1pRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					NonCC1pRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
					NonCC1pRecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);					

					NonCC1pRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					NonCC1pRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					NonCC1pRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					NonCC1pRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					NonCC1pRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { NonCC1pRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { NonCC1pRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					NonCC1pRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { NonCC1pRecoContainedMuonMomentumPlot->Fill(reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { NonCC1pRecoUncontainedMuonMomentumPlot->Fill(reco_Pmu,weight); }

					NonCC1pRecoProtonLengthPlot->Fill(l_pCandidate,weight);
//					NonCC1pRecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
					NonCC1pRecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

					NonCC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					NonCC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					NonCC1pRecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

					NonCC1pRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					if (reco_Pmu_cos_theta > 0.9) { NonCC1pRecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
					else { NonCC1pRecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
					NonCC1pRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					NonCC1pRecokMissPlot->Fill(kMiss,weight);
					NonCC1pRecoPMissMinusPlot->Fill(PMissMinus,weight);
					NonCC1pRecoPMissPlot->Fill(MissMomentum,weight);

					NonCC1pRecoDeltaPLPlot->Fill(reco_PL,weight);
					NonCC1pRecoDeltaPnPlot->Fill(reco_Pn,weight);
					NonCC1pRecoDeltaPtxPlot->Fill(reco_Ptx,weight);
					NonCC1pRecoDeltaPtyPlot->Fill(reco_Pty,weight);
					NonCC1pRecoAPlot->Fill(reco_A,weight);

					NonCC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					NonCC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					NonCC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					NonCC1pRecoECalPlot->Fill(ECal,weight);
					NonCC1pRecoEQEPlot->Fill(EQE,weight);
					NonCC1pRecoQ2Plot->Fill(reco_Q2,weight);

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					NonCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					NonCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					NonCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 2D DeltaAlphaT in DeltaPT bins
					
					NonCC1pRecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);	
					NonCC1pRecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);	
					NonCC1pRecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);						
					NonCC1pRecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
					NonCC1pRecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);													
					NonCC1pRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);	
					NonCC1pRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);					
					NonCC1pRecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);	
					NonCC1pRecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);					
					NonCC1pRecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
					NonCC1pRecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);						
					NonCC1pRecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);	
					NonCC1pRecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);	
					NonCC1pRecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);	
					NonCC1pRecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);	
					NonCC1pRecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);		

					// 2D analysis in 1D grid

					SerialNonCC1pRecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialNonCC1pRecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialNonCC1pRecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialNonCC1pRecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialNonCC1pRecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialNonCC1pRecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);	
					SerialNonCC1pRecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialNonCC1pRecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);									
					SerialNonCC1pRecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialNonCC1pRecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialNonCC1pRecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialNonCC1pRecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialNonCC1pRecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialNonCC1pRecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialNonCC1pRecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialNonCC1pRecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);		

					// 3D analysis in 1D grid

					SerialNonCC1pRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialNonCC1pRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialNonCC1pRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);																																											
					
					// 2D Analysis
				
					NonCC1pRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
					NonCC1pRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);	

					NonCC1pRecoMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight);
					if (CandidateMu_EndContainment->at(0) == 1) { NonCC1pRecoContainedMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { NonCC1pRecoUncontainedMuonMomentumVsLengthPlot->Fill(l_muCandidate,reco_Pmu,weight); }

					if (RecoCCQElike && !TrueCCQElike) {

						NonCC1pRecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						NonCC1pRecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						NonCC1pRecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						NonCC1pRecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						NonCC1pRecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						NonCC1pRecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						NonCC1pRecoCCQEECalPlot->Fill(ECal,weight);
						NonCC1pRecoCCQEQ2Plot->Fill(reco_Q2,weight);

					}

				} // End of the Non-CC1p beam related background

				// -------------------------------------------------------------------------------------------------------------------------
				// ------------------------------------------------------------------------------------------------------------------------

				// CCQE

				if (genie_mode == 0) {

					CCQERecoNuPlot->Fill(true_nu,weight);
					CCQERecoEvPlot->Fill(True_Ev,weight);
					CCQERecoNuScorePlot->Fill(NuScore,weight);
					CCQERecoFlashScorePlot->Fill(FlashScore,weight);
					CCQERecoDistancePlot->Fill(distance,weight);
					CCQERecoLengthDifferencePlot->Fill(LengthDifference,weight);

					CCQERecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCQERecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCQERecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
					CCQERecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);					

					CCQERecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					CCQERecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CCQERecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CCQERecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CCQERecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCQERecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCQERecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CCQERecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCQERecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
					CCQERecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCQERecoContainedMuonMomentumPlot->Fill(reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCQERecoUncontainedMuonMomentumPlot->Fill(reco_Pmu,weight); }

					CCQERecoProtonLengthPlot->Fill(l_pCandidate,weight);
//					CCQERecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
					CCQERecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

					CCQERecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCQERecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					CCQERecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

					CCQERecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					if (reco_Pmu_cos_theta > 0.9) { CCQERecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
					else { CCQERecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
					CCQERecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCQERecokMissPlot->Fill(kMiss,weight);
					CCQERecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCQERecoPMissPlot->Fill(MissMomentum,weight);

					CCQERecoDeltaPLPlot->Fill(reco_PL,weight);
					CCQERecoDeltaPnPlot->Fill(reco_Pn,weight);
					CCQERecoDeltaPtxPlot->Fill(reco_Ptx,weight);
					CCQERecoDeltaPtyPlot->Fill(reco_Pty,weight);
					CCQERecoAPlot->Fill(reco_A,weight);

					CCQERecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCQERecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCQERecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCQERecoECalPlot->Fill(ECal,weight);
					CCQERecoEQEPlot->Fill(EQE,weight);
					CCQERecoQ2Plot->Fill(reco_Q2,weight);

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCQERecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCQERecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCQERecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 2D DeltaAlphaT in DeltaPT bins
					
					CCQERecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);	
					CCQERecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);
					CCQERecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);					
					CCQERecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
					CCQERecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);																
					CCQERecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);
					CCQERecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);																		
					CCQERecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);	
					CCQERecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);									
					CCQERecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
					CCQERecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);					
					CCQERecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);
					CCQERecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);	
					CCQERecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);				
					CCQERecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);
					CCQERecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);	

					// 2D analysis in 1D grid

					SerialCCQERecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialCCQERecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialCCQERecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialCCQERecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialCCQERecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialCCQERecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialCCQERecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialCCQERecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);										
					SerialCCQERecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialCCQERecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialCCQERecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialCCQERecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialCCQERecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialCCQERecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialCCQERecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialCCQERecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);	

					// 3D analysis in 1D grid

					SerialCCQERecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCQERecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCQERecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);													

					// 2D Analysis
				
					CCQERecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
					CCQERecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);

					if (RecoCCQElike) {

						CCQERecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						CCQERecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						CCQERecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						CCQERecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						CCQERecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						CCQERecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						CCQERecoCCQEECalPlot->Fill(ECal,weight);
						CCQERecoCCQEQ2Plot->Fill(reco_Q2,weight);

					}								

				} // End of CCQE selection

				// ------------------------------------------------------------------------------------------------------------------------

				// CCMEC

				if (genie_mode == 10) {

					CCMECRecoNuPlot->Fill(true_nu,weight);
					CCMECRecoEvPlot->Fill(True_Ev,weight);
					CCMECRecoNuScorePlot->Fill(NuScore,weight);
					CCMECRecoFlashScorePlot->Fill(FlashScore,weight);
					CCMECRecoDistancePlot->Fill(distance,weight);
					CCMECRecoLengthDifferencePlot->Fill(LengthDifference,weight);

					CCMECRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCMECRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCMECRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
					CCMECRecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);					

					CCMECRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					CCMECRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CCMECRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CCMECRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CCMECRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCMECRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCMECRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CCMECRecoProtonLengthPlot->Fill(l_pCandidate,weight);
//					CCMECRecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
					CCMECRecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

					CCMECRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCMECRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCMECRecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
					CCMECRecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCQERecoContainedMuonMomentumPlot->Fill(reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCQERecoUncontainedMuonMomentumPlot->Fill(reco_Pmu,weight); }

					CCMECRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCMECRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					CCMECRecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

					CCMECRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					if (reco_Pmu_cos_theta > 0.9) { CCMECRecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
					else { CCMECRecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
					CCMECRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCMECRecokMissPlot->Fill(kMiss,weight);
					CCMECRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCMECRecoPMissPlot->Fill(MissMomentum,weight);

					CCMECRecoDeltaPLPlot->Fill(reco_PL,weight);
					CCMECRecoDeltaPnPlot->Fill(reco_Pn,weight);
					CCMECRecoDeltaPtxPlot->Fill(reco_Ptx,weight);
					CCMECRecoDeltaPtyPlot->Fill(reco_Pty,weight);
					CCMECRecoAPlot->Fill(reco_A,weight);

					CCMECRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCMECRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCMECRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCMECRecoECalPlot->Fill(ECal,weight);
					CCMECRecoEQEPlot->Fill(EQE,weight);
					CCMECRecoQ2Plot->Fill(reco_Q2,weight);

					//----------------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCMECRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 2D DeltaAlphaT in DeltaPT bins
					
					CCMECRecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);	
					CCMECRecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);	
					CCMECRecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);					
					CCMECRecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCMECRecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);																							
					CCMECRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);
					CCMECRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);					
					CCMECRecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);
					CCMECRecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);					
					CCMECRecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
					CCMECRecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);					
					CCMECRecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);
					CCMECRecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);
					CCMECRecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);
					CCMECRecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);
					CCMECRecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);

					// 2D analysis in 1D grid

					SerialCCMECRecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialCCMECRecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialCCMECRecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialCCMECRecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialCCMECRecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialCCMECRecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);	
					SerialCCMECRecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialCCMECRecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);									
					SerialCCMECRecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialCCMECRecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialCCMECRecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialCCMECRecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialCCMECRecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialCCMECRecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialCCMECRecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialCCMECRecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);

					// 3D analysis in 1D grid

					SerialCCMECRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCMECRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCMECRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);											

					// 2D Analysis
				
					CCMECRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
					CCMECRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);	

					if (RecoCCQElike) {

						CCMECRecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						CCMECRecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						CCMECRecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						CCMECRecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						CCMECRecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						CCMECRecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						CCMECRecoCCQEECalPlot->Fill(ECal,weight);
						CCMECRecoCCQEQ2Plot->Fill(reco_Q2,weight);

					}							

				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCRES

				if (genie_mode == 1) {

					CCRESRecoNuPlot->Fill(true_nu,weight);
					CCRESRecoEvPlot->Fill(True_Ev,weight);
					CCRESRecoNuScorePlot->Fill(NuScore,weight);
					CCRESRecoFlashScorePlot->Fill(FlashScore,weight);
					CCRESRecoDistancePlot->Fill(distance,weight);
					CCRESRecoLengthDifferencePlot->Fill(LengthDifference,weight);

					CCRESRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCRESRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCRESRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
					CCRESRecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);					

					CCRESRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					CCRESRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CCRESRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CCRESRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CCRESRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCRESRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCRESRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CCRESRecoProtonLengthPlot->Fill(l_pCandidate,weight);
//					CCRESRecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
					CCRESRecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

					CCRESRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCRESRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCRESRecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
					CCRESRecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCRESRecoContainedMuonMomentumPlot->Fill(reco_Pmu,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCRESRecoUncontainedMuonMomentumPlot->Fill(reco_Pmu,weight); }

					CCRESRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCRESRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					CCRESRecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

					CCRESRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					if (reco_Pmu_cos_theta > 0.9) { CCRESRecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
					else { CCRESRecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
					CCRESRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCRESRecokMissPlot->Fill(kMiss,weight);
					CCRESRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCRESRecoPMissPlot->Fill(MissMomentum,weight);

					CCRESRecoDeltaPLPlot->Fill(reco_PL,weight);
					CCRESRecoDeltaPnPlot->Fill(reco_Pn,weight);
					CCRESRecoDeltaPtxPlot->Fill(reco_Ptx,weight);
					CCRESRecoDeltaPtyPlot->Fill(reco_Pty,weight);
					CCRESRecoAPlot->Fill(reco_A,weight);

					CCRESRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCRESRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCRESRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCRESRecoECalPlot->Fill(ECal,weight);
					CCRESRecoEQEPlot->Fill(EQE,weight);
					CCRESRecoQ2Plot->Fill(reco_Q2,weight);

					//----------------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCRESRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCRESRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCRESRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 2D DeltaAlphaT in DeltaPT bins
					
					CCRESRecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);		
					CCRESRecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);
					CCRESRecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);					
					CCRESRecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCRESRecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);																		
					CCRESRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);
					CCRESRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);					
					CCRESRecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);	
					CCRESRecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);									
					CCRESRecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
					CCRESRecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);					
					CCRESRecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);
					CCRESRecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);
					CCRESRecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);
					CCRESRecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);
					CCRESRecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);

					// 2D analysis in 1D grid

					SerialCCRESRecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialCCRESRecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialCCRESRecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialCCRESRecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialCCRESRecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialCCRESRecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);	
					SerialCCRESRecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialCCRESRecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);									
					SerialCCRESRecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialCCRESRecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialCCRESRecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialCCRESRecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialCCRESRecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialCCRESRecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialCCRESRecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialCCRESRecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);	

					// 3D analysis in 1D grid

					SerialCCRESRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCRESRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCRESRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);									

					// 2D Analysis
				
					CCRESRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
					CCRESRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);		

					if (RecoCCQElike) {

						CCRESRecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						CCRESRecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						CCRESRecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						CCRESRecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						CCRESRecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						CCRESRecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						CCRESRecoCCQEECalPlot->Fill(ECal,weight);
						CCRESRecoCCQEQ2Plot->Fill(reco_Q2,weight);

					}						

				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCDIS

				if (genie_mode == 2) {

					CCDISRecoNuPlot->Fill(true_nu,weight);
					CCDISRecoEvPlot->Fill(True_Ev,weight);
					CCDISRecoNuScorePlot->Fill(NuScore,weight);
					CCDISRecoFlashScorePlot->Fill(FlashScore,weight);
					CCDISRecoDistancePlot->Fill(distance,weight);
					CCDISRecoLengthDifferencePlot->Fill(LengthDifference,weight);

					CCDISRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCDISRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCDISRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);
					CCDISRecoCandidateMuEndYPlot->Fill(CandidateMu_EndY->at(0),weight);					

					CCDISRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					CCDISRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CCDISRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CCDISRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CCDISRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCDISRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCDISRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CCDISRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCDISRecoMuonCosThetaSingleBinPlot->Fill(reco_Pmu_cos_theta,weight);
					CCDISRecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);

					CCDISRecoProtonLengthPlot->Fill(l_pCandidate,weight);
//					CCDISRecodProtonTracksScorePlot->Fill(ProtonTrackScore,weight);
					CCDISRecodProtonVertexDistancePlot->Fill(ProtonVertexDistance,weight);

					CCDISRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCDISRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					CCDISRecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

					CCDISRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					if (reco_Pmu_cos_theta > 0.9) { CCDISRecoDeltaForwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); } 
					else { CCDISRecoDeltaBackwardThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight); }
					CCDISRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCDISRecokMissPlot->Fill(kMiss,weight);
					CCDISRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCDISRecoPMissPlot->Fill(MissMomentum,weight);

					CCDISRecoDeltaPLPlot->Fill(reco_PL,weight);
					CCDISRecoDeltaPnPlot->Fill(reco_Pn,weight);
					CCDISRecoDeltaPtxPlot->Fill(reco_Ptx,weight);
					CCDISRecoDeltaPtyPlot->Fill(reco_Pty,weight);
					CCDISRecoAPlot->Fill(reco_A,weight);

					CCDISRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCDISRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCDISRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCDISRecoECalPlot->Fill(ECal,weight);
					CCDISRecoEQEPlot->Fill(EQE,weight);
					CCDISRecoQ2Plot->Fill(reco_Q2,weight);

					//----------------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCDISRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCDISRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCDISRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 2D DeltaAlphaT in DeltaPT bins
					
					CCDISRecoDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);		
					CCDISRecoDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(TransMissMomentum,weight);	
					CCDISRecoDeltaPn_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(reco_Pn,weight);					
					CCDISRecoECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCDISRecoECal_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(ECal,weight);																
					CCDISRecoDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);	
					CCDISRecoDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(DeltaAlphaT,weight);									
					CCDISRecoDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);
					CCDISRecoDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TransMissMomentum,weight);					
					CCDISRecoMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pmu,weight);
					CCDISRecoProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(reco_Pp_cos_theta,weight);					
					CCDISRecoProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(reco_Pp,weight);
					CCDISRecoDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(DeltaPhiT,weight);
					CCDISRecoDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(reco_Pn,weight);
					CCDISRecoDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(reco_Pty,weight);
					CCDISRecoDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(reco_Ptx,weight);

					// 2D analysis in 1D grid

					SerialCCDISRecoMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
					SerialCCDISRecoProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
					SerialCCDISRecoDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
					SerialCCDISRecoDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
					SerialCCDISRecoDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
					SerialCCDISRecoDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
					SerialCCDISRecoDeltaPn_InDeltaAlphaTPlot->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);					
					SerialCCDISRecoECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);										
					SerialCCDISRecoDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
					SerialCCDISRecoDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
					SerialCCDISRecoECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
					SerialCCDISRecoProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
					SerialCCDISRecoDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
					SerialCCDISRecoDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
					SerialCCDISRecoDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
					SerialCCDISRecoDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);	

					// 3D analysis in 1D grid

					SerialCCDISRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCDISRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCDISRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);									

					// 2D Analysis
				
					CCDISRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu,weight);
					CCDISRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);	

					if (RecoCCQElike) {

						CCDISRecoCCQEMuonMomentumPlot->Fill(reco_Pmu,weight);
						CCDISRecoCCQEMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
						CCDISRecoCCQEMuonPhiPlot->Fill(reco_Pmu_phi,weight);

						CCDISRecoCCQEProtonMomentumPlot->Fill(reco_Pp,weight);
						CCDISRecoCCQEProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
						CCDISRecoCCQEProtonPhiPlot->Fill(reco_Pp_phi,weight);

						CCDISRecoCCQEECalPlot->Fill(ECal,weight);
						CCDISRecoCCQEQ2Plot->Fill(reco_Q2,weight);

					}							

				}

				// --------------------------------------------------------------------------------------------------------------------------

				// Overlay particle breakdown using the Backtracker for PID studies

				if (CandidateMu_MCParticle_Pdg->size() > 0 && CandidateP_MCParticle_Pdg->size() > 0 ) {

					if (CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) {

						MuonRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
//						MuonRecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);		

					}

					if (CandidateP_MCParticle_Pdg->at(0) == MuonPdg) {

						MuonRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);
//						MuonRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);		

					}

					// ----------------------------------------------------------------------------------------------------------------

					if (CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg) {

						ProtonRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
//						ProtonRecoChi2Plot->Fill(reco_Pmu_chi2,weight);		

					}

					if (CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) {

						ProtonRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);
//						ProtonRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);		
						//ProtonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

					// -----------------------------------------------------------------------------------------------------------------

					if (CandidateMu_MCParticle_Pdg->at(0) == AbsChargedPionPdg) {

						PionRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
//						PionRecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);		

					}

					if (CandidateP_MCParticle_Pdg->at(0) == AbsChargedPionPdg) {

						PionRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);
//						PionRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);		

					}

				} // End of the overlay particle breakdown using the Backtracker

				else {

					CosmicRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
					CosmicRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);

					// CosmicRecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);
					// CosmicRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);

				}

			} // End of the Overlay case and the breakdown into CC1p/NonCC1p & QE,MEC,RES,DIS

			// -------------------------------------------------------------------------------------------------------------------------

			// Storing the run/subrun/event of the candidate events 

			myRunTxtFile << "Candidate #" << NEventsPassingSelectionCuts << " Run = " << Run << ", SubRun = " << SubRun << ", Event = " << Event;

			// To quickly locate the candidate events in the event displays, we need Z & X for plane 2

			myRunTxtFile << " Vertex X = " << Vertex_X->at(0) << ", Y = " << Vertex_Y->at(0) << ", Z = " << Vertex_Z->at(0) << endl; 

			// -------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the events

		std::cout << std::endl << "Created file: " << FileName << std::endl << std::endl;
//		std::cout << "Txt info file: " << TxtName << std::endl << std::endl;

		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

		// -------------------------------------------------------------------------------------------------------------------------
		
		double nentriesError = sqrt(nentries);			
		
		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed events passing the selection criteria

		double NEventsPassingSelectionCutsError = TMath::Sqrt(NEventsPassingSelectionCuts);
		double KinematicsCounterError = TMath::Sqrt(KinematicsCounter);
		
		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed CC1p events passing the selection criteria

		double CC1pEventsPassingSelectionCutsError = sqrt(CC1pEventsPassingSelectionCuts);
		TH1D* CC1pEventPlot = new TH1D("CC1pEventPlot",";CC1p events",1,0,1);
		CC1pEventPlot->SetBinContent(1,CC1pEventsPassingSelectionCuts);
		CC1pEventPlot->SetBinError(1,CC1pEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed CC2p events passing the selection criteria

		double CC2pEventsPassingSelectionCutsError = sqrt(CC2pEventsPassingSelectionCuts);
		TH1D* CC2pEventPlot = new TH1D("CC2pEventPlot",";CC2p events",1,0,1);
		CC2pEventPlot->SetBinContent(1,CC2pEventsPassingSelectionCuts);
		CC2pEventPlot->SetBinError(1,CC2pEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Multiple Vertices

		double MultipleVerticesError = sqrt(MultipleVertices);
		TH1D* MultipleVerticesEventPlot = new TH1D("MultipleVerticesEventPlot",";Multiple vertices events",1,0,1);
		MultipleVerticesEventPlot->SetBinContent(1,MultipleVertices);
		MultipleVerticesEventPlot->SetBinError(1,MultipleVerticesError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Mis-identified muon-pion events passing the selection criteria

		double MisIndetifiedMuonAsPionError = sqrt(MisIndetifiedMuonAsPion);
		TH1D* PiPEventPlot = new TH1D("PiPEventPlot",";#pi-p events",1,0,1);
		PiPEventPlot->SetBinContent(1,MisIndetifiedMuonAsPion);
		PiPEventPlot->SetBinError(1,MisIndetifiedMuonAsPionError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// mu-p events out of common range at a truth level

		double OutCommonRangeError = sqrt(OutCommonRange);
		TH1D* OutCommonRangeEventPlot = new TH1D("OutCommonRangeEventPlot",";Truth out of common range",1,0,1);
		OutCommonRangeEventPlot->SetBinContent(1,OutCommonRange);
		OutCommonRangeEventPlot->SetBinError(1,OutCommonRangeError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed CC1p1pi events passing the selection criteria

		double CC1p1piEventsPassingSelectionCutsError = sqrt(CC1p1piEventsPassingSelectionCuts);
		TH1D* CC1p1piEventPlot = new TH1D("CC1p1piEventPlot",";CC1p1#pi events",1,0,1);
		CC1p1piEventPlot->SetBinContent(1,CC1p1piEventsPassingSelectionCuts);
		CC1p1piEventPlot->SetBinError(1,CC1p1piEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Mis-identified muon-proton events passing the selection criteria

		double MisIndetifiedMuonAsProtonError = sqrt(MisIndetifiedMuonAsProton);
		TH1D* PPEventPlot = new TH1D("PPEventPlot",";p-p events",1,0,1);
		PPEventPlot->SetBinContent(1,MisIndetifiedMuonAsProton);
		PPEventPlot->SetBinError(1,MisIndetifiedMuonAsProtonError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Events with pi0's passing the selection criteria

		double pi0IncludedError = sqrt(pi0Included);
		TH1D* NeutralPiEventPlot = new TH1D("NeutralPiEventPlot",";#pi^{0} events",1,0,1);
		NeutralPiEventPlot->SetBinContent(1,pi0Included);
		NeutralPiEventPlot->SetBinError(1,pi0IncludedError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Muon Candidate MC Particle outside FV

		double CandidateMuon_MCParticle_OutFVError = sqrt(CandidateMuon_MCParticle_OutFV);
		TH1D* TrueVertexOutFVEventPlot = new TH1D("TrueVertexOutFVEventPlot",";True vertex out FV events",1,0,1);
		TrueVertexOutFVEventPlot->SetBinContent(1,CandidateMuon_MCParticle_OutFV);
		TrueVertexOutFVEventPlot->SetBinError(1,CandidateMuon_MCParticle_OutFVError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed CC3p events passing the selection criteria

		double CC3pEventsPassingSelectionCutsError = sqrt(CC3pEventsPassingSelectionCuts);
		TH1D* CC3pEventPlot = new TH1D("CC3pEventPlot",";CC3p events",1,0,1);
		CC3pEventPlot->SetBinContent(1,CC3pEventsPassingSelectionCuts);
		CC3pEventPlot->SetBinError(1,CC3pEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Events with broken muon tracks passing the selection criteria

		double BrokenMuonTrackError = sqrt(BrokenMuonTrack);
		TH1D* BrokenMuEventPlot = new TH1D("BrokenMuEventPlot",";Broken #mu events",1,0,1);
		BrokenMuEventPlot->SetBinContent(1,BrokenMuonTrack);
		BrokenMuEventPlot->SetBinError(1,BrokenMuonTrackError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Broken proton track events passing the selection criteria

		double BrokenProtonTrackError = sqrt(BrokenProtonTrack);
		TH1D* BrokenPEventPlot = new TH1D("BrokenPEventPlot",";Broken p events",1,0,1);
		BrokenPEventPlot->SetBinContent(1,BrokenProtonTrack);
		BrokenPEventPlot->SetBinError(1,BrokenProtonTrackError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed CC4p0pi events passing the selection criteria

		double CC4p0piEventsPassingSelectionCutsError = sqrt(CC4p0piEventsPassingSelectionCuts);
		TH1D* CC4pEventPlot = new TH1D("CC4pEventPlot",";CC4p events",1,0,1);
		CC4pEventPlot->SetBinContent(1,CC4p0piEventsPassingSelectionCuts);
		CC4pEventPlot->SetBinError(1,CC4p0piEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// NC events

		double NCEventsError = sqrt(NCEvents);
		TH1D* NCEventPlot = new TH1D("NCEventPlot",";NC events",1,0,1);
		NCEventPlot->SetBinContent(1,NCEvents);
		NCEventPlot->SetBinError(1,NCEventsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed CC2p1pi events passing the selection criteria

		double CC2p1piEventsPassingSelectionCutsError = sqrt(CC2p1piEventsPassingSelectionCuts);
		TH1D* CC2p1piEventPlot = new TH1D("CC2p1piEventPlot",";CC2p1pi events",1,0,1);
		CC2p1piEventPlot->SetBinContent(1,CC2p1piEventsPassingSelectionCuts);
		CC2p1piEventPlot->SetBinError(1,CC2p1piEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// In-time cosmics

		double InTimeCosmicsError = sqrt(InTimeCosmics);
		TH1D* InTimeCosmicsEventPlot = new TH1D("InTimeCosmicsEventPlot",";In-Time Cosmics events",1,0,1);
		InTimeCosmicsEventPlot->SetBinContent(1,InTimeCosmics);
		InTimeCosmicsEventPlot->SetBinError(1,InTimeCosmicsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Mis-identified muon- anti-muon events passing the selection criteria

		double MisIndetifiedMuonAsAntiMuonError = sqrt(MisIndetifiedMuonAsAntiMuon);
		TH1D* AntiMuPEventPlot = new TH1D("AntiMuPEventPlot",";#mu-p events",1,0,1);
		AntiMuPEventPlot->SetBinContent(1,MisIndetifiedMuonAsAntiMuon);
		AntiMuPEventPlot->SetBinError(1,MisIndetifiedMuonAsAntiMuonError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Mis-identified muon-pion events passing the selection criteria

		double MisIndetifiedProtonAsElectronError = sqrt(MisIndetifiedProtonAsElectron);
		TH1D* MuEEventPlot = new TH1D("MuEEventPlot",";#mu-e events",1,0,1);
		MuEEventPlot->SetBinContent(1,MisIndetifiedProtonAsElectron);
		MuEEventPlot->SetBinError(1,MisIndetifiedProtonAsElectronError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Mis-identified muon-pion events passing the selection criteria

		double MisIndetifiedProtonAsPionError = sqrt(MisIndetifiedProtonAsPion);
		TH1D* MuPiEventPlot = new TH1D("MuPiEventPlot",";#mu-#pi events",1,0,1);
		MuPiEventPlot->SetBinContent(1,MisIndetifiedProtonAsPion);
		MuPiEventPlot->SetBinError(1,MisIndetifiedProtonAsPionError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// All NonCC1p events passing the selection criteria

		double NonCC1pEventsPassingSelectionCutsError = sqrt(NonCC1pEventsPassingSelectionCuts);
		TH1D* NonCC1pEventPlot = new TH1D("NonCC1pEventPlot",";NonCC1p events",1,0,1);
		NonCC1pEventPlot->SetBinContent(1,NonCC1pEventsPassingSelectionCuts);
		NonCC1pEventPlot->SetBinError(1,NonCC1pEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		// Other MC backgrounds

		double OtherMCBkgError = sqrt(OtherMCBkg);
		TH1D* OtherMCBkgEventPlot = new TH1D("OtherMCBkgEventPlot",";Other MC Bkg events",1,0,1);
		OtherMCBkgEventPlot->SetBinContent(1,OtherMCBkg);
		OtherMCBkgEventPlot->SetBinError(1,OtherMCBkgError);

		// -------------------------------------------------------------------------------------------------------------------------	

		std::cout << std::endl << std::endl << std::scientific << "POTWeight = " << POTWeight << std::endl << std::endl; 

		// -------------------------------------------------------------------------------------------------------------------------

		// Keeping track of the number of events during preselection and copying them to the output of the event selection

		TFile* fFile = new TFile(fPathToFile,"readonly");
		TH1D* SamdefEventsPlot = (TH1D*)(fFile->Get("SamdefEventPlot"));
		TH1D* TrackLikeDaughterEventPlot = (TH1D*)(fFile->Get("TrackLikeDaughterEventPlot"));
		TH1D* MatchedTrackLikeDaughterEventPlot = (TH1D*)(fFile->Get("MatchedTrackLikeDaughterEventPlot"));
		TH1D* MomentumThresholdEventPlot = (TH1D*)(fFile->Get("MomentumThresholdEventPlot"));
		TH1D* ContainmentEventPlot = (TH1D*)(fFile->Get("ContainmentEventPlot"));

		file->cd();
		SamdefEventsPlot->Write();
		TrackLikeDaughterEventPlot->Write();
		MatchedTrackLikeDaughterEventPlot->Write();	
		MomentumThresholdEventPlot->Write();
		ContainmentEventPlot->Write();

		// -------------------------------------------------------------------------------------------------------------------------

		// Storing the event loss of the event selection 

		TH1D* ProtonEndPointContainmentEventPlot = new TH1D("ProtonEndPointContainmentEventPlot",";Proton contained end point events",1,0,1);
		TH1D* NoFlippedTrackEventPlot = new TH1D("NoFlippedTrackEventPlot",";No flipped events",1,0,1);
		TH1D* VertexContainmentEventPlot = new TH1D("VertexContainmentEventPlot",";Proton contained end point events",1,0,1);
		TH1D* MuonQualityEventPlot = new TH1D("MuonQualityEventPlot",";Muon quality events",1,0,1);
		TH1D* PidEventPlot = new TH1D("PidEventPlot",";PID events",1,0,1);
		TH1D* NuScoreEventPlot = new TH1D("NuScoreEventPlot",";#nu score events",1,0,1);
		TH1D* CommonEventPlot = new TH1D("CommonEventPlot",";Common events",1,0,1);

		ProtonEndPointContainmentEventPlot->SetBinContent(1,ContainmentCounter);
		NoFlippedTrackEventPlot->SetBinContent(1,NoFlippedTrackCounter);
		VertexContainmentEventPlot->SetBinContent(1,ContainedVertexCounter);
		MuonQualityEventPlot->SetBinContent(1,MuonQualityCutCounter);
		PidEventPlot->SetBinContent(1,PIDCounter);
		NuScoreEventPlot->SetBinContent(1,NuScoreCounter);
		CommonEventPlot->SetBinContent(1,KinematicsCounter);

		file->cd();
		ProtonEndPointContainmentEventPlot->Write();
		NoFlippedTrackEventPlot->Write();
		VertexContainmentEventPlot->Write();
		MuonQualityEventPlot->Write();
		PidEventPlot->Write();
		NuScoreEventPlot->Write();
		CommonEventPlot->Write();

		CC1pEventPlot->Write();
		CC2pEventPlot->Write();
		MultipleVerticesEventPlot->Write();
		OutCommonRangeEventPlot->Write();
		PiPEventPlot->Write();
		CC1p1piEventPlot->Write();
		PPEventPlot->Write();
		NeutralPiEventPlot->Write();
		TrueVertexOutFVEventPlot->Write();
		CC3pEventPlot->Write();
		BrokenMuEventPlot->Write();
		BrokenPEventPlot->Write();
		CC4pEventPlot->Write();
		NCEventPlot->Write();
		CC2p1piEventPlot->Write();
		InTimeCosmicsEventPlot->Write();
		AntiMuPEventPlot->Write();
		MuEEventPlot->Write();
		MuPiEventPlot->Write();

		OtherMCBkgEventPlot->Write();
		NonCC1pEventPlot->Write();

		// -------------------------------------------------------------------------------------------------------------------------

		// Keep track of the event counters and the POT normalization

		// If Runs 1/2/3, just grab the POTWeight of the last event
		POTScalePlot->SetBinContent(1,POTWeight);

		// If combined runs, take care of the overall scaling
		if ( string(fWhichSample).find("Combined") != std::string::npos ) {

			if ( string(fWhichSample).find("BeamOn") != std::string::npos ) { POTScalePlot->SetBinContent(1,1.); }
			else if ( string(fWhichSample).find("ExtBNB") != std::string::npos ) { POTScalePlot->SetBinContent(1,FullE1DCNT_wcut_Combined / FullEXT_Combined); }
			else if ( string(fWhichSample).find("Overlay") != std::string::npos ) { 

//				TString FileName = "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PeLEETuples/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
				TString FileName = "/uboone/data/users/apapadop/PeLEETuples/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
				TFile* POTFile = new TFile(FileName,"readonly");
				TH1D* POTPlot = (TH1D*)(POTFile->Get("POTCountHist"));
				double POTCount = POTPlot->GetBinContent(1);

				POTScalePlot->SetBinContent(1,Fulltor860_wcut_Combined / POTCount); 
				POTFile->Close();

			}

		}

		POTScalePlot->SetBinError(1,0);

		NEventsPlot->SetBinContent(1,nentries);
		NEventsPlot->SetBinError(1,nentriesError);

		NSelectedPlot->SetBinContent(1,NEventsPassingSelectionCuts);
		NSelectedPlot->SetBinError(1,NEventsPassingSelectionCutsError);

		NCC1pPlot->SetBinContent(1,CC1pEventsPassingSelectionCuts);
		NCC1pPlot->SetBinError(1,CC1pEventsPassingSelectionCutsError);

		NCC1p1piPlot->SetBinContent(1,CC1p1piEventsPassingSelectionCuts);
		NCC1p1piPlot->SetBinError(1,CC1p1piEventsPassingSelectionCutsError);

		NCC2pPlot->SetBinContent(1,CC2pEventsPassingSelectionCuts);
		NCC2pPlot->SetBinError(1,CC2pEventsPassingSelectionCutsError);

		NCC2p1piPlot->SetBinContent(1,CC2p1piEventsPassingSelectionCuts);
		NCC2p1piPlot->SetBinError(1,CC2p1piEventsPassingSelectionCutsError);

		NCC3pPlot->SetBinContent(1,CC3pEventsPassingSelectionCuts);
		NCC3pPlot->SetBinError(1,CC3pEventsPassingSelectionCutsError);

		//----------------------------------------//

		cout << fWhichSample << "  " << KinematicsCounter << " events passing selection criteria" << endl << endl;

		//----------------------------------------//	

		file->cd();
		file->Write();
		file->Close();
		fFile->Close();

		fFile->Close();

		//----------------------------------------//	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
