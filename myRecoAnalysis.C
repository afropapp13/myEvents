#define myRecoAnalysis_cxx
#include "myRecoAnalysis.h"
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

void myRecoAnalysis::Loop() {

	// ---------------------------------------------------------------------------------------------------------------------------------------

	int TotalCounter = 0;
	int ContainmentCounter = 0;
	int NoFlippedTrackCounter = 0;
	int SumTruncdEdxTrackCounter = 0;
	int ContainedVertexCounter = 0;
	int MuonQualityCutCounter = 0;
	int LowMCSQualityCutCounter = 0;
	int PIDCounter = 0;
	int NuScoreCounter = 0;
	int KinematicsCounter = 0;

	// ---------------------------------------------------------------------------------------------------------------------------------------

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

	// --------------------------------------------------------------------------------------------------------------------------------------------------
		
	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

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

		// For overlays only for genie, flux and reinteraction uncertainties

		if (fUniverseIndex != -1) {

			Extension = "_"+fEventWeightLabel+"_"+ToStringInt(fUniverseIndex); 

		}

		TString FileName = PathToFiles+Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root";
		TFile* file = new TFile(FileName,"recreate");
		std::cout << std::endl << "Creating a new file: " << FileName << std::endl << std::endl << std::endl;

		// ---------------------------------------------------------------------------------------------------------------------------------------

		// Txt file to keep track of the run/subrun/event of the candidate events

		TString RunTxtName = "/uboone/data/users/apapadop/myEvents/myTxtFiles/"+UBCodeVersion+"/TxtmyRunSubRunEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
		ofstream myRunTxtFile;
		myRunTxtFile.open(RunTxtName);
		myRunTxtFile << std::fixed << std::setprecision(2);
		myRunTxtFile << fWhichSample << endl << endl;

		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots

		TH1D* RecoPi0Plot = new TH1D("RecoPi0Plot",";# #pi^{0}",4,-0.5,3.5);
		TH1D* RecoNeutronPlot = new TH1D("RecoNeutronPlot",";# neutrons",6,-0.5,5.5);

		TH1D* RecoNuPlot = new TH1D("RecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* RecoEvPlot = new TH1D("RecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* RecoNuScorePlot = new TH1D("RecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* RecoFlashScorePlot = new TH1D("RecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* RecoDistancePlot = new TH1D("RecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* RecoLengthDifferencePlot = new TH1D("RecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* RecodYZPlot = new TH1D("RecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* RecoNPEPlot = new TH1D("RecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* RecoVertexActivityPlot = new TH1D("RecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);	
//		TH1D* RecoNonZeroVertexActivityPlot = new TH1D("RecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);	

		TH1D* RecoVertexXPlot = new TH1D("RecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* RecoVertexYPlot = new TH1D("RecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* RecoVertexZPlot = new TH1D("RecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* RecoMuonLLRPIDPlot = new TH1D("RecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* RecoProtonLLRPIDPlot = new TH1D("RecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* RecoMuonLengthPlot = new TH1D("RecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* RecodMuonTracksScorePlot = new TH1D("RecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* RecodMuonVertexDistancePlot = new TH1D("RecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* RecoContainedMuonLengthPlot = new TH1D("RecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* RecoUncontainedMuonLengthPlot = new TH1D("RecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
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

//		TH1D* RecokMissPlot = new TH1D("RecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* RecoPMissMinusPlot = new TH1D("RecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* RecoPMissPlot = new TH1D("RecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

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
		TH1D* CC1pTrueMuonPhiPlot = new TH1D("CC1pTrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CC1pCCQETrueMuonMomentumPlot = new TH1D("CC1pCCQETrueMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CC1pCCQETrueMuonCosThetaPlot = new TH1D("CC1pCCQETrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CC1pCCQETrueMuonPhiPlot = new TH1D("CC1pCCQETrueMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

		TH1D* CC1pTrueProtonMomentumPlot = new TH1D("CC1pTrueProtonMomentumPlot",LabelXAxisProtonMomentum, NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CC1pTrueProtonCosThetaPlot = new TH1D("CC1pTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta, NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CC1pTrueProtonPhiPlot = new TH1D("CC1pTrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CC1pCCQETrueProtonMomentumPlot = new TH1D("CC1pCCQETrueProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TH1D* CC1pCCQETrueProtonCosThetaPlot = new TH1D("CC1pCCQETrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CC1pCCQETrueProtonPhiPlot = new TH1D("CC1pCCQETrueProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CC1pTrueDeltaPTPlot = new TH1D("CC1pTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pTrueDeltaAlphaTPlot = new TH1D("CC1pTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pTrueDeltaPhiTPlot = new TH1D("CC1pTrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CC1pTrueECalPlot = new TH1D("CC1pTrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CC1pTrueEQEPlot = new TH1D("CC1pTrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CC1pTrueQ2Plot = new TH1D("CC1pTrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CC1pCCQETrueECalPlot = new TH1D("CC1pCCQETrueECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CC1pCCQETrueQ2Plot = new TH1D("CC1pCCQETrueQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for Signal CC1p

		TH1D* CC1pRecoNuPlot = new TH1D("CC1pRecoNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CC1pRecoEvPlot = new TH1D("CC1pRecoEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
		TH1D* CC1pRecoNuScorePlot = new TH1D("CC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CC1pRecoFlashScorePlot = new TH1D("CC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);
		TH1D* CC1pRecoDistancePlot = new TH1D("CC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);
		TH1D* CC1pRecoLengthDifferencePlot = new TH1D("CC1pRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* CC1pRecodYZPlot = new TH1D("CC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CC1pRecoNPEPlot = new TH1D("CC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* CC1pRecoVertexActivityPlot = new TH1D("CC1pRecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);
//		TH1D* CC1pRecoNonZeroVertexActivityPlot = new TH1D("CC1pRecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);

		TH1D* CC1pRecoVertexXPlot = new TH1D("CC1pRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CC1pRecoVertexYPlot = new TH1D("CC1pRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CC1pRecoVertexZPlot = new TH1D("CC1pRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* CC1pRecoMuonLLRPIDPlot = new TH1D("CC1pRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CC1pRecoProtonLLRPIDPlot = new TH1D("CC1pRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CC1pRecoMuonLengthPlot = new TH1D("CC1pRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CC1pRecodMuonTracksScorePlot = new TH1D("CC1pRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CC1pRecodMuonVertexDistancePlot = new TH1D("CC1pRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);	

		TH1D* CC1pRecoContainedMuonLengthPlot = new TH1D("CC1pRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CC1pRecoUncontainedMuonLengthPlot = new TH1D("CC1pRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
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

//		TH1D* CC1pRecokMissPlot = new TH1D("CC1pRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* CC1pRecoPMissMinusPlot = new TH1D("CC1pRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* CC1pRecoPMissPlot = new TH1D("CC1pRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CC1pRecoEQEPlot = new TH1D("CC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CC1pRecoQ2Plot = new TH1D("CC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		TH1D* CC1pRecoCCQEECalPlot = new TH1D("CC1pRecoCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TH1D* CC1pRecoCCQEQ2Plot = new TH1D("CC1pRecoCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);
		
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

		TH2D* CC1pRecoMuonCosThetaPlot2D = new TH2D("CC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoProtonCosThetaPlot2D = new TH2D("CC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
			ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH2D* CC1pRecoMuonPhiPlot2D = new TH2D("CC1pRecoMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,
			ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH2D* CC1pRecoProtonPhiPlot2D = new TH2D("CC1pRecoProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,
			ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH2D* CC1pRecoDeltaPTPlot2D = new TH2D("CC1pRecoDeltaPTPlot2D",LabelXAxisDeltaPT2D,NBinsDeltaPT,
			ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH2D* CC1pRecoDeltaAlphaTPlot2D = new TH2D("CC1pRecoDeltaAlphaTPlot2D",LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,
			ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH2D* CC1pRecoDeltaPhiTPlot2D = new TH2D("CC1pRecoDeltaPhiTPlot2D",LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,
			ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH2D* CC1pRecoECalPlot2D = new TH2D("CC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
		TH2D* CC1pRecoEQEPlot2D = new TH2D("CC1pRecoEQEPlot2D",LabelXAxisEQE2D,NBinsEQE,ArrayNBinsEQE,NBinsEQE,ArrayNBinsEQE);
		TH2D* CC1pRecoQ2Plot2D = new TH2D("CC1pRecoQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

		// -------------------------------------------------------------------------------------------------------------------------------------

		// 2D Reco Level Plots for Signal CC1p, POT Scaled for response matrices

		TH2D* POTScaledCC1pRecoMuonMomentumPlot2D = new TH2D("POTScaledCC1pRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
			ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH2D* POTScaledCC1pRecoProtonMomentumPlot2D = new TH2D("POTScaledCC1pRecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,NBinsProtonMomentum,
			ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH2D* POTScaledCC1pRecoCCQEMuonMomentumPlot2D = new TH2D("POTScaledCC1pRecoCCQEMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,CCQENBinsMuonMomentum,
			CCQEArrayNBinsMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH2D* POTScaledCC1pRecoCCQEProtonMomentumPlot2D = new TH2D("POTScaledCC1pCCQECCQERecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,CCQENBinsProtonMomentum,
			CCQEArrayNBinsProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);

		TH2D* POTScaledCC1pRecoMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
			ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH2D* POTScaledCC1pRecoCCQEMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoCCQEMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,CCQENBinsMuonCosTheta,
			CCQEArrayNBinsMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoCCQEProtonCosThetaPlot2D = new TH2D("POTScaledCC1pCCQECCQERecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,CCQENBinsProtonCosTheta,
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
		TH1D* NonCC1pRecodYZPlot = new TH1D("NonCC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* NonCC1pRecoNPEPlot = new TH1D("NonCC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* NonCC1pRecoVertexActivityPlot = new TH1D("NonCC1pRecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);
//		TH1D* NonCC1pRecoNonZeroVertexActivityPlot = new TH1D("NonCC1pRecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);

		TH1D* NonCC1pRecoVertexXPlot = new TH1D("NonCC1pRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* NonCC1pRecoVertexYPlot = new TH1D("NonCC1pRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* NonCC1pRecoVertexZPlot = new TH1D("NonCC1pRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* NonCC1pRecoMuonLLRPIDPlot = new TH1D("NonCC1pRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* NonCC1pRecoProtonLLRPIDPlot = new TH1D("NonCC1pRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* NonCC1pRecoMuonLengthPlot = new TH1D("NonCC1pRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* NonCC1pRecodMuonTracksScorePlot = new TH1D("NonCC1pRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* NonCC1pRecodMuonVertexDistancePlot = new TH1D("NonCC1pRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* NonCC1pRecoContainedMuonLengthPlot = new TH1D("NonCC1pRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* NonCC1pRecoUncontainedMuonLengthPlot = new TH1D("NonCC1pRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* NonCC1pRecoMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
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

//		TH1D* NonCC1pRecokMissPlot = new TH1D("NonCC1pRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* NonCC1pRecoPMissMinusPlot = new TH1D("NonCC1pRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* NonCC1pRecoPMissPlot = new TH1D("NonCC1pRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

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
		TH1D* CCQERecodYZPlot = new TH1D("CCQERecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCQERecoNPEPlot = new TH1D("CCQERecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* CCQERecoVertexActivityPlot = new TH1D("CCQERecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);
//		TH1D* CCQERecoNonZeroVertexActivityPlot = new TH1D("CCQERecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);

		TH1D* CCQERecoVertexXPlot = new TH1D("CCQERecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCQERecoVertexYPlot = new TH1D("CCQERecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCQERecoVertexZPlot = new TH1D("CCQERecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* CCQERecoMuonLLRPIDPlot = new TH1D("CCQERecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CCQERecoProtonLLRPIDPlot = new TH1D("CCQERecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCQERecoMuonLengthPlot = new TH1D("CCQERecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCQERecodMuonTracksScorePlot = new TH1D("CCQERecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);	
		TH1D* CCQERecodMuonVertexDistancePlot = new TH1D("CCQERecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCQERecoContainedMuonLengthPlot = new TH1D("CCQERecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCQERecoUncontainedMuonLengthPlot = new TH1D("CCQERecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCQERecoMuonMomentumPlot = new TH1D("CCQERecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
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

//		TH1D* CCQERecokMissPlot = new TH1D("CCQERecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* CCQERecoPMissMinusPlot = new TH1D("CCQERecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* CCQERecoPMissPlot = new TH1D("CCQERecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

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
		TH1D* CCMECRecodYZPlot = new TH1D("CCMECRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCMECRecoNPEPlot = new TH1D("CCMECRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* CCMECRecoVertexActivityPlot = new TH1D("CCMECRecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);
//		TH1D* CCMECRecoNonZeroVertexActivityPlot = new TH1D("CCMECRecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);

		TH1D* CCMECRecoVertexXPlot = new TH1D("CCMECRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCMECRecoVertexYPlot = new TH1D("CCMECRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCMECRecoVertexZPlot = new TH1D("CCMECRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* CCMECRecoMuonLLRPIDPlot = new TH1D("CCMECRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CCMECRecoProtonLLRPIDPlot = new TH1D("CCMECRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCMECRecoMuonLengthPlot = new TH1D("CCMECRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCMECRecodMuonTracksScorePlot = new TH1D("CCMECRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CCMECRecodMuonVertexDistancePlot = new TH1D("CCMECRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCMECRecoContainedMuonLengthPlot = new TH1D("CCMECRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCMECRecoUncontainedMuonLengthPlot = new TH1D("CCMECRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCMECRecoMuonMomentumPlot = new TH1D("CCMECRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoMuonPhiPlot = new TH1D("CCMECRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

		TH1D* CCMECRecoCCQEMuonMomentumPlot = new TH1D("CCMECRecoCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TH1D* CCMECRecoCCQEMuonCosThetaPlot = new TH1D("CCMERecoCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoCCQEMuonPhiPlot = new TH1D("CCMERecoCCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);

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
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCMECRecoCCQEProtonPhiPlot = new TH1D("CCMECRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CCMECRecoDeltaThetaPlot = new TH1D("CCMECRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaForwardThetaPlot = new TH1D("CCMECRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaBackwardThetaPlot = new TH1D("CCMECRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaPhiPlot = new TH1D("CCMECRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

//		TH1D* CCMECRecokMissPlot = new TH1D("CCMECRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* CCMECRecoPMissMinusPlot = new TH1D("CCMECRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* CCMECRecoPMissPlot = new TH1D("CCMECRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

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
		TH1D* CCRESRecodYZPlot = new TH1D("CCRESRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCRESRecoNPEPlot = new TH1D("CCRESRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* CCRESRecoVertexActivityPlot = new TH1D("CCRESRecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);
//		TH1D* CCRESRecoNonZeroVertexActivityPlot = new TH1D("CCRESRecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);

		TH1D* CCRESRecoVertexXPlot = new TH1D("CCRESRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCRESRecoVertexYPlot = new TH1D("CCRESRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCRESRecoVertexZPlot = new TH1D("CCRESRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* CCRESRecoMuonLLRPIDPlot = new TH1D("CCRESRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CCRESRecoProtonLLRPIDPlot = new TH1D("CCRESRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCRESRecoMuonLengthPlot = new TH1D("CCRESRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCRESRecodMuonTracksScorePlot = new TH1D("CCRESRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CCRESRecodMuonVertexDistancePlot = new TH1D("CCRESRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCRESRecoContainedMuonLengthPlot = new TH1D("CCRESRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCRESRecoUncontainedMuonLengthPlot = new TH1D("CCRESRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCRESRecoMuonMomentumPlot = new TH1D("CCRESRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
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

//		TH1D* CCRESRecokMissPlot = new TH1D("CCRESRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* CCRESRecoPMissMinusPlot = new TH1D("CCRESRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* CCRESRecoPMissPlot = new TH1D("CCRESRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

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
		TH1D* CCDISRecodYZPlot = new TH1D("CCDISRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCDISRecoNPEPlot = new TH1D("CCDISRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);
//		TH1D* CCDISRecoVertexActivityPlot = new TH1D("CCDISRecoVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);
//		TH1D* CCDISRecoNonZeroVertexActivityPlot = new TH1D("CCDISRecoNonZeroVertexActivityPlot",RecoLabelXAxisVertexActivity,NBinsVertexActivity,MinVertexActivity,MaxVertexActivity);

		TH1D* CCDISRecoVertexXPlot = new TH1D("CCDISRecoVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
		TH1D* CCDISRecoVertexYPlot = new TH1D("CCDISRecoVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
		TH1D* CCDISRecoVertexZPlot = new TH1D("CCDISRecoVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

//		TH1D* CCDISRecoMuonLLRPIDPlot = new TH1D("CCDISRecoMuonLLRPIDPlot",RecoLabelXAxisMuonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CCDISRecoProtonLLRPIDPlot = new TH1D("CCDISRecoProtonLLRPIDPlot",RecoLabelXAxisProtonLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);

		TH1D* CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCDISRecoMuonLengthPlot = new TH1D("CCDISRecoMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
//		TH1D* CCDISRecodMuonTracksScorePlot = new TH1D("CCDISRecodMuonTracksScorePlot",RecoLabelXAxisMuonTrackScore,NBinsTrackScore,MinTrackScore,MaxTrackScore);
		TH1D* CCDISRecodMuonVertexDistancePlot = new TH1D("CCDISRecodMuonVertexDistancePlot",RecoLabelXAxisMuonVertexDistanceTrackScore,NBinsMuonVertexDistance,MinMuonVertexDistance,MaxMuonVertexDistance);

		TH1D* CCDISRecoContainedMuonLengthPlot = new TH1D("CCDISRecoContainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);
		TH1D* CCDISRecoUncontainedMuonLengthPlot = new TH1D("CCDISRecoUncontainedMuonLengthPlot",RecoLabelXAxisMuonLength,NBinsMuonLength,MinMuonLength,MaxMuonLength);

		TH1D* CCDISRecoMuonMomentumPlot = new TH1D("CCDISRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
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
		TH1D* CCDISRecoCCQEProtonCosThetaPlot = new TH1D("CCDISRecoSCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TH1D* CCDISRecoCCQEProtonPhiPlot = new TH1D("CCDISRecoCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);

		TH1D* CCDISRecoDeltaThetaPlot = new TH1D("CCDISRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaForwardThetaPlot = new TH1D("CCDISRecoDeltaForwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaBackwardThetaPlot = new TH1D("CCDISRecoDeltaBackwardThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaPhiPlot = new TH1D("CCDISRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

//		TH1D* CCDISRecokMissPlot = new TH1D("CCDISRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
//		TH1D* CCDISRecoPMissMinusPlot = new TH1D("CCDISRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
//		TH1D* CCDISRecoPMissPlot = new TH1D("CCDISRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

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

//		TH1D* RecoLLRPIDPlot = new TH1D("RecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* RecoChi2Plot = new TH1D("RecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
//		TH1D* RecoThreePlaneChi2Plot = new TH1D("RecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
//			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* RecoThreePlaneChi2LogLikelihoodPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
			NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Muon 1D Reco Level Plots

//		TH1D* MuonRecoLLRPIDPlot = new TH1D("MuonRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* MuonRecoChi2Plot = new TH1D("MuonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
//		TH1D* MuonRecoThreePlaneChi2Plot = new TH1D("MuonRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
//			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* MuonRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("MuonRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Proton 1D Reco Level Plots

//		TH1D* ProtonRecoLLRPIDPlot = new TH1D("ProtonRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* ProtonRecoChi2Plot = new TH1D("ProtonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
//		TH1D* ProtonRecoThreePlaneChi2Plot = new TH1D("ProtonRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
//			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* ProtonRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("ProtonRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Pion 1D Reco Level Plots

//		TH1D* PionRecoLLRPIDPlot = new TH1D("PionRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* PionRecoChi2Plot = new TH1D("PionRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
//		TH1D* PionRecoThreePlaneChi2Plot = new TH1D("PionRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
//			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* PionRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("PionRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Cosmic 1D Reco Level Plots

//		TH1D* CosmicRecoLLRPIDPlot = new TH1D("CosmicRecoLLRPIDPlot",RecoLabelXAxisLLRPID,NBinsLLRPID,MinLLRPID,MaxLLRPID);
//		TH1D* CosmicRecoChi2Plot = new TH1D("CosmicRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
//		TH1D* CosmicRecoThreePlaneChi2Plot = new TH1D("CosmicRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
//			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* CosmicRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("CosmicRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		// Optical & calorimetric info (only Reco & CC1pReco) for selection cuts

		TH2D* RecodYZvsNPE2D = new TH2D("RecodYZvsNPE2D",RecoLabelXAxisNPE+RecoLabelXAxisdYZ,NBinsNPE,MinNPE,MaxNPE,NBinsdYZ,MindYZ,MaxdYZ);
		TH2D* CC1pRecodYZvsNPE2D = new TH2D("CC1pRecodYZvsNPE2D",RecoLabelXAxisNPE+RecoLabelXAxisdYZ,NBinsNPE,MinNPE,MaxNPE,NBinsdYZ,MindYZ,MaxdYZ);

		TH2D* RecoChi2PvsChi2Mu2D = new TH2D("RecoChi2PvsChi2Mu2D",";#chi^{2}_{#mu};#chi^{2}_{p}",NBinsChi2,MinChi2,MaxChi2,NBinsChi2,MinChi2,MaxChi2);
		TH2D* CC1pRecoChi2PvsChi2Mu2D = new TH2D("CC1pRecoChi2PvsChi2Mu2D",";#chi^{2}_{#mu};#chi^{2}_{p}",NBinsChi2,MinChi2,MaxChi2,NBinsChi2,MinChi2,MaxChi2);

		TH2D* RecoLLPvsLLMu2D = new TH2D("RecoLLPvsLLMu2D",RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood+RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,
						 2*NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood,
						 2*NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH2D* CC1pRecoLLPvsLLMu2D = new TH2D("CC1pRecoLLPvsLLMu2D",RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood+RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,
						 2*NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood,
						 2*NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

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

//		TH2D* Playground_CC1pRecoMuonMomentumResoVsMuonMomentumPlot = new TH2D("Playground_CC1pRecoMuonMomentumResoVsMuonMomentumPlot",LabelMuonMomentumResolution+LabelXAxisMuonMomentum,NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV,200,0.1,2.);

//		TH2D* Playground_CC1pRecoMuonMomentumResoVsDeltaThetaPlot = new TH2D("Playground_CC1pRecoMuonMomentumResoVsDeltaThetaPlot",LabelMuonMomentumResolution+RecoLabelXAxisDeltaTheta,NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);

		// -------------------------------------------------------------------------------------------------------------------

		TH1D* Playground_CC1pRecoMuonMomentumPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoMuonMomentumPlot_FullyContainedMuon",LabelMuonMomentumResolution,NBinsCC1pSTVReso,MinCC1pSTVMuMom,MaxCC1pSTVMuMom);
		TH1D* Playground_CC1pRecoMuonMomentumPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoMuonMomentumPlot_ExitingShortMuon",LabelMuonMomentumResolution,NBinsCC1pSTVReso,MinCC1pSTVMuMom,MaxCC1pSTVMuMom);
		TH1D* Playground_CC1pRecoMuonMomentumPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoMuonMomentumPlot_ExitingMediumMuon",LabelMuonMomentumResolution,NBinsCC1pSTVReso,MinCC1pSTVMuMom,MaxCC1pSTVMuMom);
		TH1D* Playground_CC1pRecoMuonMomentumPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoMuonMomentumPlot_ExitingLongMuon",LabelMuonMomentumResolution,NBinsCC1pSTVReso,MinCC1pSTVMuMom,MaxCC1pSTVMuMom);

		TH1D* Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);

		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);

		TH1D* Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);

		// -------------------------------------------------------------------------------------------------------------------

		TH1D* Playground_CC1pRecoDeltaPTPlot_ShortProton = new TH1D("Playground_CC1pRecoDeltaPTPlot_ShortProton",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
		TH1D* Playground_CC1pRecoDeltaPTPlot_MediumProton = new TH1D("Playground_CC1pRecoDeltaPTPlot_MediumProton",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
		TH1D* Playground_CC1pRecoDeltaPTPlot_LongProton = new TH1D("Playground_CC1pRecoDeltaPTPlot_LongProton",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);

		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ShortProton = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ShortProton",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_MediumProton = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_MediumProton",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_LongProton = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_LongProton",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);

		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ShortProton = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ShortProton",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_MediumProton = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_MediumProton",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_LongProton = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_LongProton",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);

		// -------------------------------------------------------------------------------------------------------------------

		const int NSlices = 3;

		TH1D* Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon_Slice[NSlices];

		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon_Slice[NSlices];

		TH1D* Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon_Slice[NSlices];

		TH1D* Playground_CC1pRecoDeltaPTPlot_ShortProton_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPTPlot_MediumProton_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPTPlot_LongProton_Slice[NSlices];

		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ShortProton_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_MediumProton_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_LongProton_Slice[NSlices];

		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ShortProton_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_MediumProton_Slice[NSlices];
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_LongProton_Slice[NSlices];

		for (int WhichSlice = 0; WhichSlice < NSlices; WhichSlice++) {

			Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_FullyContainedMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
			Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingShortMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
			Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingMediumMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
			Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingLongMuon",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);

			Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_FullyContainedMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
			Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingShortMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
			Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingMediumMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
			Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingLongMuon",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);

			Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_FullyContainedMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
			Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingShortMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
			Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingMediumMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
			Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ExitingLongMuon",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);

			// ----------------------------------------------------------------------------------------------------------

			Playground_CC1pRecoDeltaPTPlot_ShortProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ShortProton",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
			Playground_CC1pRecoDeltaPTPlot_MediumProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_MediumProton",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);
			Playground_CC1pRecoDeltaPTPlot_LongProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_LongProton",LabelDeltaPTResolution,NBinsCC1pSTVReso,MinCC1pSTVPt,MaxCC1pSTVPt);

			Playground_CC1pRecoDeltaAlphaTPlot_ShortProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ShortProton",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
			Playground_CC1pRecoDeltaAlphaTPlot_MediumProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_MediumProton",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);
			Playground_CC1pRecoDeltaAlphaTPlot_LongProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_LongProton",LabelDeltaAlphaTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaAlphaT,MaxCC1pSTVDeltaAlphaT);

			Playground_CC1pRecoDeltaPhiTPlot_ShortProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_ShortProton",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
			Playground_CC1pRecoDeltaPhiTPlot_MediumProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_MediumProton",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);
			Playground_CC1pRecoDeltaPhiTPlot_LongProton_Slice[WhichSlice] = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_Slice_"+ToStringInt(WhichSlice+1)+"_LongProton",LabelDeltaPhiTResolution,NBinsCC1pSTVReso,MinCC1pSTVDeltaPhiT,MaxCC1pSTVDeltaPhiT);

		}		

		// --------------------------------------------------------------------------------------------------------------------------------

		TH1D* CC1pRecoMuonTrueMomentumLongitudinalRatio = new TH1D("CC1pRecoMuonTrueMomentumLongitudinalRatio",";P^{true}_{#mu,||}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pRecoProtonTrueMomentumLongitudinalRatio = new TH1D("CC1pRecoProtonTrueMomentumLongitudinalRatio",";P^{true}_{p,||}/P^{true}_{p}",25,0.,1.);

		TH1D* CC1pRecoMuonTrueMomentumTransverseRatio = new TH1D("CC1pRecoMuonTrueMomentumTransverseRatio",";P^{true}_{#mu,T}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pRecoProtonTrueMomentumTransverseRatio = new TH1D("CC1pRecoProtonTrueMomentumTransverseRatio",";P^{true}_{p,T}/P^{true}_{p}",25,0.,1.);

		TH1D* CC1pTrueMuonTrueMomentumLongitudinalRatio = new TH1D("CC1pTrueMuonTrueMomentumLongitudinalRatio",";P^{true}_{#mu,||}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pTrueProtonTrueMomentumLongitudinalRatio = new TH1D("CC1pTrueProtonTrueMomentumLongitudinalRatio",";P^{true}_{p,||}/P^{true}_{p}",25,0.,1.);

		TH1D* CC1pTrueMuonTrueMomentumTransverseRatio = new TH1D("CC1pTrueMuonTrueMomentumTransverseRatio",";P^{true}_{#mu,T}/P^{true}_{#mu}",25,0.,1.);		
		TH1D* CC1pTrueProtonTrueMomentumTransverseRatio = new TH1D("CC1pTrueProtonTrueMomentumTransverseRatio",";P^{true}_{p,T}/P^{true}_{p}",25,0.,1.);		
		
		// ------------------------------------------------------------------------------------------------------------------------------

		TH1D* RecoMuonCosThetaPlotLLPOnlyStudy[NBinsThreePlaneChi2LogLikelihood];
		TH1D* CC1pRecoMuonCosThetaPlotLLPOnlyStudy[NBinsThreePlaneChi2LogLikelihood];

		for (int WhichBin = 0; WhichBin < NBinsThreePlaneChi2LogLikelihood; WhichBin++) {

				TString PlotLLPOnlyThres = "RecoMuonCosThetaPlot_LLPThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLLPOnlyStudy[WhichBin] = new TH1D(PlotLLPOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLLPOnlyThres = "CC1pRecoMuonCosThetaPlot_LLPThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLLPOnlyStudy[WhichBin] = new TH1D(CC1pPlotLLPOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

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

		TH1D* RecoMuonCosThetaPlotLengthOnlyStudy[NBinsMuonLength];
		TH1D* CC1pRecoMuonCosThetaPlotLengthOnlyStudy[NBinsMuonLength];

		for (int WhichBin = 0; WhichBin < NBinsMuonLength; WhichBin++) {

				TString PlotLengthOnlyThres = "RecoMuonCosThetaPlot_LengthThres_"+TString(std::to_string(WhichBin));

				RecoMuonCosThetaPlotLengthOnlyStudy[WhichBin] = new TH1D(PlotLengthOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

				TString CC1pPlotLengthOnlyThres = "CC1pRecoMuonCosThetaPlot_LengthThres_"+TString(std::to_string(WhichBin));

				CC1pRecoMuonCosThetaPlotLengthOnlyStudy[WhichBin] = new TH1D(CC1pPlotLengthOnlyThres,LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

		}

		// --------------------------------------------------------------------------------------------------------------------------------

		Tools tools;

		// --------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		// Loop over the events

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

			// Containment Demand
			// Fully contained proton candidates
			// Either fully contained or semi-contained muon candidates

			if (CandidateMu_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_EndContainment->at(0) == 0) { continue; }

			ContainmentCounter++;

			// -------------------------------------------------------------------------------------------------------------------------

			weight = POTWeight;

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 
			
				if (Weight <= 0 || Weight > 30) { continue; } // bug fix weight 
				if (T2KWeight <= 0 || T2KWeight > 30) { continue; }	// T2K tune weight			
				weight = POTWeight * Weight * T2KWeight * ROOTinoWeight; 
				
			}
			
			// --------------------------------------------------------------------------------------------------------------------------------

			// Genie, flux & reinteraction weights for systematics

			if ( fUniverseIndex != -1 && (fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || fWhichSample == "Overlay9_Run3" 
			|| fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5") ) {

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

				if (fEventWeightLabel == "expskin_FluxUnisim") { weight = weight*expskin_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "horncurrent_FluxUnisim") { weight = weight*horncurrent_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "kminus_PrimaryHadronNormalization") 
					{ weight = weight*kminus_PrimaryHadronNormalization->at(fUniverseIndex); }
				if (fEventWeightLabel == "kplus_PrimaryHadronFeynmanScaling") 
					{ weight = weight*kplus_PrimaryHadronFeynmanScaling->at(fUniverseIndex); }
				if (fEventWeightLabel == "kzero_PrimaryHadronSanfordWang") 
					{ weight = weight*kzero_PrimaryHadronSanfordWang->at(fUniverseIndex); }
				if (fEventWeightLabel == "nucleoninexsec_FluxUnisim") { weight = weight*nucleoninexsec_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "nucleonqexsec_FluxUnisim") { weight = weight*nucleonqexsec_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "nucleontotxsec_FluxUnisim") { weight = weight*nucleontotxsec_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "piminus_PrimaryHadronSWCentralSplineVariation") 
					{ weight = weight*piminus_PrimaryHadronSWCentralSplineVariation->at(fUniverseIndex); }
				if (fEventWeightLabel == "pioninexsec_FluxUnisim") { weight = weight*pioninexsec_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "pionqexsec_FluxUnisim") { weight = weight*pionqexsec_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "piontotxsec_FluxUnisim") { weight = weight*piontotxsec_FluxUnisim->at(fUniverseIndex); }
				if (fEventWeightLabel == "piplus_PrimaryHadronSWCentralSplineVariation") 
					{ weight = weight*piplus_PrimaryHadronSWCentralSplineVariation->at(fUniverseIndex); }

				// Reinteraction weights

				if (fEventWeightLabel == "reinteractions_piminus_Geant4") 
					{ weight = weight*reinteractions_piminus_Geant4->at(fUniverseIndex); }
				if (fEventWeightLabel == "reinteractions_piplus_Geant4") 
					{ weight = weight*reinteractions_piplus_Geant4->at(fUniverseIndex); }
				if (fEventWeightLabel == "reinteractions_proton_Geant4") 
					{ weight = weight*reinteractions_proton_Geant4->at(fUniverseIndex); }
//				if (fEventWeightLabel == "xsr_scc_Fa3_SCC") { weight = weight*xsr_scc_Fa3_SCC->at(fUniverseIndex); }
//				if (fEventWeightLabel == "xsr_scc_Fv3_SCC") { weight = weight*xsr_scc_Fv3_SCC->at(fUniverseIndex); }			

			}			

			// -----------------------------------------------------------------------------------------------------------------------

			if ( fabs(weight) != weight) { continue; } // Securing against infinities

			// -----------------------------------------------------------------------------------------------------------------------------

			// Contained Reconstructed Vertex

			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));
			if ( !tools.inFVVector(RecoVertex) ) { continue; }
			ContainedVertexCounter++;

			// -----------------------------------------------------------------------------------------------------------------------------

			// Jun 21 2021: Executive decision, moving on to recalibrated products

			// -----------------------------------------------------------------------------------------------------------------------------

			// Muon info

			double reco_Pmu = CandidateMu_P_Range_Recalibrate->at(0);

			// Quality cut suggested by Anne Schukraft // May 7 2021
			// For contained tracks, we have two methods for the muon momentum reconstruction
			// MCS & range, we can require that the two are withing 25%
			// to reject misreconstructed tracks 
	
			if (CandidateMu_EndContainment->at(0) == 1) { 

				double Reso =  TMath::Abs(CandidateMu_P_MCS_Recalibrate->at(0) - CandidateMu_P_Range_Recalibrate->at(0) ) / CandidateMu_P_Range_Recalibrate->at(0) ; 
				if (Reso > MuRangeMCSAgreeValue) { continue; }

			}

			MuonQualityCutCounter++;

			// -------------------------------------------------------------------------------------------------------------------------

			// Ensure that we don't have flipped tracks 
			// by demanding that muon start - vertex distance < muon end - vertex distance,
			// that proton start - vertex distance < proton end - vertex distance
			// 

			if (CandidateMuStartVertexDistance->at(0) > CandidateMuEndVertexDistance->at(0)) { continue; }
			if (CandidatePStartVertexDistance->at(0) > CandidatePEndVertexDistance->at(0)) { continue; }

			TVector3 CandidateMuonStart(CandidateMu_StartX->at(0),CandidateMu_StartY->at(0),CandidateMu_StartZ->at(0));
			TVector3 CandidateMuonEnd(CandidateMu_EndX->at(0),CandidateMu_EndY->at(0),CandidateMu_EndZ->at(0));

			TVector3 CandidateProtonStart(CandidateP_StartX->at(0),CandidateP_StartY->at(0),CandidateP_StartZ->at(0));
			TVector3 CandidateProtonEnd(CandidateP_EndX->at(0),CandidateP_EndY->at(0),CandidateP_EndZ->at(0));

			double DistanceStartPoints = (CandidateMuonStart - CandidateProtonStart).Mag();
			double DistanceEndPoints = (CandidateMuonEnd - CandidateProtonEnd).Mag(); 

			if (DistanceStartPoints > DistanceEndPoints) { continue; }

			NoFlippedTrackCounter++;

			// -------------------------------------------------------------------------------------------------------------------------

			// Jun 15 2021
			// Good reconstruction purposes on at least 2 planes
			// Require that at least 2 planes return good quality trunc dEdx

			float MuTruncdEdxEntry_Plane0 = CandidateMu_Plane0_LastEDep->at(0);
			float MuTruncdEdxEntry_Plane1 = CandidateMu_Plane1_LastEDep->at(0);
			float MuTruncdEdxEntry_Plane2 = CandidateMu_Plane2_LastEDep->at(0);

			float PTruncdEdxEntry_Plane0 = CandidateP_Plane0_LastEDep->at(0);
			float PTruncdEdxEntry_Plane1 = CandidateP_Plane1_LastEDep->at(0);
			float PTruncdEdxEntry_Plane2 = CandidateP_Plane2_LastEDep->at(0);

			if (MuTruncdEdxEntry_Plane0 < MuMinSumHitsValue && (MuTruncdEdxEntry_Plane1 < MuMinSumHitsValue || MuTruncdEdxEntry_Plane2 < MuMinSumHitsValue) ) { continue; }
			if (MuTruncdEdxEntry_Plane1 < MuMinSumHitsValue && (MuTruncdEdxEntry_Plane2 < MuMinSumHitsValue || MuTruncdEdxEntry_Plane0 < MuMinSumHitsValue) ) { continue; }
			if (MuTruncdEdxEntry_Plane2 < MuMinSumHitsValue && (MuTruncdEdxEntry_Plane0 < MuMinSumHitsValue || MuTruncdEdxEntry_Plane1 < MuMinSumHitsValue) ) { continue; }

			if (PTruncdEdxEntry_Plane0 < PMinSumHitsValue && (PTruncdEdxEntry_Plane1 < PMinSumHitsValue || PTruncdEdxEntry_Plane2 < PMinSumHitsValue) ) { continue; }
			if (PTruncdEdxEntry_Plane1 < PMinSumHitsValue && (PTruncdEdxEntry_Plane2 < PMinSumHitsValue || PTruncdEdxEntry_Plane0 < PMinSumHitsValue) ) { continue; }
			if (PTruncdEdxEntry_Plane2 < PMinSumHitsValue && (PTruncdEdxEntry_Plane0 < PMinSumHitsValue || PTruncdEdxEntry_Plane1 < PMinSumHitsValue) ) { continue; }

			SumTruncdEdxTrackCounter++;

			// -------------------------------------------------------------------------------------------------------------------------

			double l_muCandidate = CandidateMu_Length->at(0);
			double l_pCandidate = CandidateP_Length->at(0);
			double LengthDifference = l_muCandidate - l_pCandidate;
//			double MuonTrackScore = CandidateMu_TrackScore->at(0);
//			double ProtonTrackScore = CandidateP_TrackScore->at(0);
			double distance = CandidateMuP_Distance->at(0);	

			// -------------------------------------------------------------------------------------------------------------------------

			if (CandidateMu_EndContainment->at(0) == 0) { 

				// Quality cut for exiting muons
				// MCS introduces major for short exiting tracks

				// Quality cut based on PT reslution studies
				// Demand that the exiting muons are at least 1m long

				if (l_muCandidate < MuMinMCSLengthThresValue) { continue; }
				reco_Pmu = CandidateMu_P_MCS_Recalibrate->at(0); 

			}

			LowMCSQualityCutCounter++;

			double reco_Pmu_cos_theta = CandidateMu_CosTheta_Recalibrate->at(0);
			double reco_Pmu_phi = CandidateMu_Phi_Recalibrate->at(0); // deg
			double reco_Emu = TMath::Sqrt( reco_Pmu*reco_Pmu + MuonMass_GeV*MuonMass_GeV );
			
			TVector3 TVector3CandidateMuon(-1,-1,-1);
			TVector3CandidateMuon.SetMag(reco_Pmu);
			TVector3CandidateMuon.SetTheta(TMath::ACos(reco_Pmu_cos_theta));
			TVector3CandidateMuon.SetPhi(reco_Pmu_phi * TMath::Pi() / 180.);	

			// --------------------------------------------------------------------------------------------------------------------------

			// Proton info

			double reco_Pp = CandidateP_P_Range_Recalibrate->at(0);
			double reco_Pp_cos_theta = CandidateP_CosTheta_Recalibrate->at(0);
			double reco_Pp_phi = CandidateP_Phi_Recalibrate->at(0); // deg 
			double reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );
			double reco_Tp = reco_Ep - ProtonMass_GeV;
			
			TVector3 TVector3CandidateProton(-1,-1,-1);
			TVector3CandidateProton.SetMag(reco_Pp);
			TVector3CandidateProton.SetTheta(TMath::ACos(reco_Pp_cos_theta));
			TVector3CandidateProton.SetPhi(reco_Pp_phi * TMath::Pi() / 180.);	

			// --------------------------------------------------------------------------------------------------------------------------

			// Calorimetry

			// double reco_Pmu_chi2 = CandidateMu_Chi2_YPlane->at(0);
			// double reco_Pp_chi2 = CandidateP_Chi2_YPlane->at(0);

//			double reco_Pmu_ThreePlanechi2 = CandidateMu_ThreePlaneChi2->at(0);
//			double reco_Pp_ThreePlanechi2 = CandidateP_ThreePlaneChi2->at(0);

			double reco_Pmu_ThreePlaneLogLikelihood = CandidateMu_ThreePlaneLogLikelihood->at(0);
			double reco_Pp_ThreePlaneLogLikelihood = CandidateP_ThreePlaneLogLikelihood->at(0);

			// double reco_mu_LLR_Score = CandidateMu_LLR_PID->at(0);
			// double reco_p_LLR_Score = CandidateP_LLR_PID->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			// Miss quantities

//			double kMiss = Reco_kMiss->at(0);
//			double PMissMinus = Reco_PMissMinus->at(0);
//			double MissMomentum = Reco_PMiss->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			// STV quantities
			
			double TransMissMomentum = Reco_Pt_Recalibrate->at(0);

			double DeltaAlphaT = Reco_DeltaAlphaT_Recalibrate->at(0);

			double DeltaPhiT = Reco_DeltaPhiT_Recalibrate->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction

			double ECal = Reco_ECal_Recalibrate->at(0);

			// QE Energy Reconstruction

			double EQE = Reco_EQE_Recalibrate->at(0);

			// Reconstructed Q2

			double reco_Q2 = Reco_Q2_Recalibrate->at(0);

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

			// Relative angles

			double DeltaThetaProtonMuon_Deg = Reco_DeltaTheta->at(0);
			double DeltaPhiProtonMuon_Deg = Reco_DeltaPhi->at(0);		

			// ------------------------------------------------------------------------------------------------------------------------

//			double VertexActivity = VertexActivity_50x100->at(0);

			// Optical info

			double NPE = BeamFlashes_TotalPE->at(0);

			TVector3 BeamFlash(0,BeamFlashes_YCenter->at(0),BeamFlashes_ZCenter->at(0));

			double dYZ = (BeamFlash - RecoVertex).Mag();
			double MuonVertexDistance = (CandidateMuonStart-RecoVertex).Mag();	
			double ProtonVertexDistance = (CandidateProtonStart-RecoVertex).Mag();

			// Ensure that the tracks are not flipped by applying a loose max distance cut of 5 cm between the vertex and the tracks

//			if (MuonVertexDistance > 5) { continue; }
//			if (ProtonVertexDistance > 5) { continue; }  	

			// -------------------------------------------------------------------------------------------------------------------------

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

				if (VectorCuts[i] == "_PID" && !(reco_Pp_ThreePlaneLogLikelihood > ProtonThreePlaneChi2LogLikelihoodCut) ) 
					{ PassedSelection = false; }
				if (VectorCuts[i] == "_PID" && (reco_Pp_ThreePlaneLogLikelihood > ProtonThreePlaneChi2LogLikelihoodCut) ) { PIDCounter ++; }

				if (VectorCuts[i] == "_NuScore" && !( NuScore > MinimumNuScore) )  { PassedSelection = false; }
				if (VectorCuts[i] == "_NuScore" && ( NuScore > MinimumNuScore) && (reco_Pp_ThreePlaneLogLikelihood > ProtonThreePlaneChi2LogLikelihoodCut) ) { NuScoreCounter ++; }


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

			if (TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { continue; }
			if (DeltaAlphaT > ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]) { continue; }
			if (DeltaPhiT > ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]) { continue; }

//			if (ECal > ArrayNBinsECal[NBinsECal]) { continue; }
//			if (EQE > ArrayNBinsEQE[NBinsEQE]) { continue; }
//			if (reco_Q2 > ArrayNBinsQ2[NBinsQ2]) { continue; }

			KinematicsCounter++;				

			// --------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;

//			double true_kMiss = -1;
//			double true_PMissMinus = -1;
//			double true_PMiss = -1;	
		
			double true_TransMissMomentum = -1;
			double true_DeltaAlphaT = -1;
			double true_DeltaPhiT = -1;
			double true_ECal = -1;
			double true_EQE = -1;
			double true_Q2 = -1;			
			double true_nu = -1;			
			
			if (
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

//				true_kMiss = True_kMiss->at(0);
//				true_PMissMinus = True_PMissMinus->at(0);
//				true_PMiss = True_PMiss->at(0);	
			
				true_TransMissMomentum = True_Pt->at(0);
				true_DeltaAlphaT = True_DeltaAlphaT->at(0);
				true_DeltaPhiT = True_DeltaPhiT->at(0);
				true_ECal = True_ECal->at(0);
				true_EQE = True_EQE->at(0);
				true_Q2 = True_Q2->at(0);			
				true_nu = True_Ev - True_CandidateMu_P->at(0);			
				
			}

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

			for (int WhichBin = 0; WhichBin < NBinsThreePlaneChi2LogLikelihood; WhichBin++) {

				double LocalThres = MinThreePlaneChi2LogLikelihood + LLPStep * WhichBin;

				if (reco_Pp_ThreePlaneLogLikelihood > LocalThres) {
					
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

			// ----------------------------------------------------------------------------------------------------------------------

			// 1D Length Only Study

			for (int WhichBin = 0; WhichBin < NBinsMuonLength; WhichBin++) {

				double LocalThres = MinMuonLength + MuonLengthStep * WhichBin;

				if (l_muCandidate > LocalThres) {
					
					RecoMuonCosThetaPlotLengthOnlyStudy[WhichBin]->Fill(reco_Pmu_cos_theta,weight);

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			// No weight to be applied in the multiplicity plots
			RecoPi0Plot->Fill(NumberPi0); 
			RecoNeutronPlot->Fill(NumberNeutrons); 

			RecoNuPlot->Fill(true_nu,weight);
			RecoEvPlot->Fill(True_Ev,weight);
			RecoNuScorePlot->Fill(NuScore,weight);
			RecoFlashScorePlot->Fill(FlashScore,weight);
			RecoDistancePlot->Fill(distance,weight);
			RecoLengthDifferencePlot->Fill(LengthDifference,weight);
			RecodYZPlot->Fill(dYZ,weight);
			RecoNPEPlot->Fill(NPE,weight);
//			RecoVertexActivityPlot->Fill(VertexActivity,weight);
//			if (VertexActivity > 0) { RecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }	

			RecoVertexXPlot->Fill(Vertex_X->at(0),weight);
			RecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
			RecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

			RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
			RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

			// RecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
			// RecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

			RecoMuonLengthPlot->Fill(l_muCandidate,weight);
//			RecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
			RecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

			if (CandidateMu_EndContainment->at(0) == 1) { RecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
			if (CandidateMu_EndContainment->at(0) == 0) { RecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

			RecoMuonMomentumPlot->Fill(reco_Pmu,weight);
			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
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

//			RecokMissPlot->Fill(kMiss,weight);
//			RecoPMissMinusPlot->Fill(PMissMinus,weight);
//			RecoPMissPlot->Fill(MissMomentum,weight);

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

			// 2D Plot for Default Chi2 vs 3-Plane Chi2

//			if (reco_Pmu_chi2 > 0 && reco_Pmu_ThreePlanechi2 > 0) 
//				{ RecoChi2vsThreePlaneChi2TPlot->Fill(reco_Pmu_chi2,reco_Pmu_ThreePlanechi2); }


			// 2D Optical & calorimetric info

			RecodYZvsNPE2D->Fill(NPE,dYZ,weight);
//			RecoChi2PvsChi2Mu2D->Fill(reco_Pmu_chi2,reco_Pp_chi2,weight);
			RecoLLPvsLLMu2D->Fill(reco_Pmu_ThreePlaneLogLikelihood,reco_Pp_ThreePlaneLogLikelihood,weight);
			RecoLengthPvsLengthMu2D->Fill(l_muCandidate,l_pCandidate,weight);

			// -------------------------------------------------------------------------------------------------------------------------

			// Chi2 PID Studies

			// RecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);
			// RecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);

			// RecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
			// RecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);

//			RecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight/2.);
//			RecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight/2.);

			RecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);
			RecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);

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

				if ( CC1p == 1 && NumberPi0 == 0 && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg 
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
				     && true_TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				     && true_DeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]

				) {

					// --------------------------------------------------------------------------------------------------
				
					CC1pEventsPassingSelectionCuts++;
					myRunTxtFile << endl << "CC1p0pi signal event" << endl;

					// ---------------------------------------------------------------------------------------------------------------------------
					// ---------------------------------------------------------------------------------------------------------------------------

					// 1D LLP only study

					for (int WhichBin = 0; WhichBin < NBinsThreePlaneChi2LogLikelihood; WhichBin++) {

						double LocalThres = MinThreePlaneChi2LogLikelihood + LLPStep * WhichBin;

						if (reco_Pp_ThreePlaneLogLikelihood > LocalThres) {
							
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
					CC1pTrueMuonPhiPlot->Fill(True_CandidateMu_Phi->at(0),weight);

					CC1pTrueProtonMomentumPlot->Fill(True_CandidateP_P->at(0),weight);
					CC1pTrueProtonCosThetaPlot->Fill(True_CandidateP_CosTheta->at(0),weight);
					CC1pTrueProtonPhiPlot->Fill(True_CandidateP_Phi->at(0),weight);

					CC1pTrueDeltaPTPlot->Fill(true_TransMissMomentum,weight);
					CC1pTrueDeltaAlphaTPlot->Fill(true_DeltaAlphaT,weight);
					CC1pTrueDeltaPhiTPlot->Fill(true_DeltaPhiT,weight);

					CC1pTrueECalPlot->Fill(true_ECal,weight);
					CC1pTrueEQEPlot->Fill(true_EQE,weight);
					CC1pTrueQ2Plot->Fill(true_Q2,weight);

					// --------------------------------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------------------

					// 1D Reco Plots for the selected CC1p events 

					CC1pRecoNuPlot->Fill(true_nu,weight);
					CC1pRecoEvPlot->Fill(True_Ev,weight);
					CC1pRecoNuScorePlot->Fill(NuScore,weight);
					CC1pRecoFlashScorePlot->Fill(FlashScore,weight);
					CC1pRecoDistancePlot->Fill(distance,weight);
					CC1pRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CC1pRecodYZPlot->Fill(dYZ,weight);
					CC1pRecoNPEPlot->Fill(NPE,weight);
//					CC1pRecoVertexActivityPlot->Fill(VertexActivity,weight);
//					if (VertexActivity > 0) { CC1pRecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }		

					CC1pRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CC1pRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CC1pRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

					CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					// CC1pRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					// CC1pRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

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
					CC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi,weight);
					CC1pRecoProtonPhiPlot->Fill(reco_Pp_phi,weight);

//					CC1pRecokMissPlot->Fill(kMiss,weight);
//					CC1pRecoPMissMinusPlot->Fill(PMissMinus,weight);
//					CC1pRecoPMissPlot->Fill(MissMomentum,weight);

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

					// --------------------------------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------------------

					// 2D Plots Kinematic Variables

					CC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu);
					CC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp);

					CC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta);

					CC1pRecoMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi);
					CC1pRecoProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi);

					POTScaledCC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu,weight);
					POTScaledCC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp,weight);

					POTScaledCC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
					POTScaledCC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta,weight);

					POTScaledCC1pRecoMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi,weight);
					POTScaledCC1pRecoProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi,weight);

					// -----------------------------------------------------------------------------------------------------

					// True Level STV

					CC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum);
					CC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT);
					CC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT);

					POTScaledCC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum,weight);
					POTScaledCC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);
					POTScaledCC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT,weight);

					// -----------------------------------------------------------------------------------------------------

					// True level energy reconstruction & Q2

					CC1pRecoECalPlot2D->Fill(true_ECal,ECal);
					CC1pRecoEQEPlot2D->Fill(true_EQE,EQE);
					CC1pRecoQ2Plot2D->Fill(true_Q2,reco_Q2);

					POTScaledCC1pRecoECalPlot2D->Fill(true_ECal,ECal,weight);
					POTScaledCC1pRecoEQEPlot2D->Fill(true_EQE,EQE,weight);
					POTScaledCC1pRecoQ2Plot2D->Fill(true_Q2,reco_Q2,weight);

					// -----------------------------------------------------------------------------------------------------

					// 2D Optical & calorimetric info

					CC1pRecodYZvsNPE2D->Fill(NPE,dYZ,weight);
//					CC1pRecoChi2PvsChi2Mu2D->Fill(reco_Pmu_chi2,reco_Pp_chi2,weight);
					CC1pRecoLLPvsLLMu2D->Fill(reco_Pmu_ThreePlaneLogLikelihood,reco_Pp_ThreePlaneLogLikelihood,weight);
					CC1pRecoLengthPvsLengthMu2D->Fill(l_muCandidate,l_pCandidate,weight);

					// -------------------------------------------------------------------------------------------------------------------------

					// Playground for CC1p STV resolution plots
					// Warning: they are called resolution plots but they are just reco - true differences

					double CC1pMuonMomentumReso = ( reco_Pmu - True_CandidateMu_P->at(0) );
					double CC1pDeltaPTReso = ( TransMissMomentum - true_TransMissMomentum );
					double CC1pDeltaAlphaTReso = ( DeltaAlphaT - true_DeltaAlphaT );
					double CC1pDeltaPhiTReso = ( DeltaPhiT - true_DeltaPhiT );

//					Playground_CC1pRecoMuonMomentumResoVsDeltaThetaPlot->Fill(CC1pMuonMomentumReso,DeltaThetaProtonMuon_Deg,weight);
//					Playground_CC1pRecoMuonMomentumResoVsMuonMomentumPlot->Fill(CC1pMuonMomentumReso,reco_Pmu,weight);

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

					int PTSlice = -1;
					if (TransMissMomentum < 0.3) { PTSlice = 0; }  
					if (TransMissMomentum > 0.3 && TransMissMomentum < 0.6) { PTSlice = 1; }  
					if (TransMissMomentum > 0.6) { PTSlice = 2; }

					int DeltaAlphaTSlice = -1;
					if (DeltaAlphaT < 50) { DeltaAlphaTSlice = 0; }  
					if (DeltaAlphaT > 50 && DeltaAlphaT < 130) { DeltaAlphaTSlice = 1; }  
					if (DeltaAlphaT > 130) { DeltaAlphaTSlice = 2; }

					int DeltaPhiTSlice = -1;
					if (DeltaPhiT < 30) { DeltaPhiTSlice = 0; }  
					if (DeltaPhiT > 30 && DeltaPhiT < 60) { DeltaPhiTSlice = 1; }  
					if (DeltaPhiT > 60) { DeltaPhiTSlice = 2; }

					// -----------------------------------------------------------------------------------------------------------------

					// Contained Muon  

					if (CandidateMu_EndContainment->at(0) == 1) {

						Playground_CC1pRecoMuonMomentumPlot_FullyContainedMuon->Fill(CC1pMuonMomentumReso,weight);
						Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon->Fill(CC1pDeltaPhiTReso,weight);

						Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

					} else {

						// Exiting Muon

						if ( CandidateMu_Length->at(0) < 100 ) { // cm

							Playground_CC1pRecoMuonMomentumPlot_ExitingShortMuon->Fill(CC1pMuonMomentumReso,weight);
							Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon->Fill(CC1pDeltaPhiTReso,weight);

							Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

						}

						if ( CandidateMu_Length->at(0) > 100 && CandidateMu_Length->at(0) < 200 ) { // cm

							Playground_CC1pRecoMuonMomentumPlot_ExitingMediumMuon->Fill(CC1pMuonMomentumReso,weight);
							Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon->Fill(CC1pDeltaPhiTReso,weight);

							Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

						}


						if ( CandidateMu_Length->at(0) > 200) { // cm

							Playground_CC1pRecoMuonMomentumPlot_ExitingLongMuon->Fill(CC1pMuonMomentumReso,weight);
							Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon->Fill(CC1pDeltaPhiTReso,weight);

							Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

						}


					}

					// -----------------------------------------------------------------------------------------------------------------

					// By definition the proton is fully contained

					if ( CandidateP_Length->at(0) < 10 ) { // cm

						Playground_CC1pRecoDeltaPTPlot_ShortProton->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_ShortProton->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_ShortProton->Fill(CC1pDeltaPhiTReso,weight);

						Playground_CC1pRecoDeltaPTPlot_ShortProton_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_ShortProton_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_ShortProton_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

					}

					if ( CandidateP_Length->at(0) > 10 && CandidateP_Length->at(0) < 40 ) { // cm

						Playground_CC1pRecoDeltaPTPlot_MediumProton->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_MediumProton->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_MediumProton->Fill(CC1pDeltaPhiTReso,weight);

						Playground_CC1pRecoDeltaPTPlot_MediumProton_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_MediumProton_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_MediumProton_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

					}

					if ( CandidateP_Length->at(0) > 40) { // cm

						Playground_CC1pRecoDeltaPTPlot_LongProton->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_LongProton->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_LongProton->Fill(CC1pDeltaPhiTReso,weight);

						Playground_CC1pRecoDeltaPTPlot_LongProton_Slice[PTSlice]->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_LongProton_Slice[DeltaAlphaTSlice]->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_LongProton_Slice[DeltaPhiTSlice]->Fill(CC1pDeltaPhiTReso,weight);

					}

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
					) 
						{ MisIndetifiedMuonAsPion++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "pi-p event" << endl; }
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
 
							if (True_CandidateMu_EndZ->at(0) == True_CandidateP_EndZ->at(0)) { BrokenMuonTrack++; }
							else { MisIndetifiedMuPToMuMu++; OtherMCBkg++; }
							ManualNonCC1pEventsPassingSelectionCuts++;
							myRunTxtFile << endl << "Broken muon track" << endl;
					}
					else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg) 
						{ MisIndetifiedMuPToPiPi++;  ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; myRunTxtFile << endl << "pi-pi event" << endl; }
					else if ( CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) 
						{ 
							if (True_CandidateMu_EndZ->at(0) == True_CandidateP_EndZ->at(0)) { BrokenProtonTrack++; }
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
					else if ( True_CandidateMu_StartX->at(0) != True_CandidateP_StartX->at(0) 
					|| True_CandidateMu_StartY->at(0) != True_CandidateP_StartY->at(0)
					|| True_CandidateMu_StartZ->at(0) != True_CandidateP_StartZ->at(0) ) 
						{ MultipleVertices++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "Multiple vertices event" << endl; }

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
					     && true_TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
					     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
					     && true_DeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT])

					) {

						OutCommonRange++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++;
						myRunTxtFile << endl << "Out of range event" << endl;
					}

					else { 

						ManualNonCC1pEventsPassingSelectionCuts++;
						std::cout << "new MC background topology";
						std::cout << " CandidateMu_MCParticle_Pdg = " << CandidateMu_MCParticle_Pdg->at(0);
						std::cout << " CandidateP_MCParticle_Pdg = " << CandidateP_MCParticle_Pdg->at(0) << endl; 
						std::cout << " NumberProtons = " << NumberProtons << endl; 
						std::cout << " NumberChargedPions = " << NumberChargedPions << endl; 
						std::cout << " NumberMuons = " << NumberMuons << endl; 
						std::cout << " NumberPi0 = " << NumberPi0 << endl; 
						std::cout << " NumberNeutrons = " << NumberNeutrons << endl; 

					}

					NonCC1pRecoNuPlot->Fill(true_nu,weight);
					NonCC1pRecoEvPlot->Fill(True_Ev,weight);
					NonCC1pRecoNuScorePlot->Fill(NuScore,weight);
					NonCC1pRecoFlashScorePlot->Fill(FlashScore,weight);
					NonCC1pRecoDistancePlot->Fill(distance,weight);
					NonCC1pRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					NonCC1pRecodYZPlot->Fill(dYZ,weight);
					NonCC1pRecoNPEPlot->Fill(NPE,weight);
//					NonCC1pRecoVertexActivityPlot->Fill(VertexActivity,weight);	
//					if (VertexActivity > 0) { NonCC1pRecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }		

					NonCC1pRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					NonCC1pRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					NonCC1pRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

					NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					// NonCC1pRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					// NonCC1pRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					NonCC1pRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					NonCC1pRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					NonCC1pRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { NonCC1pRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { NonCC1pRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					NonCC1pRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
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

//					NonCC1pRecokMissPlot->Fill(kMiss,weight);
//					NonCC1pRecoPMissMinusPlot->Fill(PMissMinus,weight);
//					NonCC1pRecoPMissPlot->Fill(MissMomentum,weight);

					NonCC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					NonCC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					NonCC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					NonCC1pRecoECalPlot->Fill(ECal,weight);
					NonCC1pRecoEQEPlot->Fill(EQE,weight);
					NonCC1pRecoQ2Plot->Fill(reco_Q2,weight);
					
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
					CCQERecodYZPlot->Fill(dYZ,weight);
					CCQERecoNPEPlot->Fill(NPE,weight);
//					CCQERecoVertexActivityPlot->Fill(VertexActivity,weight);
//					if (VertexActivity > 0) { CCQERecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }			

					CCQERecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCQERecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCQERecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

					CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					// CCQERecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					// CCQERecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CCQERecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CCQERecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CCQERecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCQERecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCQERecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CCQERecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
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

//					CCQERecokMissPlot->Fill(kMiss,weight);
//					CCQERecoPMissMinusPlot->Fill(PMissMinus,weight);
//					CCQERecoPMissPlot->Fill(MissMomentum,weight);

					CCQERecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCQERecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCQERecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCQERecoECalPlot->Fill(ECal,weight);
					CCQERecoEQEPlot->Fill(EQE,weight);
					CCQERecoQ2Plot->Fill(reco_Q2,weight);
					
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

				}

				// ------------------------------------------------------------------------------------------------------------------------

				// CCMEC

				if (genie_mode == 10) {

					CCMECRecoNuPlot->Fill(true_nu,weight);
					CCMECRecoEvPlot->Fill(True_Ev,weight);
					CCMECRecoNuScorePlot->Fill(NuScore,weight);
					CCMECRecoFlashScorePlot->Fill(FlashScore,weight);
					CCMECRecoDistancePlot->Fill(distance,weight);
					CCMECRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CCMECRecodYZPlot->Fill(dYZ,weight);
					CCMECRecoNPEPlot->Fill(NPE,weight);
//					CCMECRecoVertexActivityPlot->Fill(VertexActivity,weight);
//					if (VertexActivity > 0) { CCMECRecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }			

					CCMECRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCMECRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCMECRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

					CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					// CCMECRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					// CCMECRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

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

//					CCMECRecokMissPlot->Fill(kMiss,weight);
//					CCMECRecoPMissMinusPlot->Fill(PMissMinus,weight);
//					CCMECRecoPMissPlot->Fill(MissMomentum,weight);

					CCMECRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCMECRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCMECRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCMECRecoECalPlot->Fill(ECal,weight);
					CCMECRecoEQEPlot->Fill(EQE,weight);
					CCMECRecoQ2Plot->Fill(reco_Q2,weight);
					
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
					CCRESRecodYZPlot->Fill(dYZ,weight);
					CCRESRecoNPEPlot->Fill(NPE,weight);
//					CCRESRecoVertexActivityPlot->Fill(VertexActivity,weight);
//					if (VertexActivity > 0) { CCRESRecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }			

					CCRESRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCRESRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCRESRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

					CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					// CCRESRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					// CCRESRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

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

//					CCRESRecokMissPlot->Fill(kMiss,weight);
//					CCRESRecoPMissMinusPlot->Fill(PMissMinus,weight);
//					CCRESRecoPMissPlot->Fill(MissMomentum,weight);

					CCRESRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCRESRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCRESRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCRESRecoECalPlot->Fill(ECal,weight);
					CCRESRecoEQEPlot->Fill(EQE,weight);
					CCRESRecoQ2Plot->Fill(reco_Q2,weight);
					
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
					CCDISRecodYZPlot->Fill(dYZ,weight);
					CCDISRecoNPEPlot->Fill(NPE,weight);
//					CCDISRecoVertexActivityPlot->Fill(VertexActivity,weight);
//					if (VertexActivity > 0) { CCDISRecoNonZeroVertexActivityPlot->Fill(VertexActivity,weight); }			

					CCDISRecoVertexXPlot->Fill(Vertex_X->at(0),weight);
					CCDISRecoVertexYPlot->Fill(Vertex_Y->at(0),weight);
					CCDISRecoVertexZPlot->Fill(Vertex_Z->at(0),weight);

					CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					// CCDISRecoMuonLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
					// CCDISRecoProtonLLRPIDPlot->Fill(reco_p_LLR_Score,weight);

					CCDISRecoMuonLengthPlot->Fill(l_muCandidate,weight);
//					CCDISRecodMuonTracksScorePlot->Fill(MuonTrackScore,weight);
					CCDISRecodMuonVertexDistancePlot->Fill(MuonVertexDistance,weight);

					if (CandidateMu_EndContainment->at(0) == 1) { CCDISRecoContainedMuonLengthPlot->Fill(l_muCandidate,weight); }
					if (CandidateMu_EndContainment->at(0) == 0) { CCDISRecoUncontainedMuonLengthPlot->Fill(l_muCandidate,weight); }

					CCDISRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
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

//					CCDISRecokMissPlot->Fill(kMiss,weight);
//					CCDISRecoPMissMinusPlot->Fill(PMissMinus,weight);
//					CCDISRecoPMissPlot->Fill(MissMomentum,weight);

					CCDISRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCDISRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCDISRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCDISRecoECalPlot->Fill(ECal,weight);
					CCDISRecoEQEPlot->Fill(EQE,weight);
					CCDISRecoQ2Plot->Fill(reco_Q2,weight);
					
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

						// MuonRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
//						MuonRecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);		
//						MuonRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight/2.);
						MuonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);	

					}

					if (CandidateP_MCParticle_Pdg->at(0) == MuonPdg) {

						// MuonRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);
//						MuonRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);		
//						MuonRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight/2.);
						MuonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

					// ----------------------------------------------------------------------------------------------------------------

					if (CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg) {

						// ProtonRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight);
//						ProtonRecoChi2Plot->Fill(reco_Pmu_chi2,weight);		
//						ProtonRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
						ProtonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);	

					}

					if (CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) {

						// ProtonRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);
//						ProtonRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);		
//						ProtonRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight/2.);
						ProtonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

					// -----------------------------------------------------------------------------------------------------------------

					if (CandidateMu_MCParticle_Pdg->at(0) == AbsChargedPionPdg) {

						// PionRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
//						PionRecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);		
//						PionRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight/2.);
						PionRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);	

					}

					if (CandidateP_MCParticle_Pdg->at(0) == AbsChargedPionPdg) {

						// PionRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);
//						PionRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);		
//						PionRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight/2.);
						PionRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

				} // End of the overlay particle breakdown using the Backtracker

				else {

					// CosmicRecoLLRPIDPlot->Fill(reco_mu_LLR_Score,weight/2.);
					// CosmicRecoLLRPIDPlot->Fill(reco_p_LLR_Score,weight/2.);

					// CosmicRecoChi2Plot->Fill(reco_Pmu_chi2,weight/2.);
					// CosmicRecoChi2Plot->Fill(reco_Pp_chi2,weight/2.);

//					CosmicRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight/2.);
//					CosmicRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight/2.);

					CosmicRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);
					CosmicRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);

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
		TH1D* SWTriggerEventPlot = (TH1D*)(fFile->Get("SWTriggerEventPlot"));
		TH1D* OneNuMuPFParticleEventPlot = (TH1D*)(fFile->Get("OneNuMuPFParticleEventPlot"));
		TH1D* TwoDaughterEventPlot = (TH1D*)(fFile->Get("TwoDaughterEventPlot"));
		TH1D* TrackLikeDaughterEventPlot = (TH1D*)(fFile->Get("TrackLikeDaughterEventPlot"));
		TH1D* MatchedTrackLikeDaughterEventPlot = (TH1D*)(fFile->Get("MatchedTrackLikeDaughterEventPlot"));
		TH1D* NuFlashScoreEventPlot = (TH1D*)(fFile->Get("NuFlashScoreEventPlot"));
		TH1D* OneBeamFlashEventPlot = (TH1D*)(fFile->Get("OneBeamFlashEventPlot"));
		TH1D* OnePairEventPlot = (TH1D*)(fFile->Get("OnePairEventPlot"));
		TH1D* MomentumThresholdEventPlot = (TH1D*)(fFile->Get("MomentumThresholdEventPlot"));
		TH1D* StartPointContainmentEventPlot = (TH1D*)(fFile->Get("StartPointContainmentEventPlot"));

		file->cd();
		SamdefEventsPlot->Write();
		SWTriggerEventPlot->Write();
		OneNuMuPFParticleEventPlot->Write();
		TwoDaughterEventPlot->Write();
		TrackLikeDaughterEventPlot->Write();
		MatchedTrackLikeDaughterEventPlot->Write();	
		NuFlashScoreEventPlot->Write();
		OneBeamFlashEventPlot->Write();
		OnePairEventPlot->Write();
		MomentumThresholdEventPlot->Write();
		StartPointContainmentEventPlot->Write();

		// -------------------------------------------------------------------------------------------------------------------------

		// Storing the event loss of the event selection 

		TH1D* ProtonEndPointContainmentEventPlot = new TH1D("ProtonEndPointContainmentEventPlot",";Proton contained end point events",1,0,1);
		TH1D* NoFlippedTrackEventPlot = new TH1D("NoFlippedTrackEventPlot",";No flipped events",1,0,1);
		TH1D* SumTruncdEdxTrackEventPlot = new TH1D("SumTruncdEdxTrackEventPlot",";Min sum trunc dEdx on 3 planes events",1,0,1);
		TH1D* VertexContainmentEventPlot = new TH1D("VertexContainmentEventPlot",";Proton contained end point events",1,0,1);
		TH1D* MuonQualityEventPlot = new TH1D("MuonQualityEventPlot",";Muon quality events",1,0,1);
		TH1D* LowMCSQualityEventPlot = new TH1D("LowMCSQualityEventPlot",";MCS muon threshold events",1,0,1);
		TH1D* PidEventPlot = new TH1D("PidEventPlot",";PID events",1,0,1);
		TH1D* NuScoreEventPlot = new TH1D("NuScoreEventPlot",";#nu score events",1,0,1);
		TH1D* CommonEventPlot = new TH1D("CommonEventPlot",";Common events",1,0,1);

		ProtonEndPointContainmentEventPlot->SetBinContent(1,ContainmentCounter);
		NoFlippedTrackEventPlot->SetBinContent(1,NoFlippedTrackCounter);
		SumTruncdEdxTrackEventPlot->SetBinContent(1,SumTruncdEdxTrackCounter);
		VertexContainmentEventPlot->SetBinContent(1,ContainedVertexCounter);
		MuonQualityEventPlot->SetBinContent(1,MuonQualityCutCounter);
		LowMCSQualityEventPlot->SetBinContent(1,LowMCSQualityCutCounter);
		PidEventPlot->SetBinContent(1,PIDCounter);
		NuScoreEventPlot->SetBinContent(1,NuScoreCounter);
		CommonEventPlot->SetBinContent(1,KinematicsCounter);

		file->cd();
		ProtonEndPointContainmentEventPlot->Write();
		NoFlippedTrackEventPlot->Write();
		SumTruncdEdxTrackEventPlot->Write();
		VertexContainmentEventPlot->Write();
		MuonQualityEventPlot->Write();
		LowMCSQualityEventPlot->Write();
		PidEventPlot->Write();
		NuScoreEventPlot->Write();
		CommonEventPlot->Write();

		CC1pEventPlot->Write();
		CC2pEventPlot->Write();
		MultipleVerticesEventPlot->Write();
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

		POTScalePlot->SetBinContent(1,POTWeight);
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

		// -------------------------------------------------------------------------------------------------------------------------	

		cout << fWhichSample << "  " << NEventsPassingSelectionCuts << " events passing selection criteria" << endl << endl;

		// -------------------------------------------------------------------------------------------------------------------------	

		file->cd();
		file->Write();
		file->Close();
		fFile->Close();

		// -------------------------------------------------------------------------------------------------------------------------	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
