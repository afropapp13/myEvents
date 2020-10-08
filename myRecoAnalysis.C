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

	int NEventsPassingSelectionCuts = 0;
	int CC1pEventsPassingSelectionCuts = 0;
	int CC1p1piEventsPassingSelectionCuts = 0;
	int CC2pEventsPassingSelectionCuts = 0;
	int CC2p1piEventsPassingSelectionCuts = 0;
	int CC3pEventsPassingSelectionCuts = 0;
	int CC3p1piEventsPassingSelectionCuts = 0;
	int CC3p2piEventsPassingSelectionCuts = 0;
	int MisIndetifiedMuonAsPion = 0;
	int MisIndetifiedMuonAsProton = 0;
	int MisIndetifiedMuonAsAntiMuon = 0;
	int MisIndetifiedProtonAsPion = 0;
	int MisIndetifiedProtonAsDeuterium = 0;
	int MisIndetifiedProtonAsElectron = 0;
	int MisIndetifiedMuPToElectronElectron = 0;
	int MisIndetifiedMuPToMuMu = 0;
	int MisIndetifiedMuPToPiPi = 0;
	int CandidateMuon_MCParticle_OutFV = 0;
	int CandidateProton_MCParticle_OutFV = 0;
	int MultipleVertices = 0;
	int InTimeCosmics = 0;
		
	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");

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

//		TString FileName = "./OutputFiles/"+UBCodeVersion+"/"+Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root";
		TString FileName = "/uboone/data/users/apapadop/myEvents/OutputFiles/"+UBCodeVersion+"/"+\
				    Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root";
		TFile* file = new TFile(FileName,"recreate");
		std::cout << std::endl << "Creating a new file: " << FileName << std::endl << std::endl << std::endl;

		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots

		TH1D* RecoNuScorePlot = new TH1D("RecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* RecoFlashScorePlot = new TH1D("RecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* RecoDistancePlot = new TH1D("RecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* RecoLengthDifferencePlot = new TH1D("RecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* RecodYZPlot = new TH1D("RecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* RecoNPEPlot = new TH1D("RecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* RecoDeltaThetaPlot = new TH1D("RecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* RecoDeltaPhiPlot = new TH1D("RecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* RecoMuonPhiPlot = new TH1D("RecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* RecoProtonPhiPlot = new TH1D("RecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* RecokMissPlot = new TH1D("RecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* RecoPMissMinusPlot = new TH1D("RecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* RecoPMissPlot = new TH1D("RecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* RecoDeltaPTPlot = new TH1D("RecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* RecoDeltaAlphaTPlot = new TH1D("RecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* RecoDeltaPhiTPlot = new TH1D("RecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* RecoEQEPlot = new TH1D("RecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* RecoQ2Plot = new TH1D("RecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
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

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for Signal CC1p

		TH1D* CC1pRecoNuScorePlot = new TH1D("CC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CC1pRecoFlashScorePlot = new TH1D("CC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CC1pRecoDistancePlot = new TH1D("CC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CC1pRecoLengthDifferencePlot = new TH1D("CC1pRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* CC1pRecodYZPlot = new TH1D("CC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CC1pRecoNPEPlot = new TH1D("CC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CC1pRecoDeltaThetaPlot = new TH1D("CC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CC1pRecoDeltaPhiPlot = new TH1D("CC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CC1pRecoMuonPhiPlot = new TH1D("CC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CC1pRecoProtonPhiPlot = new TH1D("CC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CC1pRecokMissPlot = new TH1D("CC1pRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CC1pRecoPMissMinusPlot = new TH1D("CC1pRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CC1pRecoPMissPlot = new TH1D("CC1pRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CC1pRecoEQEPlot = new TH1D("CC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CC1pRecoQ2Plot = new TH1D("CC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CC1pRecoCosThetaMuPmuPlot = new TH2D("CC1pRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CC1pRecoCosThetaPPpPlot = new TH2D("CC1pRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// 2D Reco Level Plots for Signal CC1p

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

		// 1D Reco Level Plots for non-CC1p

		TH1D* NonCC1pRecoNuScorePlot = new TH1D("NonCC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* NonCC1pRecoFlashScorePlot = new TH1D("NonCC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,
			NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* NonCC1pRecoDistancePlot = new TH1D("NonCC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* NonCC1pRecoLengthDifferencePlot = new TH1D("NonCC1pRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* NonCC1pRecodYZPlot = new TH1D("NonCC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* NonCC1pRecoNPEPlot = new TH1D("NonCC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* NonCC1pRecoDeltaThetaPlot = new TH1D("NonCC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,
			NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* NonCC1pRecoDeltaPhiPlot = new TH1D("NonCC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* NonCC1pRecoMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,
			NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* NonCC1pRecoProtonMomentumPlot = new TH1D("NonCC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoProtonCosThetaPlot = new TH1D("NonCC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* NonCC1pRecoMuonPhiPlot = new TH1D("NonCC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* NonCC1pRecoProtonPhiPlot = new TH1D("NonCC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* NonCC1pRecokMissPlot = new TH1D("NonCC1pRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* NonCC1pRecoPMissMinusPlot = new TH1D("NonCC1pRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* NonCC1pRecoPMissPlot = new TH1D("NonCC1pRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* NonCC1pRecoDeltaPTPlot = new TH1D("NonCC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* NonCC1pRecoDeltaAlphaTPlot = new TH1D("NonCC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* NonCC1pRecoDeltaPhiTPlot = new TH1D("NonCC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* NonCC1pRecoECalPlot = new TH1D("NonCC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* NonCC1pRecoEQEPlot = new TH1D("NonCC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* NonCC1pRecoQ2Plot = new TH1D("NonCC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* NonCC1pRecoCosThetaMuPmuPlot = new TH2D("NonCC1pRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* NonCC1pRecoCosThetaPPpPlot = new TH2D("NonCC1pRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBinsProtonMomentum,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ---------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCQE

		TH1D* CCQERecoNuScorePlot = new TH1D("CCQERecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCQERecoFlashScorePlot = new TH1D("CCQERecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCQERecoDistancePlot = new TH1D("CCQERecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCQERecoLengthDifferencePlot = new TH1D("CCQERecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* CCQERecodYZPlot = new TH1D("CCQERecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCQERecoNPEPlot = new TH1D("CCQERecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCQERecoDeltaThetaPlot = new TH1D("CCQERecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCQERecoDeltaPhiPlot = new TH1D("CCQERecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCQERecoMuonMomentumPlot = new TH1D("CCQERecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCQERecoProtonMomentumPlot = new TH1D("CCQERecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCQERecoProtonCosThetaPlot = new TH1D("CCQERecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCQERecoMuonPhiPlot = new TH1D("CCQERecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCQERecoProtonPhiPlot = new TH1D("CCQERecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCQERecokMissPlot = new TH1D("CCQERecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCQERecoPMissMinusPlot = new TH1D("CCQERecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCQERecoPMissPlot = new TH1D("CCQERecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCQERecoDeltaPTPlot = new TH1D("CCQERecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQERecoDeltaAlphaTPlot = new TH1D("CCQERecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQERecoDeltaPhiTPlot = new TH1D("CCQERecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCQERecoECalPlot = new TH1D("CCQERecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCQERecoEQEPlot = new TH1D("CCQERecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCQERecoQ2Plot = new TH1D("CCQERecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCQERecoCosThetaMuPmuPlot = new TH2D("CCQERecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCQERecoCosThetaPPpPlot = new TH2D("CCQERecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ---------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCMEC

		TH1D* CCMECRecoNuScorePlot = new TH1D("CCMECRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCMECRecoFlashScorePlot = new TH1D("CCMECRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCMECRecoDistancePlot = new TH1D("CCMECRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCMECRecoLengthDifferencePlot = new TH1D("CCMECRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* CCMECRecodYZPlot = new TH1D("CCMECRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCMECRecoNPEPlot = new TH1D("CCMECRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCMECRecoDeltaThetaPlot = new TH1D("CCMECRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaPhiPlot = new TH1D("CCMECRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCMECRecoMuonMomentumPlot = new TH1D("CCMECRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCMECRecoProtonMomentumPlot = new TH1D("CCMECRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoProtonCosThetaPlot = new TH1D("CCMECRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCMECRecoMuonPhiPlot = new TH1D("CCMECRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCMECRecoProtonPhiPlot = new TH1D("CCMECRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCMECRecokMissPlot = new TH1D("CCMECRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCMECRecoPMissMinusPlot = new TH1D("CCMECRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCMECRecoPMissPlot = new TH1D("CCMECRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCMECRecoDeltaPTPlot = new TH1D("CCMECRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECRecoDeltaAlphaTPlot = new TH1D("CCMECRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECRecoDeltaPhiTPlot = new TH1D("CCMECRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCMECRecoECalPlot = new TH1D("CCMECRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECRecoEQEPlot = new TH1D("CCMECRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCMECRecoQ2Plot = new TH1D("CCMECRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCMECRecoCosThetaMuPmuPlot = new TH2D("CCMECRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCMECRecoCosThetaPPpPlot = new TH2D("CCMECRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCRES

		TH1D* CCRESRecoNuScorePlot = new TH1D("CCRESRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCRESRecoFlashScorePlot = new TH1D("CCRESRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCRESRecoDistancePlot = new TH1D("CCRESRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCRESRecoLengthDifferencePlot = new TH1D("CCRESRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* CCRESRecodYZPlot = new TH1D("CCRESRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCRESRecoNPEPlot = new TH1D("CCRESRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCRESRecoDeltaThetaPlot = new TH1D("CCRESRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCRESRecoDeltaPhiPlot = new TH1D("CCRESRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCRESRecoMuonMomentumPlot = new TH1D("CCRESRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCRESRecoProtonMomentumPlot = new TH1D("CCRESRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoProtonCosThetaPlot = new TH1D("CCRESRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCRESRecoMuonPhiPlot = new TH1D("CCRESRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCRESRecoProtonPhiPlot = new TH1D("CCRESRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCRESRecokMissPlot = new TH1D("CCRESRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCRESRecoPMissMinusPlot = new TH1D("CCRESRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCRESRecoPMissPlot = new TH1D("CCRESRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCRESRecoDeltaPTPlot = new TH1D("CCRESRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESRecoDeltaAlphaTPlot = new TH1D("CCRESRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESRecoDeltaPhiTPlot = new TH1D("CCRESRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCRESRecoECalPlot = new TH1D("CCRESRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESRecoEQEPlot = new TH1D("CCRESRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCRESRecoQ2Plot = new TH1D("CCRESRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
		// 2D Analysis
		
		TH2D* CCRESRecoCosThetaMuPmuPlot = new TH2D("CCRESRecoCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
		TH2D* CCRESRecoCosThetaPPpPlot = new TH2D("CCRESRecoCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);		

		// ---------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCDIS

		TH1D* CCDISRecoNuScorePlot = new TH1D("CCDISRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCDISRecoFlashScorePlot = new TH1D("CCDISRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCDISRecoDistancePlot = new TH1D("CCDISRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCDISRecoLengthDifferencePlot = new TH1D("CCDISRecoLengthDifferencePlot",RecoLabelXAxisLengthDifference,NBinsLengthDifference,MinLengthDifference,MaxLengthDifference);
		TH1D* CCDISRecodYZPlot = new TH1D("CCDISRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCDISRecoNPEPlot = new TH1D("CCDISRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCDISRecoDeltaThetaPlot = new TH1D("CCDISRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaPhiPlot = new TH1D("CCDISRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCDISRecoMuonMomentumPlot = new TH1D("CCDISRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCDISRecoProtonMomentumPlot = new TH1D("CCDISRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoProtonCosThetaPlot = new TH1D("CCDISRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCDISRecoMuonPhiPlot = new TH1D("CCDISRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCDISRecoProtonPhiPlot = new TH1D("CCDISRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCDISRecokMissPlot = new TH1D("CCDISRecokMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TH1D* CCDISRecoPMissMinusPlot = new TH1D("CCDISRecoPMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TH1D* CCDISRecoPMissPlot = new TH1D("CCDISRecoPMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		TH1D* CCDISRecoDeltaPTPlot = new TH1D("CCDISRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISRecoDeltaAlphaTPlot = new TH1D("CCDISRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISRecoDeltaPhiTPlot = new TH1D("CCDISRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCDISRecoECalPlot = new TH1D("CCDISRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISRecoEQEPlot = new TH1D("CCDISRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCDISRecoQ2Plot = new TH1D("CCDISRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		
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

		TH1D* RecoChi2Plot = new TH1D("RecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* RecoThreePlaneChi2Plot = new TH1D("RecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* RecoThreePlaneChi2LogLikelihoodPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
			NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Muon 1D Reco Level Plots

		TH1D* MuonRecoChi2Plot = new TH1D("MuonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* MuonRecoThreePlaneChi2Plot = new TH1D("MuonRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* MuonRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("MuonRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Proton 1D Reco Level Plots

		TH1D* ProtonRecoChi2Plot = new TH1D("ProtonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* ProtonRecoThreePlaneChi2Plot = new TH1D("ProtonRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* ProtonRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("ProtonRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Pion 1D Reco Level Plots

		TH1D* PionRecoChi2Plot = new TH1D("PionRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* PionRecoThreePlaneChi2Plot = new TH1D("PionRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* PionRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("PionRecoThreePlaneChi2LogLikelihoodPlot",
			RecoLabelXAxisThreePlaneChi2LogLikelihood,NBinsThreePlaneChi2LogLikelihood,
			MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Cosmic 1D Reco Level Plots

		TH1D* CosmicRecoChi2Plot = new TH1D("CosmicRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* CosmicRecoThreePlaneChi2Plot = new TH1D("CosmicRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,
			NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
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
		TH1D* NCC3p1piPlot = new TH1D("NCC3p1piPlot","",1,0,1);
		TH1D* NCC3p2piPlot = new TH1D("NCC3p2piPlot","",1,0,1);

		// --------------------------------------------------------------------------------------------------------------------------------

		// CC1p STV Resolution Playground plots

		int NBinsCC1pSTVReso = 100; double MinCC1pSTV = -50.,MaxCC1pSTV = 50.;

		TH1D* Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon",";#deltaP_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon",";#deltaP_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon",";#deltaP_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon",";#deltaP_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);

		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon",";#delta#alpha_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon",";#delta#alpha_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon",";#delta#alpha_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon",";#delta#alpha_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);

		TH1D* Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon",";#delta#phi_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon",";#delta#phi_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon",";#delta#phi_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);
		TH1D* Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon = new TH1D("Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon",";#delta#phi_{T} resolution (%)",NBinsCC1pSTVReso,MinCC1pSTV,MaxCC1pSTV);		

		// --------------------------------------------------------------------------------------------------------------------------------

		// POT Counting

		double POTCount = -99.;

		if (string(fWhichSample).find("Overlay") != std::string::npos) {

				TString PathToPOTFile = "/uboone/data/users/apapadop/myEvents/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";

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
		
		if (string(fWhichSample).find("ExtBNB9") != std::string::npos) { POTScale = E1DCNT_wcut / EXT; }

		if (string(fWhichSample).find("Overlay") != std::string::npos) { POTScale = tor860_wcut / POTCount; }			

		// --------------------------------------------------------------------------------------------------------------------------------

		// Loop over the events

		for (Long64_t jentry=0; jentry<nentries;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

			// ---------------------------------------------------------------------------------------------------------------------

			if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

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

			weight = POTScale;

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 
			
				if (Weight <= 0 || Weight > 30) { continue; }
				if (T2KWeight <= 0 || T2KWeight > 30) { continue; }				
				weight = POTScale * Weight * T2KWeight * ROOTinoWeight; 
				
			}
			
			// --------------------------------------------------------------------------------------------------------------------------------

			// Genie, flux & reinteraction weights for systematics

			if ( fUniverseIndex != -1 && (fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || fWhichSample == "Overlay9_Run3" 
			|| fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5") ) {

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

			// -----------------------------------------------------------------------------------------------------------------------

			// hard coded limit because something still looks weird in the dirt sample at low nu-score
			// COH interaction & infinite weight
			
			//if (NuScore < 0.04) { continue; }

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
			double LengthDifference = l_muCandidate - l_pCandidate;

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

			if (reco_Pmu_mcs > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
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

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

			// Relative angles

			double DeltaThetaProtonMuon_Deg = Reco_DeltaTheta->at(0);
			double DeltaPhiProtonMuon_Deg = Reco_DeltaPhi->at(0);

			// ------------------------------------------------------------------------------------------------------------------------
			
			// CCQElike analysis
			// Use the DeltaTheta, DeltaPhi & Pt cuts
				
			if (CCQElike) {
			
				if ( !(TMath::Abs(DeltaPhiProtonMuon_Deg - 180.) < 35.) ) { continue; }
				if ( !(TMath::Abs(DeltaThetaProtonMuon_Deg - 90.) < 55.) ) { continue; }			
				if ( !(TransMissMomentum < 0.35) ) { continue; }							
			
			}			

			// ------------------------------------------------------------------------------------------------------------------------

			// Optical info

			double NPE = BeamFlashes_TotalPE->at(0);

			TVector3 BeamFlash(0,BeamFlashes_YCenter->at(0),BeamFlashes_ZCenter->at(0));
			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));
			double dYZ = (BeamFlash - RecoVertex).Mag();

			double distance = CandidateMuP_Distance->at(0);	

			// -------------------------------------------------------------------------------------------------------------------------

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

//				if (VectorCuts[i] == "_NuScore" && !( NuScore > MinimumNuScore) )  { PassedSelection = false; }

////				if (VectorCuts[i] == "_Chi2" && 
////					!(CandidateMu_Chi2_YPlane->at(0) > CandidateP_Chi2_YPlane->at(0) && 
////					CandidateP_Chi2_YPlane->at(0) < ProtonChi2Cut && CandidateMu_Chi2_YPlane->at(0) > MuonChi2Cut) )  
////					{ PassedSelection = false; }

//				if (VectorCuts[i] == "_ThreePlaneLogChi2" && 
//					!(reco_Pp_ThreePlaneLogLikelihood > ProtonThreePlaneChi2LogLikelihoodCut) ) 
//					//&& reco_Pmu_ThreePlaneLogLikelihood < MuonThreePlaneChi2LogLikelihoodCut) )  
//					{ PassedSelection = false; }

//				//if (VectorCuts[i] == "_Collinearity" && !( DeltaThetaProtonMuon_Deg < DeltaThetaCut ) )  { PassedSelection = false; }


				// WATCH OUT !!!!!!!!!!!!!!!!!!!!!!!!!!! 
				// INVERTED CUTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				if (VectorCuts[i] == "_NuScore" && !( reco_Pp_ThreePlaneLogLikelihood > ProtonThreePlaneChi2LogLikelihoodCut ) )  { PassedSelection = false; }

				if (VectorCuts[i] == "_ThreePlaneLogChi2" && 
					!( NuScore > MinimumNuScore ) ) 
					{ PassedSelection = false; }


			}

			if (PassedSelection == false) { continue; }

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

			RecoNuScorePlot->Fill(NuScore,weight);
			RecoFlashScorePlot->Fill(FlashScore,weight);

			RecoDistancePlot->Fill(distance,weight);

			RecoLengthDifferencePlot->Fill(LengthDifference,weight);
			RecodYZPlot->Fill(dYZ,weight);
			RecoNPEPlot->Fill(NPE,weight);

			RecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			RecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
			RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

			RecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			RecoProtonMomentumPlot->Fill(reco_Pp,weight);

			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			RecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			RecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			RecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

			RecokMissPlot->Fill(kMiss,weight);
			RecoPMissMinusPlot->Fill(PMissMinus,weight);
			RecoPMissPlot->Fill(MissMomentum,weight);

			RecoDeltaPTPlot->Fill(TransMissMomentum,weight);
			RecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
			RecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

			RecoECalPlot->Fill(ECal,weight);
			RecoEQEPlot->Fill(EQE,weight);
			RecoQ2Plot->Fill(reco_Q2,weight);
			
			// 2D Analysis
		
			RecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
			RecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);

			// 2D Plot for Default Chi2 vs 3-Plane Chi2

			if (reco_Pmu_chi2 > 0 && reco_Pmu_ThreePlanechi2 > 0) 
				{ RecoChi2vsThreePlaneChi2TPlot->Fill(reco_Pmu_chi2,reco_Pmu_ThreePlanechi2); }


			// 2D Optical & calorimetric info

			RecodYZvsNPE2D->Fill(NPE,dYZ,weight);
			RecoChi2PvsChi2Mu2D->Fill(reco_Pmu_chi2,reco_Pp_chi2,weight);
			RecoLLPvsLLMu2D->Fill(reco_Pmu_ThreePlaneLogLikelihood,reco_Pp_ThreePlaneLogLikelihood,weight);
			RecoLengthPvsLengthMu2D->Fill(l_muCandidate,l_pCandidate,weight);


			// -------------------------------------------------------------------------------------------------------------------------

			// Playground for Reco plots




			// ---------------------------------------------------------------------------------------------------------------------------

			// CC1p1pi Background
			
			if (CC1p1pi == 1) { CC1p1piEventsPassingSelectionCuts++; }
			
			// ---------------------------------------------------------------------------------------------------------------------------

			// CC2p Background
			
			if (CC2p == 1) { CC2pEventsPassingSelectionCuts++; }

			// ---------------------------------------------------------------------------------------------------------------------------

			// CC2p1pi Background
			
			if (CC2p1pi == 1) { CC2p1piEventsPassingSelectionCuts++; }

			// ---------------------------------------------------------------------------------------------------------------------------

			// CC3p Background
			
			if (CC3p == 1) { CC3pEventsPassingSelectionCuts++; }

			// ---------------------------------------------------------------------------------------------------------------------------

			// CC3p1pi Background
			
			if (CC3p1pi == 1) { CC3p1piEventsPassingSelectionCuts++; }

			// ---------------------------------------------------------------------------------------------------------------------------

			// CC3p2pi Background
			
			if (CC3p2pi == 1) { CC3p2piEventsPassingSelectionCuts++; }					

			// ---------------------------------------------------------------------------------------------------------------------------

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

				// CC1p Signal

				if (CC1p == 1 && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg 
				    && True_CandidateMu_StartContainment->at(0) == 1 ) {
				
					CC1pEventsPassingSelectionCuts++;

					// 1D Plots

					CC1pRecoNuScorePlot->Fill(NuScore,weight);
					CC1pRecoFlashScorePlot->Fill(FlashScore,weight);

					CC1pRecoDistancePlot->Fill(distance,weight);

					CC1pRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CC1pRecodYZPlot->Fill(dYZ,weight);
					CC1pRecoNPEPlot->Fill(NPE,weight);

					CC1pRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					CC1pRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					CC1pRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
					CC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					CC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
					CC1pRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

					CC1pRecokMissPlot->Fill(kMiss,weight);
					CC1pRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CC1pRecoPMissPlot->Fill(MissMomentum,weight);

					CC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CC1pRecoECalPlot->Fill(ECal,weight);
					CC1pRecoEQEPlot->Fill(EQE,weight);
					CC1pRecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					CC1pRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
					CC1pRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);				

					// --------------------------------------------------------------------------------------------------
					// --------------------------------------------------------------------------------------------------

					// 2D Plots Kinematic Variables

					CC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu_mcs);
					CC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp);

					CC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta);

					CC1pRecoMuonPhiPlot2D->Fill(True_CandidateMu_Phi->at(0),reco_Pmu_phi*180./TMath::Pi());
					CC1pRecoProtonPhiPlot2D->Fill(True_CandidateP_Phi->at(0),reco_Pp_phi*180./TMath::Pi());

					// -----------------------------------------------------------------------------------------------------

					// True Level STV

					CC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum);
					CC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT);
					CC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT);

					// -----------------------------------------------------------------------------------------------------

					// True level energy reconstruction & Q2

					CC1pRecoECalPlot2D->Fill(true_ECal,ECal);
					CC1pRecoEQEPlot2D->Fill(true_EQE,EQE);
					CC1pRecoQ2Plot2D->Fill(true_Q2,reco_Q2);

					// -----------------------------------------------------------------------------------------------------

					// 2D Optical & calorimetric info

					CC1pRecodYZvsNPE2D->Fill(NPE,dYZ,weight);
					CC1pRecoChi2PvsChi2Mu2D->Fill(reco_Pmu_chi2,reco_Pp_chi2,weight);
					CC1pRecoLLPvsLLMu2D->Fill(reco_Pmu_ThreePlaneLogLikelihood,reco_Pp_ThreePlaneLogLikelihood,weight);
					CC1pRecoLengthPvsLengthMu2D->Fill(l_muCandidate,l_pCandidate,weight);


					// -------------------------------------------------------------------------------------------------------------------------

					// Playground for CC1p STV resolution plots

					double CC1pDeltaPTReso = ( TransMissMomentum - true_TransMissMomentum) / true_TransMissMomentum*100.;
					double CC1pDeltaAlphaTReso = ( DeltaAlphaT - true_DeltaAlphaT) / true_DeltaAlphaT*100.;
					double CC1pDeltaPhiTReso = ( DeltaPhiT - true_DeltaPhiT) / true_DeltaPhiT*100.;

					// Contained Muon

					if (CandidateMu_EndContainment->at(0) == 1) {

						Playground_CC1pRecoDeltaPTPlot_FullyContainedMuon->Fill(CC1pDeltaPTReso,weight);
						Playground_CC1pRecoDeltaAlphaTPlot_FullyContainedMuon->Fill(CC1pDeltaAlphaTReso,weight);
						Playground_CC1pRecoDeltaPhiTPlot_FullyContainedMuon->Fill(CC1pDeltaPhiTReso,weight);

					} else {

						if ( CandidateMu_Length->at(0) < 7 ) {

							Playground_CC1pRecoDeltaPTPlot_ExitingShortMuon->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingShortMuon->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingShortMuon->Fill(CC1pDeltaPhiTReso,weight);

						}

						if ( CandidateMu_Length->at(0) > 7 && CandidateMu_Length->at(0) < 9 ) {

							Playground_CC1pRecoDeltaPTPlot_ExitingMediumMuon->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingMediumMuon->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingMediumMuon->Fill(CC1pDeltaPhiTReso,weight);

						}


						if ( CandidateMu_Length->at(0) > 9) {

							Playground_CC1pRecoDeltaPTPlot_ExitingLongMuon->Fill(CC1pDeltaPTReso,weight);
							Playground_CC1pRecoDeltaAlphaTPlot_ExitingLongMuon->Fill(CC1pDeltaAlphaTReso,weight);
							Playground_CC1pRecoDeltaPhiTPlot_ExitingLongMuon->Fill(CC1pDeltaPhiTReso,weight);

						}


					}

					// -------------------------------------------------------------------------------------------------------------------------

				}

				// -------------------------------------------------------------------------------------------------------------

				// Non-CC1p 

	//			if (CC1p != 1) {
				else {

					if (CC1p1pi != 1 && CC2p != 1 && CC2p1pi != 1 && CC3p != 1 && CC3p1pi != 1 && CC3p2pi != 1) {

						if ( TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) 
							{ MisIndetifiedMuonAsPion++;  }
						else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) 
							{ MisIndetifiedProtonAsPion++;  }
						else if ( CandidateP_MCParticle_Pdg->at(0) == DeuteriumPdg && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) { MisIndetifiedProtonAsDeuterium++;  }
						else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ElectronPdg && CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) 				
							{ MisIndetifiedProtonAsElectron++;  }
						else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == ElectronPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == ElectronPdg) 
							{ MisIndetifiedMuPToElectronElectron++;  }
						else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == MuonPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == MuonPdg) 
							{ MisIndetifiedMuPToMuMu++;  }
						else if ( TMath::Abs(CandidateP_MCParticle_Pdg->at(0)) == AbsChargedPionPdg && TMath::Abs(CandidateMu_MCParticle_Pdg->at(0)) == AbsChargedPionPdg) 
							{ MisIndetifiedMuPToPiPi++;  }
						else if ( CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) { MisIndetifiedMuonAsProton++;  }
						else if ( CandidateMu_MCParticle_Pdg->at(0) == -MuonPdg && CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) { MisIndetifiedMuonAsAntiMuon++;  }
						else if ( True_CandidateMu_StartContainment->at(0) == 0 ) { CandidateMuon_MCParticle_OutFV++; }
//						else if ( True_CandidateP_StartContainment->at(0) == 0 || True_CandidateP_EndContainment->at(0) == 0) { CandidateProton_MCParticle_OutFV++; }
						else if ( CandidateMu_MCParticle_Pdg->at(0) == -99. || CandidateP_MCParticle_Pdg->at(0) == -99.) { InTimeCosmics++; }
						else if ( True_CandidateMu_StartX->at(0) != True_CandidateP_StartX->at(0) 
						|| True_CandidateMu_StartY->at(0) != True_CandidateP_StartY->at(0)
						|| True_CandidateMu_StartZ->at(0) != True_CandidateP_StartZ->at(0) ) { MultipleVertices++; }

					}

					NonCC1pRecoNuScorePlot->Fill(NuScore,weight);
					NonCC1pRecoFlashScorePlot->Fill(FlashScore,weight);

					NonCC1pRecoDistancePlot->Fill(distance,weight);

					NonCC1pRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					NonCC1pRecodYZPlot->Fill(dYZ,weight);
					NonCC1pRecoNPEPlot->Fill(NPE,weight);

					NonCC1pRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					NonCC1pRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					NonCC1pRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
					NonCC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					NonCC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
					NonCC1pRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

					NonCC1pRecokMissPlot->Fill(kMiss,weight);
					NonCC1pRecoPMissMinusPlot->Fill(PMissMinus,weight);
					NonCC1pRecoPMissPlot->Fill(MissMomentum,weight);

					NonCC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					NonCC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					NonCC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					NonCC1pRecoECalPlot->Fill(ECal,weight);
					NonCC1pRecoEQEPlot->Fill(EQE,weight);
					NonCC1pRecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					NonCC1pRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
					NonCC1pRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);				

				}

				// -------------------------------------------------------------------------------------------------------------------------
				// ------------------------------------------------------------------------------------------------------------------------

				// CCQE

				if (genie_mode == 0) {

					CCQERecoNuScorePlot->Fill(NuScore,weight);
					CCQERecoFlashScorePlot->Fill(FlashScore,weight);

					CCQERecoDistancePlot->Fill(distance,weight);

					CCQERecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CCQERecodYZPlot->Fill(dYZ,weight);
					CCQERecoNPEPlot->Fill(NPE,weight);

					CCQERecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					CCQERecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					CCQERecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
					CCQERecoProtonMomentumPlot->Fill(reco_Pp,weight);

					CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCQERecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCQERecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
					CCQERecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

					CCQERecokMissPlot->Fill(kMiss,weight);
					CCQERecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCQERecoPMissPlot->Fill(MissMomentum,weight);

					CCQERecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCQERecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCQERecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCQERecoECalPlot->Fill(ECal,weight);
					CCQERecoEQEPlot->Fill(EQE,weight);
					CCQERecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					CCQERecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
					CCQERecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);				

				}

				// ------------------------------------------------------------------------------------------------------------------------

				// CCMEC

				if (genie_mode == 10) {

					CCMECRecoNuScorePlot->Fill(NuScore,weight);
					CCMECRecoFlashScorePlot->Fill(FlashScore,weight);

					CCMECRecoDistancePlot->Fill(distance,weight);

					CCMECRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CCMECRecodYZPlot->Fill(dYZ,weight);
					CCMECRecoNPEPlot->Fill(NPE,weight);

					CCMECRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					CCMECRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					CCMECRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
					CCMECRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					CCMECRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCMECRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCMECRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
					CCMECRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

					CCMECRecokMissPlot->Fill(kMiss,weight);
					CCMECRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCMECRecoPMissPlot->Fill(MissMomentum,weight);

					CCMECRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCMECRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCMECRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCMECRecoECalPlot->Fill(ECal,weight);
					CCMECRecoEQEPlot->Fill(EQE,weight);
					CCMECRecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					CCMECRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
					CCMECRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);				

				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCRES

				if (genie_mode == 1) {

					CCRESRecoNuScorePlot->Fill(NuScore,weight);
					CCRESRecoFlashScorePlot->Fill(FlashScore,weight);

					CCRESRecoDistancePlot->Fill(distance,weight);

					CCRESRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CCRESRecodYZPlot->Fill(dYZ,weight);
					CCRESRecoNPEPlot->Fill(NPE,weight);

					CCRESRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					CCRESRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					CCRESRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
					CCRESRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					CCRESRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCRESRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCRESRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
					CCRESRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

					CCRESRecokMissPlot->Fill(kMiss,weight);
					CCRESRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCRESRecoPMissPlot->Fill(MissMomentum,weight);

					CCRESRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCRESRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCRESRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCRESRecoECalPlot->Fill(ECal,weight);
					CCRESRecoEQEPlot->Fill(EQE,weight);
					CCRESRecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					CCRESRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
					CCRESRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);				

				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCDIS

				if (genie_mode == 2) {

					CCDISRecoNuScorePlot->Fill(NuScore,weight);
					CCDISRecoFlashScorePlot->Fill(FlashScore,weight);

					CCDISRecoDistancePlot->Fill(distance,weight);

					CCDISRecoLengthDifferencePlot->Fill(LengthDifference,weight);
					CCDISRecodYZPlot->Fill(dYZ,weight);
					CCDISRecoNPEPlot->Fill(NPE,weight);

					CCDISRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
					CCDISRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

					CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight);
					CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight);

					CCDISRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
					CCDISRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCDISRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCDISRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
					CCDISRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

					CCDISRecokMissPlot->Fill(kMiss,weight);
					CCDISRecoPMissMinusPlot->Fill(PMissMinus,weight);
					CCDISRecoPMissPlot->Fill(MissMomentum,weight);

					CCDISRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCDISRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCDISRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

					CCDISRecoECalPlot->Fill(ECal,weight);
					CCDISRecoEQEPlot->Fill(EQE,weight);
					CCDISRecoQ2Plot->Fill(reco_Q2,weight);
					
					// 2D Analysis
				
					CCDISRecoCosThetaMuPmuPlot->Fill(reco_Pmu_cos_theta,reco_Pmu_mcs,weight);
					CCDISRecoCosThetaPPpPlot->Fill(reco_Pp_cos_theta,reco_Pp,weight);				

				}

				// -------------------------------------------------------------------------------------------------------------------------
				// -------------------------------------------------------------------------------------------------------------------------

				// Chi2 PID Studies

				RecoChi2Plot->Fill(reco_Pmu_chi2,weight);
				RecoChi2Plot->Fill(reco_Pp_chi2,weight);

				RecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
				RecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight);

				RecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);
				RecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);

				// --------------------------------------------------------------------------------------------------------------------------

				// Overlay particle breakdown using the Backtracker

				if (CandidateMu_MCParticle_Pdg->size() > 0 && CandidateP_MCParticle_Pdg->size() > 0 ) {

					if (CandidateMu_MCParticle_Pdg->at(0) == MuonPdg) {

						MuonRecoChi2Plot->Fill(reco_Pmu_chi2,weight);		
						MuonRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
						MuonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);	

					}

					if (CandidateP_MCParticle_Pdg->at(0) == MuonPdg) {

						MuonRecoChi2Plot->Fill(reco_Pp_chi2,weight);		
						MuonRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight);
						MuonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

					// ----------------------------------------------------------------------------------------------------------------

					if (CandidateMu_MCParticle_Pdg->at(0) == ProtonPdg) {

						ProtonRecoChi2Plot->Fill(reco_Pmu_chi2,weight);		
						ProtonRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
						ProtonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);	

					}

					if (CandidateP_MCParticle_Pdg->at(0) == ProtonPdg) {

						ProtonRecoChi2Plot->Fill(reco_Pp_chi2,weight);		
						ProtonRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight);
						ProtonRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

					// -----------------------------------------------------------------------------------------------------------------

					if (CandidateMu_MCParticle_Pdg->at(0) == AbsChargedPionPdg) {

						PionRecoChi2Plot->Fill(reco_Pmu_chi2,weight);		
						PionRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
						PionRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);	

					}

					if (CandidateP_MCParticle_Pdg->at(0) == AbsChargedPionPdg) {

						PionRecoChi2Plot->Fill(reco_Pp_chi2,weight);		
						PionRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight);
						PionRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);	

					}

				} // End of the overlay particle breakdown using the Backtracker

				else {

					CosmicRecoChi2Plot->Fill(reco_Pmu_chi2,weight);
					CosmicRecoChi2Plot->Fill(reco_Pp_chi2,weight);

					CosmicRecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
					CosmicRecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight);

					CosmicRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);
					CosmicRecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);

				}

			} // End of the Overlay case and the breakdown into CC1p/NonCC1p & QE,MEC,RES,DIS

			// -------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the events

		std::cout << std::endl << "Created a new file: " << FileName << std::endl << std::endl << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

		// -------------------------------------------------------------------------------------------------------------------------
		
		double nentriesError = sqrt(nentries);		
		
//		std::cout << std::endl << "Number of " << fWhichSample << " initial entries = " << nentries << " +/- " << nentriesError 
//		<< " (POT normalized: " << nentries*POTScale << " +/- " << nentriesError*POTScale << ")" << std::endl;		
		
		// -------------------------------------------------------------------------------------------------------------------------	

		// All reconstructed events passing the selection criteria

		double NEventsPassingSelectionCutsError = sqrt(NEventsPassingSelectionCuts);

		std::cout << std::endl << "Number of events passing our selection criteria = " << NEventsPassingSelectionCuts << " +/- " 
		<< NEventsPassingSelectionCutsError
		<< " (POT normalized: " << NEventsPassingSelectionCuts*POTScale << " +/- " 
		<< NEventsPassingSelectionCutsError*POTScale << ")" << std::endl;
		
		// -------------------------------------------------------------------------------------------------------------------------	

		double CC1pEventsPassingSelectionCutsError = 0;
		double CC1p1piEventsPassingSelectionCutsError = 0;
		double CC2pEventsPassingSelectionCutsError = 0;
		double CC2p1piEventsPassingSelectionCutsError = 0;
		double CC3pEventsPassingSelectionCutsError = 0;
		double CC3p1piEventsPassingSelectionCutsError = 0;
		double CC3p2piEventsPassingSelectionCutsError = 0;
		double MisIndetifiedMuonAsPionError = 0;
		double MisIndetifiedProtonAsPionError = 0;
		double MisIndetifiedProtonAsDeuteriumError = 0;
		double MisIndetifiedProtonAsElectronError = 0;
		double MisIndetifiedMuPToElectronElectronError = 0;
		double MisIndetifiedMuPToMuMuError = 0;
		double MisIndetifiedMuPToPiPiError = 0;
		double MisIndetifiedMuonAsProtonError = 0;
		double MisIndetifiedMuonAsAntiMuonError = 0;
		double CandidateMuon_MCParticle_OutFVError = 0;
		double CandidateProton_MCParticle_OutFVError = 0;
		double MultipleVerticesError = 0;
		double InTimeCosmicsError = 0;

		if ( (fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || fWhichSample == "Overlay9_Run3" || fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5")
			&& fUniverseIndex == -1) {

			// All reconstructed CC1p events passing the selection criteria

			CC1pEventsPassingSelectionCutsError = sqrt(CC1pEventsPassingSelectionCuts);

			std::cout << std::endl << "CC1p & " << CC1pEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC1pEventsPassingSelectionCutsError
			<< " & " << CC1pEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC1pEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;
			
			// -------------------------------------------------------------------------------------------------------------------------	

			// All reconstructed CC1p1pi events passing the selection criteria

			CC1p1piEventsPassingSelectionCutsError = sqrt(CC1p1piEventsPassingSelectionCuts);

			std::cout << std::endl << std::endl << "CC1p1$\\pi$ & " << CC1p1piEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC1p1piEventsPassingSelectionCutsError
			<< " & " << CC1p1piEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC1p1piEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;
			
			// -------------------------------------------------------------------------------------------------------------------------	

			// All reconstructed CC2p events passing the selection criteria

			CC2pEventsPassingSelectionCutsError = sqrt(CC2pEventsPassingSelectionCuts);

			std::cout << "CC2p & " << CC2pEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC2pEventsPassingSelectionCutsError
			<< " & " << CC2pEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC2pEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// All reconstructed CC2p1pi events passing the selection criteria

			CC2p1piEventsPassingSelectionCutsError = sqrt(CC2p1piEventsPassingSelectionCuts);

			std::cout << "CC2p1$\\pi$ & " << CC2p1piEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC2p1piEventsPassingSelectionCutsError
			<< " & " << CC2p1piEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC2p1piEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;						

			// -------------------------------------------------------------------------------------------------------------------------	

			// All reconstructed CC3p events passing the selection criteria

			CC3pEventsPassingSelectionCutsError = sqrt(CC3pEventsPassingSelectionCuts);

			std::cout << "CC3p & " << CC3pEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC3pEventsPassingSelectionCutsError
			<< " & " << CC3pEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC3pEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// All reconstructed CC3p1pi events passing the selection criteria

			CC3p1piEventsPassingSelectionCutsError = sqrt(CC3p1piEventsPassingSelectionCuts);

			std::cout << "CC3p1$\\pi$ & " << CC3p1piEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC3p1piEventsPassingSelectionCutsError
			<< " & " << CC3p1piEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC3p1piEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// All reconstructed CC3p2pi events passing the selection criteria

			CC3p2piEventsPassingSelectionCutsError = sqrt(CC3p2piEventsPassingSelectionCuts);

			std::cout << "CC3p2$\\pi$ & " << CC3p2piEventsPassingSelectionCuts << " $\\pm$ " 
			<< CC3p2piEventsPassingSelectionCutsError
			<< " & " << CC3p2piEventsPassingSelectionCuts*POTScale << " $\\pm$ " 
			<< CC3p2piEventsPassingSelectionCutsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified muon-pion events passing the selection criteria

			MisIndetifiedMuonAsPionError = sqrt(MisIndetifiedMuonAsPion);

			std::cout << "$\\pi$-p & " << MisIndetifiedMuonAsPion << " $\\pm$ " 
			<< MisIndetifiedMuonAsPionError
			<< " & " << MisIndetifiedMuonAsPion*POTScale << " $\\pm$ " 
			<< MisIndetifiedMuonAsPionError*POTScale << " \\tabularnewline \n\\hline" << std::endl;	

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified muon-proton events passing the selection criteria

			MisIndetifiedMuonAsProtonError = sqrt(MisIndetifiedMuonAsProton);

			std::cout << "p-p & " << MisIndetifiedMuonAsProton << " $\\pm$ " 
			<< MisIndetifiedMuonAsProtonError
			<< " & " << MisIndetifiedMuonAsProton*POTScale << " $\\pm$ " 
			<< MisIndetifiedMuonAsProtonError*POTScale << " \\tabularnewline \n\\hline" << std::endl;	

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified muon- anti-muon events passing the selection criteria

			MisIndetifiedMuonAsAntiMuonError = sqrt(MisIndetifiedMuonAsAntiMuon);

			std::cout << "$\\mu^{+}$ - p & " << MisIndetifiedMuonAsAntiMuon << " $\\pm$ " 
			<< MisIndetifiedMuonAsAntiMuonError
			<< " & " << MisIndetifiedMuonAsAntiMuon*POTScale << " $\\pm$ " 
			<< MisIndetifiedMuonAsAntiMuonError*POTScale <<" \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified muon-pion events passing the selection criteria

			MisIndetifiedProtonAsPionError = sqrt(MisIndetifiedProtonAsPion);

			std::cout << "$\\mu$-$\\pi$ & " << MisIndetifiedProtonAsPion << " $\\pm$ " 
			<< MisIndetifiedProtonAsPionError
			<< " & " << MisIndetifiedProtonAsPion*POTScale << " $\\pm$ " 
			<< MisIndetifiedProtonAsPionError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified muon-pion events passing the selection criteria

			MisIndetifiedProtonAsDeuteriumError = sqrt(MisIndetifiedProtonAsDeuterium);

			std::cout << "$\\mu$-D & " << MisIndetifiedProtonAsDeuterium << " $\\pm$ " 
			<< MisIndetifiedProtonAsDeuteriumError
			<< " & " << MisIndetifiedProtonAsDeuterium*POTScale << " $\\pm$ " 
			<< MisIndetifiedProtonAsDeuteriumError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified muon-pion events passing the selection criteria

			MisIndetifiedProtonAsElectronError = sqrt(MisIndetifiedProtonAsElectron);

			std::cout << "$\\mu$-e & " << MisIndetifiedProtonAsElectron << " $\\pm$ " 
			<< MisIndetifiedProtonAsElectronError
			<< " & " << MisIndetifiedProtonAsElectron*POTScale << " $\\pm$ " 
			<< MisIndetifiedProtonAsElectronError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified e-e events passing the selection criteria

			MisIndetifiedMuPToElectronElectronError = sqrt(MisIndetifiedMuPToElectronElectron);

			std::cout << "e-e & " << MisIndetifiedMuPToElectronElectron << " $\\pm$ " 
			<< MisIndetifiedMuPToElectronElectronError
			<< " & " << MisIndetifiedMuPToElectronElectron*POTScale << " $\\pm$ " 
			<< MisIndetifiedMuPToElectronElectronError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified mu-mu events passing the selection criteria

			MisIndetifiedMuPToMuMuError = sqrt(MisIndetifiedMuPToMuMu);

			std::cout << "$\\mu$-$\\mu$ & " << MisIndetifiedMuPToMuMu << " $\\pm$ " 
			<< MisIndetifiedMuPToMuMuError
			<< " & " << MisIndetifiedMuPToMuMu*POTScale << " $\\pm$ " 
			<< MisIndetifiedMuPToMuMuError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Mis-identified pi-pi events passing the selection criteria

			MisIndetifiedMuPToPiPiError = sqrt(MisIndetifiedMuPToPiPi);

			std::cout << "$\\pi$-$\\pi$ & " << MisIndetifiedMuPToPiPi << " $\\pm$ " 
			<< MisIndetifiedMuPToPiPiError
			<< " & " << MisIndetifiedMuPToPiPi*POTScale << " $\\pm$ " 
			<< MisIndetifiedMuPToPiPiError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

			// -------------------------------------------------------------------------------------------------------------------------	

			// Muon Candidate MC Particle outside FV

			CandidateMuon_MCParticle_OutFVError = sqrt(CandidateMuon_MCParticle_OutFV);

			std::cout << "True vertex outside FV & " << CandidateMuon_MCParticle_OutFV << " $\\pm$ " 
			<< CandidateMuon_MCParticle_OutFVError
			<< " & " << CandidateMuon_MCParticle_OutFV*POTScale << " $\\pm$ " 
			<< CandidateMuon_MCParticle_OutFVError*POTScale << " \\tabularnewline \n\\hline" << std::endl;	

			// -------------------------------------------------------------------------------------------------------------------------	

			// In-time cosmics

			InTimeCosmicsError = sqrt(InTimeCosmics);

			std::cout << "In-time cosmics & " << InTimeCosmics << " $\\pm$ " 
			<< InTimeCosmicsError
			<< " & " << InTimeCosmics*POTScale << " $\\pm$ " 
			<< InTimeCosmicsError*POTScale << " \\tabularnewline \n\\hline" << std::endl;	

			// -------------------------------------------------------------------------------------------------------------------------	

			// Multiple Vertices

			MultipleVerticesError = sqrt(MultipleVertices);

			std::cout << "Multiple vertices & " << MultipleVertices << " $\\pm$ " 
			<< MultipleVerticesError
			<< " & " << MultipleVertices*POTScale << " $\\pm$ " 
			<< MultipleVerticesError*POTScale << " \\tabularnewline \n\\hline" << std::endl;

		}

		cout << endl;

		// -------------------------------------------------------------------------------------------------------------------------

		// Keep track of the event counters and the POT normalization

		POTScalePlot->SetBinContent(1,POTScale);
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

		NCC3p1piPlot->SetBinContent(1,CC3p1piEventsPassingSelectionCuts);
		NCC3p1piPlot->SetBinError(1,CC3p1piEventsPassingSelectionCutsError);

		NCC3p2piPlot->SetBinContent(1,CC3p2piEventsPassingSelectionCuts);
		NCC3p2piPlot->SetBinError(1,CC3p2piEventsPassingSelectionCutsError);

		// -------------------------------------------------------------------------------------------------------------------------	

		file->cd();
		file->Write();
		file->Close();

		// -------------------------------------------------------------------------------------------------------------------------	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
