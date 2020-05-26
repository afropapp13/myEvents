#define t_cxx
#include "t.h"
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

using namespace std;

void t::Loop() {

	// ---------------------------------------------------------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
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

		// ----------------------------------------------------------------------------------------------------------------------------------------------------------

		if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
		TH1D::SetDefaultSumw2();
		double weight = 1.;

		// ----------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = "./OutputFiles/"+UBCodeVersion+"/"+Cuts+"/STVStudies_"+fWhichSample+Cuts+".root";
		TFile* file = new TFile(FileName,"recreate");
		std::cout << std::endl << "Creating a new file: " << FileName << std::endl << std::endl << std::endl;

		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots

		TH1D* RecoNuScorePlot = new TH1D("RecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* RecoFlashScorePlot = new TH1D("RecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* RecoDistancePlot = new TH1D("RecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* RecodYZPlot = new TH1D("RecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* RecoNPEPlot = new TH1D("RecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* RecoDeltaThetaPlot = new TH1D("RecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* RecoDeltaPhiPlot = new TH1D("RecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* RecoMuonPhiPlot = new TH1D("RecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* RecoProtonPhiPlot = new TH1D("RecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* RecoDeltaPTPlot = new TH1D("RecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* RecoDeltaAlphaTPlot = new TH1D("RecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* RecoDeltaPhiTPlot = new TH1D("RecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* RecoEQEPlot = new TH1D("RecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* RecoQ2Plot = new TH1D("RecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// 2D Plot for Default Chi2 vs 3-Plane Chi2

		TH2D* RecoChi2vsThreePlaneChi2TPlot = new TH2D("RecoChi2vsThreePlaneChi2TPlot",RecoLabelXAxisChi2+RecoLabelXAxisThreePlaneChi2,NBinsChi2,MinChi2,MaxChi2,
											       NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);

		// ------------------------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for Signal CC1p

		TH1D* CC1pRecoNuScorePlot = new TH1D("CC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CC1pRecoFlashScorePlot = new TH1D("CC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CC1pRecoDistancePlot = new TH1D("CC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CC1pRecodYZPlot = new TH1D("CC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CC1pRecoNPEPlot = new TH1D("CC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CC1pRecoDeltaThetaPlot = new TH1D("CC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CC1pRecoDeltaPhiPlot = new TH1D("CC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CC1pRecoMuonPhiPlot = new TH1D("CC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CC1pRecoProtonPhiPlot = new TH1D("CC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CC1pRecoEQEPlot = new TH1D("CC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CC1pRecoQ2Plot = new TH1D("CC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// 2D Reco Level Plots for Signal CC1p

		TH2D* CC1pRecoMuonMomentumPlot2D = new TH2D("CC1pRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH2D* CC1pRecoProtonMomentumPlot2D = new TH2D("CC1pRecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,NBinsProtonMomentum,
												ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH2D* CC1pRecoMuonCosThetaPlot2D = new TH2D("CC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoProtonCosThetaPlot2D = new TH2D("CC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
											ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH2D* CC1pRecoMuonPhiPlot2D = new TH2D("CC1pRecoMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH2D* CC1pRecoProtonPhiPlot2D = new TH2D("CC1pRecoProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH2D* CC1pRecoDeltaPTPlot2D = new TH2D("CC1pRecoDeltaPTPlot2D",LabelXAxisDeltaPT2D,NBinsDeltaPT,ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH2D* CC1pRecoDeltaAlphaTPlot2D = new TH2D("CC1pRecoDeltaAlphaTPlot2D",LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH2D* CC1pRecoDeltaPhiTPlot2D = new TH2D("CC1pRecoDeltaPhiTPlot2D",LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH2D* CC1pRecoECalPlot2D = new TH2D("CC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
		TH2D* CC1pRecoEQEPlot2D = new TH2D("CC1pRecoEQEPlot2D",LabelXAxisEQE2D,NBinsEQE,ArrayNBinsEQE,NBinsEQE,ArrayNBinsEQE);
		TH2D* CC1pRecoQ2Plot2D = new TH2D("CC1pRecoQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

		// -------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for non-CC1p

		TH1D* NonCC1pRecoNuScorePlot = new TH1D("NonCC1pRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* NonCC1pRecoFlashScorePlot = new TH1D("NonCC1pRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* NonCC1pRecoDistancePlot = new TH1D("NonCC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* NonCC1pRecodYZPlot = new TH1D("NonCC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* NonCC1pRecoNPEPlot = new TH1D("NonCC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* NonCC1pRecoDeltaThetaPlot = new TH1D("NonCC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* NonCC1pRecoDeltaPhiPlot = new TH1D("NonCC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("NonCC1pRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* NonCC1pRecoMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* NonCC1pRecoProtonMomentumPlot = new TH1D("NonCC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoProtonCosThetaPlot = new TH1D("NonCC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* NonCC1pRecoMuonPhiPlot = new TH1D("NonCC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* NonCC1pRecoProtonPhiPlot = new TH1D("NonCC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* NonCC1pRecoDeltaPTPlot = new TH1D("NonCC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* NonCC1pRecoDeltaAlphaTPlot = new TH1D("NonCC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* NonCC1pRecoDeltaPhiTPlot = new TH1D("NonCC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* NonCC1pRecoECalPlot = new TH1D("NonCC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* NonCC1pRecoEQEPlot = new TH1D("NonCC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* NonCC1pRecoQ2Plot = new TH1D("NonCC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------
		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCQE

		TH1D* CCQERecoNuScorePlot = new TH1D("CCQERecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCQERecoFlashScorePlot = new TH1D("CCQERecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCQERecoDistancePlot = new TH1D("CCQERecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCQERecodYZPlot = new TH1D("CCQERecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCQERecoNPEPlot = new TH1D("CCQERecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCQERecoDeltaThetaPlot = new TH1D("CCQERecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCQERecoDeltaPhiPlot = new TH1D("CCQERecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCQERecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCQERecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCQERecoMuonMomentumPlot = new TH1D("CCQERecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* CCQERecoProtonMomentumPlot = new TH1D("CCQERecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCQERecoProtonCosThetaPlot = new TH1D("CCQERecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCQERecoMuonPhiPlot = new TH1D("CCQERecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCQERecoProtonPhiPlot = new TH1D("CCQERecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCQERecoDeltaPTPlot = new TH1D("CCQERecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQERecoDeltaAlphaTPlot = new TH1D("CCQERecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQERecoDeltaPhiTPlot = new TH1D("CCQERecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCQERecoECalPlot = new TH1D("CCQERecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCQERecoEQEPlot = new TH1D("CCQERecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCQERecoQ2Plot = new TH1D("CCQERecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCMEC

		TH1D* CCMECRecoNuScorePlot = new TH1D("CCMECRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCMECRecoFlashScorePlot = new TH1D("CCMECRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCMECRecoDistancePlot = new TH1D("CCMECRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCMECRecodYZPlot = new TH1D("CCMECRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCMECRecoNPEPlot = new TH1D("CCMECRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCMECRecoDeltaThetaPlot = new TH1D("CCMECRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCMECRecoDeltaPhiPlot = new TH1D("CCMECRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCMECRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCMECRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCMECRecoMuonMomentumPlot = new TH1D("CCMECRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* CCMECRecoProtonMomentumPlot = new TH1D("CCMECRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoProtonCosThetaPlot = new TH1D("CCMECRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCMECRecoMuonPhiPlot = new TH1D("CCMECRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCMECRecoProtonPhiPlot = new TH1D("CCMECRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCMECRecoDeltaPTPlot = new TH1D("CCMECRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECRecoDeltaAlphaTPlot = new TH1D("CCMECRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECRecoDeltaPhiTPlot = new TH1D("CCMECRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCMECRecoECalPlot = new TH1D("CCMECRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECRecoEQEPlot = new TH1D("CCMECRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCMECRecoQ2Plot = new TH1D("CCMECRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// --------------------------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCRES

		TH1D* CCRESRecoNuScorePlot = new TH1D("CCRESRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCRESRecoFlashScorePlot = new TH1D("CCRESRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCRESRecoDistancePlot = new TH1D("CCRESRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCRESRecodYZPlot = new TH1D("CCRESRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCRESRecoNPEPlot = new TH1D("CCRESRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCRESRecoDeltaThetaPlot = new TH1D("CCRESRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCRESRecoDeltaPhiPlot = new TH1D("CCRESRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCRESRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCRESRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCRESRecoMuonMomentumPlot = new TH1D("CCRESRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* CCRESRecoProtonMomentumPlot = new TH1D("CCRESRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoProtonCosThetaPlot = new TH1D("CCRESRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCRESRecoMuonPhiPlot = new TH1D("CCRESRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCRESRecoProtonPhiPlot = new TH1D("CCRESRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCRESRecoDeltaPTPlot = new TH1D("CCRESRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESRecoDeltaAlphaTPlot = new TH1D("CCRESRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESRecoDeltaPhiTPlot = new TH1D("CCRESRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCRESRecoECalPlot = new TH1D("CCRESRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESRecoEQEPlot = new TH1D("CCRESRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCRESRecoQ2Plot = new TH1D("CCRESRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// ---------------------------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCDIS

		TH1D* CCDISRecoNuScorePlot = new TH1D("CCDISRecoNuScorePlot",RecoLabelXAxisNuScore,NBinsNuScore,MinNuScore,MaxNuScore);
		TH1D* CCDISRecoFlashScorePlot = new TH1D("CCDISRecoFlashScorePlot",RecoLabelXAxisFlashScore,NBinsFlashScore,MinFlashScore,MaxFlashScore);

		TH1D* CCDISRecoDistancePlot = new TH1D("CCDISRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

		TH1D* CCDISRecodYZPlot = new TH1D("CCDISRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
		TH1D* CCDISRecoNPEPlot = new TH1D("CCDISRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

		TH1D* CCDISRecoDeltaThetaPlot = new TH1D("CCDISRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
		TH1D* CCDISRecoDeltaPhiPlot = new TH1D("CCDISRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

		TH1D* CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot = new TH1D("CCDISRecoThreePlaneChi2LogLikelihoodCandidateMuonPlot",
			RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);
		TH1D* CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot = new TH1D("CCDISRecoThreePlaneChi2LogLikelihoodCandidateProtonPlot",
			RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood,NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		TH1D* CCDISRecoMuonMomentumPlot = new TH1D("CCDISRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* CCDISRecoProtonMomentumPlot = new TH1D("CCDISRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c

		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoProtonCosThetaPlot = new TH1D("CCDISRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCDISRecoMuonPhiPlot = new TH1D("CCDISRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TH1D* CCDISRecoProtonPhiPlot = new TH1D("CCDISRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

		TH1D* CCDISRecoDeltaPTPlot = new TH1D("CCDISRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISRecoDeltaAlphaTPlot = new TH1D("CCDISRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISRecoDeltaPhiTPlot = new TH1D("CCDISRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

		TH1D* CCDISRecoECalPlot = new TH1D("CCDISRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISRecoEQEPlot = new TH1D("CCDISRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TH1D* CCDISRecoQ2Plot = new TH1D("CCDISRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

		// --------------------------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------------------------

		// Chi2 PID Studies

		// 1D Reco Level Plots

		TH1D* RecoChi2Plot = new TH1D("RecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* RecoThreePlaneChi2Plot = new TH1D("RecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* RecoThreePlaneChi2LogLikelihoodPlot = new TH1D("RecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
													NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Muon 1D Reco Level Plots

		TH1D* MuonRecoChi2Plot = new TH1D("MuonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* MuonRecoThreePlaneChi2Plot = new TH1D("MuonRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* MuonRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("MuonRecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
													NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Proton 1D Reco Level Plots

		TH1D* ProtonRecoChi2Plot = new TH1D("ProtonRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* ProtonRecoThreePlaneChi2Plot = new TH1D("ProtonRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* ProtonRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("ProtonRecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
													NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Pion 1D Reco Level Plots

		TH1D* PionRecoChi2Plot = new TH1D("PionRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* PionRecoThreePlaneChi2Plot = new TH1D("PionRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* PionRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("PionRecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
													NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// Cosmic 1D Reco Level Plots

		TH1D* CosmicRecoChi2Plot = new TH1D("CosmicRecoChi2Plot",RecoLabelXAxisChi2,NBinsChi2,MinChi2,MaxChi2);
		TH1D* CosmicRecoThreePlaneChi2Plot = new TH1D("CosmicRecoThreePlaneChi2Plot",RecoLabelXAxisThreePlaneChi2,NBinsThreePlaneChi2,MinThreePlaneChi2,MaxThreePlaneChi2);
		TH1D* CosmicRecoThreePlaneChi2LogLikelihoodPlot = new TH1D("CosmicRecoThreePlaneChi2LogLikelihoodPlot",RecoLabelXAxisThreePlaneChi2LogLikelihood,
													NBinsThreePlaneChi2LogLikelihood,MinThreePlaneChi2LogLikelihood,MaxThreePlaneChi2LogLikelihood);

		// ---------------------------------------------------------------------------------------------------------------------------------------------------------

		// POT Counting

		double POTCount = -99.;

		if (string(fWhichSample).find("Overlay") != std::string::npos) {

				// Locally
				//TString PathToPOTFile = "mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
				// gpvm's
				TString PathToPOTFile = "../PreSelection/OutputFiles/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";

				TFile* POTFile = TFile::Open(PathToPOTFile,"readonly");
				TH1D* POTCountHist = (TH1D*)(POTFile->Get("POTCountHist"));
				POTCount = POTCountHist->GetBinContent(1);
				POTFile->Close();
		}

		// ---------------------------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the events

		for (Long64_t jentry=0; jentry<nentries;jentry++) {
	//	for (Long64_t jentry=0; jentry<2000;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

			// ----------------------------------------------------------------------------------------------------------------------------------------------------

			if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// Demand for exactly one muon/proton candidate

			if (CandidateMu_P->size() != 1) { continue; }
			if (CandidateP_P->size() != 1) { continue; }

			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// Containment Demand
			// Fully contained proton candidates
			// Either fully contained or semi-contained muon candidates

			if (CandidateMu_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_StartContainment->at(0) == 0) { continue; }
			if (CandidateP_EndContainment->at(0) == 0) { continue; }

			// ------------------------------------------------------------------------------------------------------------------------------------------------

			// POT Scaling


			double tor860_wcut = 1;
			double E1DCNT_wcut = 1.;
			double EXT = 1.;

			if (string(fWhichSample).find("Run1") != std::string::npos) {

				tor860_wcut = tor860_wcut_Run1;
				E1DCNT_wcut = E1DCNT_wcut_Run1;
				EXT = EXT_Run1;

			}

			// ------------------------------------------------------------------------------------------------------------------------------------------------

			if (string(fWhichSample).find("ExtBNB9") != std::string::npos) { weight = E1DCNT_wcut / EXT;}

			else if (string(fWhichSample).find("Overlay") != std::string::npos) { weight = ( tor860_wcut / POTCount) * Weight * T2KWeight; }

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// EventWeight weights

// Fix it !!!!

			// Genie
//			if (string(fWhichSample).find("Genie_All") != std::string::npos) { weight = weight * EventWeightValues[fWhichSample][fUniverse];}

			// Flux
//			if (string(fWhichSample).find("FluxUnisim") != std::string::npos || string(fWhichSample).find("Primary") != std::string::npos) 
//				{ weight = weight * EventWeightValues[fWhichSample][fUniverse];}

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			if ( fabs(weight) != weight || weight <= 0.) { continue; } // Securing against infinities & negative weights

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			if (NuScore < 0.04) { continue; } // hard coded limit because something still looks weird in the dirt sample at low nu-score // COH interaction & infinite weight

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Muon & proton kinetic variables

			double reco_Pmu_mcs = CandidateMu_P->at(0);
			double reco_Pmu_cos_theta = CandidateMu_CosTheta->at(0);
			double reco_Pmu_phi = CandidateMu_Phi->at(0) * TMath::Pi() / 180.;
			double reco_Emu = TMath::Sqrt( reco_Pmu_mcs*reco_Pmu_mcs + MuonMass_GeV*MuonMass_GeV );

			double reco_Pp = CandidateP_P->at(0);
			double reco_Pp_cos_theta = CandidateP_CosTheta->at(0);
			double reco_Pp_phi = CandidateP_Phi->at(0) * TMath::Pi() / 180.;
			double reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );
			double reco_Tp = reco_Ep - ProtonMass_GeV;

			// ---------------------------------------------------------------------------------------------------------------------------------------------------

			// Calorimetry

			double reco_Pmu_chi2 = CandidateMu_Chi2_YPlane->at(0);
			double reco_Pp_chi2 = CandidateP_Chi2_YPlane->at(0);

			double reco_Pmu_ThreePlanechi2 = CandidateMu_ThreePlaneChi2->at(0);
			double reco_Pp_ThreePlanechi2 = CandidateP_ThreePlaneChi2->at(0);

			double reco_Pmu_ThreePlaneLogLikelihood = log(CandidateMu_ThreePlaneLogLikelihood->at(0));
			double reco_Pp_ThreePlaneLogLikelihood = log(CandidateP_ThreePlaneLogLikelihood->at(0));

			// -----------------------------------------------------------------------------------------------------------------------------------------------------------

			// Transverse variables

			TVector3 TVector3CandidateMuon;
			TVector3CandidateMuon.SetMagThetaPhi(CandidateMu_P->at(0),TMath::ACos(CandidateMu_CosTheta->at(0)),CandidateMu_Phi->at(0)*TMath::Pi()/180.);

			TVector3 TVector3CandidateMuonTrans;
			TVector3CandidateMuonTrans.SetXYZ(TVector3CandidateMuon.X(),TVector3CandidateMuon.Y(),0.);

			TVector3 TVector3CandidateProton;
			TVector3CandidateProton.SetMagThetaPhi(CandidateP_P->at(0),TMath::ACos(CandidateP_CosTheta->at(0)),CandidateP_Phi->at(0)*TMath::Pi()/180.);

			TVector3 TVector3CandidateProtonTrans;
			TVector3CandidateProtonTrans.SetXYZ(TVector3CandidateProton.X(),TVector3CandidateProton.Y(),0.);

			TVector3 TransMissMomentumV3 = TVector3CandidateMuonTrans + TVector3CandidateProtonTrans;

			double TransMissMomentum = TransMissMomentumV3.Mag();

			double DeltaAlphaT = TMath::ACos( (- TVector3CandidateMuonTrans*TransMissMomentumV3) /
					( TVector3CandidateMuonTrans.Mag()*TransMissMomentum ) ) * 180./TMath::Pi();
			if (DeltaAlphaT >= 180.) { DeltaAlphaT -= 180.; }
			if (DeltaAlphaT < 0.) { DeltaAlphaT += 180.; }

			double DeltaPhiT = TMath::ACos( (- TVector3CandidateMuonTrans*TVector3CandidateProtonTrans) /
					( TVector3CandidateMuonTrans.Mag()*TVector3CandidateProtonTrans.Mag() ) ) * 180./TMath::Pi();
			if (DeltaPhiT >= 180.) { DeltaPhiT -= 180.; }
			if (DeltaPhiT < 0.) { DeltaPhiT += 180.; }

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction

			double ECal = reco_Emu + reco_Tp +BE;

			// QE Energy Reconstruction

			double EQE = ( ProtonMass_GeV * BE + ProtonMass_GeV * reco_Emu ) / ( ProtonMass_GeV - reco_Emu + reco_Pmu_mcs * reco_Pmu_cos_theta);

			// Reconstructed Q2

			TLorentzVector reco_nu(0.,0.,ECal,ECal);
			TLorentzVector TLorentzVector4CandidateMuon(TVector3CandidateMuon,reco_Emu);
			TLorentzVector reco_q = reco_nu - TLorentzVector4CandidateMuon;
			double reco_Q2 = - reco_q.Mag2();

			// -------------------------------------------------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Make sure that the same events fill the same plots

			if (reco_Pmu_mcs < ArrayNBinsMuonMomentum[0]) { continue; }
			if (reco_Pmu_cos_theta < ArrayNBinsMuonCosTheta[0]) { continue; }
			if (reco_Pmu_phi < ArrayNBinsMuonPhi[0]) { continue; }

			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }
			if (reco_Pp_cos_theta < ArrayNBinsProtonCosTheta[0]) { continue; }
			if (reco_Pp_phi < ArrayNBinsProtonPhi[0]) { continue; }

			if (TransMissMomentum < ArrayNBinsDeltaPT[0]) { continue; }
			if (DeltaAlphaT < ArrayNBinsDeltaAlphaT[0]) { continue; }
			if (DeltaPhiT < ArrayNBinsDeltaPhiT[0]) { continue; }

			if (ECal < ArrayNBinsECal[0]) { continue; }
			if (EQE < ArrayNBinsEQE[0]) { continue; }
			if (reco_Q2 < ArrayNBinsQ2[0]) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			if (reco_Pmu_mcs > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
			if (reco_Pmu_cos_theta > ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]) { continue; }
			if (reco_Pmu_phi > ArrayNBinsMuonPhi[NBinsMuonPhi]) { continue; }

			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }
			if (reco_Pp_cos_theta > ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]) { continue; }
			if (reco_Pp_phi > ArrayNBinsProtonPhi[NBinsProtonPhi]) { continue; }

			if (TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { continue; }
			if (DeltaAlphaT > ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]) { continue; }
			if (DeltaPhiT > ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]) { continue; }

			if (ECal > ArrayNBinsECal[NBinsECal]) { continue; }
			if (EQE > ArrayNBinsEQE[NBinsEQE]) { continue; }
			if (reco_Q2 > ArrayNBinsQ2[NBinsQ2]) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Relative angles

			double DeltaThetaProtonMuon = TVector3CandidateMuon.Angle(TVector3CandidateProton);
			double DeltaThetaProtonMuon_Deg = DeltaThetaProtonMuon * 180. / TMath::Pi();
			if (DeltaThetaProtonMuon_Deg >= 180.) { DeltaThetaProtonMuon_Deg -= 180.; }
			if (DeltaThetaProtonMuon_Deg < 0.) { DeltaThetaProtonMuon_Deg += 180.; }
			double DeltaPhiProtonMuon = TVector3CandidateMuon.DeltaPhi(TVector3CandidateProton);
			double DeltaPhiProtonMuon_Deg = DeltaPhiProtonMuon * 180. / TMath::Pi();
			if (DeltaPhiProtonMuon_Deg >= 360.) { DeltaPhiProtonMuon_Deg -= 360.; }
			if (DeltaPhiProtonMuon_Deg < 0.) { DeltaPhiProtonMuon_Deg += 360.; }

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Optical info

			double NPE = BeamFlashes_TotalPE->at(0);

			TVector3 BeamFlash(0,BeamFlashes_YCenter->at(0),BeamFlashes_ZCenter->at(0));
			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));
			double dYZ = (BeamFlash - RecoVertex).Mag();

			double distance = CandidateMuP_Distance->at(0);	

			// ------------------------------------------------------------------------------------------------------------------------------------------

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

				if (VectorCuts[i] == "_NuScore" && !( NuScore > MinimumNuScore) )  { PassedSelection = false; }

				if (VectorCuts[i] == "_Chi2" && 
					!(CandidateMu_Chi2_YPlane->at(0) > CandidateP_Chi2_YPlane->at(0) && 
					CandidateP_Chi2_YPlane->at(0) < ProtonChi2Cut && CandidateMu_Chi2_YPlane->at(0) > MuonChi2Cut) )  
					{ PassedSelection = false; }

				if (VectorCuts[i] == "_ThreePlaneLogChi2" && 
					!(reco_Pp_ThreePlaneLogLikelihood > ProtonThreePlaneChi2LogLikelihoodCut) ) 
					//&& reco_Pmu_ThreePlaneLogLikelihood < MuonThreePlaneChi2LogLikelihoodCut) )  
					{ PassedSelection = false; }

				if (VectorCuts[i] == "_Collinearity" && !( DeltaThetaProtonMuon_Deg < DeltaThetaCut ) )  { PassedSelection = false; }


			}

			if (PassedSelection == false) { continue; }

			NEventsPassingSelectionCuts++;

			// -------------------------------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;
			if (
	//			fWhichSample == "Overlay9"
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { genie_mode = MCParticle_Mode; }

			// ---------------------------------------------------------------------------------------------------------------------------------------------------------

			int fCC1p = 0;

			double True_TransMissMomentum = -99., True_DeltaAlphaT = -99., True_DeltaPhiT = -99.;
			double true_ECal = -99., true_EQE = -99., true_Q2 = -99.;

			if (    
	//			fWhichSample == "Overlay9"			
				string(fWhichSample).find("Overlay") != std::string::npos
				&& CandidateMu_MCParticle_Pdg->at(0) == MuonPdg
				&& CandidateP_MCParticle_Pdg->at(0) == ProtonPdg

				&& True_CandidateMu_StartContainment->at(0) == 1
				&& True_CandidateP_StartContainment->at(0) == 1
				&& True_CandidateP_EndContainment->at(0) == 1

				&& True_CandidateMu_P->at(0) > ArrayNBinsMuonMomentum[0]
				&& True_CandidateP_P->at(0) > ArrayNBinsProtonMomentum[0]
				&& True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
				&& True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]

				&& True_CandidateMu_CosTheta->at(0) > ArrayNBinsMuonCosTheta[0]
				&& True_CandidateP_CosTheta->at(0) > ArrayNBinsProtonCosTheta[0]
				&& True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
				&& True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
				&& CC1p == 1 // Must be CC1p at truth level for it to be CC1p at reco level
			) { 

				// ----------------------------------------------------------------------------------------------------------------------------------

				// True Level STV

				TVector3 True_TVector3CandidateMuon;
				TVector3 True_TVector3CandidateMuonTrans;
				TVector3 True_TVector3CandidateProton;
				TVector3 True_TVector3CandidateProtonTrans;

				True_TVector3CandidateMuon.SetMagThetaPhi(True_CandidateMu_P->at(0),TMath::ACos(True_CandidateMu_CosTheta->at(0)),True_CandidateMu_Phi->at(0)*TMath::Pi()/180.);
				True_TVector3CandidateMuonTrans.SetXYZ(True_TVector3CandidateMuon.X(),True_TVector3CandidateMuon.Y(),0.);

				True_TVector3CandidateProton.SetMagThetaPhi(True_CandidateP_P->at(0),TMath::ACos(True_CandidateP_CosTheta->at(0)),True_CandidateP_Phi->at(0)*TMath::Pi()/180.);

				True_TVector3CandidateProtonTrans.SetXYZ(True_TVector3CandidateProton.X(),True_TVector3CandidateProton.Y(),0.);

				TVector3 True_TransMissMomentumV3 = True_TVector3CandidateMuonTrans + True_TVector3CandidateProtonTrans;

				True_TransMissMomentum = True_TransMissMomentumV3.Mag();

				True_DeltaAlphaT = TMath::ACos( (- True_TVector3CandidateMuonTrans*True_TransMissMomentumV3) /
						( True_TVector3CandidateMuonTrans.Mag()*True_TransMissMomentum ) ) * 180./TMath::Pi();

				True_DeltaPhiT = TMath::ACos( (- True_TVector3CandidateMuonTrans*True_TVector3CandidateProtonTrans) /
						( True_TVector3CandidateMuonTrans.Mag()*True_TVector3CandidateProtonTrans.Mag() ) ) * 180./TMath::Pi();

				// ----------------------------------------------------------------------------------------------------------------------------------

				double true_Emu = TMath::Sqrt( True_CandidateMu_P->at(0)*True_CandidateMu_P->at(0) + MuonMass_GeV*MuonMass_GeV );
				double true_Ep = TMath::Sqrt( True_CandidateP_P->at(0)*True_CandidateP_P->at(0) + ProtonMass_GeV*ProtonMass_GeV );
				double true_Tp = true_Ep - ProtonMass_GeV;

				true_ECal = true_Emu + true_Tp + BE;
				true_EQE = ( ProtonMass_GeV * BE + ProtonMass_GeV * true_Emu ) / ( ProtonMass_GeV - true_Emu + True_CandidateMu_P->at(0) * True_CandidateMu_CosTheta->at(0));

				TLorentzVector true_nu(0.,0.,true_ECal,true_ECal);
				TVector3 TVector3TrueCandidateMuon;
				TVector3TrueCandidateMuon.SetMagThetaPhi(True_CandidateMu_P->at(0),TMath::ACos(True_CandidateMu_CosTheta->at(0)),True_CandidateMu_Phi->at(0));

				TLorentzVector TLorentzVector4TrueCandidateMuon(TVector3TrueCandidateMuon,true_Emu);
				TLorentzVector true_q = true_nu - TLorentzVector4TrueCandidateMuon;
				true_Q2 = - true_q.Mag2();

				// ----------------------------------------------------------------------------------------------------------------------------------

				if ( 
				    True_TransMissMomentum > ArrayNBinsDeltaPT[0]
				 && True_TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				 && True_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]
				 && True_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				 && True_DeltaPhiT > ArrayNBinsDeltaPhiT[0]
				 && True_DeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
				 && true_ECal > ArrayNBinsECal[0]
				 && true_ECal < ArrayNBinsECal[NBinsECal]
				 && true_EQE > ArrayNBinsEQE[0]
				 && true_EQE < ArrayNBinsEQE[NBinsEQE]
				 && true_Q2 > ArrayNBinsQ2[0]
				 && true_Q2 < ArrayNBinsQ2[NBinsQ2]
				) {

					fCC1p = 1;

				}

			}

			// ----------------------------------------------------------------------------------------------------------------------

			RecoNuScorePlot->Fill(NuScore,weight);
			RecoFlashScorePlot->Fill(FlashScore,weight);

			RecoDistancePlot->Fill(distance,weight);

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

			RecoDeltaPTPlot->Fill(TransMissMomentum,weight);
			RecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
			RecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

			RecoECalPlot->Fill(ECal,weight);
			RecoEQEPlot->Fill(EQE,weight);
			RecoQ2Plot->Fill(reco_Q2,weight);

			// 2D Plot for Default Chi2 vs 3-Plane Chi2

			if (reco_Pmu_chi2 > 0 && reco_Pmu_ThreePlanechi2 > 0) 
				{ RecoChi2vsThreePlaneChi2TPlot->Fill(reco_Pmu_chi2,reco_Pmu_ThreePlanechi2); }

			// ---------------------------------------------------------------------------------------------------------------------------

			// CC1p Signal

			if (fCC1p == 1) {

				// 1D Plots

				CC1pRecoNuScorePlot->Fill(NuScore,weight);
				CC1pRecoFlashScorePlot->Fill(FlashScore,weight);

				CC1pRecoDistancePlot->Fill(distance,weight);

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

				CC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				CC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
				CC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

				CC1pRecoECalPlot->Fill(ECal,weight);
				CC1pRecoEQEPlot->Fill(EQE,weight);
				CC1pRecoQ2Plot->Fill(reco_Q2,weight);

				// ---------------------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------------------

				// 2D Plots Kinematic Variables

				CC1pRecoMuonMomentumPlot2D->Fill(reco_Pmu_mcs,True_CandidateMu_P->at(0));
				CC1pRecoProtonMomentumPlot2D->Fill(reco_Pp,True_CandidateP_P->at(0));

				CC1pRecoMuonCosThetaPlot2D->Fill(reco_Pmu_cos_theta,True_CandidateMu_CosTheta->at(0));
				CC1pRecoProtonCosThetaPlot2D->Fill(reco_Pp_cos_theta,True_CandidateP_CosTheta->at(0));

				CC1pRecoMuonPhiPlot2D->Fill(reco_Pmu_phi*180./TMath::Pi(),True_CandidateMu_Phi->at(0));
				CC1pRecoProtonPhiPlot2D->Fill(reco_Pp_phi*180./TMath::Pi(),True_CandidateP_Phi->at(0));

				// ---------------------------------------------------------------------------------------------------------------------

				// True Level STV

				CC1pRecoDeltaPTPlot2D->Fill(True_TransMissMomentum,TransMissMomentum);
				CC1pRecoDeltaAlphaTPlot2D->Fill(True_DeltaAlphaT,DeltaAlphaT);
				CC1pRecoDeltaPhiTPlot2D->Fill(True_DeltaPhiT,DeltaPhiT);

				// ---------------------------------------------------------------------------------------------------------------------

				// True level energy reconstruction & Q2

				CC1pRecoECalPlot2D->Fill(true_ECal,ECal);
				CC1pRecoEQEPlot2D->Fill(true_EQE,EQE);
				CC1pRecoQ2Plot2D->Fill(true_Q2,reco_Q2);

			}

			// ------------------------------------------------------------------------------------------------------------------------------------------

			// Non-CC1p 

			if (fCC1p != 1) {

				NonCC1pRecoNuScorePlot->Fill(NuScore,weight);
				NonCC1pRecoFlashScorePlot->Fill(FlashScore,weight);

				NonCC1pRecoDistancePlot->Fill(distance,weight);

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

				NonCC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				NonCC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
				NonCC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

				NonCC1pRecoECalPlot->Fill(ECal,weight);
				NonCC1pRecoEQEPlot->Fill(EQE,weight);
				NonCC1pRecoQ2Plot->Fill(reco_Q2,weight);

			}

			// ---------------------------------------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------------------------------------

			// CCQE

			if (genie_mode == 0) {

				CCQERecoNuScorePlot->Fill(NuScore,weight);
				CCQERecoFlashScorePlot->Fill(FlashScore,weight);

				CCQERecoDistancePlot->Fill(distance,weight);

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

				CCQERecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				CCQERecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
				CCQERecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

				CCQERecoECalPlot->Fill(ECal,weight);
				CCQERecoEQEPlot->Fill(EQE,weight);
				CCQERecoQ2Plot->Fill(reco_Q2,weight);

			}

			// ---------------------------------------------------------------------------------------------------------------------------------------------------

			// CCMEC

			if (genie_mode == 10) {

				CCMECRecoNuScorePlot->Fill(NuScore,weight);
				CCMECRecoFlashScorePlot->Fill(FlashScore,weight);

				CCMECRecoDistancePlot->Fill(distance,weight);

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

				CCMECRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				CCMECRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
				CCMECRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

				CCMECRecoECalPlot->Fill(ECal,weight);
				CCMECRecoEQEPlot->Fill(EQE,weight);
				CCMECRecoQ2Plot->Fill(reco_Q2,weight);

			}

			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// CCRES

			if (genie_mode == 1) {

				CCRESRecoNuScorePlot->Fill(NuScore,weight);
				CCRESRecoFlashScorePlot->Fill(FlashScore,weight);

				CCRESRecoDistancePlot->Fill(distance,weight);

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

				CCRESRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				CCRESRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
				CCRESRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

				CCRESRecoECalPlot->Fill(ECal,weight);
				CCRESRecoEQEPlot->Fill(EQE,weight);
				CCRESRecoQ2Plot->Fill(reco_Q2,weight);

			}

			// ------------------------------------------------------------------------------------------------------------------------------------------------------------

			// CCDIS

			if (genie_mode == 2) {

				CCDISRecoNuScorePlot->Fill(NuScore,weight);
				CCDISRecoFlashScorePlot->Fill(FlashScore,weight);

				CCDISRecoDistancePlot->Fill(distance,weight);

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

				CCDISRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				CCDISRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
				CCDISRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);

				CCDISRecoECalPlot->Fill(ECal,weight);
				CCDISRecoEQEPlot->Fill(EQE,weight);
				CCDISRecoQ2Plot->Fill(reco_Q2,weight);

			}

			// ----------------------------------------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// Chi2 PID Studies

			RecoChi2Plot->Fill(reco_Pmu_chi2,weight);
			RecoChi2Plot->Fill(reco_Pp_chi2,weight);

			RecoThreePlaneChi2Plot->Fill(reco_Pmu_ThreePlanechi2,weight);
			RecoThreePlaneChi2Plot->Fill(reco_Pp_ThreePlanechi2,weight);

			RecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pmu_ThreePlaneLogLikelihood,weight/2.);
			RecoThreePlaneChi2LogLikelihoodPlot->Fill(reco_Pp_ThreePlaneLogLikelihood,weight/2.);

			// --------------------------------------------------------------------------------------------------------------------------------------------------------

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

				// -----------------------------------------------------------------------------------------------------------------------------------------------------

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

				// --------------------------------------------------------------------------------------------------------------------------------------------------

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

			// --------------------------------------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the events

		file->cd();
		file->Write();
		file->Close();

		double nentriesError = sqrt(nentries);
		double NEventsPassingSelectionCutsError = sqrt(NEventsPassingSelectionCuts);

		std::cout << std::endl << "Number of " << fWhichSample << " initial entries = " << nentries << " +/- " << nentriesError 
		<< " (POT normalized: " << nentries*weight << " +/- " << nentriesError*weight << ")" << std::endl;

		std::cout << std::endl << "Number of events passing our selection criteria = " << NEventsPassingSelectionCuts << " +/- " << NEventsPassingSelectionCutsError
		<< " (POT normalized: " << NEventsPassingSelectionCuts*weight << " +/- " << NEventsPassingSelectionCutsError*weight << ")" << std::endl;

		std::cout << std::endl << "Created a new file: " << FileName << std::endl << std::endl << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
