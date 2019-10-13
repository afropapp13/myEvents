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

#include "../myCCQEAnalysis/Constants.h"

using namespace std;
using namespace Constants;

void t::Loop() {

	if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
	TH1D::SetDefaultSumw2();
	double weight = 1.;

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
//	VectorCuts.push_back("_Chi2");
//	VectorCuts.push_back("_Collinearity");
//	VectorCuts.push_back("_MatchedFlash");
//	VectorCuts.push_back("_Distance");
//	VectorCuts.push_back("_Coplanarity");
//	VectorCuts.push_back("_TransImb");

//	VectorCuts.push_back("_Distance");
//	VectorCuts.push_back("_Chi2");
//	VectorCuts.push_back("_Collinearity");
//	VectorCuts.push_back("_MatchedFlash");
//	VectorCuts.push_back("_Coplanarity");
//	VectorCuts.push_back("_TransImb");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		Cuts = Cuts + VectorCuts[i];

	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	TString FileName = "./OutputFiles/"+UBCodeVersion+"/CCQEStudies_"+WhichSample+Cuts+".root";
	TFile* file = new TFile(FileName,"recreate");

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Binning & Ranges 

	static const int NBinsMuonMomentum = 7; TString RecoLabelXAxisMuonMomentum = ";Reco P_{#mu} (GeV/c)"; 
	static const double BinsMuonMomentum[NBinsMuonMomentum+1] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};

	static const int NBinsProtonMomentum = 7; TString RecoLabelXAxisProtonMomentum = ";Reco P_{p} (GeV/c)";
	static const double BinsProtonMomentum[NBinsProtonMomentum+1] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

	static const int NBinsMuonCosTheta = 7; TString RecoLabelXAxisMuonCosTheta = ";Reco cos(#theta_{#mu})"; 
	static const double BinsMuonCosTheta[NBinsMuonCosTheta+1] = {-0.65,-0.41,-0.17,0.07,0.32,0.56,0.80,0.95};

	static const int NBinsProtonCosTheta = 7; TString RecoLabelXAxisProtonCosTheta = ";Reco cos(#theta_{p})"; 
	static const double BinsProtonCosTheta[NBinsProtonCosTheta+1] = {0.15,0.26,0.37,0.47,0.58,0.69,0.80,0.95};

	static const int NBinsMuonPhi = 7; TString RecoLabelXAxisMuonPhi = ";Reco #phi_{#mu} (deg)"; 
	static const double BinsMuonPhi[NBinsMuonPhi+1] = {-180.0,-128.6,-77.1,-25.7,25.7,77.1,128.6,180.0};

	static const int NBinsProtonPhi = 7; TString RecoLabelXAxisProtonPhi = ";Reco #phi_{p} (deg)"; 
	static const double BinsProtonPhi[NBinsMuonPhi+1] = {-180.0,-128.6,-77.1,-25.7,25.7,77.1,128.6,180.0};

	static const int NBinsChi2CandidateMuon = 25; TString RecoLabelXAxisChi2CandidateMuon = ";(#chi^{2}_{p})^{#mu}";
	static const double MinChi2CandidateMuon = 0., MaxChi2CandidateMuon = 500.;

	static const int NBinsChi2CandidateProton = 25; TString RecoLabelXAxisChi2CandidateProton = ";(#chi^{2}_{p})^{p}";
	static const double MinChi2CandidateProton = 0., MaxChi2CandidateProton = 500.;

	static const int NBinsThreePlaneChi2CandidateMuon = 15; TString RecoLabelXAxisThreePlaneChi2CandidateMuon = ";3-Plane (#chi^{2}_{p})^{#mu}";
	static const double MinThreePlaneChi2CandidateMuon = 0., MaxThreePlaneChi2CandidateMuon = 1.;

	static const int NBinsThreePlaneChi2CandidateProton = 15; TString RecoLabelXAxisThreePlaneChi2CandidateProton = ";3-Plane (#chi^{2}_{p})^{p}";
	static const double MinThreePlaneChi2CandidateProton = 0., MaxThreePlaneChi2CandidateProton = 1.;

	static const int NBinsDeltaTheta = 18; TString RecoLabelXAxisDeltaTheta = ";#Delta#theta";
	static const double MinDeltaTheta = 0., MaxDeltaTheta = 180.;

	static const int NBinsDeltaPhi = 15; TString RecoLabelXAxisDeltaPhi = ";#Delta#phi";
	static const double MinDeltaPhi = 0., MaxDeltaPhi = 360.;

	static const int NBinsdYZ = 10; TString RecoLabelXAxisdYZ = ";d_{YZ} (cm)";
	static const double MindYZ = 0., MaxdYZ = 500.;

	static const int NBinsNPE = 10; TString RecoLabelXAxisNPE = ";# PE";
	static const double MinNPE = 0., MaxNPE = 1000.;

	static const int NBinsDistance = 15; TString RecoLabelXAxisDistance = ";#mu - p distance (cm)";
	static const double MinDistance = 0., MaxDistance = 5.;

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots

	TH1D* RecoDistancePlot = new TH1D("RecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* RecodYZPlot = new TH1D("RecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* RecoNPEPlot = new TH1D("RecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* RecoDeltaThetaPlot = new TH1D("RecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* RecoDeltaPhiPlot = new TH1D("RecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* RecoChi2CandidateMuonPlot = new TH1D("RecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* RecoChi2CandidateProtonPlot = new TH1D("RecoChi2CandidateProtonPlot",RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* RecoThreePlaneChi2CandidateMuonPlot = new TH1D("RecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* RecoThreePlaneChi2CandidateProtonPlot = new TH1D("RecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);

	TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* RecoMuonPhiPlot = new TH1D("RecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* RecoProtonPhiPlot = new TH1D("RecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots for Signal CC1p

	TH1D* CC1pRecoDistancePlot = new TH1D("CC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* CC1pRecodYZPlot = new TH1D("CC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* CC1pRecoNPEPlot = new TH1D("CC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* CC1pRecoDeltaThetaPlot = new TH1D("CC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* CC1pRecoDeltaPhiPlot = new TH1D("CC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* CC1pRecoChi2CandidateMuonPlot = new TH1D("CC1pRecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* CC1pRecoChi2CandidateProtonPlot = new TH1D("CC1pRecoChi2CandidateProtonPlot",RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* CC1pRecoThreePlaneChi2CandidateMuonPlot = new TH1D("CC1pRecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* CC1pRecoThreePlaneChi2CandidateProtonPlot = new TH1D("CC1pRecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);


	TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* CC1pRecoMuonPhiPlot = new TH1D("CC1pRecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* CC1pRecoProtonPhiPlot = new TH1D("CC1pRecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots for non-CC1p

	TH1D* NonCC1pRecoDistancePlot = new TH1D("NonCC1pRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* NonCC1pRecodYZPlot = new TH1D("NonCC1pRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* NonCC1pRecoNPEPlot = new TH1D("NonCC1pRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* NonCC1pRecoDeltaThetaPlot = new TH1D("NonCC1pRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* NonCC1pRecoDeltaPhiPlot = new TH1D("NonCC1pRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* NonCC1pRecoChi2CandidateMuonPlot = new TH1D("NonCC1pRecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* NonCC1pRecoChi2CandidateProtonPlot = new TH1D("NonCC1pRecoChi2CandidateProtonPlot"
		,RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* NonCC1pRecoThreePlaneChi2CandidateMuonPlot = new TH1D("NonCC1pRecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* NonCC1pRecoThreePlaneChi2CandidateProtonPlot = new TH1D("NonCC1pRecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);

	TH1D* NonCC1pRecoMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* NonCC1pRecoProtonMomentumPlot = new TH1D("NonCC1pRecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* NonCC1pRecoProtonCosThetaPlot = new TH1D("NonCC1pRecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* NonCC1pRecoMuonPhiPlot = new TH1D("NonCC1pRecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* NonCC1pRecoProtonPhiPlot = new TH1D("NonCC1pRecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);


	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots for CCQE

	TH1D* CCQERecoDistancePlot = new TH1D("CCQERecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* CCQERecodYZPlot = new TH1D("CCQERecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* CCQERecoNPEPlot = new TH1D("CCQERecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* CCQERecoDeltaThetaPlot = new TH1D("CCQERecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* CCQERecoDeltaPhiPlot = new TH1D("CCQERecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* CCQERecoChi2CandidateMuonPlot = new TH1D("CCQERecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* CCQERecoChi2CandidateProtonPlot = new TH1D("CCQERecoChi2CandidateProtonPlot",RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* CCQERecoThreePlaneChi2CandidateMuonPlot = new TH1D("CCQERecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* CCQERecoThreePlaneChi2CandidateProtonPlot = new TH1D("CCQERecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);

	TH1D* CCQERecoMuonMomentumPlot = new TH1D("CCQERecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* CCQERecoProtonMomentumPlot = new TH1D("CCQERecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* CCQERecoProtonCosThetaPlot = new TH1D("CCQERecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* CCQERecoMuonPhiPlot = new TH1D("CCQERecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* CCQERecoProtonPhiPlot = new TH1D("CCQERecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots for CCMEC

	TH1D* CCMECRecoDistancePlot = new TH1D("CCMECRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* CCMECRecodYZPlot = new TH1D("CCMECRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* CCMECRecoNPEPlot = new TH1D("CCMECRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* CCMECRecoDeltaThetaPlot = new TH1D("CCMECRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* CCMECRecoDeltaPhiPlot = new TH1D("CCMECRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* CCMECRecoChi2CandidateMuonPlot = new TH1D("CCMECRecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* CCMECRecoChi2CandidateProtonPlot = new TH1D("CCMECRecoChi2CandidateProtonPlot",RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* CCMECRecoThreePlaneChi2CandidateMuonPlot = new TH1D("CCMECRecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* CCMECRecoThreePlaneChi2CandidateProtonPlot = new TH1D("CCMECRecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);


	TH1D* CCMECRecoMuonMomentumPlot = new TH1D("CCMECRecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* CCMECRecoProtonMomentumPlot = new TH1D("CCMECRecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* CCMECRecoProtonCosThetaPlot = new TH1D("CCMECRecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* CCMECRecoMuonPhiPlot = new TH1D("CCMECRecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* CCMECRecoProtonPhiPlot = new TH1D("CCMECRecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots for CCRES

	TH1D* CCRESRecoDistancePlot = new TH1D("CCRESRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* CCRESRecodYZPlot = new TH1D("CCRESRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* CCRESRecoNPEPlot = new TH1D("CCRESRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* CCRESRecoDeltaThetaPlot = new TH1D("CCRESRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* CCRESRecoDeltaPhiPlot = new TH1D("CCRESRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* CCRESRecoChi2CandidateMuonPlot = new TH1D("CCRESRecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* CCRESRecoChi2CandidateProtonPlot = new TH1D("CCRESRecoChi2CandidateProtonPlot",RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* CCRESRecoThreePlaneChi2CandidateMuonPlot = new TH1D("CCRESRecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* CCRESRecoThreePlaneChi2CandidateProtonPlot = new TH1D("CCRESRecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);

	TH1D* CCRESRecoMuonMomentumPlot = new TH1D("CCRESRecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* CCRESRecoProtonMomentumPlot = new TH1D("CCRESRecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* CCRESRecoProtonCosThetaPlot = new TH1D("CCRESRecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* CCRESRecoMuonPhiPlot = new TH1D("CCRESRecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* CCRESRecoProtonPhiPlot = new TH1D("CCRESRecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Level Plots for CCDIS

	TH1D* CCDISRecoDistancePlot = new TH1D("CCDISRecoDistancePlot",RecoLabelXAxisDistance,NBinsDistance,MinDistance,MaxDistance);

	TH1D* CCDISRecodYZPlot = new TH1D("CCDISRecodYZPlot",RecoLabelXAxisdYZ,NBinsdYZ,MindYZ,MaxdYZ);
	TH1D* CCDISRecoNPEPlot = new TH1D("CCDISRecoNPEPlot",RecoLabelXAxisNPE,NBinsNPE,MinNPE,MaxNPE);

	TH1D* CCDISRecoDeltaThetaPlot = new TH1D("CCDISRecoDeltaThetaPlot",RecoLabelXAxisDeltaTheta,NBinsDeltaTheta,MinDeltaTheta,MaxDeltaTheta);
	TH1D* CCDISRecoDeltaPhiPlot = new TH1D("CCDISRecoDeltaPhiPlot",RecoLabelXAxisDeltaPhi,NBinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);

	TH1D* CCDISRecoChi2CandidateMuonPlot = new TH1D("CCDISRecoChi2CandidateMuonPlot",RecoLabelXAxisChi2CandidateMuon,NBinsChi2CandidateMuon,MinChi2CandidateMuon,MaxChi2CandidateMuon);
	TH1D* CCDISRecoChi2CandidateProtonPlot = new TH1D("CCDISRecoChi2CandidateProtonPlot",RecoLabelXAxisChi2CandidateProton,NBinsChi2CandidateProton,MinChi2CandidateProton,MaxChi2CandidateProton);

	TH1D* CCDISRecoThreePlaneChi2CandidateMuonPlot = new TH1D("CCDISRecoThreePlaneChi2CandidateMuonPlot",RecoLabelXAxisThreePlaneChi2CandidateMuon,
							NBinsThreePlaneChi2CandidateMuon,MinThreePlaneChi2CandidateMuon,MaxThreePlaneChi2CandidateMuon);
	TH1D* CCDISRecoThreePlaneChi2CandidateProtonPlot = new TH1D("CCDISRecoThreePlaneChi2CandidateProtonPlot",RecoLabelXAxisThreePlaneChi2CandidateProton,
							NBinsThreePlaneChi2CandidateProton,MinThreePlaneChi2CandidateProton,MaxThreePlaneChi2CandidateProton);

	TH1D* CCDISRecoMuonMomentumPlot = new TH1D("CCDISRecoMuonMomentumPlot",RecoLabelXAxisMuonMomentum,NBinsMuonMomentum,BinsMuonMomentum); // GeV/c
	TH1D* CCDISRecoProtonMomentumPlot = new TH1D("CCDISRecoProtonMomentumPlot",RecoLabelXAxisProtonMomentum,NBinsProtonMomentum,BinsProtonMomentum); // GeV/c

	TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",RecoLabelXAxisMuonCosTheta,NBinsMuonCosTheta,BinsMuonCosTheta);
	TH1D* CCDISRecoProtonCosThetaPlot = new TH1D("CCDISRecoProtonCosThetaPlot",RecoLabelXAxisProtonCosTheta,NBinsProtonCosTheta,BinsProtonCosTheta);

	TH1D* CCDISRecoMuonPhiPlot = new TH1D("CCDISRecoMuonPhiPlot",RecoLabelXAxisMuonPhi,NBinsMuonPhi,BinsMuonPhi);
	TH1D* CCDISRecoProtonPhiPlot = new TH1D("CCDISRecoProtonPhiPlot",RecoLabelXAxisProtonPhi,NBinsProtonPhi,BinsProtonPhi);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Loop over the events

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//	for (Long64_t jentry=0; jentry<2000;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);	nbytes += nb;

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//std::cout << "CandidateMu_P->size() = " << CandidateMu_P->size() << std::endl;
//		if (CandidateMuP_Distance->size() != 1) { continue; }
		if (CandidateMu_P->size() != 1) { continue; }

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// POT Scaling

		if (WhichSample == "ExtBNB9") { weight = E1DCNT_wcut / EXT;}
		if (WhichSample == "Overlay9") { weight = ( tor860_wcut / OverlayPOT) * Weight;}
		if (WhichSample == "Overlay9_DLdown") { weight = ( tor860_wcut / OverlayPOT_DLdown) * Weight;}
		if (WhichSample == "Overlay9_SCE") { weight = ( tor860_wcut / OverlayPOT_SCE) * Weight;}
		if (WhichSample == "OverlayDirt9") { weight = ( tor860_wcut / DirtPOT) * Weight;}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		double reco_Pmu_mcs = CandidateMu_P->at(0);
		double reco_Pp = CandidateP_P->at(0);

		double reco_Pmu_cos_theta = CandidateMu_CosTheta->at(0);
		double reco_Pp_cos_theta = CandidateP_CosTheta->at(0);

		double reco_Pmu_phi = CandidateMu_Phi->at(0) * TMath::Pi() / 180.;
		double reco_Pp_phi = CandidateP_Phi->at(0) * TMath::Pi() / 180.;

		double reco_Pmu_chi2 = CandidateMu_Chi2_YPlane->at(0);
		double reco_Pp_chi2 = CandidateP_Chi2_YPlane->at(0);

		double reco_Pmu_ThreePlanechi2 = CandidateMu_ThreePlaneChi2->at(0);
		double reco_Pp_ThreePlanechi2 = CandidateP_ThreePlaneChi2->at(0);

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		TVector3 TVector3CandidateMuon;
		TVector3CandidateMuon.SetMagThetaPhi(CandidateMu_P->at(0),TMath::ACos(CandidateMu_CosTheta->at(0)),CandidateMu_Phi->at(0)*TMath::Pi()/180.);

//		TVector3 TVector3CandidateMuonTrans;
//		TVector3CandidateMuonTrans.SetXYZ(TVector3CandidateMuon.X(),TVector3CandidateMuon.Y(),0.);

		TVector3 TVector3CandidateProton;
		TVector3CandidateProton.SetMagThetaPhi(CandidateP_P->at(0),TMath::ACos(CandidateP_CosTheta->at(0)),CandidateP_Phi->at(0)*TMath::Pi()/180.);

//		TVector3 TVector3CandidateProtonTrans;
//		TVector3CandidateProtonTrans.SetXYZ(TVector3CandidateProton.X(),TVector3CandidateProton.Y(),0.);

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		double DeltaThetaProtonMuon = TVector3CandidateMuon.Angle(TVector3CandidateProton);
		double DeltaThetaProtonMuon_Deg = DeltaThetaProtonMuon * 180. / TMath::Pi();
		if (DeltaThetaProtonMuon_Deg >= 180.) { DeltaThetaProtonMuon_Deg -= 180.; }
		if (DeltaThetaProtonMuon_Deg < 0.) { DeltaThetaProtonMuon_Deg += 180.; }
		double DeltaPhiProtonMuon = TVector3CandidateMuon.DeltaPhi(TVector3CandidateProton);
		double DeltaPhiProtonMuon_Deg = DeltaPhiProtonMuon * 180. / TMath::Pi();
		if (DeltaPhiProtonMuon_Deg >= 360.) { DeltaPhiProtonMuon_Deg -= 360.; }
		if (DeltaPhiProtonMuon_Deg < 0.) { DeltaPhiProtonMuon_Deg += 360.; }

		double NPE = BeamFlashes_TotalPE->at(0);

		TVector3 BeamFlash(0,BeamFlashes_YCenter->at(0),BeamFlashes_ZCenter->at(0));
		TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));
		double dYZ = (BeamFlash - RecoVertex).Mag();

		double distance = CandidateMuP_Distance->at(0);	
		double chi2ratio = CandidateP_Chi2_YPlane->at(0) / CandidateMu_Chi2_YPlane->at(0);

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Selection Cuts

		bool PassedSelection = true;

		for (int i = 0; i < NCuts; i++) {

			if (VectorCuts[i] == "_Chi2" && 
				!(CandidateMu_Chi2_YPlane->at(0) > CandidateP_Chi2_YPlane->at(0) && CandidateP_Chi2_YPlane->at(0) < ProtonChi2Cut && CandidateMu_Chi2_YPlane->at(0) > MuonChi2Cut) )  
				{ PassedSelection = false; }

			if (VectorCuts[i] == "_MatchedFlash" && !(NPE > BeamFlashPEThreshold && dYZ < YZBeamFlashVertexMaxDistance) )
				{ PassedSelection = false; }

			if (VectorCuts[i] == "_Collinearity" && !( fabs(DeltaThetaProtonMuon_Deg - DeltaThetaCentralValue) < DeltaThetaOpeningAngle ) )  { PassedSelection = false; }

			if (VectorCuts[i] == "_Distance" && !( distance < MaxMuPDistance) )  { PassedSelection = false; }

//			if (VectorCuts[i] == "_Coplanarity" && !( fabs(delta_phi - 180.) < 35. ) )  { PassedSelection = false; }

//			if (VectorCuts[i] == "_TransImb" && !( reco_Pt_mcs < 0.35 ) )  { PassedSelection = false; }


		}

		if (PassedSelection == false) { continue; }

		NEventsPassingSelectionCuts++;

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		int genie_mode = -1;
		if (
			WhichSample == "Overlay9"
//			string(WhichSample).find("Overlay") != std::string::npos 
			&& MCParticle_Mode != -1 ) { genie_mode = MCParticle_Mode; }

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		int fCC1p = 0;

		if (    
			WhichSample == "Overlay9"			
//			string(WhichSample).find("Overlay") != std::string::npos
//			&& CC1p == 1
			&& CandidateMu_MCParticle_Pdg->at(0) == MuonPdg
			&& CandidateP_MCParticle_Pdg->at(0) == ProtonPdg
			&& CandidateMu_P->at(0) > ArrayNBinsMuonMomentum[0]
			&& CandidateP_P->at(0) > ArrayNBinsProtonMomentum[0]
			&& True_CandidateMu_P->at(0) > ArrayNBinsMuonMomentum[0]
			&& True_CandidateP_P->at(0) > ArrayNBinsProtonMomentum[0]
		) { fCC1p = 1;}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		RecoDistancePlot->Fill(distance,weight);

		RecodYZPlot->Fill(dYZ,weight);
		RecoNPEPlot->Fill(NPE,weight);

		RecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
		RecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

		RecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2,weight);
		RecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2,weight);

		RecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2,weight);
		RecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2,weight);

		RecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
		RecoProtonMomentumPlot->Fill(reco_Pp,weight);

		RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
		RecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

		RecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
		RecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// CC1p Signal

		if (fCC1p == 1) {

			CC1pRecoDistancePlot->Fill(distance,weight);

			CC1pRecodYZPlot->Fill(dYZ,weight);
			CC1pRecoNPEPlot->Fill(NPE,weight);

			CC1pRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			CC1pRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			CC1pRecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2,weight);
			CC1pRecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2,weight);

			CC1pRecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2,weight);
			CC1pRecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2,weight);

			CC1pRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			CC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);

			CC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			CC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			CC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			CC1pRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Non-CC1p 

		if (fCC1p != 1) {

			NonCC1pRecoDistancePlot->Fill(distance,weight);

			NonCC1pRecodYZPlot->Fill(dYZ,weight);
			NonCC1pRecoNPEPlot->Fill(NPE,weight);

			NonCC1pRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			NonCC1pRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			NonCC1pRecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2,weight);
			NonCC1pRecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2,weight);

			NonCC1pRecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2,weight);
			NonCC1pRecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2,weight);

			NonCC1pRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			NonCC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);

			NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			NonCC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			NonCC1pRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			NonCC1pRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// CCQE

		if (genie_mode == 0) {

			CCQERecoDistancePlot->Fill(distance,weight);

			CCQERecodYZPlot->Fill(dYZ,weight);
			CCQERecoNPEPlot->Fill(NPE,weight);

			CCQERecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			CCQERecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			CCQERecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2,weight);
			CCQERecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2,weight);

			CCQERecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2,weight);
			CCQERecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2,weight);

			CCQERecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			CCQERecoProtonMomentumPlot->Fill(reco_Pp,weight);

			CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			CCQERecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			CCQERecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			CCQERecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// CCMEC

		if (genie_mode == 10) {

			CCMECRecoDistancePlot->Fill(distance,weight);

			CCMECRecodYZPlot->Fill(dYZ,weight);
			CCMECRecoNPEPlot->Fill(NPE,weight);

			CCMECRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			CCMECRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			CCMECRecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2,weight);
			CCMECRecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2,weight);

			CCMECRecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2,weight);
			CCMECRecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2,weight);

			CCMECRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			CCMECRecoProtonMomentumPlot->Fill(reco_Pp,weight);

			CCMECRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			CCMECRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			CCMECRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			CCMECRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// CCRES

		if (genie_mode == 1) {

			CCRESRecoDistancePlot->Fill(distance,weight);

			CCRESRecodYZPlot->Fill(dYZ,weight);
			CCRESRecoNPEPlot->Fill(NPE,weight);

			CCRESRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			CCRESRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			CCRESRecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2,weight);
			CCRESRecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2,weight);

			CCRESRecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2,weight);
			CCRESRecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2,weight);

			CCRESRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			CCRESRecoProtonMomentumPlot->Fill(reco_Pp,weight);

			CCRESRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			CCRESRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			CCRESRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			CCRESRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// CCDIS

		if (genie_mode == 2) {

			CCDISRecoDistancePlot->Fill(distance,weight);

			CCDISRecodYZPlot->Fill(dYZ,weight);
			CCDISRecoNPEPlot->Fill(NPE,weight);

			CCDISRecoDeltaThetaPlot->Fill(DeltaThetaProtonMuon_Deg,weight);
			CCDISRecoDeltaPhiPlot->Fill(DeltaPhiProtonMuon_Deg,weight);

			CCDISRecoChi2CandidateMuonPlot->Fill(reco_Pmu_chi2);
			CCDISRecoChi2CandidateProtonPlot->Fill(reco_Pp_chi2);

			CCDISRecoThreePlaneChi2CandidateMuonPlot->Fill(reco_Pmu_ThreePlanechi2);
			CCDISRecoThreePlaneChi2CandidateProtonPlot->Fill(reco_Pp_ThreePlanechi2);

			CCDISRecoMuonMomentumPlot->Fill(reco_Pmu_mcs,weight);
			CCDISRecoProtonMomentumPlot->Fill(reco_Pp,weight);

			CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			CCDISRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			CCDISRecoMuonPhiPlot->Fill(reco_Pmu_phi*180./TMath::Pi(),weight);
			CCDISRecoProtonPhiPlot->Fill(reco_Pp_phi*180./TMath::Pi(),weight);

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	} // End of the loop over the events

	file->cd();
	file->Write();
	file->Close();

	double nentriesError = sqrt(nentries);
	double NEventsPassingSelectionCutsError = sqrt(NEventsPassingSelectionCuts);

	std::cout << std::endl << "Number of " << WhichSample << " initial entries = " << nentries << " +/- " << nentriesError 
	<< " (POT normalized: " << nentries*weight << " +/- " << nentriesError*weight << ")" << std::endl;

	std::cout << std::endl << "Number of events passing our selection criteria = " << NEventsPassingSelectionCuts << " +/- " << NEventsPassingSelectionCutsError
	<< " (POT normalized: " << NEventsPassingSelectionCuts*weight << " +/- " << NEventsPassingSelectionCutsError*weight << ")" << std::endl;

	std::cout << std::endl << "Created a new file: " << FileName << std::endl << std::endl;


} // End of the program
