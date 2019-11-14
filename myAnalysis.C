#define myAnalysis_cxx
#include "myAnalysis.h"
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
#include <fstream>
#include <iomanip>
#include <vector>

#include "../myCCQEAnalysis/Constants.h"
#include "../../MyClasses/Tools.h"

using namespace std;
using namespace Constants;

void myAnalysis::Loop()
{

	if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
	TH1D::SetDefaultSumw2();

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Output Files

	TString FileName = "./OutputFiles/"+UBCodeVersion+"/CCQEAnalysis_"+WhichSample+"_"+WhichRun+"_"+UBCodeVersion+".root";
	TFile* OutputFile = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

	ofstream TxtFile;
	TxtFile.open ("./TxtFiles/"+UBCodeVersion+"/CCQEAnalysis_"+WhichSample+"_"+WhichRun+"_"+UBCodeVersion+".txt");

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D Reco Transverse Variables

	TH1D* RecoDeltaPTPlot = new TH1D("RecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* RecoDeltaAlphaTPlot = new TH1D("RecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* RecoDeltaPhiTPlot = new TH1D("RecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* RecoMuonPhiPlot = new TH1D("RecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* RecoProtonPhiPlot = new TH1D("RecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* RecoQ2Plot = new TH1D("RecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// CC1p 1D Reco True Transverse Variables

	TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* CC1pRecoMuonPhiPlot = new TH1D("CC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* CC1pRecoProtonPhiPlot = new TH1D("CC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* CC1pRecoQ2Plot = new TH1D("CC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	TH1D* CC1pTrueDeltaPTPlot = new TH1D("CC1pTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* CC1pTrueDeltaAlphaTPlot = new TH1D("CC1pTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* CC1pTrueDeltaPhiTPlot = new TH1D("CC1pTrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* CC1pTrueMuonMomentumPlot = new TH1D("CC1pTrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* CC1pTrueMuonPhiPlot = new TH1D("CC1pTrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* CC1pTrueProtonMomentumPlot = new TH1D("CC1pTrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* CC1pTrueProtonPhiPlot = new TH1D("CC1pTrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* CC1pTrueProtonCosThetaPlot = new TH1D("CC1pTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* CC1pTrueECalPlot = new TH1D("CC1pTrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* CC1pTrueQ2Plot = new TH1D("CC1pTrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// 2D Reco True Transverse Variables

	TH2D* CC1pRecoTrueDeltaPTPlot2D = new TH2D("CC1pRecoTrueDeltaPTPlot2D",
		LabelXAxisDeltaPT2D,NBinsDeltaPT,ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH2D* CC1pRecoTrueDeltaAlphaTPlot2D = new TH2D("CC1pRecoTrueDeltaAlphaTPlot2D",
		LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH2D* CC1pRecoTrueDeltaPhiTPlot2D = new TH2D("CC1pRecoTrueDeltaPhiTPlot2D",
		LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH2D* CC1pRecoTrueMuonMomentumPlot2D = new TH2D("CC1pRecoTrueMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH2D* CC1pRecoTrueMuonPhiPlot2D = new TH2D("CC1pRecoTrueMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH2D* CC1pRecoTrueMuonCosThetaPlot2D = new TH2D("CC1pRecoTrueMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH2D* CC1pRecoTrueProtonMomentumPlot2D = new TH2D("CC1pRecoTrueProtonMomentumPlot2D",
		LabelXAxisProtonMomentum2D,NBinsProtonMomentum,ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH2D* CC1pRecoTrueProtonPhiPlot2D = new TH2D("CC1pRecoTrueProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH2D* CC1pRecoTrueProtonCosThetaPlot2D = new TH2D("CC1pRecoTrueProtonCosThetaPlot2D",
		LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH2D* CC1pRecoTrueECalPlot2D = new TH2D("CC1pRecoTrueECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
	TH2D* CC1pRecoTrueQ2Plot2D = new TH2D("CC1pRecoTrueQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Bkg 1D Reco True Transverse Variables

	TH1D* BkgRecoTrueDeltaPTPlot = new TH1D("BkgRecoTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* BkgRecoTrueDeltaAlphaTPlot = new TH1D("BkgRecoTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* BkgRecoTrueDeltaPhiTPlot = new TH1D("BkgRecoTrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* BkgRecoTrueMuonMomentumPlot = new TH1D("BkgRecoTrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* BkgRecoTrueMuonPhiPlot = new TH1D("BkgRecoTrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* BkgRecoTrueMuonCosThetaPlot = new TH1D("BkgRecoTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* BkgRecoTrueProtonMomentumPlot = new TH1D("BkgRecoTrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* BkgRecoTrueProtonPhiPlot = new TH1D("BkgRecoTrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* BkgRecoTrueProtonCosThetaPlot = new TH1D("BkgRecoTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* BkgRecoTrueECalPlot = new TH1D("BkgRecoTrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* BkgRecoTrueQ2Plot = new TH1D("BkgRecoTrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// 2D Reco True Transverse Variables

	TH2D* BkgRecoTrueDeltaPTPlot2D = new TH2D("BkgRecoTrueDeltaPTPlot2D",
		LabelXAxisDeltaPT2D,NBinsDeltaPT,ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH2D* BkgRecoTrueDeltaAlphaTPlot2D = new TH2D("BkgRecoTrueDeltaAlphaTPlot2D",
		LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH2D* BkgRecoTrueDeltaPhiTPlot2D = new TH2D("BkgRecoTrueDeltaPhiTPlot2D",
		LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH2D* BkgRecoMuonMomentumPlot2D = new TH2D("BkgRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH2D* BkgRecoMuonPhiPlot2D = new TH2D("BkgRecoMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH2D* BkgRecoMuonCosThetaPlot2D = new TH2D("BkgRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH2D* BkgRecoProtonMomentumPlot2D = new TH2D("BkgRecoProtonMomentumPlot2D",
		LabelXAxisProtonMomentum2D,NBinsProtonMomentum,ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH2D* BkgRecoProtonPhiPlot2D = new TH2D("BkgRecoProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH2D* BkgRecoProtonCosThetaPlot2D = new TH2D("BkgRecoProtonCosThetaPlot2D",
		LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH2D* BkgRecoTrueECalPlot2D = new TH2D("BkgRecoTrueECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
	TH2D* BkgRecoTrueQ2Plot2D = new TH2D("BkgRecoTrueQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);


	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// 1D True Transverse Variables

	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	Tools tools;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	int TrueCC1pCounter = 0;
	int TrueCCQElikeCounter = 0;

	int CandidateRecoCC1pCounter = 0;
	int CandidateRecoCCQElikeCounter = 0;

	int CandidateTrueCC1pCounter = 0;
	int CandidateTrueCCQElikeCounter = 0;

	int BkgNonCC1pCounter = 0;
	int BkgNonCCQElikeCounter = 0;

	int InTimeCosmicsNonCC1pCounter = 0;
	int InTimeCosmicsNonCCQElikeCounter = 0;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//	for (Long64_t jentry=0; jentry<50000;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//		double weight = 1.;
		double weight = Weight;

		// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Tracks

		int NCandidatePairs = CandidateMuP_Distance->size();
		if(NCandidatePairs != 1) { continue; }

		for (int WhichCandidatePair = 0; WhichCandidatePair < NCandidatePairs; WhichCandidatePair++) {

			if (CandidateMu_StartContainment->at(WhichCandidatePair) != 1) { continue; }
			if (CandidateP_StartContainment->at(WhichCandidatePair) != 1 || CandidateP_EndContainment->at(WhichCandidatePair) != 1) { continue; }
			if (CandidateMu_P->at(WhichCandidatePair) < ArrayNBinsMuonMomentum[0]) { continue; }
			if (CandidateP_P->at(WhichCandidatePair) < ArrayNBinsProtonMomentum[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Candidate muon

			TVector3 TVector3CandidateMuon;
			double CandidateMuonTrackCosTheta = CandidateMu_CosTheta->at(WhichCandidatePair);
			double CandidateMuonTrackTheta = TMath::ACos(CandidateMuonTrackCosTheta);
			double CandidateMuonTrackPhi = CandidateMu_Phi->at(WhichCandidatePair) * TMath::Pi() / 180.;
			double CandidateMuonTrackPhi_Deg = CandidateMu_Phi->at(WhichCandidatePair);
			double CandidateMu_Momentum_GeV = CandidateMu_P->at(WhichCandidatePair); // GeV/c
			double CandidateMu_Momentum_MeV = 1000. * CandidateMu_Momentum_GeV; // MeV/c
			double CandidateMuonTrack_KE_MeV = tools.PToKE(MuonPdg,CandidateMu_Momentum_MeV); // MeV/c
			double CandidateMuonTrack_KE_GeV = CandidateMuonTrack_KE_MeV / 1000.; // GeV/c	
			double CandidateMuonTrack_E_GeV = CandidateMuonTrack_KE_GeV + MuonMass_GeV; // GeV
			TVector3CandidateMuon.SetMagThetaPhi(CandidateMu_Momentum_GeV,CandidateMuonTrackTheta,CandidateMuonTrackPhi);
			TLorentzVector CandidateMuon4V(TVector3CandidateMuon,CandidateMuonTrack_E_GeV);
			TVector3 TVector3CandidateMuonTransverse(TVector3CandidateMuon.X(),TVector3CandidateMuon.Y(),0);
			double TVector3CandidateMuonTransverseMag = TVector3CandidateMuonTransverse.Mag();

			// Candidate proton

			TVector3 TVector3CandidateProton;
			double CandidateProtonTrackCosTheta = CandidateP_CosTheta->at(WhichCandidatePair);
			double CandidateProtonTrackTheta = TMath::ACos(CandidateProtonTrackCosTheta);
			double CandidateProtonTrackPhi = CandidateP_Phi->at(WhichCandidatePair) * TMath::Pi() / 180.;
			double CandidateProtonTrackPhi_Deg = CandidateP_Phi->at(WhichCandidatePair);
			double CandidateP_Momentum_GeV = CandidateP_P->at(WhichCandidatePair); // GeV/c
			double CandidateP_Momentum_MeV = 1000. * CandidateP_Momentum_GeV; // MeV/c
			double CandidateProtonTrack_KE_MeV = tools.PToKE(ProtonPdg,CandidateP_Momentum_MeV); // MeV/c
			double CandidateProtonTrack_KE_GeV = CandidateProtonTrack_KE_MeV / 1000.; // GeV/c	
			double CandidateProtonTrack_E_GeV = CandidateProtonTrack_KE_GeV + ProtonMass_GeV; // GeV
			TVector3CandidateProton.SetMagThetaPhi(CandidateP_Momentum_GeV,CandidateProtonTrackTheta,CandidateProtonTrackPhi);
			TLorentzVector CandidateProton4V(TVector3CandidateProton,CandidateProtonTrack_E_GeV);
			TVector3 TVector3CandidateProtonTransverse(TVector3CandidateProton.X(),TVector3CandidateProton.Y(),0);
			double TVector3CandidateProtonTransverseMag = TVector3CandidateProtonTransverse.Mag();

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Relative angles

			double DeltaPhiProtonMuon = TVector3CandidateMuon.DeltaPhi(TVector3CandidateProton);
			double DeltaPhiProtonMuon_Deg = DeltaPhiProtonMuon * 180. / TMath::Pi();
			if (DeltaPhiProtonMuon_Deg < 0.) { DeltaPhiProtonMuon_Deg += 360.; }
			if (DeltaPhiProtonMuon_Deg >= 360.) { DeltaPhiProtonMuon_Deg -= 360.; }

			double DeltaThetaProtonMuon = TVector3CandidateMuon.Angle(TVector3CandidateProton);
			double DeltaThetaProtonMuon_Deg = DeltaThetaProtonMuon * 180. / TMath::Pi();
			if (DeltaThetaProtonMuon_Deg >= 180.) { DeltaThetaProtonMuon_Deg -= 180.; }
			if (DeltaThetaProtonMuon_Deg < 0.) { DeltaThetaProtonMuon_Deg += 180.; }

			// Selection Cuts

			if (CandidateMu_Chi2_YPlane->at(WhichCandidatePair)  < CandidateP_Chi2_YPlane->at(WhichCandidatePair) 
				|| CandidateP_Chi2_YPlane->at(WhichCandidatePair) > ProtonChi2Cut ) { continue; }
			TVector3 BeamFlash(0,BeamFlashes_YCenter->at(0),BeamFlashes_ZCenter->at(0));
			TVector3 RecoVertex(Vertex_X->at(WhichCandidatePair),Vertex_Y->at(WhichCandidatePair),Vertex_Z->at(WhichCandidatePair));
			double dYZ = (BeamFlash - RecoVertex).Mag();	
			if (dYZ > YZBeamFlashVertexMaxDistance || BeamFlashes_TotalPE->at(0) < NPE) { continue; }
			if ( fabs(DeltaThetaProtonMuon_Deg - DeltaThetaCentralValue) > DeltaThetaOpeningAngle ) { continue; }
//			if ( CandidateMuP_Distance->at(WhichCandidatePair) > MaxMuPDistance ) { continue; } // 6 cm at PreSelection stage
			if ( CandidateMuP_Distance->at(WhichCandidatePair) > 2 ) { continue; } // 2 cm to improve purity & reject cosmics

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Reconstruction

			double RecoECal = CandidateProtonTrack_KE_GeV + CandidateMuonTrack_KE_GeV + MuonMass_GeV + BE; // GeV

			// Reconstructed Neutrino

			TVector3 RecoNeutrinoZaxis(0,0,RecoECal);
			TLorentzVector RecoNeutrinoZaxis4V(0,0,RecoECal,RecoECal);

			TVector3 RecoNeutrino = TVector3CandidateMuon + TVector3CandidateProton;
			double InitialRecoNeutrinoMag = RecoNeutrino.Mag();

			TVector3 TVector3TransMissMomentum(RecoNeutrino.X(),RecoNeutrino.Y(),0);
			double TransMissMomentum = TVector3TransMissMomentum.Mag();

			// Reconstructed q vector using reco level info and calorimetry

			TLorentzVector Recoq4V = RecoNeutrinoZaxis4V - CandidateMuon4V;
			double RecoQ2 = -Recoq4V.Mag2();

			// Define the reco level Transverse Variables

			double RecoDeltaAlphaT = TMath::ACos( (- TVector3CandidateMuonTransverse*TVector3TransMissMomentum) /
				( TVector3CandidateMuonTransverseMag*TransMissMomentum ) ) * 180./TMath::Pi();

			double RecoDeltaPhiT = TMath::ACos( (- TVector3CandidateMuonTransverse*TVector3CandidateProtonTransverse) /
				( TVector3CandidateMuonTransverseMag*TVector3CandidateProtonTransverseMag ) ) * 180./TMath::Pi();

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------
			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// STV // CC1p reco analysis

			if ( 
				TransMissMomentum > ArrayNBinsDeltaPT[0] && TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				&& RecoDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] && RecoDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				&& RecoDeltaPhiT > ArrayNBinsDeltaPhiT[0] && RecoDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]	
			) {
				CandidateRecoCC1pCounter++;
				
				RecoDeltaPTPlot->Fill(TransMissMomentum,weight);
				RecoDeltaAlphaTPlot->Fill(RecoDeltaAlphaT,weight);
				RecoDeltaPhiTPlot->Fill(RecoDeltaPhiT,weight);

				if (RecoECal > ArrayNBinsECal[0] && RecoECal < ArrayNBinsECal[NBinsECal]) { RecoECalPlot->Fill(RecoECal,weight); }
				if (RecoQ2 > ArrayNBinsQ2[0] && RecoQ2 < ArrayNBinsQ2[NBinsQ2]) { RecoQ2Plot->Fill(RecoQ2,weight); }

/*			}

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// CCQElike reco analysis

			if ( 
				fabs(DeltaPhiProtonMuon_Deg - DeltaPhiCentralValue) < DeltaPhiOpeningAngle 
				&& TransMissMomentum < MaxTransMissMomentum
				&& CandidateMu_Momentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
				&& CandidateP_Momentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				&& CandidateMuonTrackCosTheta > ArrayNBinsMuonCosTheta[0] && CandidateMuonTrackCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
				&& CandidateProtonTrackCosTheta > ArrayNBinsProtonCosTheta[0] && CandidateProtonTrackCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			) {*/

				CandidateRecoCCQElikeCounter++;

				RecoMuonMomentumPlot->Fill(CandidateMu_Momentum_GeV,weight);
				RecoMuonPhiPlot->Fill(CandidateMuonTrackPhi_Deg,weight);
				RecoMuonCosThetaPlot->Fill(CandidateMuonTrackCosTheta,weight);

				RecoProtonMomentumPlot->Fill(CandidateP_Momentum_GeV,weight);
				RecoProtonPhiPlot->Fill(CandidateProtonTrackPhi_Deg,weight);
				RecoProtonCosThetaPlot->Fill(CandidateProtonTrackCosTheta,weight);

			}

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------------------------------------------------------

			// Backtracker only for the Overlay sample

			if ( 
//				string(WhichSample).find("Overlay9") != std::string::npos
				WhichSample == "Overlay9"
				&& CandidateMu_MCParticle_Pdg->at(WhichCandidatePair) == MuonPdg
				&& CandidateP_MCParticle_Pdg->at(WhichCandidatePair) == ProtonPdg
				&& CandidateMu_MCParticle_Purity->at(WhichCandidatePair) > PurityThreshold
				&& CandidateP_MCParticle_Purity->at(WhichCandidatePair) > PurityThreshold
			) {


				// Candidate Muon Associated True Info From the BackTracker

				bool TrueCandidateMuonTrackStartContainment = True_CandidateMu_StartContainment->at(WhichCandidatePair);
				bool TrueCandidateMuonTrackEndContainment = True_CandidateMu_EndContainment->at(WhichCandidatePair);
				double TrueCandidateMuonTrackPhi_Deg = True_CandidateMu_Phi->at(WhichCandidatePair); // deg
				double TrueCandidateMuonTrackPhi = TrueCandidateMuonTrackPhi_Deg * TMath::Pi()/180.; // rad
				double TrueCandidateMuonTrackCosTheta = True_CandidateMu_CosTheta->at(WhichCandidatePair);
				double TrueCandidateMuonTrackTheta = TMath::ACos(TrueCandidateMuonTrackCosTheta); // rad
//				double TrueCandidateMuonTrackLength = True_CandidateMu_Length->at(WhichCandidatePair); // cm
				double TrueCandidateMuonTrackMomentum_GeV = True_CandidateMu_P->at(WhichCandidatePair); // GeV
				double TrueCandidateMuonTrackMomentum_MeV = 1000.* TrueCandidateMuonTrackMomentum_GeV; // MeV
				double TrueCandidateMuonTrack_KE_MeV = tools.PToKE(MuonPdg,TrueCandidateMuonTrackMomentum_MeV); // MeV
				double TrueCandidateMuonTrack_KE_GeV = TrueCandidateMuonTrack_KE_MeV / 1000.; // GeV
				double TrueCandidateMuonTrack_E_GeV = TrueCandidateMuonTrack_KE_GeV + MuonMass_GeV; // GeV

				TVector3 TrueCandidateMuonTrackMomentumVector(-99.,-99.,-99.);
				TrueCandidateMuonTrackMomentumVector.SetMagThetaPhi(TrueCandidateMuonTrackMomentum_GeV,TrueCandidateMuonTrackTheta, TrueCandidateMuonTrackPhi);
				TLorentzVector TrueCandidateMuonTrackMomentumVector4V(TrueCandidateMuonTrackMomentumVector,TrueCandidateMuonTrack_E_GeV);
				TVector3 TrueCandidateMuonTrackMomentumVectorTransverse(TrueCandidateMuonTrackMomentumVector.X(),TrueCandidateMuonTrackMomentumVector.Y(),0.);
				double TrueCandidateMuonTrackMomentumVectorTransverseMag = TrueCandidateMuonTrackMomentumVectorTransverse.Mag();

				// Candidate Proton Associated True Info From the BackTracker

				bool TrueCandidateProtonTrackStartContainment = True_CandidateP_StartContainment->at(WhichCandidatePair);
				bool TrueCandidateProtonTrackEndContainment = True_CandidateP_EndContainment->at(WhichCandidatePair);
				double TrueCandidateProtonTrackPhi_Deg = True_CandidateP_Phi->at(WhichCandidatePair); // deg
				double TrueCandidateProtonTrackPhi = TrueCandidateProtonTrackPhi_Deg * TMath::Pi()/180.; // rad
				double TrueCandidateProtonTrackCosTheta = True_CandidateP_CosTheta->at(WhichCandidatePair);
				double TrueCandidateProtonTrackTheta = TMath::ACos(TrueCandidateProtonTrackCosTheta); // rad
//				double TrueCandidateProtonTrackLength = True_CandidateP_Length->at(WhichCandidatePair); // cm
				double TrueCandidateProtonTrackMomentum_GeV = True_CandidateP_P->at(WhichCandidatePair); // GeV
				double TrueCandidateProtonTrackMomentum_MeV = 1000.* TrueCandidateProtonTrackMomentum_GeV; // MeV
				double TrueCandidateProtonTrack_KE_MeV = tools.PToKE(ProtonPdg,TrueCandidateProtonTrackMomentum_MeV); // MeV
				double TrueCandidateProtonTrack_KE_GeV = TrueCandidateProtonTrack_KE_MeV / 1000.; // GeV
				double TrueCandidateProtonTrack_E_GeV = TrueCandidateProtonTrack_KE_GeV + ProtonMass_GeV; // GeV

				TVector3 TrueCandidateProtonTrackMomentumVector(-99.,-99.,-99.);
				TrueCandidateProtonTrackMomentumVector.SetMagThetaPhi(TrueCandidateProtonTrackMomentum_GeV,TrueCandidateProtonTrackTheta, TrueCandidateProtonTrackPhi);
				TLorentzVector TrueCandidateProtonTrackMomentumVector4V(TrueCandidateProtonTrackMomentumVector,TrueCandidateProtonTrack_E_GeV);
				TVector3 TrueCandidateProtonTrackMomentumVectorTransverse(TrueCandidateProtonTrackMomentumVector.X(),TrueCandidateProtonTrackMomentumVector.Y(),0.);
				double TrueCandidateProtonTrackMomentumVectorTransverseMag = TrueCandidateProtonTrackMomentumVectorTransverse.Mag();

				// ---------------------------------------------------------------------------------------------------------------------------------------

				double RecoTrueDeltaPhiProtonMuon = TrueCandidateMuonTrackMomentumVector.DeltaPhi(TrueCandidateProtonTrackMomentumVector);
				double RecoTrueDeltaThetaProtonMuon = TrueCandidateMuonTrackMomentumVector.Angle(TrueCandidateProtonTrackMomentumVector);

				double RecoTrueDeltaPhiProtonMuon_Deg = DeltaPhiProtonMuon * 180. / TMath::Pi();
				double RecoTrueDeltaThetaProtonMuon_Deg = DeltaThetaProtonMuon * 180. / TMath::Pi();

				if (RecoTrueDeltaThetaProtonMuon_Deg < 0.) { RecoTrueDeltaThetaProtonMuon_Deg += 180.; }
				if (RecoTrueDeltaThetaProtonMuon_Deg > 180.) { RecoTrueDeltaThetaProtonMuon_Deg -= 180.; }

				if (RecoTrueDeltaPhiProtonMuon_Deg < 0.) { RecoTrueDeltaPhiProtonMuon_Deg += 360.; }
				if (RecoTrueDeltaPhiProtonMuon_Deg > 360.) { RecoTrueDeltaPhiProtonMuon_Deg -= 360.; }

				// Calorimetric Reconstruction

				double RecoTrueECal = TrueCandidateProtonTrack_KE_GeV + TrueCandidateMuonTrack_E_GeV + BE; // GeV

				// Reconstructed Neutrino

				TVector3 RecoTrueNeutrinoZaxis(0,0,RecoTrueECal);
				TLorentzVector RecoTrueNeutrinoZaxis4V(0,0,RecoTrueECal,RecoTrueECal);

				TVector3 RecoTrueNeutrino = TrueCandidateMuonTrackMomentumVector + TrueCandidateProtonTrackMomentumVector;
				double InitialRecoTrueNeutrinoMag = RecoTrueNeutrino.Mag();

				TVector3 TVector3RecoTrueTransMissMomentum(RecoTrueNeutrino.X(),RecoTrueNeutrino.Y(),0);
				double RecoTrueTransMissMomentum = TVector3RecoTrueTransMissMomentum.Mag();

				// Reconstructed q vector using reco true level info and calorimetry

				TLorentzVector RecoTrueq4V = RecoTrueNeutrinoZaxis4V - TrueCandidateMuonTrackMomentumVector4V;
				double RecoTrueQ2 = -RecoTrueq4V.Mag2();

				// Define the Transverse Variables

				double RecoTrueDeltaAlphaT = 
					TMath::ACos( (- TrueCandidateMuonTrackMomentumVectorTransverse*TVector3RecoTrueTransMissMomentum) /
					( TrueCandidateMuonTrackMomentumVectorTransverseMag*RecoTrueTransMissMomentum ) ) * 180./TMath::Pi();

				double RecoTrueDeltaPhiT = 
					TMath::ACos( (- TrueCandidateMuonTrackMomentumVectorTransverse*TVector3CandidateProtonTransverse) /
					( TrueCandidateMuonTrackMomentumVectorTransverseMag*TVector3CandidateProtonTransverseMag ) ) 
					* 180./TMath::Pi();

				// -----------------------------------------------------------------------------------------------------------------------------------------------
				// -----------------------------------------------------------------------------------------------------------------------------------------------

				// STV analysis // CC1p analysis with Backtracker

				if ( 
					TransMissMomentum > ArrayNBinsDeltaPT[0] && TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
					&& RecoDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] && RecoDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
					&& RecoDeltaPhiT > ArrayNBinsDeltaPhiT[0] && RecoDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
				) {

					// True candidate CC1p events

					if ( 
						RecoTrueTransMissMomentum > ArrayNBinsDeltaPT[0] 
						&& RecoTrueTransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
						&& RecoTrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] 
						&& RecoTrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
						&& RecoTrueDeltaPhiT > ArrayNBinsDeltaPhiT[0] 
						&& RecoTrueDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
						&& TrueCandidateMuonTrackMomentum_GeV > ArrayNBinsMuonMomentum[0] 
						&& TrueCandidateProtonTrackMomentum_GeV > ArrayNBinsProtonMomentum[0]
						&& fabs(RecoTrueDeltaThetaProtonMuon_Deg - DeltaThetaCentralValue) < DeltaThetaOpeningAngle
//						&& TrueCC1pEvent == true
						&& TrueCandidateMuonTrackStartContainment == true
						&& TrueCandidateProtonTrackStartContainment == true
						&& TrueCandidateProtonTrackEndContainment == true
					) {

						CandidateTrueCC1pCounter++;

						CC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CC1pRecoDeltaAlphaTPlot->Fill(RecoDeltaAlphaT,weight);
						CC1pRecoDeltaPhiTPlot->Fill(RecoDeltaPhiT,weight);

						CC1pTrueDeltaPTPlot->Fill(RecoTrueTransMissMomentum,weight);
						CC1pTrueDeltaAlphaTPlot->Fill(RecoTrueDeltaAlphaT,weight);
						CC1pTrueDeltaPhiTPlot->Fill(RecoTrueDeltaPhiT,weight);

						CC1pRecoTrueDeltaPTPlot2D->Fill(RecoTrueTransMissMomentum,TransMissMomentum,weight);
						CC1pRecoTrueDeltaAlphaTPlot2D->Fill(RecoTrueDeltaAlphaT,RecoDeltaAlphaT,weight);
						CC1pRecoTrueDeltaPhiTPlot2D->Fill(RecoTrueDeltaPhiT,RecoDeltaPhiT,weight);

						CC1pRecoMuonMomentumPlot->Fill(CandidateMu_Momentum_GeV,weight);
						CC1pRecoMuonPhiPlot->Fill(CandidateMuonTrackPhi_Deg,weight);
						CC1pRecoMuonCosThetaPlot->Fill(CandidateMuonTrackCosTheta,weight);

						CC1pTrueMuonMomentumPlot->Fill(TrueCandidateMuonTrackMomentum_GeV,weight);
						CC1pTrueMuonPhiPlot->Fill(TrueCandidateMuonTrackPhi_Deg,weight);
						CC1pTrueMuonCosThetaPlot->Fill(TrueCandidateMuonTrackCosTheta,weight);

						CC1pRecoTrueMuonMomentumPlot2D->Fill(TrueCandidateMuonTrackMomentum_GeV,CandidateMu_Momentum_GeV,weight);
						CC1pRecoTrueMuonPhiPlot2D->Fill(TrueCandidateMuonTrackPhi_Deg,CandidateMuonTrackPhi_Deg,weight);
						CC1pRecoTrueMuonCosThetaPlot2D->Fill(TrueCandidateMuonTrackCosTheta,CandidateMuonTrackCosTheta,weight);

						CC1pRecoProtonMomentumPlot->Fill(CandidateP_Momentum_GeV,weight);
						CC1pRecoProtonPhiPlot->Fill(CandidateProtonTrackPhi_Deg,weight);
						CC1pRecoProtonCosThetaPlot->Fill(CandidateProtonTrackCosTheta,weight);

						CC1pTrueProtonMomentumPlot->Fill(TrueCandidateProtonTrackMomentum_GeV,weight);
						CC1pTrueProtonPhiPlot->Fill(TrueCandidateProtonTrackPhi_Deg,weight);
						CC1pTrueProtonCosThetaPlot->Fill(TrueCandidateProtonTrackCosTheta,weight);

						CC1pRecoTrueProtonMomentumPlot2D->Fill(TrueCandidateProtonTrackMomentum_GeV,CandidateP_Momentum_GeV,weight);
						CC1pRecoTrueProtonPhiPlot2D->Fill(TrueCandidateProtonTrackPhi_Deg,CandidateProtonTrackPhi_Deg,weight);
						CC1pRecoTrueProtonCosThetaPlot2D->Fill(TrueCandidateProtonTrackCosTheta,CandidateProtonTrackCosTheta,weight);

						if (RecoECal > ArrayNBinsECal[0] && RecoECal < ArrayNBinsECal[NBinsECal] &&
						    RecoTrueECal > ArrayNBinsECal[0] && RecoTrueECal < ArrayNBinsECal[NBinsECal]) { 
			
							CC1pRecoECalPlot->Fill(RecoECal,weight); 
							CC1pTrueECalPlot->Fill(RecoTrueECal,weight); 
							CC1pRecoTrueECalPlot2D->Fill(RecoTrueECal,RecoECal,weight); 

						}

						if (RecoQ2 > ArrayNBinsQ2[0] && RecoQ2 < ArrayNBinsQ2[NBinsQ2] &&
						    RecoTrueQ2 > ArrayNBinsQ2[0] && RecoTrueQ2 < ArrayNBinsQ2[NBinsQ2]) { 

							CC1pRecoQ2Plot->Fill(RecoQ2,weight); 
							CC1pTrueQ2Plot->Fill(RecoTrueQ2,weight); 
							CC1pRecoTrueQ2Plot2D->Fill(RecoTrueQ2,RecoQ2,weight); 

						}

					} // End of the true candidate CC1p events

					// -----------------------------------------------------------------------------------------------------------------------------------------------
				
					else { // BNB event but NonCC1p

						BkgNonCC1pCounter++;

						BkgRecoTrueDeltaPTPlot->Fill(TransMissMomentum,weight);
						BkgRecoTrueDeltaAlphaTPlot->Fill(RecoDeltaAlphaT,weight);
						BkgRecoTrueDeltaPhiTPlot->Fill(RecoDeltaPhiT,weight);

						BkgRecoTrueDeltaPTPlot2D->Fill(RecoTrueTransMissMomentum,TransMissMomentum,weight);
						BkgRecoTrueDeltaAlphaTPlot2D->Fill(RecoTrueDeltaAlphaT,RecoDeltaAlphaT,weight);
						BkgRecoTrueDeltaPhiTPlot2D->Fill(RecoTrueDeltaPhiT,RecoDeltaPhiT,weight);

						BkgRecoTrueMuonMomentumPlot->Fill(CandidateMu_Momentum_GeV,weight);
						BkgRecoTrueMuonPhiPlot->Fill(CandidateMuonTrackPhi_Deg,weight);
						BkgRecoTrueMuonCosThetaPlot->Fill(CandidateMuonTrackCosTheta,weight);

						BkgRecoMuonMomentumPlot2D->Fill(TrueCandidateMuonTrackMomentum_GeV,CandidateMu_Momentum_GeV,weight);
						BkgRecoMuonPhiPlot2D->Fill(TrueCandidateMuonTrackPhi_Deg,CandidateMuonTrackPhi_Deg,weight);
						BkgRecoMuonCosThetaPlot2D->Fill(TrueCandidateMuonTrackCosTheta,CandidateMuonTrackCosTheta,weight);

						BkgRecoTrueProtonMomentumPlot->Fill(CandidateP_Momentum_GeV,weight);
						BkgRecoTrueProtonPhiPlot->Fill(CandidateProtonTrackPhi_Deg,weight);
						BkgRecoTrueProtonCosThetaPlot->Fill(CandidateProtonTrackCosTheta,weight);

						BkgRecoProtonMomentumPlot2D->Fill(TrueCandidateProtonTrackMomentum_GeV,CandidateP_Momentum_GeV,weight);
						BkgRecoProtonPhiPlot2D->Fill(TrueCandidateProtonTrackPhi_Deg,CandidateProtonTrackPhi_Deg,weight);
						BkgRecoProtonCosThetaPlot2D->Fill(TrueCandidateProtonTrackCosTheta,CandidateProtonTrackCosTheta,weight);

						if (RecoECal > ArrayNBinsECal[0] && RecoECal < ArrayNBinsECal[NBinsECal] &&
						    RecoTrueECal > ArrayNBinsECal[0] && RecoTrueECal < ArrayNBinsECal[NBinsECal]) { 
			
							BkgRecoTrueECalPlot->Fill(RecoECal,weight); 
							BkgRecoTrueECalPlot2D->Fill(RecoTrueECal,RecoECal,weight); 

						}

						if (RecoQ2 > ArrayNBinsQ2[0] && RecoQ2 < ArrayNBinsQ2[NBinsQ2] &&
						    RecoTrueQ2 > ArrayNBinsQ2[0] && RecoTrueQ2 < ArrayNBinsQ2[NBinsQ2]) { 

							BkgRecoTrueQ2Plot->Fill(RecoQ2,weight); 
							BkgRecoTrueQ2Plot2D->Fill(RecoTrueQ2,RecoQ2,weight); 

						}

					} // End of BNB event but NonCC1p

				} // End of the case within reco CC1p range

				// ----------------------------------------------------------------------------------------------------------------------------------------------
				// ----------------------------------------------------------------------------------------------------------------------------------------------

				// CCQElike analysis with Backtracker

/*				if ( 
					CandidateMu_Momentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
					&& CandidateP_Momentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
					&& CandidateMuonTrackCosTheta > ArrayNBinsMuonCosTheta[0] 
					&& CandidateMuonTrackCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
					&& CandidateProtonTrackCosTheta > ArrayNBinsProtonCosTheta[0] 
					&& CandidateProtonTrackCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
					&& fabs(DeltaPhiProtonMuon_Deg - DeltaPhiCentralValue) < DeltaPhiOpeningAngle 
					&& TransMissMomentum < MaxTransMissMomentum
				) {

					if ( 
						TrueCandidateMuonTrackMomentum_GeV > ArrayNBinsMuonMomentum[0] 
						&& TrueCandidateProtonTrackMomentum_GeV > ArrayNBinsProtonMomentum[0]
						&& TrueCandidateMuonTrackMomentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
						&& TrueCandidateProtonTrackMomentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
						&& TrueCandidateMuonTrackCosTheta > ArrayNBinsMuonCosTheta[0] 
						&& TrueCandidateMuonTrackCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
						&& TrueCandidateProtonTrackCosTheta > ArrayNBinsProtonCosTheta[0] 
						&& TrueCandidateProtonTrackCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
						&& fabs(RecoTrueDeltaThetaProtonMuon_Deg - DeltaThetaCentralValue) < DeltaThetaOpeningAngle 
						&& fabs(RecoTrueDeltaPhiProtonMuon_Deg - DeltaPhiCentralValue) < DeltaPhiOpeningAngle 
						&& RecoTrueTransMissMomentum < MaxTransMissMomentum
//						&& TrueCCQElikeEvent == true
					) {

						CandidateTrueCCQElikeCounter++;

						CC1pRecoMuonMomentumPlot->Fill(CandidateMu_Momentum_GeV,weight);
						CC1pRecoMuonPhiPlot->Fill(CandidateMuonTrackPhi_Deg,weight);
						CC1pRecoMuonCosThetaPlot->Fill(CandidateMuonTrackCosTheta,weight);

						CC1pTrueMuonMomentumPlot->Fill(TrueCandidateMuonTrackMomentum_GeV,weight);
						CC1pTrueMuonPhiPlot->Fill(TrueCandidateMuonTrackPhi_Deg,weight);
						CC1pTrueMuonCosThetaPlot->Fill(TrueCandidateMuonTrackCosTheta,weight);

						CC1pRecoTrueMuonMomentumPlot2D->Fill(TrueCandidateMuonTrackMomentum_GeV,CandidateMu_Momentum_GeV,weight);
						CC1pRecoTrueMuonPhiPlot2D->Fill(TrueCandidateMuonTrackPhi_Deg,CandidateMuonTrackPhi_Deg,weight);
						CC1pRecoTrueMuonCosThetaPlot2D->Fill(TrueCandidateMuonTrackCosTheta,CandidateMuonTrackCosTheta,weight);

						CC1pRecoProtonMomentumPlot->Fill(CandidateP_Momentum_GeV,weight);
						CC1pRecoProtonPhiPlot->Fill(CandidateProtonTrackPhi_Deg,weight);
						CC1pRecoProtonCosThetaPlot->Fill(CandidateProtonTrackCosTheta,weight);

						CC1pTrueProtonMomentumPlot->Fill(TrueCandidateProtonTrackMomentum_GeV,weight);
						CC1pTrueProtonPhiPlot->Fill(TrueCandidateProtonTrackPhi_Deg,weight);
						CC1pTrueProtonCosThetaPlot->Fill(TrueCandidateProtonTrackCosTheta,weight);

						CC1pRecoTrueProtonMomentumPlot2D->Fill(TrueCandidateProtonTrackMomentum_GeV,CandidateP_Momentum_GeV,weight);
						CC1pRecoTrueProtonPhiPlot2D->Fill(TrueCandidateProtonTrackPhi_Deg,CandidateProtonTrackPhi_Deg,weight);
						CC1pRecoTrueProtonCosThetaPlot2D->Fill(TrueCandidateProtonTrackCosTheta,CandidateProtonTrackCosTheta,weight);



					} // End of the true candidate CC1p events

					// -----------------------------------------------------------------------------------------------------------------------------------------------
					
					else { // BNB event but NonCCQElike

						BkgNonCCQElikeCounter++;

						BkgRecoTrueMuonMomentumPlot->Fill(CandidateMu_Momentum_GeV,weight);
						BkgRecoTrueMuonPhiPlot->Fill(CandidateMuonTrackPhi_Deg,weight);
						BkgRecoTrueMuonCosThetaPlot->Fill(CandidateMuonTrackCosTheta,weight);

						BkgRecoMuonMomentumPlot2D->Fill(TrueCandidateMuonTrackMomentum_GeV,CandidateMu_Momentum_GeV,weight);
						BkgRecoMuonPhiPlot2D->Fill(TrueCandidateMuonTrackPhi_Deg,CandidateMuonTrackPhi_Deg,weight);
						BkgRecoMuonCosThetaPlot2D->Fill(TrueCandidateMuonTrackCosTheta,CandidateMuonTrackCosTheta,weight);

						BkgRecoTrueProtonMomentumPlot->Fill(CandidateP_Momentum_GeV,weight);
						BkgRecoTrueProtonPhiPlot->Fill(CandidateProtonTrackPhi_Deg,weight);
						BkgRecoTrueProtonCosThetaPlot->Fill(CandidateProtonTrackCosTheta,weight);

						BkgRecoProtonMomentumPlot2D->Fill(TrueCandidateProtonTrackMomentum_GeV,CandidateP_Momentum_GeV,weight);
						BkgRecoProtonPhiPlot2D->Fill(TrueCandidateProtonTrackPhi_Deg,CandidateProtonTrackPhi_Deg,weight);
						BkgRecoProtonCosThetaPlot2D->Fill(TrueCandidateProtonTrackCosTheta,CandidateProtonTrackCosTheta,weight);

					}

				} // End of the CCQElike analysis with Backtracker
*/
				// ----------------------------------------------------------------------------------------------------------------------------------------------
				// ----------------------------------------------------------------------------------------------------------------------------------------------

			} // End of the Backtracker

			// -------------------------------------------------------------------------------------------------------------------------------------------------------
				
			else { // In-time Cosmics

				// STV-CC1p analysis

				if ( 
					TransMissMomentum > ArrayNBinsDeltaPT[0] && TransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
					&& RecoDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] && RecoDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
					&& RecoDeltaPhiT > ArrayNBinsDeltaPhiT[0] && RecoDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
				) {

					InTimeCosmicsNonCC1pCounter++;

					BkgRecoTrueDeltaPTPlot->Fill(TransMissMomentum,weight);
					BkgRecoTrueDeltaAlphaTPlot->Fill(RecoDeltaAlphaT,weight);
					BkgRecoTrueDeltaPhiTPlot->Fill(RecoDeltaPhiT,weight);

					if (RecoECal > ArrayNBinsECal[0] && RecoECal < ArrayNBinsECal[NBinsECal]) { BkgRecoTrueECalPlot->Fill(RecoECal,weight); }
					if (RecoQ2 > ArrayNBinsQ2[0] && RecoQ2 < ArrayNBinsQ2[NBinsQ2]) { BkgRecoTrueQ2Plot->Fill(RecoQ2,weight); }

/*				} // End of the STV-CC1p analysis

				// -----------------------------------------------------------------------------------------------------------------------------------------------------

				// CCQElike analysis

				if ( 
					CandidateMu_Momentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
					&& CandidateP_Momentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
					&& CandidateMuonTrackCosTheta > ArrayNBinsMuonCosTheta[0] 
					&& CandidateMuonTrackCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
					&& CandidateProtonTrackCosTheta > ArrayNBinsProtonCosTheta[0] 
					&& CandidateProtonTrackCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
					&& fabs(DeltaPhiProtonMuon_Deg - DeltaPhiCentralValue) < DeltaPhiOpeningAngle 
					&& TransMissMomentum < MaxTransMissMomentum
				) {

					InTimeCosmicsNonCCQElikeCounter++;
*/
					BkgRecoTrueMuonMomentumPlot->Fill(CandidateMu_Momentum_GeV,weight);
					BkgRecoTrueMuonPhiPlot->Fill(CandidateMuonTrackPhi_Deg,weight);
					BkgRecoTrueMuonCosThetaPlot->Fill(CandidateMuonTrackCosTheta,weight);

					BkgRecoTrueProtonMomentumPlot->Fill(CandidateP_Momentum_GeV,weight);
					BkgRecoTrueProtonPhiPlot->Fill(CandidateProtonTrackPhi_Deg,weight);
					BkgRecoTrueProtonCosThetaPlot->Fill(CandidateProtonTrackCosTheta,weight);

				} // End of the CCQElike analysis

			} // End of in-time Cosmics

		} // End of the loop over the candidate pairs

	} // End of the loop over the events

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// STV-CC1p analysis summary

	std::cout << std::endl;

/*	if ( string(WhichSample).find("Overlay9") != std::string::npos) {

		std::cout << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;
		TxtFile << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

	}*/

	std::cout << "Candidate Reco CC1p events = " << CandidateRecoCC1pCounter << std::endl;
	TxtFile << "Candidate Reco CC1p events = " << CandidateRecoCC1pCounter << std::endl;

	if ( string(WhichSample).find("Overlay9") != std::string::npos) {

		std::cout << "Candidate True CC1p events = " << CandidateTrueCC1pCounter << std::endl;
		TxtFile << "Candidate True CC1p events = " << CandidateTrueCC1pCounter << std::endl;

		std::cout << "Candidate True NonCC1p events = " << BkgNonCC1pCounter << std::endl;
		TxtFile << "Candidate True NonCC1p events = " << BkgNonCC1pCounter << std::endl;

		std::cout << "Candidate In Time Cosmics NonCC1p events = " << InTimeCosmicsNonCC1pCounter << std::endl;
		TxtFile << "Candidate In Time Cosmics NonCC1p events = " << InTimeCosmicsNonCC1pCounter << std::endl;

		double CC1pEfficiency = double(CandidateTrueCC1pCounter) / double(TrueCC1pCounter) * 100.;
		double CC1pPurity = double(CandidateTrueCC1pCounter) / double(CandidateRecoCC1pCounter) * 100.;

/*		std::cout << "STV Efficiency = " << CC1pEfficiency << " %" << std::endl;
		TxtFile << "STV Efficiency = " << CC1pEfficiency << " %" << std::endl;*/

		std::cout << "STV Purity = " << CC1pPurity << " %" << std::endl;
		TxtFile << "STV Purity = " << CC1pPurity << " %" << std::endl;

	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// CCQElike analysis summary

/*	std::cout << std::endl;

	if ( string(WhichSample).find("Overlay9") != std::string::npos) {

		std::cout << std::endl << "True CCQElike events = " << TrueCCQElikeCounter << std::endl;
		TxtFile << std::endl << "True CCQElike events = " << TrueCCQElikeCounter << std::endl;

	}

	std::cout << "Candidate Reco CCQElike events = " << CandidateRecoCCQElikeCounter << std::endl;
	TxtFile << "Candidate Reco CCQElike events = " << CandidateRecoCCQElikeCounter << std::endl;

	if ( string(WhichSample).find("Overlay9") != std::string::npos) {

		std::cout << "Candidate True CCQElike events = " << CandidateTrueCCQElikeCounter << std::endl;
		TxtFile << "Candidate True CCQElike events = " << CandidateTrueCCQElikeCounter << std::endl;

		std::cout << "Candidate True NonCCQElike events = " << BkgNonCCQElikeCounter << std::endl;
		TxtFile << "Candidate True NonCCQElike events = " << BkgNonCCQElikeCounter << std::endl;

		std::cout << "Candidate In Time Cosmics NonCCQElike events = " << InTimeCosmicsNonCCQElikeCounter << std::endl;
		TxtFile << "Candidate In Time Cosmics NonCCQElike events = " << InTimeCosmicsNonCCQElikeCounter << std::endl;

		double CCQElikeEfficiency = double(CandidateTrueCCQElikeCounter) / double(TrueCCQElikeCounter) * 100.;
		double CCQElikePurity = double(CandidateTrueCCQElikeCounter) / double(CandidateRecoCCQElikeCounter) * 100.;

		std::cout << "CCQElike Efficiency = " << CCQElikeEfficiency << " %" << std::endl;
		TxtFile << "CCQElike Efficiency = " << CCQElikeEfficiency << " %" << std::endl;

		std::cout << "CCQElike Purity = " << CCQElikePurity << " %" << std::endl;
		TxtFile << "CCQElike Purity = " << CCQElikePurity << " %" << std::endl;

	}
*/
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	OutputFile->Write();
	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;

} // End of the program
