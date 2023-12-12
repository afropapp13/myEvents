#define reco_selection_cxx
#include "reco_selection.h"
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

#include "../myClasses/Tools.h"
#include "../myClasses/STV_Tools.h"

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

	int TotalCounter = 0;
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

		// Txt file to keep track of the run/subrun/event of the candidate events

		TString RunTxtName = "/uboone/data/users/apapadop/my3DEvents/myTxtFiles/"+UBCodeVersion+"/TxtmyRunSubRunEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
		ofstream myRunTxtFile;
		myRunTxtFile.open(RunTxtName);
		myRunTxtFile << std::fixed << std::setprecision(2);
		myRunTxtFile << fWhichSample << endl << endl;

		//----------------------------------------//

		// 1D Reco Level Plots

		TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* RecoMuonCosThetaSingleBinPlot = new TH1D("RecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* RecoDeltaPnPlot = new TH1D("RecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* RecoDeltaPTPlot = new TH1D("RecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* RecoDeltaAlphaTPlot = new TH1D("RecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* RecoDeltaAlpha3DqPlot = new TH1D("RecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);	
		TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D True Level Plots for Signal CC1p reconstructed candidates

		TH1D* CC1pTrueEvPlot = new TH1D("CC1pTrueEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);

		TH1D* CC1pTrueMuonMomentumPlot = new TH1D("CC1pTrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pTrueMuonCosThetaSingleBinPlot = new TH1D("CC1pTrueMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* CC1pTrueProtonMomentumPlot = new TH1D("CC1pTrueProtonMomentumPlot",LabelXAxisProtonMomentum, NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CC1pTrueProtonCosThetaPlot = new TH1D("CC1pTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta, NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CC1pTrueDeltaPTPlot = new TH1D("CC1pTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pTrueDeltaAlphaTPlot = new TH1D("CC1pTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pTrueDeltaAlpha3DqPlot = new TH1D("CC1pTrueDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* CC1pTrueDeltaPnPlot = new TH1D("CC1pTrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		TH1D* CC1pTrueECalPlot = new TH1D("CC1pTrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for Signal CC1p

		TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoMuonCosThetaSingleBinPlot = new TH1D("CC1pRecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CC1pRecoDeltaPnPlot = new TH1D("CC1pRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pRecoDeltaAlpha3DqPlot = new TH1D("CC1pRecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		// 2D Reco Level Plots for Signal CC1p, unweighted

		TH2D* CC1pRecoMuonMomentumPlot2D = new TH2D("CC1pRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
			ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH2D* CC1pRecoProtonMomentumPlot2D = new TH2D("CC1pRecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,NBinsProtonMomentum,
			ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH2D* CC1pRecoMuonCosThetaPlot2D = new TH2D("CC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("CC1pRecoMuonCosThetaSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH2D* CC1pRecoProtonCosThetaPlot2D = new TH2D("CC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
			ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH2D* CC1pRecoDeltaPTPlot2D = new TH2D("CC1pRecoDeltaPTPlot2D",LabelXAxisDeltaPT2D,NBinsDeltaPT,
			ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH2D* CC1pRecoDeltaAlphaTPlot2D = new TH2D("CC1pRecoDeltaAlphaTPlot2D",LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,
			ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH2D* CC1pRecoDeltaAlpha3DqPlot2D = new TH2D("CC1pRecoDeltaAlpha3DqPlot2D",LabelXAxisDeltaAlpha3Dq2D,NBinsDeltaAlpha3Dq,
			ArrayNBinsDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH2D* CC1pRecoECalPlot2D = new TH2D("CC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
		TH2D* CC1pRecoDeltaPnPlot2D = new TH2D("CC1pRecoDeltaPnPlot2D",LabelXAxisDeltaPn2D,NBinsDeltaPn,ArrayNBinsDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		// -------------------------------------------------------------------------------------------------------------------------------------

		// 2D Reco Level Plots for Signal CC1p, POT Scaled for response matrices

		TH2D* POTScaledCC1pRecoMuonMomentumPlot2D = new TH2D("POTScaledCC1pRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
			ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH2D* POTScaledCC1pRecoProtonMomentumPlot2D = new TH2D("POTScaledCC1pRecoProtonMomentumPlot2D",LabelXAxisProtonMomentum2D,NBinsProtonMomentum,
			ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);

		TH2D* POTScaledCC1pRecoMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH2D* POTScaledCC1pRecoProtonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoProtonCosThetaPlot2D",LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,
			ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);


		TH2D* POTScaledCC1pRecoDeltaPTPlot2D = new TH2D("POTScaledCC1pRecoDeltaPTPlot2D",LabelXAxisDeltaPT2D,NBinsDeltaPT,
			ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH2D* POTScaledCC1pRecoDeltaAlphaTPlot2D = new TH2D("POTScaledCC1pRecoDeltaAlphaTPlot2D",LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,
			ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH2D* POTScaledCC1pRecoDeltaAlpha3DqPlot2D = new TH2D("POTScaledCC1pRecoDeltaAlpha3DqPlot2D",LabelXAxisDeltaAlpha3Dq2D,NBinsDeltaAlpha3Dq,
			ArrayNBinsDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH2D* POTScaledCC1pRecoECalPlot2D = new TH2D("POTScaledCC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
		TH2D* POTScaledCC1pRecoDeltaPnPlot2D = new TH2D("POTScaledCC1pRecoDeltaPnPlot2D",LabelXAxisDeltaPn2D,NBinsDeltaPn,ArrayNBinsDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		// -------------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for non-CC1p

		TH1D* NonCC1pRecoMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum); // GeV/c
		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoMuonCosThetaSingleBinPlot = new TH1D("NonCC1pRecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* NonCC1pRecoProtonMomentumPlot = new TH1D("NonCC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum); // GeV/c
		TH1D* NonCC1pRecoProtonCosThetaPlot = new TH1D("NonCC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* NonCC1pRecoDeltaPnPlot = new TH1D("NonCC1pRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* NonCC1pRecoDeltaPTPlot = new TH1D("NonCC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* NonCC1pRecoDeltaAlphaTPlot = new TH1D("NonCC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* NonCC1pRecoDeltaAlpha3DqPlot = new TH1D("NonCC1pRecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* NonCC1pRecoECalPlot = new TH1D("NonCC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		// ---------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCQE

		TH1D* CCQERecoMuonMomentumPlot = new TH1D("CCQERecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCQERecoMuonCosThetaSingleBinPlot = new TH1D("CCQERecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* CCQERecoProtonMomentumPlot = new TH1D("CCQERecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCQERecoProtonCosThetaPlot = new TH1D("CCQERecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCQERecoDeltaPnPlot = new TH1D("CCQERecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		TH1D* CCQERecoDeltaPTPlot = new TH1D("CCQERecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQERecoDeltaAlphaTPlot = new TH1D("CCQERecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQERecoDeltaAlpha3DqPlot = new TH1D("CCQERecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* CCQERecoECalPlot = new TH1D("CCQERecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		// ---------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCMEC

		TH1D* CCMECRecoMuonMomentumPlot = new TH1D("CCMECRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECRecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* CCMECRecoProtonMomentumPlot = new TH1D("CCMECRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCMECRecoProtonCosThetaPlot = new TH1D("CCMECRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CCMECRecoDeltaPnPlot = new TH1D("CCMECRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		TH1D* CCMECRecoDeltaPTPlot = new TH1D("CCMECRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECRecoDeltaAlphaTPlot = new TH1D("CCMECRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECRecoDeltaAlpha3DqPlot = new TH1D("CCMECRecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* CCMECRecoECalPlot = new TH1D("CCMECRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		// ------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCRES
		TH1D* CCRESRecoMuonMomentumPlot = new TH1D("CCRESRecoMuonMomentumPlot",LabelXAxisMuonMomentum,
			NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,
			NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESRecoMuonCosThetaSingleBinPlot","",
			1,0.,1.);
	
		TH1D* CCRESRecoProtonMomentumPlot = new TH1D("CCRESRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCRESRecoProtonCosThetaPlot = new TH1D("CCRESRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCRESRecoDeltaPnPlot = new TH1D("CCRESRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);

		TH1D* CCRESRecoDeltaPTPlot = new TH1D("CCRESRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESRecoDeltaAlphaTPlot = new TH1D("CCRESRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESRecoDeltaAlpha3DqPlot = new TH1D("CCRESRecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* CCRESRecoECalPlot = new TH1D("CCRESRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

		// ---------------------------------------------------------------------------------------------------------------------------------

		// 1D Reco Level Plots for CCDIS

		TH1D* CCDISRecoMuonMomentumPlot = new TH1D("CCDISRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISRecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* CCDISRecoProtonMomentumPlot = new TH1D("CCDISRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
			NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CCDISRecoProtonCosThetaPlot = new TH1D("CCDISRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
			NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TH1D* CCDISRecoDeltaPnPlot = new TH1D("CCDISRecoDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TH1D* CCDISRecoDeltaPTPlot = new TH1D("CCDISRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISRecoDeltaAlphaTPlot = new TH1D("CCDISRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISRecoDeltaAlpha3DqPlot = new TH1D("CCDISRecoDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);

		TH1D* CCDISRecoECalPlot = new TH1D("CCDISRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);

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

		Tools tools;		

		//----------------------------------------//

		// Cosmic rejection plots (crt, min distance)

		TString CRTVetoPlotName = "RecoCRTVetoPlot";
		TString CRTVetoXLabel = ";crt veto;";
		int CRTVetoNBins = 2;
		double CRTVetoMin = -0.5;
		double CRTVetoMax = 1.5;

                TH1D* RecoCRTVeto = new TH1D(CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);
                TH1D* CC1pRecoCRTVeto = new TH1D("CC1p"+CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);
                TH1D* NonCC1pRecoCRTVeto = new TH1D("NonCC1p"+CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);
                TH1D* CCQERecoCRTVeto = new TH1D("CCQE"+CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);
                TH1D* CCMECRecoCRTVeto = new TH1D("CCMEC"+CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);
                TH1D* CCRESRecoCRTVeto = new TH1D("CCRES"+CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);
                TH1D* CCDISRecoCRTVeto = new TH1D("CCDIS"+CRTVetoPlotName,CRTVetoXLabel,CRTVetoNBins,CRTVetoMin,CRTVetoMax);

		//----------------------------------------//

		TString CRTHitPEPlotName = "RecoCRTHitPEPlot";
		TString CRTHitPEXLabel = ";crt hit pe;";
		int CRTHitPENBins = 30;
		double CRTHitPEMin = 0.;
		double CRTHitPEMax = 300.;

                TH1D* RecoCRTHitPE = new TH1D(CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);
                TH1D* CC1pRecoCRTHitPE = new TH1D("CC1p"+CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);
                TH1D* NonCC1pRecoCRTHitPE = new TH1D("NonCC1p"+CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);
                TH1D* CCQERecoCRTHitPE = new TH1D("CCQE"+CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);
                TH1D* CCMECRecoCRTHitPE = new TH1D("CCMEC"+CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);
                TH1D* CCRESRecoCRTHitPE = new TH1D("CCRES"+CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);
                TH1D* CCDISRecoCRTHitPE = new TH1D("CCDIS"+CRTHitPEPlotName,CRTHitPEXLabel,CRTHitPENBins,CRTHitPEMin,CRTHitPEMax);

		//----------------------------------------//

		TString CosmicIPAll3DPlotName = "RecoCosmicIPAll3DPlot";
		TString CosmicIPAll3DXLabel = ";cosmic IP 3D [cm];";
		int CosmicIPAll3DNBins = 30;
		double CosmicIPAll3DMin = 0.;
		double CosmicIPAll3DMax = 300.;

                TH1D* RecoCosmicIPAll3D = new TH1D(CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);
                TH1D* CC1pRecoCosmicIPAll3D = new TH1D("CC1p"+CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);
                TH1D* NonCC1pRecoCosmicIPAll3D = new TH1D("NonCC1p"+CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);
                TH1D* CCQERecoCosmicIPAll3D = new TH1D("CCQE"+CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);
                TH1D* CCMECRecoCosmicIPAll3D = new TH1D("CCMEC"+CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);
                TH1D* CCRESRecoCosmicIPAll3D = new TH1D("CCRES"+CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);
                TH1D* CCDISRecoCosmicIPAll3D = new TH1D("CCDIS"+CosmicIPAll3DPlotName,CosmicIPAll3DXLabel,CosmicIPAll3DNBins,CosmicIPAll3DMin,CosmicIPAll3DMax);

		//----------------------------------------//

		TString CosmicDirAll3DPlotName = "RecoCosmicDirAll3DPlot";
		TString CosmicDirAll3DXLabel = ";cosmic dir 3D [rad];";
		int CosmicDirAll3DNBins = 30;
		double CosmicDirAll3DMin = -1.;
		double CosmicDirAll3DMax = 1.;

                TH1D* RecoCosmicDirAll3D = new TH1D(CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);
                TH1D* CC1pRecoCosmicDirAll3D = new TH1D("CC1p"+CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);
                TH1D* NonCC1pRecoCosmicDirAll3D = new TH1D("NonCC1p"+CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);
                TH1D* CCQERecoCosmicDirAll3D = new TH1D("CCQE"+CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);
                TH1D* CCMECRecoCosmicDirAll3D = new TH1D("CCMEC"+CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);
                TH1D* CCRESRecoCosmicDirAll3D = new TH1D("CCRES"+CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);
                TH1D* CCDISRecoCosmicDirAll3D = new TH1D("CCDIS"+CosmicDirAll3DPlotName,CosmicDirAll3DXLabel,CosmicDirAll3DNBins,CosmicDirAll3DMin,CosmicDirAll3DMax);

		//----------------------------------------//

		// 3D analysis
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

		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

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

		// 3D analysis
		// Ecal in DeltaPn and DeltaAlpha3Dq bins		

		TH1D* RecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
		TH1D* CC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
		TH1D* CC1pTrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];				
		TH1D* NonCC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];		
		TH1D* CCQERecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
		TH1D* CCMECRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
		TH1D* CCRESRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
		TH1D* CCDISRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];
		TH2D* CC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot2D[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];	
		TH2D* POTScaledCC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot2D[TwoDNBinsDeltaPn][TwoDNBinsDeltaAlpha3Dq];

		for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

			//------------------------------//

			for (int WhichDeltaAlpha3Dq = 0; WhichDeltaAlpha3Dq < TwoDNBinsDeltaAlpha3Dq; WhichDeltaAlpha3Dq++) {	

				TString ECalTwoDInDeltaPnDeltaAlpha3DqLabel = "ECal_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"_DeltaAlpha3Dq_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq+1])+"Plot";
		
				RecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("Reco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);
				CC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("CC1pReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);
				CC1pTrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("CC1pTrue"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);					
				NonCC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("NonCC1pReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);
				CCQERecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("CCQEReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);
				CCMECRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("CCMECReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);
				CCRESRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("CCRESReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);
				CCDISRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH1D("CCDISReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);	
				CC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot2D[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH2D("CC1pReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0],TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);	
				POTScaledCC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot2D[WhichDeltaPn][WhichDeltaAlpha3Dq] = new TH2D("POTScaledCC1pReco"+ECalTwoDInDeltaPnDeltaAlpha3DqLabel+"2D",LabelXAxisECal2D,TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0],TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices[WhichDeltaPn][WhichDeltaAlpha3Dq][0]);																									

			} // End of the loop over the 2D DeltaAlpha3Dq bins

		} // End of the loop over the 2D DeltaPn bins
	
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

		for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

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

		for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {	

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

		//---------------------//

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

		//---------------------//

		// ECal in DeltaPT & DeltaAlphaT bins
		TH1D* SerialRecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("RecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH1D* SerialCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("CC1pRecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH1D* SerialCC1pTrueECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("CC1pTrueSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);		
		TH1D* SerialNonCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("NonCC1pRecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH1D* SerialCCQERecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("CCQERecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH1D* SerialCCMECRecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("CCMECRecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);								
		TH1D* SerialCCRESRecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("CCRESRecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH1D* SerialCCDISRecoECal_InDeltaPnDeltaAlpha3DqPlot = new TH1D("CCDISRecoSerialECal_DeltaPnDeltaAlpha3DqPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH2D* SerialCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot2D = new TH2D("CC1pRecoSerialECal_DeltaPnDeltaAlpha3DqPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot2D = new TH2D("POTScaledCC1pRecoSerialECal_DeltaPnDeltaAlpha3DqPlot2D",LabelXAxisECal2D,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0],tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices)[0]);

	
		//---------------------//

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

		// Signal / Bkg for STV PRD

		TH1D* CCQESignalRecoDeltaPTPlot = new TH1D("CCQESignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQESignalRecoDeltaAlphaTPlot = new TH1D("CCQESignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQESignalRecoECalPlot = new TH1D("CCQESignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCQESignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCQESignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCQEBkgRecoDeltaPTPlot = new TH1D("CCQEBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQEBkgRecoDeltaAlphaTPlot = new TH1D("CCQEBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQEBkgRecoECalPlot = new TH1D("CCQEBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);	
		TH1D* CCQEBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCQEBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);		

		TH1D* CCMECSignalRecoDeltaPTPlot = new TH1D("CCMECSignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECSignalRecoDeltaAlphaTPlot = new TH1D("CCMECSignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECSignalRecoECalPlot = new TH1D("CCMECSignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECSignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECSignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCMECBkgRecoDeltaPTPlot = new TH1D("CCMECBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECBkgRecoDeltaAlphaTPlot = new TH1D("CCMECBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECBkgRecoECalPlot = new TH1D("CCMECBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCRESSignalRecoDeltaPTPlot = new TH1D("CCRESSignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESSignalRecoDeltaAlphaTPlot = new TH1D("CCRESSignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESSignalRecoECalPlot = new TH1D("CCRESSignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESSignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESSignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCRESBkgRecoDeltaPTPlot = new TH1D("CCRESBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESBkgRecoDeltaAlphaTPlot = new TH1D("CCRESBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESBkgRecoECalPlot = new TH1D("CCRESBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);		

		TH1D* CCDISSignalRecoDeltaPTPlot = new TH1D("CCDISSignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISSignalRecoDeltaAlphaTPlot = new TH1D("CCDISSignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISSignalRecoECalPlot = new TH1D("CCDISSignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISSignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISSignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCDISBkgRecoDeltaPTPlot = new TH1D("CCDISBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISBkgRecoDeltaAlphaTPlot = new TH1D("CCDISBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISBkgRecoECalPlot = new TH1D("CCDISBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);											
		//----------------------------------------//
		//----------------------------------------//

		// Atmospherics
		
		TH1D* CC1pRecoDeltaThetaPlot[TwoDNBinsDeltaPT][NInte];
		TH1D* CC1pRecoDeltaThetaPercPlot[TwoDNBinsDeltaPT][NInte];
		TH1D* CC1pRecoThetaPlot[TwoDNBinsDeltaPT][NInte];	
		TH1D* CC1pTrueThetaPlot[TwoDNBinsDeltaPT][NInte];	
		TH2D* POTScaledCC1pRecoDeltaThetaPlot2D[TwoDNBinsDeltaPT][NInte];

		double delta_theta_min = -60;
		double delta_theta_max = 60;
		int delta_theta_bins = 30;
	
		// Loop over the DeltaPT slices
		for (int ideltapt = 0; ideltapt < TwoDNBinsDeltaPT; ideltapt++) {

			// Loop over the interaction processes
			for (int iinte = 0; iinte < NInte; iinte++) {

				CC1pRecoDeltaThetaPlot[ideltapt][iinte] = new TH1D("CC1p" + InteractionLabels[iinte] + "RecoDeltaTheta_" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt+1]) +"Plot",";#delta#theta [deg]",30,-15,15);
		
				CC1pRecoDeltaThetaPercPlot[ideltapt][iinte] = new TH1D("CC1p" + InteractionLabels[iinte] + "RecoDeltaThetaPerc_" +tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt+1]) + "Plot",";#delta#theta [%]",25,-100,100);

				CC1pRecoThetaPlot[ideltapt][iinte] = new TH1D("CC1p" + InteractionLabels[iinte] + "RecoTheta_" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt+1]) +"Plot",";#theta_{reco} [deg]",15,0.,30.);

				CC1pTrueThetaPlot[ideltapt][iinte] = new TH1D("CC1p" + InteractionLabels[iinte] + "TrueTheta_" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt+1]) +"Plot",";#theta_{true} [deg]",15,0.,30.);

				POTScaledCC1pRecoDeltaThetaPlot2D[ideltapt][iinte] = new TH2D("POTScaledCC1p" + InteractionLabels[iinte] + "RecoDeltaTheta_" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[ideltapt+1]) +"2DPlot",";#theta_{true} [deg];#theta_{reco} [deg]",15,0.,30.,15,0.,30);
	
			} // End of the loop over the interaction processes

		} // End of the loop over DeltaPT slices
		

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

			double l_muCandidate = CandidateMu_Length->at(0);
			double l_pCandidate = CandidateP_Length->at(0);
			double LengthDifference = l_muCandidate - l_pCandidate;
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
			double DeltaAlpha3Dq = Reco_DeltaAlpha3Dq->at(0);
			double DeltaPhiT = Reco_DeltaPhiT->at(0);
			double DeltaPhi3D = Reco_DeltaPhi3D->at(0);			

			// Minerva variables
			double reco_Pn = Reco_Pn->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction
			double ECal = Reco_ECal->at(0);

			// -------------------------------------------------------------------------------------------------------------------------
			
			STV_Tools reco_stv_tool(TVector3CandidateMuon,TVector3CandidateProton,reco_Emu,reco_Ep);
	
			// STV redefinition if P_p < 0.5 where the biases have been observed
			if (reco_Pp < 0.5) { 

				reco_Pp = ( 1.-0.01*fPP->Eval(reco_Pp) ) * reco_Pp ;
				TVector3CandidateProton.SetMag(reco_Pp);
				reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );

				TransMissMomentum = reco_stv_tool.ReturnPt();
				DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
				DeltaAlpha3Dq = reco_stv_tool.ReturnDeltaAlpha3Dq();
				DeltaPhiT = reco_stv_tool.ReturnDeltaPhiT();
				DeltaPhi3D = reco_stv_tool.ReturnDeltaPhi3D();				
				ECal = reco_stv_tool.ReturnECal();

				reco_Pn = reco_stv_tool.ReturnPn();

			}

			// ----------------------------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------------------------

			// Overflow bins
			// Affects

			// DeltaPT
			// DeltaPn
			// ECal

			if (TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { TransMissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
			if (reco_Pn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { reco_Pn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

			if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }

			// ---------------------------------------------------------------------------------------------------------------------------

			// Underflow bins
			// Affects

			// ECal
			// EQE

			if (ECal < ArrayNBinsECal[0]) { ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			

			// -------------------------------------------------------------------------------------------------------------------------

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

				if (VectorCuts[i] == "_PID_NuScore" && !(reco_p_LLR_Score < ProtonLLRPIDScore) ) 
					{ PassedSelection = false; }

				if (VectorCuts[i] == "_PID_NuScore_CRT" && !(crtveto == 0) ) 
					{ PassedSelection = false; }


			}

			if (PassedSelection == false) { continue; }

			NEventsPassingSelectionCuts++;

			// -------------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Make sure that the same events fill the same plots

			if (reco_Pmu < ArrayNBinsMuonMomentum[0]) { continue; }
			if (reco_Pmu_cos_theta < ArrayNBinsMuonCosTheta[0]) { continue; }

			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }
			if (reco_Pp_cos_theta < ArrayNBinsProtonCosTheta[0]) { continue; }

			if (TransMissMomentum < ArrayNBinsDeltaPT[0]) { continue; }
			if (DeltaAlphaT < ArrayNBinsDeltaAlphaT[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

			if (reco_Pmu > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
			if (reco_Pmu_cos_theta > ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]) { continue; }

			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }
			if (reco_Pp_cos_theta > ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]) { continue; }
			if (DeltaAlphaT > ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]) { continue; }

			KinematicsCounter++;				

			// --------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;

			double true_Pn = -1;
		
			double true_TransMissMomentum = -1;
			double true_DeltaAlphaT = -1;
			double true_DeltaAlpha3Dq = -1;
			double true_DeltaPhiT = -1;
			double true_DeltaPhi3D = -1;			
			double true_ECal = -1;
			double true_nu = -1;

			// Only for MC to obtain true vales			
			
			if (
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

				true_Pn = True_Pn->at(0);
			
				true_TransMissMomentum = True_Pt->at(0);
				true_DeltaAlphaT = True_DeltaAlphaT->at(0);
				true_DeltaAlpha3Dq = True_DeltaAlpha3Dq->at(0);
				true_DeltaPhiT = True_DeltaPhiT->at(0);
				true_DeltaPhi3D = True_DeltaPhi3D->at(0);				
				true_ECal = True_ECal->at(0);
				true_nu = True_Ev - True_CandidateMu_P->at(0);	

				// ----------------------------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------------------------

				// Overflow bins
				// Affects

				// DeltaPT
				// DeltaPn
				// ECal

				if (true_TransMissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { true_TransMissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
				if (true_Pn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { true_Pn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

				if (true_ECal > ArrayNBinsECal[NBinsECal]) { true_ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }

				// ---------------------------------------------------------------------------------------------------------------------------

				// Underflow bins
				// Affects

				// ECal
				
				if (true_ECal < ArrayNBinsECal[0]) { true_ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			

				// ----------------------------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------------------------						
				
			} // End of if statement: Only for MC to obtain true vales

			//----------------------------------------//

			// Indices for 2D analysis

			int DeltaPTTwoDIndex = tools.ReturnIndex(TransMissMomentum, TwoDArrayNBinsDeltaPT);
			int DeltaPnTwoDIndex = tools.ReturnIndex(reco_Pn, TwoDArrayNBinsDeltaPn);			
			int DeltaAlphaTTwoDIndex = tools.ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
			int DeltaAlpha3DqTwoDIndex = tools.ReturnIndex(DeltaAlpha3Dq, TwoDArrayNBinsDeltaAlpha3Dq);
			int MuonCosThetaTwoDIndex = tools.ReturnIndex(reco_Pmu_cos_theta, TwoDArrayNBinsMuonCosTheta);
			int ProtonCosThetaTwoDIndex = tools.ReturnIndex(reco_Pp_cos_theta, TwoDArrayNBinsProtonCosTheta);
			int MuonMomentumTwoDIndex = tools.ReturnIndex(reco_Pmu, TwoDArrayNBinsMuonMomentum);
			int ProtonMomentumTwoDIndex = tools.ReturnIndex(reco_Pp, TwoDArrayNBinsProtonMomentum);			

			// 3D analysis treated in 1D grid

			int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,ECal);
			int SerialECalInDeltaPnDeltaAlpha3DqIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices,DeltaPnTwoDIndex,DeltaAlpha3DqTwoDIndex,ECal);
			int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,ECal);
			int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,ECal);								

			//----------------------------------------//

			// No weight to be applied in the multiplicity plots

			RecoMuonMomentumPlot->Fill(reco_Pmu,weight);
			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			RecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

			RecoProtonMomentumPlot->Fill(reco_Pp,weight);
			RecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

			RecoDeltaPnPlot->Fill(reco_Pn,weight);

			RecoDeltaPTPlot->Fill(TransMissMomentum,weight);
			RecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
			RecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

			RecoECalPlot->Fill(ECal,weight);

			//------------------------------//

			// Cosmic rejection (crt, min distance)

			RecoCRTVeto->Fill(crtveto,weight);
			RecoCRTHitPE->Fill(crthitpe,weight);
			RecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
			RecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);

			// 3D ECal uncorrelated

			RecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
			RecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(ECal,weight);
			RecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
			RecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);							

			// 3D analysis in 1D grid / correlated

			SerialRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
			SerialRecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
			SerialRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
			SerialRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);						

			// ---------------------------------------------------------------------------------------------------------------------------

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
				     && True_CandidateMu_CosTheta->at(0) > ArrayNBinsMuonCosTheta[0] 
                     		     && True_CandidateP_CosTheta->at(0) > ArrayNBinsProtonCosTheta[0]
				     && true_TransMissMomentum > ArrayNBinsDeltaPT[0]
				     && true_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]

				     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
				     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				     && True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta] 
				     && True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
				     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]

				) {

					// --------------------------------------------------------------------------------------------------
				
					CC1pEventsPassingSelectionCuts++;
					isSignal = true;
					myRunTxtFile << endl << "CC1p0pi signal event" << endl;

					// ---------------------------------------------------------------------------------------------------------------------------
					// ---------------------------------------------------------------------------------------------------------------------------
					// 1D Plots using True level info for selected CC1p events  

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(True_CandidateMu_P->at(0),2.) );
					double true_Nu = True_Ev - true_MuonEnergy;

					CC1pTrueEvPlot->Fill(True_Ev,weight);

					CC1pTrueMuonMomentumPlot->Fill(True_CandidateMu_P->at(0),weight);
					CC1pTrueMuonCosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
					CC1pTrueMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					CC1pTrueProtonMomentumPlot->Fill(True_CandidateP_P->at(0),weight);
					CC1pTrueProtonCosThetaPlot->Fill(True_CandidateP_CosTheta->at(0),weight);

					CC1pTrueDeltaPTPlot->Fill(true_TransMissMomentum,weight);
					CC1pTrueDeltaAlphaTPlot->Fill(true_DeltaAlphaT,weight);
					CC1pTrueDeltaAlpha3DqPlot->Fill(true_DeltaAlpha3Dq,weight);

					CC1pTrueDeltaPnPlot->Fill(true_Pn,weight);

					CC1pTrueECalPlot->Fill(true_ECal,weight);

					//------------------------------//

					int TrueDeltaPTTwoDIndex = tools.ReturnIndex(true_TransMissMomentum, TwoDArrayNBinsDeltaPT);
					int TrueDeltaPnTwoDIndex = tools.ReturnIndex(true_Pn, TwoDArrayNBinsDeltaPn);					
					int TrueDeltaAlphaTTwoDIndex = tools.ReturnIndex(true_DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
					int TrueDeltaAlpha3DqTwoDIndex = tools.ReturnIndex(true_DeltaAlpha3Dq, TwoDArrayNBinsDeltaAlpha3Dq);
					int TrueMuonCosThetaTwoDIndex = tools.ReturnIndex(True_CandidateMu_CosTheta->at(0), TwoDArrayNBinsMuonCosTheta);
					int TrueProtonCosThetaTwoDIndex = tools.ReturnIndex(True_CandidateP_CosTheta->at(0), TwoDArrayNBinsProtonCosTheta);
					int TrueMuonMomentumTwoDIndex = tools.ReturnIndex(True_CandidateMu_P->at(0), TwoDArrayNBinsMuonMomentum);
					int TrueProtonMomentumTwoDIndex = tools.ReturnIndex(True_CandidateP_P->at(0), TwoDArrayNBinsProtonMomentum);						

					// 3D analysis treated in 1D grid

					int TrueSerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,TrueDeltaPTTwoDIndex,TrueDeltaAlphaTTwoDIndex,true_ECal);
					int TrueSerialECalInDeltaPnDeltaAlpha3DqIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPnDeltaAlpha3DqSlices,TrueDeltaPnTwoDIndex,TrueDeltaAlpha3DqTwoDIndex,true_ECal);
					int TrueSerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,TrueMuonCosThetaTwoDIndex,TrueMuonMomentumTwoDIndex,true_ECal);
					int TrueSerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,TrueProtonCosThetaTwoDIndex,TrueProtonMomentumTwoDIndex,true_ECal);						

					// 3D analysis in 1D grid

					SerialCC1pTrueECal_InDeltaPTDeltaAlphaTPlot->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCC1pTrueECal_InDeltaPnDeltaAlpha3DqPlot->Fill(TrueSerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialCC1pTrueECal_InMuonCosThetaMuonMomentumPlot->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCC1pTrueECal_InProtonCosThetaProtonMomentumPlot->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,weight);																															
					//----------------------------------------//

					// 3D analysis (uncorrelated)

					CC1pTrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TrueDeltaPTTwoDIndex][TrueDeltaAlphaTTwoDIndex]->Fill(true_ECal,weight);	
					CC1pTrueECal_InDeltaPnDeltaAlpha3DqTwoDPlot[TrueDeltaPnTwoDIndex][TrueDeltaAlpha3DqTwoDIndex]->Fill(true_ECal,weight);	
					CC1pTrueECal_InMuonCosThetaMuonMomentumTwoDPlot[TrueMuonCosThetaTwoDIndex][TrueMuonMomentumTwoDIndex]->Fill(ECal,weight);
					CC1pTrueECal_InProtonCosThetaProtonMomentumTwoDPlot[TrueProtonCosThetaTwoDIndex][TrueProtonMomentumTwoDIndex]->Fill(ECal,weight);										
					//----------------------------------------//

					// 1D Reco Plots for the selected CC1p events 

					CC1pRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);

					CC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CC1pRecoDeltaPnPlot->Fill(reco_Pn,weight);

					CC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CC1pRecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

					CC1pRecoECalPlot->Fill(ECal,weight);
					
					//------------------------------//

                        		// Cosmic rejection (crt, min distance)
                        
                                        CC1pRecoCRTVeto->Fill(crtveto,weight);
                        		CC1pRecoCRTHitPE->Fill(crthitpe,weight);
                		        CC1pRecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
		                        CC1pRecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);

					CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(ECal,weight);	
					CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						
					// 3D analysis in 1D grid

					SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialCC1pRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCC1pRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);																																		
					// --------------------------------------------------------------------------------------------------

					// 2D Plots Kinematic Variables

					CC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu);
					CC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp);

					CC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoMuonCosThetaSingleBinPlot2D->Fill(0.5,0.5);
					CC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta);

					POTScaledCC1pRecoMuonMomentumPlot2D->Fill(True_CandidateMu_P->at(0),reco_Pmu,weight);
					POTScaledCC1pRecoProtonMomentumPlot2D->Fill(True_CandidateP_P->at(0),reco_Pp,weight);

					POTScaledCC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
					POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D->Fill(0.5,0.5,weight);
					POTScaledCC1pRecoProtonCosThetaPlot2D->Fill(True_CandidateP_CosTheta->at(0),reco_Pp_cos_theta,weight);


					// 3D analysis in 1D grid

					SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot2D->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,SerialECalInDeltaPTDeltaAlphaTIndex);
					SerialCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot2D->Fill(TrueSerialECalInDeltaPnDeltaAlpha3DqIndex,SerialECalInDeltaPnDeltaAlpha3DqIndex);
					SerialCC1pRecoECal_InMuonCosThetaMuonMomentumPlot2D->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,SerialECalInMuonCosThetaMuonMomentumIndex);
					SerialCC1pRecoECal_InProtonCosThetaProtonMomentumPlot2D->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,SerialECalInProtonCosThetaProtonMomentumIndex);

					// 3D analysis in 1D grid

					SerialPOTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTPlot2D->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialPOTScaledCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot2D->Fill(TrueSerialECalInDeltaPnDeltaAlpha3DqIndex,SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialPOTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumPlot2D->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialPOTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumPlot2D->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,SerialECalInProtonCosThetaProtonMomentumIndex,weight);

					// -----------------------------------------------------------------------------------------------------

					// True Level STV

					CC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum);
					CC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT);
					CC1pRecoDeltaAlpha3DqPlot2D->Fill(true_DeltaAlpha3Dq,DeltaAlpha3Dq);

					CC1pRecoDeltaPnPlot2D->Fill(true_Pn,reco_Pn);

					POTScaledCC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum,weight);
					POTScaledCC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);
					POTScaledCC1pRecoDeltaAlpha3DqPlot2D->Fill(true_DeltaAlpha3Dq,DeltaAlpha3Dq,weight);

					POTScaledCC1pRecoDeltaPnPlot2D->Fill(true_Pn,reco_Pn,weight);

					// -----------------------------------------------------------------------------------------------------

					// True level energy reconstruction & Q2

					CC1pRecoECalPlot2D->Fill(true_ECal,ECal);

					POTScaledCC1pRecoECalPlot2D->Fill(true_ECal,ECal,weight);

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal,weight);

					CC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot2D[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot2D[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(true_ECal,ECal,weight);

					CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(true_ECal,ECal,weight);
										
					CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(true_ECal,ECal);	
					POTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(true_ECal,ECal,weight);														
					// -------------------------------------------------------------------------------------------------------------------------

					// Playground for CC1p true momenta (longitudinal & perpendicular) ratios

					TVector3 TrueCandidateMuon(1,1,1);
					TrueCandidateMuon.SetMag(True_CandidateMu_P->at(0));
					TrueCandidateMuon.SetPhi(True_CandidateMu_Phi->at(0) * TMath::Pi() / 180.); // rad
					TrueCandidateMuon.SetTheta(TMath::ACos(True_CandidateMu_CosTheta->at(0))); // rad
					double TrueCandidateMuon_E = TMath::Sqrt( TMath::Power(True_CandidateMu_P->at(0),2.) + TMath::Power(MuonMass_GeV,2.));

					TVector3 TrueCandidateProton(1,1,1);
					TrueCandidateProton.SetMag(True_CandidateP_P->at(0));
					TrueCandidateProton.SetPhi(True_CandidateP_Phi->at(0) * TMath::Pi() / 180.); // rad
					TrueCandidateProton.SetTheta(TMath::ACos(True_CandidateP_CosTheta->at(0))); // rad
					double TrueCandidateProton_E = TMath::Sqrt( TMath::Power(True_CandidateP_P->at(0),2.) + TMath::Power(ProtonMass_GeV,2.));
					//-----------------------------------------------------------------------------

					// Atmospherics
				
					STV_Tools true_stv_tool(TrueCandidateMuon,TrueCandidateProton,TrueCandidateMuon_E,TrueCandidateProton_E);

					double reco_delta_theta = reco_stv_tool.ReturnDeltaTheta(); // deg
					double true_delta_theta = true_stv_tool.ReturnDeltaTheta(); // deg
					double diff = reco_delta_theta - true_delta_theta; // deg
					double perc = diff / true_delta_theta * 100.; // deg

					// DeltaPT slices

					CC1pRecoDeltaThetaPlot[DeltaPTTwoDIndex][0]->Fill(diff,weight);	
					CC1pRecoDeltaThetaPlot[DeltaPTTwoDIndex][inte_mode]->Fill(diff,weight);

					CC1pRecoThetaPlot[DeltaPTTwoDIndex][0]->Fill(reco_delta_theta,weight);	
					CC1pRecoThetaPlot[DeltaPTTwoDIndex][inte_mode]->Fill(reco_delta_theta,weight);
	
					CC1pTrueThetaPlot[DeltaPTTwoDIndex][0]->Fill(true_delta_theta,weight);	
					CC1pTrueThetaPlot[DeltaPTTwoDIndex][inte_mode]->Fill(true_delta_theta,weight);
				
					CC1pRecoDeltaThetaPercPlot[DeltaPTTwoDIndex][0]->Fill(perc,weight);	
					CC1pRecoDeltaThetaPercPlot[DeltaPTTwoDIndex][inte_mode]->Fill(perc,weight);
	
					POTScaledCC1pRecoDeltaThetaPlot2D[DeltaPTTwoDIndex][0]->Fill(true_delta_theta,reco_delta_theta,weight);	
					POTScaledCC1pRecoDeltaThetaPlot2D[DeltaPTTwoDIndex][inte_mode]->Fill(true_delta_theta,reco_delta_theta,weight);
	
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
					}
					else if (CC1p1pi == 1) { CC1p1piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC1p1pi event" << endl; }
					else if (CC2p1pi == 1) { CC2p1piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC2p1pi event" << endl;}
					else if (CC2p == 1) { CC2pEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC2p event" << endl; }
					else if (CC3p == 1) { CC3pEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; myRunTxtFile << endl << "CC3p event" << endl; }
					else if (CC3p1pi == 1) { CC3p1piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								 myRunTxtFile << endl << "CC3p1pi event" << endl; }
					else if (CC3p2pi == 1) { CC3p2piEventsPassingSelectionCuts++; ManualNonCC1pEventsPassingSelectionCuts++; OtherMCBkg++; 
								 myRunTxtFile << endl << "CC3p2pi event" << endl;}
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
					     && true_TransMissMomentum > ArrayNBinsDeltaPT[0]
					     && true_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]

					     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
					     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
					     && True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta] 
					     && True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
					     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
						)
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

					}

					NonCC1pRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					NonCC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					NonCC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					NonCC1pRecoDeltaPnPlot->Fill(reco_Pn,weight);

					NonCC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					NonCC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					NonCC1pRecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

					NonCC1pRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

                        		// Cosmic rejection (crt, min distance)
                        
                        		NonCC1pRecoCRTVeto->Fill(crtveto,weight);
                                        NonCC1pRecoCRTHitPE->Fill(crthitpe,weight);
                                        NonCC1pRecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
                                        NonCC1pRecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);

					NonCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					NonCC1pRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(ECal,weight);	
					NonCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					NonCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 3D analysis in 1D grid

					SerialNonCC1pRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialNonCC1pRecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialNonCC1pRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialNonCC1pRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);																																											
				} // End of the Non-CC1p beam related background

				// -------------------------------------------------------------------------------------------------------------------------
				// ------------------------------------------------------------------------------------------------------------------------

				// CCQE

				if (genie_mode == 0) {

					CCQERecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCQERecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					CCQERecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCQERecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCQERecoDeltaPnPlot->Fill(reco_Pn,weight);

					CCQERecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCQERecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCQERecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

					CCQERecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCQE Signal / Bkg	

					if (isSignal) {				

						CCQESignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCQESignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCQESignalRecoECalPlot->Fill(ECal,weight);	
						CCQESignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					} else {

						CCQEBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCQEBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCQEBkgRecoECalPlot->Fill(ECal,weight);
						CCQEBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					}				

					//------------------------------//

                                        // Cosmic rejection (crt, min distance)
                                        
                                        CCQERecoCRTVeto->Fill(crtveto,weight);
                                        CCQERecoCRTHitPE->Fill(crthitpe,weight);
                                        CCQERecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
                                        CCQERecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);

					CCQERecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCQERecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(ECal,weight);	
					CCQERecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCQERecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 3D analysis in 1D grid

					SerialCCQERecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCQERecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialCCQERecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCQERecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);													

				} // End of CCQE selection

				// ------------------------------------------------------------------------------------------------------------------------

				// CCMEC

				if (genie_mode == 10) {

					CCMECRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCMECRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCMECRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					CCMECRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCMECRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCMECRecoDeltaPnPlot->Fill(reco_Pn,weight);

					CCMECRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCMECRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCMECRecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

					CCMECRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCMEC Signal / Bkg	

					if (isSignal) {				

						CCMECSignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCMECSignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCMECSignalRecoECalPlot->Fill(ECal,weight);
						CCMECSignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);							

					} else {

						CCMECBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCMECBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCMECBkgRecoECalPlot->Fill(ECal,weight);
						CCMECBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					}						

					//----------------------------------------//

                                        // Cosmic rejection (crt, min distance)
                                        
                                        CCMECRecoCRTVeto->Fill(crtveto,weight);
                                        CCMECRecoCRTHitPE->Fill(crthitpe,weight);
                                        CCMECRecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
                                        CCMECRecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);
 
					// 3D uncorrelated

					CCMECRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 3D correlated

					SerialCCMECRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCMECRecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialCCMECRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCMECRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);											

				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCRES

				if (genie_mode == 1) {

					CCRESRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCRESRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCRESRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					CCRESRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCRESRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);

					CCRESRecoDeltaPnPlot->Fill(reco_Pn,weight);

					CCRESRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCRESRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCRESRecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

					CCRESRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCRES Signal / Bkg	

					if (isSignal) {				

						CCRESSignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCRESSignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCRESSignalRecoECalPlot->Fill(ECal,weight);
						CCRESSignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);							

					} else {

						CCRESBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCRESBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCRESBkgRecoECalPlot->Fill(ECal,weight);
						CCRESBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					}					

					//----------------------------------------//

                                        // Cosmic rejection (crt, min distance)
                                        
                                        CCRESRecoCRTVeto->Fill(crtveto,weight);
                                        CCRESRecoCRTHitPE->Fill(crthitpe,weight);
                                        CCRESRecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
                                        CCRESRecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);

					// 3D uncorrelated

					CCRESRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCRESRecoECal_InDeltaPnDeltaAlpha3DqTwoDPlot[DeltaPnTwoDIndex][DeltaAlpha3DqTwoDIndex]->Fill(ECal,weight);	
					CCRESRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCRESRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 3D correlated

					SerialCCRESRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCRESRecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialCCRESRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCRESRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);									

				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCDIS

				if (genie_mode == 2) {

					CCDISRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCDISRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					CCDISRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					CCDISRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					CCDISRecoDeltaPnPlot->Fill(reco_Pn,weight);

					CCDISRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					CCDISRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					CCDISRecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);

					CCDISRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCDIS Signal / Bkg	

					if (isSignal) {				

						CCDISSignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCDISSignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCDISSignalRecoECalPlot->Fill(ECal,weight);
						CCDISSignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);							

					} else {

						CCDISBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCDISBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCDISBkgRecoECalPlot->Fill(ECal,weight);
						CCDISBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					}					

					//----------------------------------------//

                                        // Cosmic rejection (crt, min distance)
                                        
                                        CCDISRecoCRTVeto->Fill(crtveto,weight);
                                        CCDISRecoCRTHitPE->Fill(crthitpe,weight);
                                        CCDISRecoCosmicIPAll3D->Fill(CosmicIPAll3D,weight);
                                        CCDISRecoCosmicDirAll3D->Fill(CosmicDirAll3D,weight);

					// 3D uncorrelated

					CCDISRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCDISRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCDISRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 3D correlated

					SerialCCDISRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCDISRecoECal_InDeltaPnDeltaAlpha3DqPlot->Fill(SerialECalInDeltaPnDeltaAlpha3DqIndex,weight);
					SerialCCDISRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCDISRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);									

				}

				// --------------------------------------------------------------------------------------------------------------------------

				// Overlay particle breakdown using the Backtracker for PID studies

				if ( !(CandidateMu_MCParticle_Pdg->size() > 0 && CandidateP_MCParticle_Pdg->size() > 0) ) { 

					cout << "muon candidate true pdg = " << CandidateMu_MCParticle_Pdg->at(0) << endl;
					cout << "proton candidate true pdg = " << CandidateP_MCParticle_Pdg->at(0) << endl;

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

		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

		//----------------------------------------//

		cout << fWhichSample << "  " << KinematicsCounter << " events passing selection criteria" << endl << endl;

		//----------------------------------------//	

		file->cd();
		file->Write();
		file->Close();
		fFile->Close();

		//----------------------------------------//	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
