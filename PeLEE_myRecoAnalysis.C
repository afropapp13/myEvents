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
		TH1D* RecoDeltaPhiTPlot = new TH1D("RecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* RecoDeltaPhi3DPlot = new TH1D("RecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

		TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		

		// -------------------------------------------------------------------------------------------------------------------------------

		// 1D True Level Plots for Signal CC1p reconstructed candidates

		TH1D* CC1pTrueNuPlot = new TH1D("CC1pTrueNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
		TH1D* CC1pTrueEvPlot = new TH1D("CC1pTrueEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);

		TH1D* CC1pTrueMuonMomentumPlot = new TH1D("CC1pTrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pTrueMuonCosThetaSingleBinPlot = new TH1D("CC1pTrueMuonCosThetaSingleBinPlot","",1,0.,1.);

		TH1D* CC1pTrueProtonMomentumPlot = new TH1D("CC1pTrueProtonMomentumPlot",LabelXAxisProtonMomentum, NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TH1D* CC1pTrueProtonCosThetaPlot = new TH1D("CC1pTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta, NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

		TH1D* CC1pTrueDeltaPTPlot = new TH1D("CC1pTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CC1pTrueDeltaAlphaTPlot = new TH1D("CC1pTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CC1pTrueDeltaAlpha3DqPlot = new TH1D("CC1pTrueDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);
		TH1D* CC1pTrueDeltaPhiTPlot = new TH1D("CC1pTrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CC1pTrueDeltaPhi3DPlot = new TH1D("CC1pTrueDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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
		TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CC1pRecoDeltaPhi3DPlot = new TH1D("CC1pRecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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
		TH2D* CC1pRecoDeltaPhiTPlot2D = new TH2D("CC1pRecoDeltaPhiTPlot2D",LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,
			ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH2D* CC1pRecoDeltaPhi3DPlot2D = new TH2D("CC1pRecoDeltaPhi3DPlot2D",LabelXAxisDeltaPhi3D2D,NBinsDeltaPhi3D,
			ArrayNBinsDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);			

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
		TH2D* POTScaledCC1pRecoDeltaPhiTPlot2D = new TH2D("POTScaledCC1pRecoDeltaPhiTPlot2D",LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,
			ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH2D* POTScaledCC1pRecoDeltaPhi3DPlot2D = new TH2D("POTScaledCC1pRecoDeltaPhi3DPlot2D",LabelXAxisDeltaPhi3D2D,NBinsDeltaPhi3D,
			ArrayNBinsDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);			

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
		TH1D* NonCC1pRecoDeltaPhiTPlot = new TH1D("NonCC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* NonCC1pRecoDeltaPhi3DPlot = new TH1D("NonCC1pRecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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
		TH1D* CCQERecoDeltaPhiTPlot = new TH1D("CCQERecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCQERecoDeltaPhi3DPlot = new TH1D("CCQERecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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
		TH1D* CCMECRecoDeltaPhiTPlot = new TH1D("CCMECRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCMECRecoDeltaPhi3DPlot = new TH1D("CCMECRecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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
		TH1D* CCRESRecoDeltaPhiTPlot = new TH1D("CCRESRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCRESRecoDeltaPhi3DPlot = new TH1D("CCRESRecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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
		TH1D* CCDISRecoDeltaPhiTPlot = new TH1D("CCDISRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCDISRecoDeltaPhi3DPlot = new TH1D("CCDISRecoDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);		

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

		//----------------------------------------//			

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

		// Signal / Bkg for STV PRD

		TH1D* CCQESignalRecoDeltaPTPlot = new TH1D("CCQESignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQESignalRecoDeltaAlphaTPlot = new TH1D("CCQESignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQESignalRecoDeltaPhiTPlot = new TH1D("CCQESignalRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCQESignalRecoECalPlot = new TH1D("CCQESignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCQESignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCQESignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCQEBkgRecoDeltaPTPlot = new TH1D("CCQEBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCQEBkgRecoDeltaAlphaTPlot = new TH1D("CCQEBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCQEBkgRecoDeltaPhiTPlot = new TH1D("CCQEBkgRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCQEBkgRecoECalPlot = new TH1D("CCQEBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);	
		TH1D* CCQEBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCQEBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);		

		TH1D* CCMECSignalRecoDeltaPTPlot = new TH1D("CCMECSignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECSignalRecoDeltaAlphaTPlot = new TH1D("CCMECSignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECSignalRecoDeltaPhiTPlot = new TH1D("CCMECSignalRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCMECSignalRecoECalPlot = new TH1D("CCMECSignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECSignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECSignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCMECBkgRecoDeltaPTPlot = new TH1D("CCMECBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCMECBkgRecoDeltaAlphaTPlot = new TH1D("CCMECBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCMECBkgRecoDeltaPhiTPlot = new TH1D("CCMECBkgRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCMECBkgRecoECalPlot = new TH1D("CCMECBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCMECBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCRESSignalRecoDeltaPTPlot = new TH1D("CCRESSignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESSignalRecoDeltaAlphaTPlot = new TH1D("CCRESSignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESSignalRecoDeltaPhiTPlot = new TH1D("CCRESSignalRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCRESSignalRecoECalPlot = new TH1D("CCRESSignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESSignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESSignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCRESBkgRecoDeltaPTPlot = new TH1D("CCRESBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCRESBkgRecoDeltaAlphaTPlot = new TH1D("CCRESBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCRESBkgRecoDeltaPhiTPlot = new TH1D("CCRESBkgRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCRESBkgRecoECalPlot = new TH1D("CCRESBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCRESBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);		

		TH1D* CCDISSignalRecoDeltaPTPlot = new TH1D("CCDISSignalRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISSignalRecoDeltaAlphaTPlot = new TH1D("CCDISSignalRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISSignalRecoDeltaPhiTPlot = new TH1D("CCDISSignalRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCDISSignalRecoECalPlot = new TH1D("CCDISSignalRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISSignalRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISSignalRecoMuonCosThetaSingleBinPlot","",1,0.,1.);			

		TH1D* CCDISBkgRecoDeltaPTPlot = new TH1D("CCDISBkgRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TH1D* CCDISBkgRecoDeltaAlphaTPlot = new TH1D("CCDISBkgRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TH1D* CCDISBkgRecoDeltaPhiTPlot = new TH1D("CCDISBkgRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TH1D* CCDISBkgRecoECalPlot = new TH1D("CCDISBkgRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TH1D* CCDISBkgRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISBkgRecoMuonCosThetaSingleBinPlot","",1,0.,1.);											

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
				if (fTune == "GENIEv2") { weight = POTWeight; }
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
				fWhichSample == "Overlay9_Run3" || fWhichSample == "Overlay9_Run4a" || fWhichSample == "Overlay9_Run4aRutgers"
			|| fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5" || fWhichSample == "Overlay9_Combined" 
			|| fWhichSample == "OverlayDirt9_Run1" || fWhichSample == "OverlayDirt9_Run2" || 
				fWhichSample == "OverlayDirt9_Run3" || fWhichSample == "OverlayDirt9_Run4a" || fWhichSample == "OverlayDirt9_Run4aRutgers"
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

					int concat = tools.ConcatRunSubRunEvent(Run,SubRun,Event,fUniverseIndex);
					weight = weight*tools.PoissonRandomNumber(concat); 
				
				}				

			}	

			//----------------------------------------//

			// detailed xsec genie uncertainties	

			if ( 
				fUniverseIndex != -1 && fWhichSample == "Overlay9_Run1_DecompXSecUnc"
			) {

				// Watch out: The EventWeight weights already include the weight for the tune

				double SF = 1.;

				// Genie weights
				if (fEventWeightLabel == "AGKYpT1pi_UBGenie") { SF = double(AGKYpT1pi_UBGenie->at(fUniverseIndex)) / 1000. ; }
				if (fEventWeightLabel == "AGKYxF1pi_UBGenie") { SF = AGKYxF1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "AhtBY_UBGenie") { SF =  AhtBY_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "BhtBY_UBGenie") { SF =  BhtBY_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "CV1uBY_UBGenie") { SF =  CV1uBY_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "CV2uBY_UBGenie") { SF =  CV2uBY_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "EtaNCEL_UBGenie") { SF =  EtaNCEL_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FrAbs_N_UBGenie") { SF =  FrAbs_N_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FrAbs_pi_UBGenie") { SF =  FrAbs_pi_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "FrCEx_N_UBGenie") { SF =  FrCEx_N_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FrCEx_pi_UBGenie") { SF =  FrCEx_pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FrInel_N_UBGenie") { SF =  FrInel_N_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FrInel_pi_UBGenie") { SF =  FrInel_pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FrPiProd_N_UBGenie") { SF =  FrPiProd_N_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "FrPiProd_pi_UBGenie") { SF =  FrPiProd_pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FracDelta_CCMEC_UBGenie") { SF =  FracDelta_CCMEC_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "FracPN_CCMEC_UBGenie") { SF =  FracPN_CCMEC_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "MFP_N_UBGenie") { SF =  MFP_N_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "MFP_pi_UBGenie") { SF =  MFP_pi_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "MaCCQE_UBGenie") { SF =  double(MaCCQE_UBGenie->at(fUniverseIndex)) / 1000. ; }	
				if (fEventWeightLabel == "MaCCRES_UBGenie") { SF =  MaCCRES_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "MaNCEL_UBGenie") { SF =  MaNCEL_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "MaNCRES_UBGenie") { SF =  MaNCRES_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "MvCCRES_UBGenie") { SF =  MvCCRES_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "MvNCRES_UBGenie") { SF =  MvNCRES_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarnCC1pi_UBGenie") { SF =  NonRESBGvbarnCC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarnCC2pi_UBGenie") { SF =  NonRESBGvbarnCC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarnNC1pi_UBGenie") { SF =  NonRESBGvbarnNC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarnNC2pi_UBGenie") { SF =  NonRESBGvbarnNC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "NonRESBGvbarpCC1pi_UBGenie") { SF =  NonRESBGvbarpCC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarpCC2pi_UBGenie") { SF =  NonRESBGvbarpCC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarpNC1pi_UBGenie") { SF =  NonRESBGvbarpNC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvbarpNC2pi_UBGenie") { SF =  NonRESBGvbarpNC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvnCC1pi_UBGenie") { SF =  NonRESBGvnCC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "NonRESBGvnCC2pi_UBGenie") { SF =  NonRESBGvnCC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvnNC1pi_UBGenie") { SF =  NonRESBGvnNC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvnNC2pi_UBGenie") { SF =  NonRESBGvnNC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvpCC1pi_UBGenie") { SF =  NonRESBGvpCC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "NonRESBGvpCC2pi_UBGenie") { SF =  NonRESBGvpCC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvpNC1pi_UBGenie") { SF =  NonRESBGvpNC1pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NonRESBGvpNC2pi_UBGenie") { SF =  NonRESBGvpNC2pi_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NormCCMEC_UBGenie") { SF =  NormCCMEC_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "NormNCMEC_UBGenie") { SF =  NormNCMEC_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "RDecBR1eta_UBGenie") { SF =  RDecBR1eta_UBGenie->at(fUniverseIndex) / 1000. ; }	
				if (fEventWeightLabel == "RDecBR1gamma_UBGenie") { SF =  RDecBR1gamma_UBGenie->at(fUniverseIndex) / 1000. ; }

				// Unisims

				if (fEventWeightLabel == "UnShortAxFFCCQEshape_UBGenie") { SF =  UnShortAxFFCCQEshape_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "UnShortDecayAngMEC_UBGenie") { SF =  UnShortDecayAngMEC_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "UnShortRPA_CCQE_UBGenie") { SF =  UnShortRPA_CCQE_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "UnShortTheta_Delta2Npi_UBGenie") { SF =  UnShortTheta_Delta2Npi_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "UnShortVecFFCCQEshape_UBGenie") { SF =  UnShortVecFFCCQEshape_UBGenie->at(fUniverseIndex) / 1000. ; }
				if (fEventWeightLabel == "UnShortXSecShape_CCMEC_UBGenie") { SF =  UnShortXSecShape_CCMEC_UBGenie->at(fUniverseIndex) / 1000. ; }																																

				weight = weight * SF / T2KWeight;		

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
			double DeltaAlpha3Dq = Reco_DeltaAlpha3Dq->at(0);
			double DeltaPhiT = Reco_DeltaPhiT->at(0);
			double DeltaPhi3D = Reco_DeltaPhi3D->at(0);			

			// Minerva variables
			double reco_Pn = Reco_Pn->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// Calorimetric Energy Reconstruction
			double ECal = Reco_ECal->at(0);

			// -------------------------------------------------------------------------------------------------------------------------

			// STV redefinition if P_p < 0.5 where the biases have been observed
			if (reco_Pp < 0.5) { 

				reco_Pp = ( 1.-0.01*fPP->Eval(reco_Pp) ) * reco_Pp ;
				TVector3CandidateProton.SetMag(reco_Pp);
				reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );

				STV_Tools reco_stv_tool(TVector3CandidateMuon,TVector3CandidateProton,reco_Emu,reco_Ep);
	
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

			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }
			if (reco_Pp_cos_theta < ArrayNBinsProtonCosTheta[0]) { continue; }

			if (TransMissMomentum < ArrayNBinsDeltaPT[0]) { continue; }
			if (DeltaAlphaT < ArrayNBinsDeltaAlphaT[0]) { continue; }
			if (DeltaPhiT < ArrayNBinsDeltaPhiT[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

			if (reco_Pmu > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
			if (reco_Pmu_cos_theta > ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]) { continue; }

			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }
			if (reco_Pp_cos_theta > ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]) { continue; }
			if (DeltaAlphaT > ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]) { continue; }
			if (DeltaPhiT > ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]) { continue; }

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
			RecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
			RecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);			

			RecoECalPlot->Fill(ECal,weight);

			//------------------------------//

			// 3D ECal uncorrelated

			RecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
			RecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
			RecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);							

			// 3D analysis in 1D grid / correlated

			SerialRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
			SerialRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
			SerialRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);						

			// ---------------------------------------------------------------------------------------------------------------------------

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

				bool isSignal = false;

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
				     && true_TransMissMomentum > ArrayNBinsDeltaPT[0]
				     && true_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]
				     && true_DeltaPhiT > ArrayNBinsDeltaPhiT[0]

				     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
				     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				     && True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta] 
				     && True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
				     && true_DeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				     && true_DeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]

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

					CC1pTrueNuPlot->Fill(true_Nu,weight);
					CC1pTrueEvPlot->Fill(True_Ev,weight);

					CC1pTrueMuonMomentumPlot->Fill(True_CandidateMu_P->at(0),weight);
					CC1pTrueMuonCosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
					CC1pTrueMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					CC1pTrueProtonMomentumPlot->Fill(True_CandidateP_P->at(0),weight);
					CC1pTrueProtonCosThetaPlot->Fill(True_CandidateP_CosTheta->at(0),weight);

					CC1pTrueDeltaPTPlot->Fill(true_TransMissMomentum,weight);
					CC1pTrueDeltaAlphaTPlot->Fill(true_DeltaAlphaT,weight);
					CC1pTrueDeltaAlpha3DqPlot->Fill(true_DeltaAlpha3Dq,weight);
					CC1pTrueDeltaPhiTPlot->Fill(true_DeltaPhiT,weight);
					CC1pTrueDeltaPhi3DPlot->Fill(true_DeltaPhi3D,weight);					

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
					int TrueSerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,TrueMuonCosThetaTwoDIndex,TrueMuonMomentumTwoDIndex,true_ECal);
					int TrueSerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,TrueProtonCosThetaTwoDIndex,TrueProtonMomentumTwoDIndex,true_ECal);						

					// 3D analysis in 1D grid

					SerialCC1pTrueECal_InDeltaPTDeltaAlphaTPlot->Fill(TrueSerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCC1pTrueECal_InMuonCosThetaMuonMomentumPlot->Fill(TrueSerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCC1pTrueECal_InProtonCosThetaProtonMomentumPlot->Fill(TrueSerialECalInProtonCosThetaProtonMomentumIndex,weight);																																						

					//----------------------------------------//

					// 3D analysis (uncorrelated)

					CC1pTrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TrueDeltaPTTwoDIndex][TrueDeltaAlphaTTwoDIndex]->Fill(true_ECal,weight);	
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
					CC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
					CC1pRecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);					

					CC1pRecoECalPlot->Fill(ECal,weight);
					
					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 3D analysis in 1D grid

					SerialCC1pRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
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
					CC1pRecoDeltaAlpha3DqPlot2D->Fill(true_DeltaAlpha3Dq,DeltaAlpha3Dq);
					CC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT);
					CC1pRecoDeltaPhi3DPlot2D->Fill(true_DeltaPhi3D,DeltaPhi3D);					

					CC1pRecoDeltaPnPlot2D->Fill(true_Pn,reco_Pn);

					POTScaledCC1pRecoDeltaPTPlot2D->Fill(true_TransMissMomentum,TransMissMomentum,weight);
					POTScaledCC1pRecoDeltaAlphaTPlot2D->Fill(true_DeltaAlphaT,DeltaAlphaT,weight);
					POTScaledCC1pRecoDeltaAlpha3DqPlot2D->Fill(true_DeltaAlpha3Dq,DeltaAlpha3Dq,weight);
					POTScaledCC1pRecoDeltaPhiTPlot2D->Fill(true_DeltaPhiT,DeltaPhiT,weight);
					POTScaledCC1pRecoDeltaPhi3DPlot2D->Fill(true_DeltaPhi3D,DeltaPhi3D,weight);					

					POTScaledCC1pRecoDeltaPnPlot2D->Fill(true_Pn,reco_Pn,weight);

					// -----------------------------------------------------------------------------------------------------

					// True level energy reconstruction & Q2

					CC1pRecoECalPlot2D->Fill(true_ECal,ECal);

					POTScaledCC1pRecoECalPlot2D->Fill(true_ECal,ECal,weight);

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot2D[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(true_ECal,ECal,weight);


					CC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(true_ECal,ECal);
					POTScaledCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot2D[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(true_ECal,ECal,weight);
										
					CC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(true_ECal,ECal);	
					POTScaledCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot2D[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(true_ECal,ECal,weight);														
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

					// -------------------------------------------------------------------------------------------------------------------------

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
					     && true_TransMissMomentum > ArrayNBinsDeltaPT[0]
					     && true_DeltaAlphaT > ArrayNBinsDeltaAlphaT[0]
					     && true_DeltaPhiT > ArrayNBinsDeltaPhiT[0]

					     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
					     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
					     && True_CandidateMu_CosTheta->at(0) < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta] 
					     && True_CandidateP_CosTheta->at(0) < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
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

					NonCC1pRecoMuonMomentumPlot->Fill(reco_Pmu,weight);
					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					NonCC1pRecoProtonMomentumPlot->Fill(reco_Pp,weight);
					NonCC1pRecoProtonCosThetaPlot->Fill(reco_Pp_cos_theta,weight);
					NonCC1pRecoDeltaPnPlot->Fill(reco_Pn,weight);

					NonCC1pRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
					NonCC1pRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);
					NonCC1pRecoDeltaAlpha3DqPlot->Fill(DeltaAlpha3Dq,weight);
					NonCC1pRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
					NonCC1pRecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);					

					NonCC1pRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					NonCC1pRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					NonCC1pRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					NonCC1pRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 3D analysis in 1D grid

					SerialNonCC1pRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
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
					CCQERecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
					CCQERecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);					

					CCQERecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCQE Signal / Bkg	

					if (isSignal) {				

						CCQESignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCQESignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCQESignalRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCQESignalRecoECalPlot->Fill(ECal,weight);	
						CCQESignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					} else {

						CCQEBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCQEBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCQEBkgRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCQEBkgRecoECalPlot->Fill(ECal,weight);
						CCQEBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					}				

					//------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCQERecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCQERecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCQERecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 3D analysis in 1D grid

					SerialCCQERecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
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
					CCMECRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
					CCMECRecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);					

					CCMECRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCMEC Signal / Bkg	

					if (isSignal) {				

						CCMECSignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCMECSignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCMECSignalRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCMECSignalRecoECalPlot->Fill(ECal,weight);
						CCMECSignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);							

					} else {

						CCMECBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCMECBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCMECBkgRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCMECBkgRecoECalPlot->Fill(ECal,weight);
						CCMECBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);

					}						

					//----------------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCMECRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCMECRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 3D analysis in 1D grid

					SerialCCMECRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
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
					CCRESRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
					CCRESRecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);					

					CCRESRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCRES Signal / Bkg	

					if (isSignal) {				

						CCRESSignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCRESSignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCRESSignalRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCRESSignalRecoECalPlot->Fill(ECal,weight);
						CCRESSignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);							

					} else {

						CCRESBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCRESBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCRESBkgRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCRESBkgRecoECalPlot->Fill(ECal,weight);
						CCRESBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					}					

					//----------------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCRESRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCRESRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCRESRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);					

					// 3D analysis in 1D grid

					SerialCCRESRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
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
					CCDISRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
					CCDISRecoDeltaPhi3DPlot->Fill(DeltaPhi3D,weight);					

					CCDISRecoECalPlot->Fill(ECal,weight);

					//------------------------------//

					// CCDIS Signal / Bkg	

					if (isSignal) {				

						CCDISSignalRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCDISSignalRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCDISSignalRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCDISSignalRecoECalPlot->Fill(ECal,weight);
						CCDISSignalRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);							

					} else {

						CCDISBkgRecoDeltaPTPlot->Fill(TransMissMomentum,weight);
						CCDISBkgRecoDeltaAlphaTPlot->Fill(DeltaAlphaT,weight);									
						CCDISBkgRecoDeltaPhiTPlot->Fill(DeltaPhiT,weight);
						CCDISBkgRecoECalPlot->Fill(ECal,weight);
						CCDISBkgRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);						

					}					

					//----------------------------------------//

					// Ecal in DeltaPT & DeltaAlphaT bins

					CCDISRecoECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);	
					CCDISRecoECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
					CCDISRecoECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);						

					// 3D analysis in 1D grid

					SerialCCDISRecoECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
					SerialCCDISRecoECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
					SerialCCDISRecoECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);									

				}

				// --------------------------------------------------------------------------------------------------------------------------

				// Overlay particle breakdown using the Backtracker for PID studies

				// Sanity check

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
//		std::cout << "Txt info file: " << TxtName << std::endl << std::endl;

		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

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
