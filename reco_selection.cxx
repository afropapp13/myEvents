#define reco_selection_cxx
#include "reco_selection.h"
#include <TH2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TF1.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include "../myClasses/Tools.h"
#include "../myClasses/STV_Tools.h"

using namespace std;

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

	int NEventsPassingSelectionCuts = 0;
	int CC1pEventsPassingSelectionCuts = 0;
	int NonCC1pEventsPassingSelectionCuts = 0;
	
	//----------------------------------------//

	Tools tools;			

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

		// Muon CosTheta

		TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);	
		TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* CC1pRecoMuonCosThetaPlot2D = new TH2D("CC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH2D* POTScaledCC1pRecoMuonCosThetaPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
			ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* NonCC1pRecoMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCQERecoMuonCosThetaPlot = new TH1D("CCQERecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCMECRecoMuonCosThetaPlot = new TH1D("CCMECRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCRESRecoMuonCosThetaPlot = new TH1D("CCRESRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TH1D* CCDISRecoMuonCosThetaPlot = new TH1D("CCDISRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		
		//----------------------------------------//

		// Muon CosTheta Single Bin

		TH1D* RecoMuonCosThetaSingleBinPlot = new TH1D("RecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CC1pRecoMuonCosThetaSingleBinPlot = new TH1D("CC1pRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CC1pTrueMuonCosThetaSingleBinPlot = new TH1D("CC1pTrueMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH2D* CC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("CC1pRecoMuonCosThetaSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH2D* POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D = new TH2D("POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH1D* NonCC1pRecoMuonCosThetaSingleBinPlot = new TH1D("NonCC1pRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCQERecoMuonCosThetaSingleBinPlot = new TH1D("CCQERecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCMECRecoMuonCosThetaSingleBinPlot = new TH1D("CCMECRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCRESRecoMuonCosThetaSingleBinPlot = new TH1D("CCRESRecoMuonCosThetaSingleBinPlot","",1,0.,1.);
		TH1D* CCDISRecoMuonCosThetaSingleBinPlot = new TH1D("CCDISRecoMuonCosThetaSingleBinPlot","",1,0.,1.);

		//----------------------------------------//

		// ThetaZ

		TH1D* RecoThetaZPlot = new TH1D("RecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CC1pRecoThetaZPlot = new TH1D("CC1pRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);	
		TH1D* CC1pTrueThetaZPlot = new TH1D("CC1pTrueThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH2D* CC1pRecoThetaZPlot2D = new TH2D("CC1pRecoThetaZPlot2D",LabelXAxisThetaZ2D,NBinsThetaZ,
			ArrayNBinsThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH2D* POTScaledCC1pRecoThetaZPlot2D = new TH2D("POTScaledCC1pRecoThetaZPlot2D",LabelXAxisThetaZ2D,NBinsThetaZ,
			ArrayNBinsThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* NonCC1pRecoThetaZPlot = new TH1D("NonCC1pRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCQERecoThetaZPlot = new TH1D("CCQERecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCMECRecoThetaZPlot = new TH1D("CCMECRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCRESRecoThetaZPlot = new TH1D("CCRESRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
		TH1D* CCDISRecoThetaZPlot = new TH1D("CCDISRecoThetaZPlot",LabelXAxisThetaZ,NBinsThetaZ,ArrayNBinsThetaZ);
	
		//----------------------------------------//

		//CosThetaZ

		TH1D* RecoCosThetaZPlot = new TH1D("RecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH1D* CC1pRecoCosThetaZPlot = new TH1D("CC1pRecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);	
		TH1D* CC1pTrueCosThetaZPlot = new TH1D("CC1pTrueCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH2D* CC1pRecoCosThetaZPlot2D = new TH2D("CC1pRecoCosThetaZPlot2D",LabelXAxisCosThetaZ2D,NBinsCosThetaZ,
			ArrayNBinsCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH2D* POTScaledCC1pRecoCosThetaZPlot2D = new TH2D("POTScaledCC1pRecoCosThetaZPlot2D",LabelXAxisCosThetaZ2D,NBinsCosThetaZ,
			ArrayNBinsCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH1D* NonCC1pRecoCosThetaZPlot = new TH1D("NonCC1pRecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH1D* CCQERecoCosThetaZPlot = new TH1D("CCQERecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH1D* CCMECRecoCosThetaZPlot = new TH1D("CCMECRecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH1D* CCRESRecoCosThetaZPlot = new TH1D("CCRESRecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);
		TH1D* CCDISRecoCosThetaZPlot = new TH1D("CCDISRecoCosThetaZPlot",LabelXAxisCosThetaZ,NBinsCosThetaZ,ArrayNBinsCosThetaZ);

		//----------------------------------------//

		double ecal_diff_min = -0.2;
		double ecal_diff_max = 0.2;
		double ecal_reso_min = -50;
		double ecal_reso_max = 50;

		double thetaz_diff_min = -30;
		double thetaz_diff_max = 30;
		double thetaz_reso_min = -100;
		double thetaz_reso_max = 100;

		//----------------------------------------//

		TH1D* CC1pThetaZDiff_ECalSlicesPlot[TwoDNBinsECal];	
		TH1D* CC1pThetaZReso_ECalSlicesPlot[TwoDNBinsECal];	
		TH1D* CC1pTrueThetaZ_TrueECalSlicesPlot[TwoDNBinsECal];	
		TH1D* CC1pTrueThetaZ_TrueEnuSlicesPlot[TwoDNBinsECal];	
		TH1D* CC1pThetaZ_ECalSlicesPlot[TwoDNBinsECal];	
	
		//----------------------------------------//

		// ThetaZ in ECal slices 
		// Uncorrelated

		TH1D* RecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];
		TH1D* CC1pRecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];	
		TH1D* CC1pTrueThetaZ_ECalSlicesPlot[TwoDNBinsECal];
		TH2D* CC1pRecoThetaZ_ECalSlicesPlot2D[TwoDNBinsECal];
		TH2D* POTScaledCC1pRecoThetaZ_ECalSlicesPlot2D[TwoDNBinsECal];
		TH1D* NonCC1pRecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];
		TH1D* CCQERecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];
		TH1D* CCMECRecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];
		TH1D* CCRESRecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];
		TH1D* CCDISRecoThetaZ_ECalSlicesPlot[TwoDNBinsECal];

		// ThetaZ in ECal slices
		// Correlated

		TH1D* SerialRecoThetaZ_InECalPlot = new TH1D("RecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH1D* SerialCC1pRecoThetaZ_InECalPlot = new TH1D("CC1pRecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH1D* SerialCC1pTrueThetaZ_InECalPlot = new TH1D("CC1pTrueSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);		
		TH1D* SerialNonCC1pRecoThetaZ_InECalPlot = new TH1D("NonCC1pRecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH1D* SerialCCQERecoThetaZ_InECalPlot = new TH1D("CCQERecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH1D* SerialCCMECRecoThetaZ_InECalPlot = new TH1D("CCMECRecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);								
		TH1D* SerialCCRESRecoThetaZ_InECalPlot = new TH1D("CCRESRecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH1D* SerialCCDISRecoThetaZ_InECalPlot = new TH1D("CCDISRecoSerialThetaZ_ECalPlot",LabelXAxisThetaZ,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH2D* SerialCC1pRecoThetaZ_InECalPlot2D = new TH2D("CC1pRecoSerialThetaZ_ECalPlot2D",LabelXAxisThetaZ2D,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0],tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);
		TH2D* SerialPOTScaledCC1pRecoThetaZ_InECalPlot2D = new TH2D("POTScaledCC1pRecoSerialThetaZ_ECalPlot2D",LabelXAxisThetaZ2D,tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0],tools.Return2DNBins(TwoDArrayNBinsThetaZInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaZInECalSlices)[0]);				
		//
	
		// Loop over the ECal slices
		
		for (int iecal = 0; iecal < TwoDNBinsECal; iecal++) {

			CC1pThetaZDiff_ECalSlicesPlot[iecal] = new TH1D("CC1pThetaZDiff_ECalSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_ECalSlicesPlot[iecal] = new TH1D("CC1pThetaZReso_ECalSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
			CC1pThetaZ_ECalSlicesPlot[iecal] = new TH1D("CC1pThetaZ_ECalSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";#theta_{z}^{reco} [deg]",30,0,180);
			CC1pTrueThetaZ_TrueECalSlicesPlot[iecal] = new TH1D("CC1pTrueThetaZ_TrueECalSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";#theta_{z}^{true} [deg]",30,0,180);
			CC1pTrueThetaZ_TrueEnuSlicesPlot[iecal] = new TH1D("CC1pTrueThetaZ_TrueEnuSlices" + tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1]) +"Plot",";#theta_{z}^{true} [deg]",30,0,180);

			TString ThetaZTwoDInECalLabel = "ThetaZ_ECal_"+tools.ConvertToString(TwoDArrayNBinsECal[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[iecal+1])+"Plot";			
			RecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("Reco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			CC1pRecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("CC1pReco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);	
			CC1pTrueThetaZ_ECalSlicesPlot[iecal] = new TH1D("CC1pTrue" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			CC1pRecoThetaZ_ECalSlicesPlot2D[iecal] = new TH2D("CC1pReco" + ThetaZTwoDInECalLabel + "2D",LabelXAxisThetaZ2D,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0],TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			POTScaledCC1pRecoThetaZ_ECalSlicesPlot2D[iecal] = new TH2D("POTScaledCC1pReco" + ThetaZTwoDInECalLabel + "2D",LabelXAxisThetaZ2D,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0],TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			NonCC1pRecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("NonCC1pReco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			CCQERecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("CCQEReco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			CCMECRecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("CCMECReco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			CCRESRecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("CCRESReco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
			CCDISRecoThetaZ_ECalSlicesPlot[iecal] = new TH1D("CCDISReco" + ThetaZTwoDInECalLabel,LabelXAxisThetaZ,TwoDArrayNBinsThetaZInECalSlices[iecal].size()-1,&TwoDArrayNBinsThetaZInECalSlices[iecal][0]);
		
		} // End of the loop over ECal slices

		//----------------------------------------//

		TH1D* CC1pECalDiff_MuonMomentumSlicesPlot[TwoDNBinsMuonMomentum];	
		TH1D* CC1pECalReso_MuonMomentumSlicesPlot[TwoDNBinsMuonMomentum];	

		TH1D* CC1pThetaZDiff_MuonMomentumSlicesPlot[TwoDNBinsMuonMomentum];	
		TH1D* CC1pThetaZReso_MuonMomentumSlicesPlot[TwoDNBinsMuonMomentum];	

		// Loop over the MuonMomentum slices
		
		for (int iecal = 0; iecal < TwoDNBinsMuonMomentum; iecal++) {

			CC1pThetaZDiff_MuonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZDiff_MuonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_MuonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZReso_MuonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
	
			CC1pECalDiff_MuonMomentumSlicesPlot[iecal] = new TH1D("CC1pECalDiff_MuonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_MuonMomentumSlicesPlot[iecal] = new TH1D("CC1pECalReso_MuonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[iecal+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over MuonMomentum slices
	
		//----------------------------------------//

		TH1D* CC1pECalDiff_MuonCosThetaSlicesPlot[TwoDNBinsMuonCosTheta];	
		TH1D* CC1pECalReso_MuonCosThetaSlicesPlot[TwoDNBinsMuonCosTheta];	

		// Loop over the MuonCosTheta slices
		
		for (int i = 0; i < TwoDNBinsMuonCosTheta; i++) {

			CC1pECalDiff_MuonCosThetaSlicesPlot[i] = new TH1D("CC1pECalDiff_MuonCosThetaSlices" + tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[i])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[i+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_MuonCosThetaSlicesPlot[i] = new TH1D("CC1pECalReso_MuonCosThetaSlices" + tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[i])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[i+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over MuonMomentum slices
	
		//----------------------------------------//

		TH1D* CC1pECalDiff_ProtonMomentumSlicesPlot[TwoDNBinsProtonMomentum];	
		TH1D* CC1pECalReso_ProtonMomentumSlicesPlot[TwoDNBinsProtonMomentum];	

		TH1D* CC1pThetaZDiff_ProtonMomentumSlicesPlot[TwoDNBinsProtonMomentum];	
		TH1D* CC1pThetaZReso_ProtonMomentumSlicesPlot[TwoDNBinsProtonMomentum];	

		// Loop over the ProtonMomentum slices
		
		for (int iecal = 0; iecal < TwoDNBinsProtonMomentum; iecal++) {

			CC1pThetaZDiff_ProtonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZDiff_ProtonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_ProtonMomentumSlicesPlot[iecal] = new TH1D("CC1pThetaZReso_ProtonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
	
			CC1pECalDiff_ProtonMomentumSlicesPlot[iecal] = new TH1D("CC1pECalDiff_ProtonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_ProtonMomentumSlicesPlot[iecal] = new TH1D("CC1pECalReso_ProtonMomentumSlices" + tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[iecal+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over ProtonMomentum slices
	
		//----------------------------------------//

		TH1D* CC1pECalDiff_ProtonCosThetaSlicesPlot[TwoDNBinsProtonCosTheta];	
		TH1D* CC1pECalReso_ProtonCosThetaSlicesPlot[TwoDNBinsProtonCosTheta];	

		// Loop over the MuonCosTheta slices
		
		for (int i = 0; i < TwoDNBinsProtonCosTheta; i++) {

			CC1pECalDiff_ProtonCosThetaSlicesPlot[i] = new TH1D("CC1pECalDiff_ProtonCosThetaSlices" + tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[i])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[i+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_ProtonCosThetaSlicesPlot[i] = new TH1D("CC1pECalReso_ProtonCosThetaSlices" + tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[i])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[i+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over MuonMomentum slices
	
		//----------------------------------------//

		TH1D* CC1pECalDiff_DeltaPTSlicesPlot[TwoDNBinsDeltaPT];	
		TH1D* CC1pECalReso_DeltaPTSlicesPlot[TwoDNBinsDeltaPT];	

		TH1D* CC1pThetaZDiff_DeltaPTSlicesPlot[TwoDNBinsDeltaPT];	
		TH1D* CC1pThetaZReso_DeltaPTSlicesPlot[TwoDNBinsDeltaPT];	

		// Loop over the DeltaPT slices
		
		for (int i = 0; i < TwoDNBinsDeltaPT; i++) {

			CC1pThetaZDiff_DeltaPTSlicesPlot[i] = new TH1D("CC1pThetaZDiff_DeltaPTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[i+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_DeltaPTSlicesPlot[i] = new TH1D("CC1pThetaZReso_DeltaPTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[i+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
	
			CC1pECalDiff_DeltaPTSlicesPlot[i] = new TH1D("CC1pECalDiff_DeltaPTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[i+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_DeltaPTSlicesPlot[i] = new TH1D("CC1pECalReso_DeltaPTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[i+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over DeltaPT slices
	
		//----------------------------------------//

		TH1D* CC1pECalDiff_DeltaAlphaTSlicesPlot[TwoDNBinsDeltaAlphaT];	
		TH1D* CC1pECalReso_DeltaAlphaTSlicesPlot[TwoDNBinsDeltaAlphaT];	

		TH1D* CC1pThetaZDiff_DeltaAlphaTSlicesPlot[TwoDNBinsDeltaAlphaT];	
		TH1D* CC1pThetaZReso_DeltaAlphaTSlicesPlot[TwoDNBinsDeltaAlphaT];	

		// Loop over the DeltaAlphaT slices
		
		for (int i = 0; i < TwoDNBinsDeltaAlphaT; i++) {

			CC1pThetaZDiff_DeltaAlphaTSlicesPlot[i] = new TH1D("CC1pThetaZDiff_DeltaAlphaTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_DeltaAlphaTSlicesPlot[i] = new TH1D("CC1pThetaZReso_DeltaAlphaTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
	
			CC1pECalDiff_DeltaAlphaTSlicesPlot[i] = new TH1D("CC1pECalDiff_DeltaAlphaTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_DeltaAlphaTSlicesPlot[i] = new TH1D("CC1pECalReso_DeltaAlphaTSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[i+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over DeltaAlphaT slices

		//----------------------------------------//

		TH1D* CC1pECalDiff_DeltaPnSlicesPlot[TwoDNBinsDeltaPn];	
		TH1D* CC1pECalReso_DeltaPnSlicesPlot[TwoDNBinsDeltaPn];	

		TH1D* CC1pThetaZDiff_DeltaPnSlicesPlot[TwoDNBinsDeltaPn];	
		TH1D* CC1pThetaZReso_DeltaPnSlicesPlot[TwoDNBinsDeltaPn];	

		// Loop over the DeltaPn slices
		
		for (int i = 0; i < TwoDNBinsDeltaPn; i++) {

			CC1pThetaZDiff_DeltaPnSlicesPlot[i] = new TH1D("CC1pThetaZDiff_DeltaPnSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPn[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[i+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_DeltaPnSlicesPlot[i] = new TH1D("CC1pThetaZReso_DeltaPnSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPn[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[i+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
	
			CC1pECalDiff_DeltaPnSlicesPlot[i] = new TH1D("CC1pECalDiff_DeltaPnSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPn[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[i+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_DeltaPnSlicesPlot[i] = new TH1D("CC1pECalReso_DeltaPnSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaPn[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[i+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over DeltaPn slices
	
		//----------------------------------------//

		TH1D* CC1pECalDiff_DeltaAlpha3DSlicesPlot[TwoDNBinsDeltaAlpha3D];	
		TH1D* CC1pECalReso_DeltaAlpha3DSlicesPlot[TwoDNBinsDeltaAlpha3D];	

		TH1D* CC1pThetaZDiff_DeltaAlpha3DSlicesPlot[TwoDNBinsDeltaAlpha3D];	
		TH1D* CC1pThetaZReso_DeltaAlpha3DSlicesPlot[TwoDNBinsDeltaAlpha3D];	

		// Loop over the DeltaAlpha3D slices
		
		for (int i = 0; i < TwoDNBinsDeltaAlpha3D; i++) {

			CC1pThetaZDiff_DeltaAlpha3DSlicesPlot[i] = new TH1D("CC1pThetaZDiff_DeltaAlpha3DSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i+1]) +"Plot",";#theta_{z}^{reco} - #theta_{z}^{true} [deg]",31,thetaz_diff_min,thetaz_diff_max);
			CC1pThetaZReso_DeltaAlpha3DSlicesPlot[i] = new TH1D("CC1pThetaZReso_DeltaAlpha3DSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i+1]) +"Plot",";(#theta_{z}^{reco} - #theta_{z}^{true})/#theta_{z}^{true} [%]",51,thetaz_reso_min,thetaz_reso_max);
	
			CC1pECalDiff_DeltaAlpha3DSlicesPlot[i] = new TH1D("CC1pECalDiff_DeltaAlpha3DSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i+1]) +"Plot",";E_{Cal}^{reco} - E_{Cal}^{true} [GeV]",31,ecal_diff_min,ecal_diff_max);
			CC1pECalReso_DeltaAlpha3DSlicesPlot[i] = new TH1D("CC1pECalReso_DeltaAlpha3DSlices" + tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3D[i+1]) +"Plot",";(E_{Cal}^{reco} - E_{Cal}^{true})/E_{Cal}^{true} [%]",51,ecal_reso_min,ecal_reso_max);
		
		} // End of the loop over DeltaAlphaT slices


		//----------------------------------------//

		// 2D plots

		// Reco

                TH2D* POTScaledCC1pRecoThetaZRecoECalPlot2D = new TH2D("POTScaledCC1pRecoThetaZRecoECalPlot2D",";#theta_{z}^{reco} [deg];E_{Cal}^{reco} [GeV]",
                        NBinsThetaZ,ArrayNBinsThetaZ,20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal]);

                TH2D* POTScaledCC1pRecoECalTrueECalPlot2D = new TH2D("POTScaledCC1pRecoECalTrueECalPlot2D",";E_{Cal}^{true} [GeV];E_{Cal}^{reco} [GeV]",
                        20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal],20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal]);

                TH2D* POTScaledCC1pRecoECalTrueEnuPlot2D = new TH2D("POTScaledCC1pRecoECalTrueEnuPlot2D",";E_{#nu}^{true} [GeV];E_{Cal}^{reco} [GeV]",
                        20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal],20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal]);

		// True

                TH2D* POTScaledCC1pTrueThetaZTrueECalPlot2D = new TH2D("POTScaledCC1pTrueThetaZTrueECalPlot2D",";#theta_{z}^{true} [deg];E_{Cal}^{true} [GeV]",
                        NBinsThetaZ,ArrayNBinsThetaZ,20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal]);

                TH2D* POTScaledCC1pTrueThetaZTrueEnuPlot2D = new TH2D("POTScaledCC1pTrueThetaZTrueEnuPlot2D",";#theta_{z}^{true} [deg];E_{#nu}^{true} [GeV]",
                        NBinsThetaZ,ArrayNBinsThetaZ,20,ArrayNBinsECal[0],ArrayNBinsECal[NBinsECal]);

		//----------------------------------------//
		//----------------------------------------//

		// Loop over the events

		cout << nentries << " events included in the file" << endl;

		for (Long64_t jentry=0; jentry<nentries;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

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

				// For detector variations runs 1-3, the eventweight weights are -1., set them back to 1.
				if (Weight == -1.) { Weight = 1.; }
				if (T2KWeight == -1.) { T2KWeight = 1.; }
			
				// For detector variations runs 4-5, the event weights are NOT -1
				// However setting the weights to 1 for consistency
				if (
					string(fWhichSample).find("CV") != std::string::npos || 
					string(fWhichSample).find("CVextra") != std::string::npos || 
					string(fWhichSample).find("LYDown") != std::string::npos || 
					string(fWhichSample).find("LYRayleigh") != std::string::npos || 
					string(fWhichSample).find("LYAttenuation") != std::string::npos || 
					string(fWhichSample).find("SCE") != std::string::npos || 
					string(fWhichSample).find("Recombination2") != std::string::npos || 
					string(fWhichSample).find("X") != std::string::npos || 
					string(fWhichSample).find("YZ") != std::string::npos || 
					string(fWhichSample).find("ThetaXZ") != std::string::npos || 
					string(fWhichSample).find("ThetaYZ") != std::string::npos
				) {

					Weight = 1.;
					T2KWeight = 1.;
				
				}

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

			double reco_p_LLR_Score = CandidateP_LLR_PID->at(0);

			// -----------------------------------------------------------------------------------------------------------------------------

			// Kinematic variables
	
			double ECal = Reco_ECal->at(0);
			double ThetaZ = Reco_ThetaZ->at(0); // deg
			double CosThetaZ = TMath::Cos(ThetaZ * TMath::Pi() / 180.);

			double DeltaPT = Reco_Pt->at(0);
			double DeltaAlphaT = Reco_DeltaAlphaT->at(0);
			double DeltaPn = Reco_Pn->at(0);
			double DeltaAlpha3D = Reco_DeltaAlpha3Dq->at(0);

			// -------------------------------------------------------------------------------------------------------------------------
			
			STV_Tools reco_stv_tool(TVector3CandidateMuon,TVector3CandidateProton,reco_Emu,reco_Ep);
	
			// Redefinition if P_p < 0.5 where the biases have been observed

			if (reco_Pp < 0.5) { 

				reco_Pp = ( 1.-0.01*fPP->Eval(reco_Pp) ) * reco_Pp ;
				TVector3CandidateProton.SetMag(reco_Pp);
				reco_Ep = TMath::Sqrt( reco_Pp*reco_Pp + ProtonMass_GeV*ProtonMass_GeV );

				ECal = reco_stv_tool.ReturnECal();
				ThetaZ = reco_stv_tool.ReturnThetaZ(); // deg
				CosThetaZ = TMath::Cos(ThetaZ * TMath::Pi() / 180.);

				DeltaPT = reco_stv_tool.ReturnPt();
				DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
				DeltaPn = reco_stv_tool.ReturnPn();
				DeltaAlpha3D = reco_stv_tool.ReturnDeltaAlpha3Dq();

			}

			// Underflow / overflow
			if (ThetaZ < ArrayNBinsThetaZ[0]) { ThetaZ = (ArrayNBinsThetaZ[0] + ArrayNBinsThetaZ[1])/2.; }
			if (ThetaZ > ArrayNBinsThetaZ[NBinsThetaZ]) { ThetaZ = (ArrayNBinsThetaZ[NBinsThetaZ] + ArrayNBinsThetaZ[NBinsThetaZ-1])/2.; }

			if (DeltaPT > ArrayNBinsDeltaPT[NBinsDeltaPT]) { DeltaPT = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }

			if (DeltaPn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { DeltaPn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

			if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
			if (ECal < ArrayNBinsECal[0]) { ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }	

			//----------------------------------------//

			//Reco  2D indices

			int ECalTwoDIndex = tools.ReturnIndex(ECal, TwoDArrayNBinsECal);
			int SerialThetaZInECalIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaZInECalSlices,ECalTwoDIndex,ThetaZ);

			//----------------------------------------//

			// Selection Cuts

			bool PassedSelection = true;

			for (int i = 0; i < NCuts; i++) {

				if (VectorCuts[i] == "_PID_NuScore" && !(reco_p_LLR_Score < ProtonLLRPIDScore) ) 
					{ PassedSelection = false; }

				if (VectorCuts[i] == "_PID_NuScore_CRT" && !(crtveto == 0) ) 
					{ PassedSelection = false; }


			}

			if (PassedSelection == false) { continue; }

			// -------------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Make sure that the same events fill the same plots

			if (reco_Pmu < ArrayNBinsMuonMomentum[0]) { continue; }
			if (reco_Pp < ArrayNBinsProtonMomentum[0]) { continue; }

			// --------------------------------------------------------------------------------------------------------------------

			if (reco_Pmu > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }
			if (reco_Pp > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }

			NEventsPassingSelectionCuts++;

			// --------------------------------------------------------------------------------------------------------------------------------

			int genie_mode = -1;

			double true_ECal = -1;
			double true_ThetaZ = -1;
			double true_CosThetaZ = -1;

			double true_DeltaPT = -1;
			double true_DeltaAlphaT = -1;
			double true_DeltaPn = -1;
			double true_DeltaAlpha3D = -1;

			//----------------------------------------//

			//True 2D indices

			int TrueECalTwoDIndex = -1;
			int TrueSerialThetaZInECalIndex = -1;

			//----------------------------------------//

			// Only for MC to obtain true vales			
			
			if (
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

				true_ECal = True_ECal->at(0);
				true_ThetaZ = True_ThetaZ->at(0); // deg
				true_CosThetaZ = TMath::Cos(true_ThetaZ * TMath::Pi() / 180.);
	
				true_DeltaPT = True_Pt->at(0);
				true_DeltaAlphaT = True_DeltaAlphaT->at(0);
				true_DeltaPn = True_Pn->at(0);
				true_DeltaAlpha3D = True_DeltaAlpha3Dq->at(0);

                        	// Underflow / overflow
                                if (true_ThetaZ < ArrayNBinsThetaZ[0]) { true_ThetaZ = (ArrayNBinsThetaZ[0] + ArrayNBinsThetaZ[1])/2.; }
                                if (true_ThetaZ > ArrayNBinsThetaZ[NBinsThetaZ]) { true_ThetaZ = (ArrayNBinsThetaZ[NBinsThetaZ] + ArrayNBinsThetaZ[NBinsThetaZ-1])/2.; }
                        
				if (true_DeltaPT > ArrayNBinsDeltaPT[NBinsDeltaPT]) { true_DeltaPT = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }

				if (true_DeltaPn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { true_DeltaPn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

	                        if (true_ECal > ArrayNBinsECal[NBinsECal]) { true_ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
        	                if (true_ECal < ArrayNBinsECal[0]) { true_ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }
	
                                if (True_Ev > ArrayNBinsECal[NBinsECal]) { True_Ev = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
                                if (True_Ev < ArrayNBinsECal[0]) { True_Ev = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }

				TrueECalTwoDIndex = tools.ReturnIndex(true_ECal, TwoDArrayNBinsECal);
				TrueSerialThetaZInECalIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaZInECalSlices,TrueECalTwoDIndex,true_ThetaZ);

			} // End of if statement: Only for MC to obtain true vales

			//----------------------------------------//

			RecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
			RecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
			RecoThetaZPlot->Fill(ThetaZ,weight); // deg
			RecoCosThetaZPlot->Fill(CosThetaZ,weight); // deg

			// 2D analysis
			RecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
			SerialRecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);

			//------------------------------//

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

				     && True_CandidateMu_P->at(0) < ArrayNBinsMuonMomentum[NBinsMuonMomentum] 
				     && True_CandidateP_P->at(0) < ArrayNBinsProtonMomentum[NBinsProtonMomentum]

				) {

					// --------------------------------------------------------------------------------------------------
				
					CC1pEventsPassingSelectionCuts++;
					isSignal = true;

					// ---------------------------------------------------------------------------------------------------------------------------
					// ---------------------------------------------------------------------------------------------------------------------------
					// 1D Plots using True level info for selected CC1p events  

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(True_CandidateMu_P->at(0),2.) );
					double true_Nu = True_Ev - true_MuonEnergy;

					CC1pTrueMuonCosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
					CC1pTrueMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CC1pTrueThetaZPlot->Fill(true_ThetaZ,weight);
					CC1pTrueCosThetaZPlot->Fill(true_CosThetaZ,weight);

					// 2D analysis
					CC1pTrueThetaZ_ECalSlicesPlot[TrueECalTwoDIndex]->Fill(true_ThetaZ,weight);
					SerialCC1pTrueThetaZ_InECalPlot->Fill(TrueSerialThetaZInECalIndex,weight);

					//----------------------------------------//

					// 1D Reco Plots for the selected CC1p events 

					CC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CC1pRecoThetaZPlot->Fill(ThetaZ,weight);
					CC1pRecoCosThetaZPlot->Fill(CosThetaZ,weight);
	
					// 2D analysis
					CC1pRecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
					SerialCC1pRecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);

					//------------------------------//

					CC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
					CC1pRecoMuonCosThetaSingleBinPlot2D->Fill(0.5,0.5);
					CC1pRecoThetaZPlot2D->Fill(true_ThetaZ,ThetaZ);
					CC1pRecoCosThetaZPlot2D->Fill(true_CosThetaZ,CosThetaZ);

					// 2D analysis
					CC1pRecoThetaZ_ECalSlicesPlot2D[ECalTwoDIndex]->Fill(true_ThetaZ,ThetaZ,weight);
					SerialCC1pRecoThetaZ_InECalPlot2D->Fill(TrueSerialThetaZInECalIndex,SerialThetaZInECalIndex,weight);

					POTScaledCC1pRecoMuonCosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta,weight);
					POTScaledCC1pRecoMuonCosThetaSingleBinPlot2D->Fill(0.5,0.5,weight);
					POTScaledCC1pRecoThetaZPlot2D->Fill(true_ThetaZ,ThetaZ,weight);
					POTScaledCC1pRecoCosThetaZPlot2D->Fill(true_CosThetaZ,CosThetaZ,weight);

					// 2D analysis
					POTScaledCC1pRecoThetaZ_ECalSlicesPlot2D[ECalTwoDIndex]->Fill(true_ThetaZ,ThetaZ,weight);
					SerialPOTScaledCC1pRecoThetaZ_InECalPlot2D->Fill(TrueSerialThetaZInECalIndex,SerialThetaZInECalIndex,weight);
					
					//------------------------------//

					// Atmospherics
				
					double diff = ThetaZ - true_ThetaZ; // deg
					double reso = diff / true_ThetaZ * 100.; // %

					double ECal_diff = ECal - true_ECal; // GeV
					double ECal_reso = ECal_diff / true_ECal * 100.; // %

					//------------------------------//

					// ECal slices

					int ECalTwoDIndex = tools.ReturnIndex(ECal, TwoDArrayNBinsECal);
					int TrueECalTwoDIndex = tools.ReturnIndex(true_ECal, TwoDArrayNBinsECal);
					int TrueEnuTwoDIndex = tools.ReturnIndex(True_Ev, TwoDArrayNBinsECal);

					CC1pThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);	
					CC1pTrueThetaZ_TrueECalSlicesPlot[TrueECalTwoDIndex]->Fill(true_ThetaZ,weight);	
					CC1pTrueThetaZ_TrueEnuSlicesPlot[TrueEnuTwoDIndex]->Fill(true_ThetaZ,weight);	
					CC1pThetaZReso_ECalSlicesPlot[ECalTwoDIndex]->Fill(reso,weight);
					CC1pThetaZDiff_ECalSlicesPlot[ECalTwoDIndex]->Fill(diff,weight);	
	
					//------------------------------//

					// MuonMomentum slices

					int MuonMomentumTwoDIndex = tools.ReturnIndex(reco_Pmu, TwoDArrayNBinsMuonMomentum);

					CC1pThetaZDiff_MuonMomentumSlicesPlot[MuonMomentumTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_MuonMomentumSlicesPlot[MuonMomentumTwoDIndex]->Fill(reso,weight);

					CC1pECalDiff_MuonMomentumSlicesPlot[MuonMomentumTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_MuonMomentumSlicesPlot[MuonMomentumTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// MuonCosTheta slices

					int MuonCosThetaTwoDIndex = tools.ReturnIndex(reco_Pmu_cos_theta,TwoDArrayNBinsMuonCosTheta);

					CC1pECalDiff_MuonCosThetaSlicesPlot[MuonCosThetaTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_MuonCosThetaSlicesPlot[MuonCosThetaTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// ProtonMomentum slices

					int ProtonMomentumTwoDIndex = tools.ReturnIndex(reco_Pp, TwoDArrayNBinsProtonMomentum);

					CC1pThetaZDiff_ProtonMomentumSlicesPlot[ProtonMomentumTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_ProtonMomentumSlicesPlot[ProtonMomentumTwoDIndex]->Fill(reso,weight);

					CC1pECalDiff_ProtonMomentumSlicesPlot[ProtonMomentumTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_ProtonMomentumSlicesPlot[ProtonMomentumTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// ProtonCosTheta slices

					int ProtonCosThetaTwoDIndex = tools.ReturnIndex(reco_Pp_cos_theta,TwoDArrayNBinsProtonCosTheta);

					CC1pECalDiff_ProtonCosThetaSlicesPlot[ProtonCosThetaTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_ProtonCosThetaSlicesPlot[ProtonCosThetaTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// DeltaPT slices

					int DeltaPTTwoDIndex = tools.ReturnIndex(DeltaPT, TwoDArrayNBinsDeltaPT);

					CC1pThetaZDiff_DeltaPTSlicesPlot[DeltaPTTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_DeltaPTSlicesPlot[DeltaPTTwoDIndex]->Fill(reso,weight);

					CC1pECalDiff_DeltaPTSlicesPlot[DeltaPTTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_DeltaPTSlicesPlot[DeltaPTTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// DeltaAlphaT slices

					int DeltaAlphaTTwoDIndex = tools.ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);

					CC1pThetaZDiff_DeltaAlphaTSlicesPlot[DeltaAlphaTTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_DeltaAlphaTSlicesPlot[DeltaAlphaTTwoDIndex]->Fill(reso,weight);
	
					CC1pECalDiff_DeltaAlphaTSlicesPlot[DeltaAlphaTTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_DeltaAlphaTSlicesPlot[DeltaAlphaTTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// DeltaPn slices

					int DeltaPnTwoDIndex = tools.ReturnIndex(DeltaPn, TwoDArrayNBinsDeltaPn);

					CC1pThetaZDiff_DeltaPnSlicesPlot[DeltaPnTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_DeltaPnSlicesPlot[DeltaPnTwoDIndex]->Fill(reso,weight);

					CC1pECalDiff_DeltaPnSlicesPlot[DeltaPnTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_DeltaPnSlicesPlot[DeltaPnTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// DeltaAlpha3D slices

					int DeltaAlpha3DTwoDIndex = tools.ReturnIndex(DeltaAlpha3D, TwoDArrayNBinsDeltaAlpha3D);

					CC1pThetaZDiff_DeltaAlpha3DSlicesPlot[DeltaAlpha3DTwoDIndex]->Fill(diff,weight);	
					CC1pThetaZReso_DeltaAlpha3DSlicesPlot[DeltaAlpha3DTwoDIndex]->Fill(reso,weight);

					CC1pECalDiff_DeltaAlpha3DSlicesPlot[DeltaAlpha3DTwoDIndex]->Fill(ECal_diff,weight);	
					CC1pECalReso_DeltaAlpha3DSlicesPlot[DeltaAlpha3DTwoDIndex]->Fill(ECal_reso,weight);

					//------------------------------//

					// 2D plots

					// Reco

					POTScaledCC1pRecoThetaZRecoECalPlot2D->Fill(ThetaZ,ECal,weight);
					POTScaledCC1pRecoECalTrueECalPlot2D->Fill(true_ECal,ECal,weight);
					POTScaledCC1pRecoECalTrueEnuPlot2D->Fill(True_Ev,ECal,weight);
				
					// True
					
					POTScaledCC1pTrueThetaZTrueECalPlot2D->Fill(true_ThetaZ,true_ECal,weight);
					POTScaledCC1pTrueThetaZTrueEnuPlot2D->Fill(true_ThetaZ,True_Ev,weight);

					//------------------------------//


				} // End of the CC1p signal

				// -------------------------------------------------------------------------------------------------------------

				// Non-CC1p beam related background or EXT BNB

				else {

					NonCC1pEventsPassingSelectionCuts++;

					NonCC1pRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					NonCC1pRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					NonCC1pRecoThetaZPlot->Fill(ThetaZ,weight);
					NonCC1pRecoCosThetaZPlot->Fill(CosThetaZ,weight);

					// 2D analysis
					NonCC1pRecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
					SerialNonCC1pRecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);

					//------------------------------//

				} // End of the Non-CC1p beam related background

				// -------------------------------------------------------------------------------------------------------------------------
				// ------------------------------------------------------------------------------------------------------------------------

				// CCQE

				if (genie_mode == 0) {

					CCQERecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCQERecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCQERecoThetaZPlot->Fill(ThetaZ,weight);
					CCQERecoCosThetaZPlot->Fill(CosThetaZ,weight);

					// 2D analysis
					CCQERecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
					SerialCCQERecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);

				} // End of CCQE selection

				// ------------------------------------------------------------------------------------------------------------------------

				// CCMEC

				if (genie_mode == 10) {

					CCMECRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCMECRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCMECRecoThetaZPlot->Fill(ThetaZ,weight);
					CCMECRecoCosThetaZPlot->Fill(CosThetaZ,weight);

					// 2D analysis
					CCMECRecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
					SerialCCMECRecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);
		
				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCRES

				if (genie_mode == 1) {

					CCRESRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCRESRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCRESRecoThetaZPlot->Fill(ThetaZ,weight);
					CCRESRecoCosThetaZPlot->Fill(CosThetaZ,weight);
	
					// 2D analysis
					CCRESRecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
					SerialCCRESRecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);
	
				}

				// -------------------------------------------------------------------------------------------------------------------------

				// CCDIS

				if (genie_mode == 2) {

					CCDISRecoMuonCosThetaPlot->Fill(reco_Pmu_cos_theta,weight);
					CCDISRecoMuonCosThetaSingleBinPlot->Fill(0.5,weight);
					CCDISRecoThetaZPlot->Fill(ThetaZ,weight);
					CCDISRecoCosThetaZPlot->Fill(CosThetaZ,weight);

					// 2D analysis
					CCRESRecoThetaZ_ECalSlicesPlot[ECalTwoDIndex]->Fill(ThetaZ,weight);
					SerialCCRESRecoThetaZ_InECalPlot->Fill(SerialThetaZInECalIndex,weight);

				}

				// --------------------------------------------------------------------------------------------------------------------------

				// Overlay particle breakdown using the Backtracker for PID studies

				if ( !(CandidateMu_MCParticle_Pdg->size() > 0 && CandidateP_MCParticle_Pdg->size() > 0) ) { 

					cout << "muon candidate true pdg = " << CandidateMu_MCParticle_Pdg->at(0) << endl;
					cout << "proton candidate true pdg = " << CandidateP_MCParticle_Pdg->at(0) << endl;

				}

			} // End of the Overlay case and the breakdown into CC1p/NonCC1p & QE,MEC,RES,DIS

			// -------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the events

		std::cout << std::endl << "Created file: " << FileName << std::endl << std::endl;

		std::cout << "---------------------------------------------------------------------" << std::endl << std::endl;

		//----------------------------------------//

		cout << fWhichSample << "  " << NEventsPassingSelectionCuts << " events passing selection criteria" << endl << endl;

		//----------------------------------------//	

		file->cd();
		file->Write();
		file->Close();
		fFile->Close();

		//----------------------------------------//	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program
