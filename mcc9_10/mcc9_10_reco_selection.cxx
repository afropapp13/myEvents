#define mcc9_10_reco_selection_cxx
#include "mcc9_10_reco_selection.h"
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

#include "../../../generators/Tools.h"

using namespace std;

//----------------------------------------//

TString ToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

void mcc9_10_reco_selection::Loop() {

	//----------------------------------------//

	Tools tools;			

	//----------------------------------------//
		
	TString Cuts = "_nocuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	VectorCuts.push_back("");
	//VectorCuts.push_back("_PID_NuScore");
	//VectorCuts.push_back("_CRT");

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

		TString FileName = event_selection_file_path+Cuts+"/"+fTune+"ncpi0_"+fWhichSample+Extension+Cuts+".root";
		TFile* file = new TFile(FileName,"recreate");
		std::cout << std::endl << "Creating a new file: " << FileName << std::endl << std::endl << std::endl;

		//----------------------------------------//

		// Txt file to keep track of the run/subrun/event of the candidate events

		TString RunTxtName = event_selection_file_path+Cuts+"/"+fTune+"ncpi0_"+fWhichSample+Extension+Cuts+".txt";
		ofstream myRunTxtFile;
		myRunTxtFile.open(RunTxtName);
		myRunTxtFile << std::fixed << std::setprecision(2);
		myRunTxtFile << fWhichSample << endl << endl;			

		//----------------------------------------//

		// kine_pio_vtx_dis

		TH1D* Recokine_pio_vtx_disPlot = new TH1D("Recokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* NCCOHRecokine_pio_vtx_disPlot = new TH1D("NCCOHRecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);	
		TH1D* NCCOHTruekine_pio_vtx_disPlot = new TH1D("NCCOHTruekine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH2D* NCCOHRecokine_pio_vtx_disPlot2D = new TH2D("NCCOHRecokine_pio_vtx_disPlot2D",LabelXAxiskine_pio_vtx_dis2D,NBinskine_pio_vtx_dis,
			  min_kine_pio_vtx_dis,max_kine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH2D* POTScaledNCCOHRecokine_pio_vtx_disPlot2D = new TH2D("POTScaledNCCOHRecokine_pio_vtx_disPlot2D",LabelXAxiskine_pio_vtx_dis2D,NBinskine_pio_vtx_dis,
			  min_kine_pio_vtx_dis,max_kine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* NonNCCOHRecokine_pio_vtx_disPlot = new TH1D("NonNCCOHRecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* QERecokine_pio_vtx_disPlot = new TH1D("QERecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* MECRecokine_pio_vtx_disPlot = new TH1D("MECRecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* RESRecokine_pio_vtx_disPlot = new TH1D("RESRecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* DISRecokine_pio_vtx_disPlot = new TH1D("DISRecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);
		TH1D* COHRecokine_pio_vtx_disPlot = new TH1D("COHRecokine_pio_vtx_disPlot",LabelXAxiskine_pio_vtx_dis,NBinskine_pio_vtx_dis,min_kine_pio_vtx_dis,max_kine_pio_vtx_dis);	

		//----------------------------------------//

		// kine_pio_flag

		TH1D* Recokine_pio_flagPlot = new TH1D("Recokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* NCCOHRecokine_pio_flagPlot = new TH1D("NCCOHRecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);	
		TH1D* NCCOHTruekine_pio_flagPlot = new TH1D("NCCOHTruekine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH2D* NCCOHRecokine_pio_flagPlot2D = new TH2D("NCCOHRecokine_pio_flagPlot2D",LabelXAxiskine_pio_flag2D,NBinskine_pio_flag,
			  min_kine_pio_flag,max_kine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH2D* POTScaledNCCOHRecokine_pio_flagPlot2D = new TH2D("POTScaledNCCOHRecokine_pio_flagPlot2D",LabelXAxiskine_pio_flag2D,NBinskine_pio_flag,
			  min_kine_pio_flag,max_kine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* NonNCCOHRecokine_pio_flagPlot = new TH1D("NonNCCOHRecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* QERecokine_pio_flagPlot = new TH1D("QERecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* MECRecokine_pio_flagPlot = new TH1D("MECRecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* RESRecokine_pio_flagPlot = new TH1D("RESRecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* DISRecokine_pio_flagPlot = new TH1D("DISRecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);
		TH1D* COHRecokine_pio_flagPlot = new TH1D("COHRecokine_pio_flagPlot",LabelXAxiskine_pio_flag,NBinskine_pio_flag,min_kine_pio_flag,max_kine_pio_flag);	

		//----------------------------------------//

		// numu_score

		TH1D* Reconumu_scorePlot = new TH1D("Reconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* NCCOHReconumu_scorePlot = new TH1D("NCCOHReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);	
		TH1D* NCCOHTruenumu_scorePlot = new TH1D("NCCOHTruenumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH2D* NCCOHReconumu_scorePlot2D = new TH2D("NCCOHReconumu_scorePlot2D",LabelXAxisnumu_score2D,NBinsnumu_score,
			  min_numu_score,max_numu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH2D* POTScaledNCCOHReconumu_scorePlot2D = new TH2D("POTScaledNCCOHReconumu_scorePlot2D",LabelXAxisnumu_score2D,NBinsnumu_score,
			  min_numu_score,max_numu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* NonNCCOHReconumu_scorePlot = new TH1D("NonNCCOHReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* QEReconumu_scorePlot = new TH1D("QEReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* MECReconumu_scorePlot = new TH1D("MECReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* RESReconumu_scorePlot = new TH1D("RESReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* DISReconumu_scorePlot = new TH1D("DISReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);
		TH1D* COHReconumu_scorePlot = new TH1D("COHReconumu_scorePlot",LabelXAxisnumu_score,NBinsnumu_score,min_numu_score,max_numu_score);	

		//----------------------------------------//

		// nc_pio_score

		TH1D* Reconc_pio_scorePlot = new TH1D("Reconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* NCCOHReconc_pio_scorePlot = new TH1D("NCCOHReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);	
		TH1D* NCCOHTruenc_pio_scorePlot = new TH1D("NCCOHTruenc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH2D* NCCOHReconc_pio_scorePlot2D = new TH2D("NCCOHReconc_pio_scorePlot2D",LabelXAxisnc_pio_score2D,NBinsnc_pio_score,
			  min_nc_pio_score,max_nc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH2D* POTScaledNCCOHReconc_pio_scorePlot2D = new TH2D("POTScaledNCCOHReconc_pio_scorePlot2D",LabelXAxisnc_pio_score2D,NBinsnc_pio_score,
			  min_nc_pio_score,max_nc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* NonNCCOHReconc_pio_scorePlot = new TH1D("NonNCCOHReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* QEReconc_pio_scorePlot = new TH1D("QEReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* MECReconc_pio_scorePlot = new TH1D("MECReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* RESReconc_pio_scorePlot = new TH1D("RESReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* DISReconc_pio_scorePlot = new TH1D("DISReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);
		TH1D* COHReconc_pio_scorePlot = new TH1D("COHReconc_pio_scorePlot",LabelXAxisnc_pio_score,NBinsnc_pio_score,min_nc_pio_score,max_nc_pio_score);	

		//----------------------------------------//

		// single_photon_numu_score

		TH1D* Recosingle_photon_numu_scorePlot = new TH1D("Recosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* NCCOHRecosingle_photon_numu_scorePlot = new TH1D("NCCOHRecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);	
		TH1D* NCCOHTruesingle_photon_numu_scorePlot = new TH1D("NCCOHTruesingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH2D* NCCOHRecosingle_photon_numu_scorePlot2D = new TH2D("NCCOHRecosingle_photon_numu_scorePlot2D",LabelXAxissingle_photon_numu_score2D,NBinssingle_photon_numu_score,
			  min_single_photon_numu_score,max_single_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH2D* POTScaledNCCOHRecosingle_photon_numu_scorePlot2D = new TH2D("POTScaledNCCOHRecosingle_photon_numu_scorePlot2D",LabelXAxissingle_photon_numu_score2D,NBinssingle_photon_numu_score,
			  min_single_photon_numu_score,max_single_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* NonNCCOHRecosingle_photon_numu_scorePlot = new TH1D("NonNCCOHRecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* QERecosingle_photon_numu_scorePlot = new TH1D("QERecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* MECRecosingle_photon_numu_scorePlot = new TH1D("MECRecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* RESRecosingle_photon_numu_scorePlot = new TH1D("RESRecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* DISRecosingle_photon_numu_scorePlot = new TH1D("DISRecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		TH1D* COHRecosingle_photon_numu_scorePlot = new TH1D("COHRecosingle_photon_numu_scorePlot",LabelXAxissingle_photon_numu_score,NBinssingle_photon_numu_score,min_single_photon_numu_score,max_single_photon_numu_score);
		
		//----------------------------------------//

		// single_photon_other_score

		TH1D* Recosingle_photon_other_scorePlot = new TH1D("Recosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* NCCOHRecosingle_photon_other_scorePlot = new TH1D("NCCOHRecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);	
		TH1D* NCCOHTruesingle_photon_other_scorePlot = new TH1D("NCCOHTruesingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH2D* NCCOHRecosingle_photon_other_scorePlot2D = new TH2D("NCCOHRecosingle_photon_other_scorePlot2D",LabelXAxissingle_photon_other_score2D,NBinssingle_photon_other_score,
			  min_single_photon_other_score,max_single_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH2D* POTScaledNCCOHRecosingle_photon_other_scorePlot2D = new TH2D("POTScaledNCCOHRecosingle_photon_other_scorePlot2D",LabelXAxissingle_photon_other_score2D,NBinssingle_photon_other_score,
			  min_single_photon_other_score,max_single_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* NonNCCOHRecosingle_photon_other_scorePlot = new TH1D("NonNCCOHRecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* QERecosingle_photon_other_scorePlot = new TH1D("QERecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* MECRecosingle_photon_other_scorePlot = new TH1D("MECRecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* RESRecosingle_photon_other_scorePlot = new TH1D("RESRecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* DISRecosingle_photon_other_scorePlot = new TH1D("DISRecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);
		TH1D* COHRecosingle_photon_other_scorePlot = new TH1D("COHRecosingle_photon_other_scorePlot",LabelXAxissingle_photon_other_score,NBinssingle_photon_other_score,min_single_photon_other_score,max_single_photon_other_score);

		//----------------------------------------//

		// single_photon_ncpi0_score

		TH1D* Recosingle_photon_ncpi0_scorePlot = new TH1D("Recosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* NCCOHRecosingle_photon_ncpi0_scorePlot = new TH1D("NCCOHRecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);	
		TH1D* NCCOHTruesingle_photon_ncpi0_scorePlot = new TH1D("NCCOHTruesingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH2D* NCCOHRecosingle_photon_ncpi0_scorePlot2D = new TH2D("NCCOHRecosingle_photon_ncpi0_scorePlot2D",LabelXAxissingle_photon_ncpi0_score2D,NBinssingle_photon_ncpi0_score,
			  min_single_photon_ncpi0_score,max_single_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH2D* POTScaledNCCOHRecosingle_photon_ncpi0_scorePlot2D = new TH2D("POTScaledNCCOHRecosingle_photon_ncpi0_scorePlot2D",LabelXAxissingle_photon_ncpi0_score2D,NBinssingle_photon_ncpi0_score,
			  min_single_photon_ncpi0_score,max_single_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* NonNCCOHRecosingle_photon_ncpi0_scorePlot = new TH1D("NonNCCOHRecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* QERecosingle_photon_ncpi0_scorePlot = new TH1D("QERecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* MECRecosingle_photon_ncpi0_scorePlot = new TH1D("MECRecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* RESRecosingle_photon_ncpi0_scorePlot = new TH1D("RESRecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* DISRecosingle_photon_ncpi0_scorePlot = new TH1D("DISRecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);
		TH1D* COHRecosingle_photon_ncpi0_scorePlot = new TH1D("COHRecosingle_photon_ncpi0_scorePlot",LabelXAxissingle_photon_ncpi0_score,NBinssingle_photon_ncpi0_score,min_single_photon_ncpi0_score,max_single_photon_ncpi0_score);

		//----------------------------------------//

		// single_photon_nue_score

		TH1D* Recosingle_photon_nue_scorePlot = new TH1D("Recosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* NCCOHRecosingle_photon_nue_scorePlot = new TH1D("NCCOHRecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);	
		TH1D* NCCOHTruesingle_photon_nue_scorePlot = new TH1D("NCCOHTruesingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH2D* NCCOHRecosingle_photon_nue_scorePlot2D = new TH2D("NCCOHRecosingle_photon_nue_scorePlot2D",LabelXAxissingle_photon_nue_score2D,NBinssingle_photon_nue_score,
			  min_single_photon_nue_score,max_single_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH2D* POTScaledNCCOHRecosingle_photon_nue_scorePlot2D = new TH2D("POTScaledNCCOHRecosingle_photon_nue_scorePlot2D",LabelXAxissingle_photon_nue_score2D,NBinssingle_photon_nue_score,
			  min_single_photon_nue_score,max_single_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* NonNCCOHRecosingle_photon_nue_scorePlot = new TH1D("NonNCCOHRecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* QERecosingle_photon_nue_scorePlot = new TH1D("QERecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* MECRecosingle_photon_nue_scorePlot = new TH1D("MECRecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* RESRecosingle_photon_nue_scorePlot = new TH1D("RESRecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* DISRecosingle_photon_nue_scorePlot = new TH1D("DISRecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);
		TH1D* COHRecosingle_photon_nue_scorePlot = new TH1D("COHRecosingle_photon_nue_scorePlot",LabelXAxissingle_photon_nue_score,NBinssingle_photon_nue_score,min_single_photon_nue_score,max_single_photon_nue_score);

		//----------------------------------------//

		// Blips

		TH1D* RecoBlip_xPlot = new TH1D("RecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* NCCOHRecoBlip_xPlot = new TH1D("NCCOHRecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);	
		TH1D* NCCOHTrueBlip_xPlot = new TH1D("NCCOHTrueBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH2D* NCCOHRecoBlip_xPlot2D = new TH2D("NCCOHRecoBlip_xPlot2D",LabelXAxisBlip_x2D,NBinsBlip_x,
			  borderx,FVx-borderx,NBinsBlip_x,borderx,FVx-borderx);
		TH2D* POTScaledNCCOHRecoBlip_xPlot2D = new TH2D("POTScaledNCCOHRecoBlip_xPlot2D",LabelXAxisBlip_x2D,NBinsBlip_x,
			  borderx,FVx-borderx,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* NonNCCOHRecoBlip_xPlot = new TH1D("NonNCCOHRecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* QERecoBlip_xPlot = new TH1D("QERecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* MECRecoBlip_xPlot = new TH1D("MECRecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* RESRecoBlip_xPlot = new TH1D("RESRecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* DISRecoBlip_xPlot = new TH1D("DISRecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);
		TH1D* COHRecoBlip_xPlot = new TH1D("COHRecoBlip_xPlot",LabelXAxisBlip_x,NBinsBlip_x,borderx,FVx-borderx);	

		TH1D* RecoBlip_yPlot = new TH1D("RecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* NCCOHRecoBlip_yPlot = new TH1D("NCCOHRecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);	
		TH1D* NCCOHTrueBlip_yPlot = new TH1D("NCCOHTrueBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH2D* NCCOHRecoBlip_yPlot2D = new TH2D("NCCOHRecoBlip_yPlot2D",LabelXAxisBlip_y2D,NBinsBlip_y,
			  -FVy/2. + bordery,FVy/2. - bordery,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH2D* POTScaledNCCOHRecoBlip_yPlot2D = new TH2D("POTScaledNCCOHRecoBlip_yPlot2D",LabelXAxisBlip_y2D,NBinsBlip_y,
			  -FVy/2. + bordery,FVy/2. - bordery,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* NonNCCOHRecoBlip_yPlot = new TH1D("NonNCCOHRecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* QERecoBlip_yPlot = new TH1D("QERecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* MECRecoBlip_yPlot = new TH1D("MECRecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* RESRecoBlip_yPlot = new TH1D("RESRecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* DISRecoBlip_yPlot = new TH1D("DISRecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		TH1D* COHRecoBlip_yPlot = new TH1D("COHRecoBlip_yPlot",LabelXAxisBlip_y,NBinsBlip_y,-FVy/2. + bordery,FVy/2. - bordery);
		
		TH1D* RecoBlip_zPlot = new TH1D("RecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* NCCOHRecoBlip_zPlot = new TH1D("NCCOHRecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);	
		TH1D* NCCOHTrueBlip_zPlot = new TH1D("NCCOHTrueBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH2D* NCCOHRecoBlip_zPlot2D = new TH2D("NCCOHRecoBlip_zPlot2D",LabelXAxisBlip_z2D,NBinsBlip_z,
			  borderz,FVz-borderz,NBinsBlip_z,borderz,FVz-borderz);
		TH2D* POTScaledNCCOHRecoBlip_zPlot2D = new TH2D("POTScaledNCCOHRecoBlip_zPlot2D",LabelXAxisBlip_z2D,NBinsBlip_z,
			  borderz,FVz-borderz,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* NonNCCOHRecoBlip_zPlot = new TH1D("NonNCCOHRecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* QERecoBlip_zPlot = new TH1D("QERecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* MECRecoBlip_zPlot = new TH1D("MECRecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* RESRecoBlip_zPlot = new TH1D("RESRecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* DISRecoBlip_zPlot = new TH1D("DISRecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		TH1D* COHRecoBlip_zPlot = new TH1D("COHRecoBlip_zPlot",LabelXAxisBlip_z,NBinsBlip_z,borderz,FVz-borderz);
		
		TH1D* ReconBlips_savedPlot = new TH1D("ReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* NCCOHReconBlips_savedPlot = new TH1D("NCCOHReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);	
		TH1D* NCCOHTruenBlips_savedPlot = new TH1D("NCCOHTruenBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH2D* NCCOHReconBlips_savedPlot2D = new TH2D("NCCOHReconBlips_savedPlot2D",LabelXAxisnBlips_saved2D,NBinsnBlips_saved,
			  nBlips_saved_min,nBlips_saved_max,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH2D* POTScaledNCCOHReconBlips_savedPlot2D = new TH2D("POTScaledNCCOHReconBlips_savedPlot2D",LabelXAxisnBlips_saved2D,NBinsnBlips_saved,
			  nBlips_saved_min,nBlips_saved_max,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* NonNCCOHReconBlips_savedPlot = new TH1D("NonNCCOHReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* QEReconBlips_savedPlot = new TH1D("QEReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* MECReconBlips_savedPlot = new TH1D("MECReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* RESReconBlips_savedPlot = new TH1D("RESReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* DISReconBlips_savedPlot = new TH1D("DISReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);
		TH1D* COHReconBlips_savedPlot = new TH1D("COHReconBlips_savedPlot",LabelXAxisnBlips_saved,NBinsnBlips_saved,nBlips_saved_min,nBlips_saved_max);		

		TH1D* ReconBlips_radiusPlot = new TH1D("ReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* NCCOHReconBlips_radiusPlot = new TH1D("NCCOHReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);	
		TH1D* NCCOHTruenBlips_radiusPlot = new TH1D("NCCOHTruenBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH2D* NCCOHReconBlips_radiusPlot2D = new TH2D("NCCOHReconBlips_radiusPlot2D",LabelXAxisnBlips_radius2D,NBinsnBlips_radius,
			  nBlips_radius_min,nBlips_radius_max,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH2D* POTScaledNCCOHReconBlips_radiusPlot2D = new TH2D("POTScaledNCCOHReconBlips_radiusPlot2D",LabelXAxisnBlips_radius2D,NBinsnBlips_radius,
			  nBlips_radius_min,nBlips_radius_max,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* NonNCCOHReconBlips_radiusPlot = new TH1D("NonNCCOHReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* QEReconBlips_radiusPlot = new TH1D("QEReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* MECReconBlips_radiusPlot = new TH1D("MECReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* RESReconBlips_radiusPlot = new TH1D("RESReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* DISReconBlips_radiusPlot = new TH1D("DISReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);
		TH1D* COHReconBlips_radiusPlot = new TH1D("COHReconBlips_radiusPlot",LabelXAxisnBlips_radius,NBinsnBlips_radius,nBlips_radius_min,nBlips_radius_max);		
		
		TH1D* RecoBlip_energyPlot = new TH1D("RecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* NCCOHRecoBlip_energyPlot = new TH1D("NCCOHRecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);	
		TH1D* NCCOHTrueBlip_energyPlot = new TH1D("NCCOHTrueBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH2D* NCCOHRecoBlip_energyPlot2D = new TH2D("NCCOHRecoBlip_energyPlot2D",LabelXAxisBlip_energy2D,NBinsBlip_energy,
			  Blip_energy_min,Blip_energy_max,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH2D* POTScaledNCCOHRecoBlip_energyPlot2D = new TH2D("POTScaledNCCOHRecoBlip_energyPlot2D",LabelXAxisBlip_energy2D,NBinsBlip_energy,
			  Blip_energy_min,Blip_energy_max,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* NonNCCOHRecoBlip_energyPlot = new TH1D("NonNCCOHRecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* QERecoBlip_energyPlot = new TH1D("QERecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* MECRecoBlip_energyPlot = new TH1D("MECRecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* RESRecoBlip_energyPlot = new TH1D("RESRecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* DISRecoBlip_energyPlot = new TH1D("DISRecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);
		TH1D* COHRecoBlip_energyPlot = new TH1D("COHRecoBlip_energyPlot",LabelXAxisBlip_energy,NBinsBlip_energy,Blip_energy_min,Blip_energy_max);		

		TH1D* RecoBlip_proxtrkdistPlot = new TH1D("RecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* NCCOHRecoBlip_proxtrkdistPlot = new TH1D("NCCOHRecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);	
		TH1D* NCCOHTrueBlip_proxtrkdistPlot = new TH1D("NCCOHTrueBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH2D* NCCOHRecoBlip_proxtrkdistPlot2D = new TH2D("NCCOHRecoBlip_proxtrkdistPlot2D",LabelXAxisBlip_proxtrkdist2D,NBinsBlip_proxtrkdist,
			  Blip_proxtrkdist_min,Blip_proxtrkdist_max,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH2D* POTScaledNCCOHRecoBlip_proxtrkdistPlot2D = new TH2D("POTScaledNCCOHRecoBlip_proxtrkdistPlot2D",LabelXAxisBlip_proxtrkdist2D,NBinsBlip_proxtrkdist,
			  Blip_proxtrkdist_min,Blip_proxtrkdist_max,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* NonNCCOHRecoBlip_proxtrkdistPlot = new TH1D("NonNCCOHRecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* QERecoBlip_proxtrkdistPlot = new TH1D("QERecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* MECRecoBlip_proxtrkdistPlot = new TH1D("MECRecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* RESRecoBlip_proxtrkdistPlot = new TH1D("RESRecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* DISRecoBlip_proxtrkdistPlot = new TH1D("DISRecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);
		TH1D* COHRecoBlip_proxtrkdistPlot = new TH1D("COHRecoBlip_proxtrkdistPlot",LabelXAxisBlip_proxtrkdist,NBinsBlip_proxtrkdist,Blip_proxtrkdist_min,Blip_proxtrkdist_max);	

		TH1D* Recoblip_vrtPlot = new TH1D("Recoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* NCCOHRecoblip_vrtPlot = new TH1D("NCCOHRecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);	
		TH1D* NCCOHTrueblip_vrtPlot = new TH1D("NCCOHTrueblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH2D* NCCOHRecoblip_vrtPlot2D = new TH2D("NCCOHRecoblip_vrtPlot2D",LabelXAxisblip_vrt2D,NBinsblip_vrt,
			  blip_vrt_min,blip_vrt_max,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH2D* POTScaledNCCOHRecoblip_vrtPlot2D = new TH2D("POTScaledNCCOHRecoblip_vrtPlot2D",LabelXAxisblip_vrt2D,NBinsblip_vrt,
			  blip_vrt_min,blip_vrt_max,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* NonNCCOHRecoblip_vrtPlot = new TH1D("NonNCCOHRecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* QERecoblip_vrtPlot = new TH1D("QERecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* MECRecoblip_vrtPlot = new TH1D("MECRecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* RESRecoblip_vrtPlot = new TH1D("RESRecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* DISRecoblip_vrtPlot = new TH1D("DISRecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		TH1D* COHRecoblip_vrtPlot = new TH1D("COHRecoblip_vrtPlot",LabelXAxisblip_vrt,NBinsblip_vrt,blip_vrt_min,blip_vrt_max);
		
		TH1D* Recoblip_cos_alphapi0Plot = new TH1D("Recoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* NCCOHRecoblip_cos_alphapi0Plot = new TH1D("NCCOHRecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);	
		TH1D* NCCOHTrueblip_cos_alphapi0Plot = new TH1D("NCCOHTrueblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH2D* NCCOHRecoblip_cos_alphapi0Plot2D = new TH2D("NCCOHRecoblip_cos_alphapi0Plot2D",LabelXAxisblip_cos_alphapi02D,NBinsblip_cos_alphapi0,
			  blip_cos_alphapi0_min,blip_cos_alphapi0_max,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH2D* POTScaledNCCOHRecoblip_cos_alphapi0Plot2D = new TH2D("POTScaledNCCOHRecoblip_cos_alphapi0Plot2D",LabelXAxisblip_cos_alphapi02D,NBinsblip_cos_alphapi0,
			  blip_cos_alphapi0_min,blip_cos_alphapi0_max,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* NonNCCOHRecoblip_cos_alphapi0Plot = new TH1D("NonNCCOHRecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* QERecoblip_cos_alphapi0Plot = new TH1D("QERecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* MECRecoblip_cos_alphapi0Plot = new TH1D("MECRecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* RESRecoblip_cos_alphapi0Plot = new TH1D("RESRecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* DISRecoblip_cos_alphapi0Plot = new TH1D("DISRecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);
		TH1D* COHRecoblip_cos_alphapi0Plot = new TH1D("COHRecoblip_cos_alphapi0Plot",LabelXAxisblip_cos_alphapi0,NBinsblip_cos_alphapi0,blip_cos_alphapi0_min,blip_cos_alphapi0_max);		

		TH1D* Recoblip_cos_alphag1Plot = new TH1D("Recoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* NCCOHRecoblip_cos_alphag1Plot = new TH1D("NCCOHRecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);	
		TH1D* NCCOHTrueblip_cos_alphag1Plot = new TH1D("NCCOHTrueblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH2D* NCCOHRecoblip_cos_alphag1Plot2D = new TH2D("NCCOHRecoblip_cos_alphag1Plot2D",LabelXAxisblip_cos_alphag12D,NBinsblip_cos_alphag1,
			  blip_cos_alphag1_min,blip_cos_alphag1_max,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH2D* POTScaledNCCOHRecoblip_cos_alphag1Plot2D = new TH2D("POTScaledNCCOHRecoblip_cos_alphag1Plot2D",LabelXAxisblip_cos_alphag12D,NBinsblip_cos_alphag1,
			  blip_cos_alphag1_min,blip_cos_alphag1_max,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* NonNCCOHRecoblip_cos_alphag1Plot = new TH1D("NonNCCOHRecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* QERecoblip_cos_alphag1Plot = new TH1D("QERecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* MECRecoblip_cos_alphag1Plot = new TH1D("MECRecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* RESRecoblip_cos_alphag1Plot = new TH1D("RESRecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* DISRecoblip_cos_alphag1Plot = new TH1D("DISRecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);
		TH1D* COHRecoblip_cos_alphag1Plot = new TH1D("COHRecoblip_cos_alphag1Plot",LabelXAxisblip_cos_alphag1,NBinsblip_cos_alphag1,blip_cos_alphag1_min,blip_cos_alphag1_max);		

		TH1D* Recoblip_cos_alphag2Plot = new TH1D("Recoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* NCCOHRecoblip_cos_alphag2Plot = new TH1D("NCCOHRecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);	
		TH1D* NCCOHTrueblip_cos_alphag2Plot = new TH1D("NCCOHTrueblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH2D* NCCOHRecoblip_cos_alphag2Plot2D = new TH2D("NCCOHRecoblip_cos_alphag2Plot2D",LabelXAxisblip_cos_alphag22D,NBinsblip_cos_alphag2,
			  blip_cos_alphag2_min,blip_cos_alphag2_max,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH2D* POTScaledNCCOHRecoblip_cos_alphag2Plot2D = new TH2D("POTScaledNCCOHRecoblip_cos_alphag2Plot2D",LabelXAxisblip_cos_alphag22D,NBinsblip_cos_alphag2,
			  blip_cos_alphag2_min,blip_cos_alphag2_max,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* NonNCCOHRecoblip_cos_alphag2Plot = new TH1D("NonNCCOHRecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* QERecoblip_cos_alphag2Plot = new TH1D("QERecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* MECRecoblip_cos_alphag2Plot = new TH1D("MECRecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* RESRecoblip_cos_alphag2Plot = new TH1D("RESRecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* DISRecoblip_cos_alphag2Plot = new TH1D("DISRecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);
		TH1D* COHRecoblip_cos_alphag2Plot = new TH1D("COHRecoblip_cos_alphag2Plot",LabelXAxisblip_cos_alphag2,NBinsblip_cos_alphag2,blip_cos_alphag2_min,blip_cos_alphag2_max);		


		//----------------------------------------//		

		// g1 Momentum

		TH1D* Recog1MomentumPlot = new TH1D("Recog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* NCCOHRecog1MomentumPlot = new TH1D("NCCOHRecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);	
		TH1D* NCCOHTrueg1MomentumPlot = new TH1D("NCCOHTrueg1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH2D* NCCOHRecog1MomentumPlot2D = new TH2D("NCCOHRecog1MomentumPlot2D",LabelXAxisg1Momentum2D,NBinsg1Momentum,
			  ArrayNBinsg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH2D* POTScaledNCCOHRecog1MomentumPlot2D = new TH2D("POTScaledNCCOHRecog1MomentumPlot2D",LabelXAxisg1Momentum2D,NBinsg1Momentum,
			  ArrayNBinsg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* NonNCCOHRecog1MomentumPlot = new TH1D("NonNCCOHRecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* QERecog1MomentumPlot = new TH1D("QERecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* MECRecog1MomentumPlot = new TH1D("MECRecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* RESRecog1MomentumPlot = new TH1D("RESRecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* DISRecog1MomentumPlot = new TH1D("DISRecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);
		TH1D* COHRecog1MomentumPlot = new TH1D("COHRecog1MomentumPlot",LabelXAxisg1Momentum,NBinsg1Momentum,ArrayNBinsg1Momentum);

		//----------------------------------------//		

		// g2 Momentum

		TH1D* Recog2MomentumPlot = new TH1D("Recog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* NCCOHRecog2MomentumPlot = new TH1D("NCCOHRecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);	
		TH1D* NCCOHTrueg2MomentumPlot = new TH1D("NCCOHTrueg2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH2D* NCCOHRecog2MomentumPlot2D = new TH2D("NCCOHRecog2MomentumPlot2D",LabelXAxisg2Momentum2D,NBinsg2Momentum,
			  ArrayNBinsg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH2D* POTScaledNCCOHRecog2MomentumPlot2D = new TH2D("POTScaledNCCOHRecog2MomentumPlot2D",LabelXAxisg2Momentum2D,NBinsg2Momentum,
			  ArrayNBinsg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* NonNCCOHRecog2MomentumPlot = new TH1D("NonNCCOHRecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* QERecog2MomentumPlot = new TH1D("QERecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* MECRecog2MomentumPlot = new TH1D("MECRecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* RESRecog2MomentumPlot = new TH1D("RESRecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* DISRecog2MomentumPlot = new TH1D("DISRecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);
		TH1D* COHRecog2MomentumPlot = new TH1D("COHRecog2MomentumPlot",LabelXAxisg2Momentum,NBinsg2Momentum,ArrayNBinsg2Momentum);

		//----------------------------------------//		

		// Pi0 Momentum

		TH1D* RecoPi0MomentumPlot = new TH1D("RecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* NCCOHRecoPi0MomentumPlot = new TH1D("NCCOHRecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);	
		TH1D* NCCOHTruePi0MomentumPlot = new TH1D("NCCOHTruePi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH2D* NCCOHRecoPi0MomentumPlot2D = new TH2D("NCCOHRecoPi0MomentumPlot2D",LabelXAxisPi0Momentum2D,NBinsPi0Momentum,
			  ArrayNBinsPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH2D* POTScaledNCCOHRecoPi0MomentumPlot2D = new TH2D("POTScaledNCCOHRecoPi0MomentumPlot2D",LabelXAxisPi0Momentum2D,NBinsPi0Momentum,
			  ArrayNBinsPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* NonNCCOHRecoPi0MomentumPlot = new TH1D("NonNCCOHRecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* QERecoPi0MomentumPlot = new TH1D("QERecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* MECRecoPi0MomentumPlot = new TH1D("MECRecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* RESRecoPi0MomentumPlot = new TH1D("RESRecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* DISRecoPi0MomentumPlot = new TH1D("DISRecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);
		TH1D* COHRecoPi0MomentumPlot = new TH1D("COHRecoPi0MomentumPlot",LabelXAxisPi0Momentum,NBinsPi0Momentum,ArrayNBinsPi0Momentum);

		//----------------------------------------//		

		// two-shower angle

		TH1D* Recotwo_shower_anglePlot = new TH1D("Recotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* NCCOHRecotwo_shower_anglePlot = new TH1D("NCCOHRecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);	
		TH1D* NCCOHTruetwo_shower_anglePlot = new TH1D("NCCOHTruetwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH2D* NCCOHRecotwo_shower_anglePlot2D = new TH2D("NCCOHRecotwo_shower_anglePlot2D",LabelXAxistwo_shower_angle2D,NBinstwo_shower_angle,
			  ArrayNBinstwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH2D* POTScaledNCCOHRecotwo_shower_anglePlot2D = new TH2D("POTScaledNCCOHRecotwo_shower_anglePlot2D",LabelXAxistwo_shower_angle2D,NBinstwo_shower_angle,
			  ArrayNBinstwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* NonNCCOHRecotwo_shower_anglePlot = new TH1D("NonNCCOHRecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* QERecotwo_shower_anglePlot = new TH1D("QERecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* MECRecotwo_shower_anglePlot = new TH1D("MECRecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* RESRecotwo_shower_anglePlot = new TH1D("RESRecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* DISRecotwo_shower_anglePlot = new TH1D("DISRecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);
		TH1D* COHRecotwo_shower_anglePlot = new TH1D("COHRecotwo_shower_anglePlot",LabelXAxistwo_shower_angle,NBinstwo_shower_angle,ArrayNBinstwo_shower_angle);

		//----------------------------------------//		

		// g1 CosTheta

		TH1D* Recog1CosThetaPlot = new TH1D("Recog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* NCCOHRecog1CosThetaPlot = new TH1D("NCCOHRecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);	
		TH1D* NCCOHTrueg1CosThetaPlot = new TH1D("NCCOHTrueg1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH2D* NCCOHRecog1CosThetaPlot2D = new TH2D("NCCOHRecog1CosThetaPlot2D",LabelXAxisg1CosTheta2D,NBinsg1CosTheta,
			  ArrayNBinsg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH2D* POTScaledNCCOHRecog1CosThetaPlot2D = new TH2D("POTScaledNCCOHRecog1CosThetaPlot2D",LabelXAxisg1CosTheta2D,NBinsg1CosTheta,
			  ArrayNBinsg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* NonNCCOHRecog1CosThetaPlot = new TH1D("NonNCCOHRecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* QERecog1CosThetaPlot = new TH1D("QERecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* MECRecog1CosThetaPlot = new TH1D("MECRecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* RESRecog1CosThetaPlot = new TH1D("RESRecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* DISRecog1CosThetaPlot = new TH1D("DISRecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);
		TH1D* COHRecog1CosThetaPlot = new TH1D("COHRecog1CosThetaPlot",LabelXAxisg1CosTheta,NBinsg1CosTheta,ArrayNBinsg1CosTheta);

		//----------------------------------------//		

		// g2 CosTheta

		TH1D* Recog2CosThetaPlot = new TH1D("Recog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* NCCOHRecog2CosThetaPlot = new TH1D("NCCOHRecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);	
		TH1D* NCCOHTrueg2CosThetaPlot = new TH1D("NCCOHTrueg2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH2D* NCCOHRecog2CosThetaPlot2D = new TH2D("NCCOHRecog2CosThetaPlot2D",LabelXAxisg2CosTheta2D,NBinsg2CosTheta,
			  ArrayNBinsg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH2D* POTScaledNCCOHRecog2CosThetaPlot2D = new TH2D("POTScaledNCCOHRecog2CosThetaPlot2D",LabelXAxisg2CosTheta2D,NBinsg2CosTheta,
			  ArrayNBinsg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* NonNCCOHRecog2CosThetaPlot = new TH1D("NonNCCOHRecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* QERecog2CosThetaPlot = new TH1D("QERecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* MECRecog2CosThetaPlot = new TH1D("MECRecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* RESRecog2CosThetaPlot = new TH1D("RESRecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* DISRecog2CosThetaPlot = new TH1D("DISRecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);
		TH1D* COHRecog2CosThetaPlot = new TH1D("COHRecog2CosThetaPlot",LabelXAxisg2CosTheta,NBinsg2CosTheta,ArrayNBinsg2CosTheta);

		//----------------------------------------//		

		// Pi0 CosTheta

		TH1D* RecoPi0CosThetaPlot = new TH1D("RecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* NCCOHRecoPi0CosThetaPlot = new TH1D("NCCOHRecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);	
		TH1D* NCCOHTruePi0CosThetaPlot = new TH1D("NCCOHTruePi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH2D* NCCOHRecoPi0CosThetaPlot2D = new TH2D("NCCOHRecoPi0CosThetaPlot2D",LabelXAxisPi0CosTheta2D,NBinsPi0CosTheta,
			  ArrayNBinsPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH2D* POTScaledNCCOHRecoPi0CosThetaPlot2D = new TH2D("POTScaledNCCOHRecoPi0CosThetaPlot2D",LabelXAxisPi0CosTheta2D,NBinsPi0CosTheta,
			  ArrayNBinsPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* NonNCCOHRecoPi0CosThetaPlot = new TH1D("NonNCCOHRecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* QERecoPi0CosThetaPlot = new TH1D("QERecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* MECRecoPi0CosThetaPlot = new TH1D("MECRecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* RESRecoPi0CosThetaPlot = new TH1D("RESRecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* DISRecoPi0CosThetaPlot = new TH1D("DISRecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		TH1D* COHRecoPi0CosThetaPlot = new TH1D("COHRecoPi0CosThetaPlot",LabelXAxisPi0CosTheta,NBinsPi0CosTheta,ArrayNBinsPi0CosTheta);
		
		//----------------------------------------//

		// Single Bin

		TH1D* RecoSingleBinPlot = new TH1D("RecoSingleBinPlot","",1,0.,1.);
		TH1D* NCCOHRecoSingleBinPlot = new TH1D("NCCOHRecoSingleBinPlot","",1,0.,1.);
		TH1D* NCCOHTrueSingleBinPlot = new TH1D("NCCOHTrueSingleBinPlot","",1,0.,1.);
		TH2D* NCCOHRecoSingleBinPlot2D = new TH2D("NCCOHRecoSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH2D* POTScaledNCCOHRecoSingleBinPlot2D = new TH2D("POTScaledNCCOHRecoSingleBinPlot2D","; ;",1,0.,1.,1,0.,1.);
		TH1D* NonNCCOHRecoSingleBinPlot = new TH1D("NonNCCOHRecoSingleBinPlot","",1,0.,1.);
		TH1D* QERecoSingleBinPlot = new TH1D("QERecoSingleBinPlot","",1,0.,1.);
		TH1D* MECRecoSingleBinPlot = new TH1D("MECRecoSingleBinPlot","",1,0.,1.);
		TH1D* RESRecoSingleBinPlot = new TH1D("RESRecoSingleBinPlot","",1,0.,1.);
		TH1D* DISRecoSingleBinPlot = new TH1D("DISRecoSingleBinPlot","",1,0.,1.);
		TH1D* COHRecoSingleBinPlot = new TH1D("COHRecoSingleBinPlot","",1,0.,1.);

		//----------------------------------------//

		// counters

		int signal_counter = 0;
		int coh_counter = 0;
		int pass_selection_counter = 0;

		int bkg_0pi0_X_counter = 0;
		int bkg_Mpi0_X_counter = 0; 
		int bkg_bwds_1pi0_X_counter = 0;
		int bkg_1n_0p_1pi0_X_counter = 0;
		int bkg_Nn_0p_1pi0_X_counter = 0;
		int bkg_1p_0n_1pi0_X_counter = 0;
		int bkg_Np_0n_1pi0_X_counter = 0;
		int bkg_1pi0_Npipm_X_counter = 0;
		int bkg_1pi0_Np_Nn_0pipm_X_counter = 0;
		int bkg_1pi0_Np_Nn_Npipm_X_counter = 0;
		int bkg_1pi0_Nmh_X_counter = 0;
		int bkg_1pi0_Nl_X_counter = 0;								
		int bkg_other_counter = 0;		

		//----------------------------------------//

		// Loop over the events

		cout << nentries << " events included in the file" << endl;

		for (Long64_t jentry=0; jentry<nentries;jentry++) {

			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);	nbytes += nb;

			if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;

			//--------------------//

			weight = POTWeight;
			if (jentry == nentries -1) { cout << "pot scale = " << POTWeight << endl; }

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
			
			//--------------------//

			// Genie, flux & reinteraction weights for systematics

			if ( 
				   fUniverseIndex != -1 && (
				   fWhichSample == "mcc9_10_Overlay9_Run1" 
				|| fWhichSample == "mcc9_10_Overlay9_Run2" 
				|| fWhichSample == "mcc9_10_Overlay9_Run3" 
				|| fWhichSample == "mcc9_10_Overlay9_Run4a" 
				|| fWhichSample == "mcc9_10_Overlay9_Run4b"
				|| fWhichSample == "mcc9_10_Overlay9_Run4c" 
				|| fWhichSample == "mcc9_10_Overlay9_Run4d" 
				|| fWhichSample == "mcc9_10_Overlay9_Run5" 
				|| fWhichSample == "mcc9_10_Overlay9_Combined" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run1" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run2" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run3" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run4a" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run4b"
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run4b_standalone"
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run4b_unified"
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run4c" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run4d" 
				|| fWhichSample == "mcc9_10_OverlayDirt9_Run5"
				|| fWhichSample == "mcc9_10_OverlayDirt9_Combined"
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

				// Flux weights
				if (fEventWeightLabel == "fluxes") { 

					if ( fUniverseIndex < (int)(fluxes->size()) ) {

						weight = weight*fluxes->at(fUniverseIndex); 

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

			//--------------------//

			if ( fabs(weight) != weight) { continue; } // Securing against infinities

			//--------------------//

			// Contained Reconstructed Vertex

			TVector3 RecoVertex(Vertex_X->at(0),Vertex_Y->at(0),Vertex_Z->at(0));

			if ( !tools.inFVVector(RecoVertex) ) { continue; }
			if ( Vertex_Z->at(0) < 250 ) { continue; }
			if ( Vertex_Z->at(0) > 660 && Vertex_Z->at(0) < 760 ) { continue; }			

			int cut_reco_blip_counter_radius = 0;

			for (int iblip = 0; iblip < nBlips_saved; iblip++) {

				TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
				double blip_vrt_dist = (RecoVertex - blip_3d).Mag();

				Recoblip_vrtPlot->Fill(blip_vrt_dist,weight);	
				
				if (blip_vrt_dist < radius) { cut_reco_blip_counter_radius++; }

			}

			//if (cut_reco_blip_counter_radius != 0) { continue; }

			//--------------------//

// need to check containment
			//TVector3 CandidateMuonStart(CandidateMu_StartX->at(0),CandidateMu_StartY->at(0),CandidateMu_StartZ->at(0));
			//TVector3 CandidateMuonEnd(CandidateMu_EndX->at(0),CandidateMu_EndY->at(0),CandidateMu_EndZ->at(0));

			//TVector3 CandidateProtonStart(CandidateP_StartX->at(0),CandidateP_StartY->at(0),CandidateP_StartZ->at(0));
			//TVector3 CandidateProtonEnd(CandidateP_EndX->at(0),CandidateP_EndY->at(0),CandidateP_EndZ->at(0));

			//--------------------//

			// Pi0 kinematics

			double pi0_p = reco_pi0_p->at(0); // GeV
			double pi0_costheta = reco_pi0_costheta->at(0);	
			double pi0_phi = reco_pi0_phi->at(0); // rad	

			TVector3 pi0_v(-1.,-1.,-1.);
			pi0_v.SetMag(pi0_p);
			pi0_v.SetPhi(pi0_phi);
			pi0_v.SetTheta(TMath::ACos(pi0_costheta));

			if (pi0_costheta < pi0_costheta_thres) { continue; }

			//--------------------//
			
			// g1 kinematics

			double g1_p = reco_g1_p->at(0); // GeV
			double g1_costheta = reco_g1_costheta->at(0);	
			double g1_phi = reco_g1_phi->at(0); // rad	

			TVector3 g1_v(-1.,-1.,-1.);
			g1_v.SetMag(g1_p);
			g1_v.SetPhi(g1_phi);
			g1_v.SetTheta(TMath::ACos(g1_costheta));

			if (g1_costheta < gamma1_costheta_thres) { continue; }			

			TVector3 g1_start(g1_start_x,g1_start_y,g1_start_z);
			TVector3 g1_end(g1_end_x,g1_end_y,g1_end_z);	
			
			if ( FVz - g1_end_z < 10 ) { continue; } // cm
			//if ( !tools.loose_inFVVector(g1_end) ) { continue; }			

			//if ( !tools.inFVVector(g1_start) ) { continue; }
			//if ( !tools.inFVVector(g1_end) ) { continue; }			
			
			//--------------------//			

			// g2 kinematics

			double g2_p = reco_g2_p->at(0); // GeV
			double g2_costheta = reco_g2_costheta->at(0);	
			double g2_phi = reco_g2_phi->at(0); // rad	

			TVector3 g2_v(-1.,-1.,-1.);
			g2_v.SetMag(g2_p);
			g2_v.SetPhi(g2_phi);
			g2_v.SetTheta(TMath::ACos(g2_costheta));

			if (g2_costheta < gamma2_costheta_thres) { continue; }		
			
			TVector3 g2_start(g2_start_x,g2_start_y,g2_start_z);
			TVector3 g2_end(g2_end_x,g2_end_y,g2_end_z);

			//if ( FVz - g2_end_z < 10 ) { continue; } // cm			
			
			//if ( !tools.inFVVector(g2_start) ) { continue; }
			//if ( !tools.inFVVector(g2_end) ) { continue; }				

			// two-shower opening angle

			double two_shower_angle = reco_shower_opening_angle->at(0); // deg

			if ( two_shower_angle > 70 ) { continue; }

			//--------------------//

			// Event selection

			// Needs to be fine-tuned
			if (pi0_costheta < pi0_costheta_thres) { continue; }
			//if (wc_nc_pio_score < -0.5) {continue;}
			//if (wc_kine_pio_flag!=2) {continue;}
			//if (wc_kine_pio_flag==0) {continue;}

			//--------------------//

			// Loop over primary pfparticles
			// Reject those events with protons above kinetic energy threshold
			// The proton threshold can be located under NCpi0/generators/constants.h 

			int primary_proton_counter = 0;
			int primary_muon_counter = 0;
			int primary_charged_pion_counter = 0;
			int primary_electron_counter = 0;
			int primary_photon_counter = 0;
			int primary_neutron_counter = 0;
			int primary_neutral_pion_counter = 0;

			int secondary_proton_counter = 0;
			int secondary_muon_counter = 0;
			int secondary_charged_pion_counter = 0;
			int secondary_electron_counter = 0;
			int secondary_photon_counter = 0;
			int secondary_neutron_counter = 0;
			int secondary_neutral_pion_counter = 0;

			int pfps = wc_reco_pdg->size();

			for (int ipfp = 0; ipfp < pfps; ipfp++ ) {

				// only primaries (mother = 0) 
				if (wc_reco_mother->at(ipfp) == 0) {

					TVector3 v_mom(wc_reco_p->at(ipfp).at(0), wc_reco_p->at(ipfp).at(1), wc_reco_p->at(ipfp).at(2)); 
					double mom = v_mom.Mag();

					// Only proton candidates
					if (wc_reco_pdg->at(ipfp) == ProtonPdg) {

						double e = TMath::Sqrt( mom*mom + ProtonMass_GeV * ProtonMass_GeV);
						double ke = e - ProtonMass_GeV;

						if (ke > proton_ke_thres) { primary_proton_counter++; }

					} // end of the primary protons

					// Only muon candidates
					else if (wc_reco_pdg->at(ipfp) == MuonPdg) {

						if (mom > 0.) { primary_muon_counter++; }
					
					} // end of the primary muons

					// Only charged pion candidates
					else if (wc_reco_pdg->at(ipfp) == AbsChargedPionPdg) {

						if (mom > 0.) { primary_charged_pion_counter++; }
					
					} // end of the primary charged pions

					// Only neutral pion candidates
					else if (wc_reco_pdg->at(ipfp) == NeutralPionPdg) {

						if (mom > 0.) { primary_neutral_pion_counter++; }
					
					} // end of the primary neutral pions

					// Only electron candidates
					else if (wc_reco_pdg->at(ipfp) == ElectronPdg) {

						if (mom > 0.07) { primary_electron_counter++; }
										
					} // end of the primary charged pions

					// Only photon candidates
					else if (wc_reco_pdg->at(ipfp) == PhotonPdg) {

						if (mom > 0.07) { primary_photon_counter++; }
										
					} // end of the primary photons

					// Only neutron candidates
					else if (wc_reco_pdg->at(ipfp) == NeutronPdg) {

						double e = TMath::Sqrt( mom*mom + NeutronMass_GeV * NeutronMass_GeV);
						double ke = e - NeutronMass_GeV;						
						//cout << "neutron mom = " << mom << " ke = " << ke << " bkg_1n_0p_1pi0_X = " << bkg_1n_0p_1pi0_X << endl;
						//if (mom > 0.) { primary_neutron_counter++; }
						if (ke > 0.01) { primary_neutron_counter++; }						
										
					} // end of the primary neutrons

					else { 
						
						cout << "primary non proton/muon/charged pion/electron candidate with pdg = " << wc_reco_pdg->at(ipfp) << " run = " << Run << "  subrun = " << SubRun << " event = " << Event << endl; 
					
					}

				} else {

					// secondary particles

					//cout << "wc_reco_mother->at(ipfp) = " << wc_reco_mother->at(ipfp) << endl;
					int mother = wc_reco_mother->at(ipfp);

					// loop over the secondary particles
					for (int ipfp_s = 0; ipfp_s < pfps; ipfp_s++ ) {

						if ( mother == wc_reco_id->at(ipfp_s) ) {

							int secondary_pdg = wc_reco_pdg->at(ipfp_s); 

							//secondary protons
							if ( TMath::Abs(secondary_pdg) == ProtonPdg) {

								TVector3 v_mom(wc_reco_p->at(ipfp_s).at(0), wc_reco_p->at(ipfp_s).at(1), wc_reco_p->at(ipfp_s).at(2)); 
								double mom = v_mom.Mag();
								double e = TMath::Sqrt( mom*mom + ProtonMass_GeV * ProtonMass_GeV);
								double ke = e - ProtonMass_GeV;

								if (ke > proton_ke_thres) { secondary_proton_counter++; }

							}

							//secondary charged pions
							else if ( TMath::Abs(secondary_pdg) == AbsChargedPionPdg) {

								secondary_charged_pion_counter++;
							
							}

							//secondary neutral pions
							else if ( TMath::Abs(secondary_pdg) == NeutralPionPdg) {

								secondary_neutral_pion_counter++;
							
							}

							//secondary muons
							else if ( TMath::Abs(secondary_pdg) == MuonPdg) {

								secondary_muon_counter++;
							
							}

							//secondary electrons
							else if ( TMath::Abs(secondary_pdg) == ElectronPdg) {

								secondary_electron_counter++;
							
							}

							//secondary photons
							else if ( TMath::Abs(secondary_pdg) == PhotonPdg) {

								secondary_photon_counter++;
							
							}

							//secondary neutron
							else if ( TMath::Abs(secondary_pdg) == NeutronPdg) {

								secondary_neutron_counter++;
							
							}

							else { 
						
								cout << "secondary non proton/muon/charged pion/electron candidate with pdg = " << wc_reco_pdg->at(ipfp_s) << " run = " << Run << "  subrun = " << SubRun << " event = " << Event << endl; 
							
							}

						} // end of grabbing the correct secondary particle

					} // end of the loop over the secondary particles

				}

			}

			// wc

			//if (primary_neutron_counter != 0) { continue; }			
			if (primary_proton_counter != 0) { continue; }
			if (primary_muon_counter != 0) { continue; }
			if (primary_charged_pion_counter != 0) { continue; }

			if (secondary_proton_counter != 0) { continue; }
			if (secondary_muon_counter != 0) { continue; }
			if (secondary_charged_pion_counter != 0) { continue; }

			if (wc_single_photon_other_score < 0 || wc_single_photon_other_score > 2.) { continue; }
			if (wc_single_photon_numu_score < -1) { continue; }
			if (wc_single_photon_ncpi0_score > 0.6) { continue; }
			if (wc_single_photon_nue_score < -2.5) { continue; }

			// cout << "primary muon counter = " << primary_muon_counter << endl;
			// cout << "primary proton counter = " << primary_proton_counter << endl;
			// cout << "primary charged pion counter = " << primary_charged_pion_counter << endl;
			//cout << "primary neutral pion counter = " << primary_neutral_pion_counter << endl;
			//cout << "primary electron counter = " << primary_electron_counter << endl;
			//cout << "primary photon counter = " << primary_photon_counter << endl;
			//cout << "primary neutron counter = " << primary_neutron_counter << endl;

			// cout << "secondary muon counter = " << secondary_muon_counter << endl;
			// cout << "secondary proton counter = " << secondary_proton_counter << endl;
			// cout << "secondary charged pion counter = " << secondary_charged_pion_counter << endl;
			// cout << "secondary electron counter = " << secondary_electron_counter << endl << endl;
			// cout << "secondary photon counter = " << secondary_photon_counter << endl << endl;
			//cout << "secondary neutron counter = " << secondary_neutron_counter << endl;

			//cout << endl;

			//--------------------//

			// Reject events that do not have two showers

			int nshowers = 0;
			int nmuontracks = 0;
			int nprotontracks = 0;
			int npiontracks = 0;

			for(int i=0; i < (int)wc_kine_energy_particle->size(); i++)
			{
				int pdgcode = wc_kine_particle_type->at(i);

				if( TMath::Abs(pdgcode) == ElectronPdg && wc_kine_energy_particle->at(i) > 10) { // KE in MeV
					
					nshowers++;

				}

				if( TMath::Abs(pdgcode) == ProtonPdg && wc_kine_energy_particle->at(i) > 0) { // KE in MeV
					
					nprotontracks++;

				}			
				
				if( TMath::Abs(pdgcode) == MuonPdg && wc_kine_energy_particle->at(i) > 10) { // KE in MeV
					
					nmuontracks++;

				}	
				
				if( TMath::Abs(pdgcode) == AbsChargedPionPdg && wc_kine_energy_particle->at(i) > 10) { // KE in MeV
					
					npiontracks++;

				}	

			}		
			
			//if (nshowers != 2) { continue; }
			//if (nprotontracks != 0) { continue; }
			if (nmuontracks != 0) { continue; }
			if (npiontracks != 0) { continue; }

			//if (reco_pi0_invmass->at(0) > 0.4) { continue; } // GeV	

			//--------------------//

			// Quality cut
			// pi0_p uses the expression with alpha/ the opening angle between the two showers
			// reco_pi0_p_gammas uses teh vector sum of the two gammas

			if ( TMath::Abs(pi0_p - reco_pi0_p_gammas->at(0)) / pi0_p * 100. > 200) { continue; }
			if (two_shower_start_dist > 100) { continue; } // cm

			pass_selection_counter++;			

			//--------------------//
	
			// Underflow / overflow
			if (pi0_p < ArrayNBinsPi0Momentum[0]) { pi0_p = (ArrayNBinsPi0Momentum[0] + ArrayNBinsPi0Momentum[1])/2.; }
			if (pi0_p > ArrayNBinsPi0Momentum[NBinsPi0Momentum]) { pi0_p = (ArrayNBinsPi0Momentum[NBinsPi0Momentum] + ArrayNBinsPi0Momentum[NBinsPi0Momentum-1])/2.; }

			if (g1_p < ArrayNBinsg1Momentum[0]) { g1_p = (ArrayNBinsg1Momentum[0] + ArrayNBinsg1Momentum[1])/2.; }
			if (g1_p > ArrayNBinsg1Momentum[NBinsg1Momentum]) { g1_p = (ArrayNBinsg1Momentum[NBinsg1Momentum] + ArrayNBinsg1Momentum[NBinsg1Momentum-1])/2.; }

			if (g2_p < ArrayNBinsg2Momentum[0]) { g2_p = (ArrayNBinsg2Momentum[0] + ArrayNBinsg2Momentum[1])/2.; }
			if (g2_p > ArrayNBinsg2Momentum[NBinsg2Momentum]) { g2_p = (ArrayNBinsg2Momentum[NBinsg2Momentum] + ArrayNBinsg2Momentum[NBinsg2Momentum-1])/2.; }

			//--------------------//

			// Selection Cuts

			//bool PassedSelection = true;

			// for (int i = 0; i < NCuts; i++) {

			// 	if (VectorCuts[i] == "_PID_NuScore" && !(reco_p_LLR_Score < ProtonLLRPIDScore) ) 
			// 		{ PassedSelection = false; }

			// 	if (VectorCuts[i] == "_PID_NuScore_CRT" && !( reco_p_LLR_Score < ProtonLLRPIDScore) ) 
			// 		{ PassedSelection = false; }


			// }

			// if (PassedSelection == false) { continue; }

			myRunTxtFile << "run = " << Run << ",  subrun = " << SubRun << ", event = " << Event << ", coh = " << coh << ", signal = " << signal << endl;
			myRunTxtFile << "vertex x = " << RecoVertex.X() << ",  y = " << RecoVertex.Y() << ", z = " << RecoVertex.Z() << endl;
			myRunTxtFile << "g1_start_x = " << g1_start_x << ",  g1_start_y = " << g1_start_y << ", g1_start_z = " << g1_start_z << endl;
			myRunTxtFile << "g1_end_x = " << g1_end_x << ",  g1_end_y = " << g1_end_y << ", g1_end_z = " << g1_end_z << endl;
			myRunTxtFile << "g2_start_x = " << g2_start_x << ",  g2_start_y = " << g2_start_y << ", g2_start_z = " << g2_start_z << endl;			
			myRunTxtFile << "g2_end_x = " << g2_end_x << ",  g2_end_y = " << g2_end_y << ", g2_end_z = " << g2_end_z << endl;						
			myRunTxtFile << "bkg_1n_0p_1pi0_X = " << bkg_1n_0p_1pi0_X << endl << endl;			
			
			//--------------------//

			int genie_mode = -1;

			// backtracked vars
			double true_pi0_p = -1;
			double true_pi0_costheta = -1;
			double true_g1_p = -1;
			double true_g1_costheta = -1;
			double true_g2_p = -1;
			double true_g2_costheta = -1;
			double true_two_shower_angle = -1;

			//----------------------------------------//

			// Only for MC to obtain true vales			
			
			if (
				string(fWhichSample).find("Overlay") != std::string::npos 
				&& MCParticle_Mode != -1 ) { 
				
				genie_mode = MCParticle_Mode; 

				//true_ECal = True_ECal->at(0);

                // // Underflow / overflow
                // if (true_ThetaVis < ArrayNBinsThetaVis[0]) { true_ThetaVis = (ArrayNBinsThetaVis[0] + ArrayNBinsThetaVis[1])/2.; }
                // if (true_ThetaVis > ArrayNBinsThetaVis[NBinsThetaVis]) { true_ThetaVis = (ArrayNBinsThetaVis[NBinsThetaVis] + ArrayNBinsThetaVis[NBinsThetaVis-1])/2.; }

			} // End of if statement: Only for MC to obtain true vales

			//----------------------------------------//

			RecoPi0MomentumPlot->Fill(pi0_p,weight);
			RecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
			RecoSingleBinPlot->Fill(0.5,weight);
			Reconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
			Reconumu_scorePlot->Fill(wc_numu_score,weight);
			Recokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
			Recog1MomentumPlot->Fill(g1_p,weight);
			Recog1CosThetaPlot->Fill(g1_costheta,weight);
			Recog2MomentumPlot->Fill(g2_p,weight);
			Recog2CosThetaPlot->Fill(g2_costheta,weight);
			Recotwo_shower_anglePlot->Fill(two_shower_angle,weight);
			Recokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
			Recosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
			Recosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
			Recosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
			Recosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

			// Blips
			ReconBlips_savedPlot->Fill(nBlips_saved,weight);
			int reco_blip_counter_radius = 0;

			for (int iblip = 0; iblip < nBlips_saved; iblip++) {

				RecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
				RecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
				RecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);												
				RecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

				TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
				double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
				Recoblip_vrtPlot->Fill(blip_vrt_dist,weight);	
				
				if (blip_vrt_dist < radius) { 
					
					RecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);
					reco_blip_counter_radius++; 
					Recoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v),weight);	
					Recoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v),weight);	
					Recoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v),weight);					
				
				}

			}

			ReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);

			//------------------------------//

			if (string(fWhichSample).find("Overlay") != std::string::npos) { 

				int inte_mode = -1;
				if (genie_mode == 0) { inte_mode = 1; } // QE 
				else if (genie_mode == 10) { inte_mode = 2; } // MEC 
				else if (genie_mode == 1) { inte_mode = 3; } // RES 
				else if (genie_mode == 2) { inte_mode = 4; } // DIS 
				else if (genie_mode == 3) { inte_mode = 5; coh_counter++; } // COH 
				else { inte_mode = 6; } // other 

				// signal

				if (signal) {

					signal_counter++;

			// 		NCCOHTruePi0CosThetaPlot->Fill(True_CandidateMu_CosTheta->at(0),weight);
			// 		NCCOHTrueSingleBinPlot->Fill(0.5,weight);
			// 		NCCOHTruePi0MomentumPlot->Fill(true_pmiss,weight);
			// Reconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
			// Reconumu_scorePlot->Fill(wc_numu_score,weight);
			// 			Recokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
			// 		Recog1MomentumPlot->Fill(g1_p,weight);
			// Recog1CosThetaPlot->Fill(g1_costheta,weight);
			// Recog2MomentumPlot->Fill(g2_p,weight);
			// Recog2CosThetaPlot->Fill(g2_costheta,weight);
			// Recotwo_shower_anglePlot->Fill(two_shower_angle,weight);
			// Recokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
			// Recosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
			// Recosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
			// Recosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
			// Recosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

			// 		// Blips
			// 		NCCOHTruenBlips_savedPlot->Fill(nBlips_saved,weight);
			// 		int nchohtrue_blip_counter_radius = 0;

			// 		for (int iblip = 0; iblip < nBlips_saved; iblip++) {

			// 			NCCOHTrueBlip_xPlot->Fill(Blip_x->at(iblip),weight);
			// 			NCCOHTrueBlip_yPlot->Fill(Blip_y->at(iblip),weight);
			// 			NCCOHTrueBlip_zPlot->Fill(Blip_z->at(iblip),weight);
			// 			NCCOHTrueBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

			// 			TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
			// 			double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
			// 			NCCOHTrueblip_vrtPlot->Fill(blip_vrt_dist,weight);	
						
			// 			if (blip_vrt_dist < radius) { 
			// 
						// 			NCCOHTrueBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);
			// nchohtrue_blip_counter_radius++; 
						//				Recoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v),weight);	
							//Recoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v),weight);	
				//Recoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v),weight);
		
		//}

			// 		}

			// 		NCCOHTruenBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);		

					//----------------------------------------//

					// 1D Reco Plots for the selected NCCOH events 

					NCCOHReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
					NCCOHReconumu_scorePlot->Fill(wc_numu_score,weight);
					NCCOHRecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
					NCCOHRecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					NCCOHRecoSingleBinPlot->Fill(0.5,weight);
					NCCOHRecoPi0MomentumPlot->Fill(pi0_p,weight);		
					NCCOHRecog1MomentumPlot->Fill(g1_p,weight);
					NCCOHRecog1CosThetaPlot->Fill(g1_costheta,weight);
					NCCOHRecog2MomentumPlot->Fill(g2_p,weight);
					NCCOHRecog2CosThetaPlot->Fill(g2_costheta,weight);
					NCCOHRecotwo_shower_anglePlot->Fill(two_shower_angle,weight);	
					NCCOHRecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);	
					NCCOHRecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					NCCOHRecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					NCCOHRecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					NCCOHRecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					NCCOHReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int nccohreco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						NCCOHRecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						NCCOHRecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						NCCOHRecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						NCCOHRecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						NCCOHRecoblip_vrtPlot->Fill(blip_vrt_dist,weight);	
						
						if (blip_vrt_dist < radius) { 
							
							NCCOHRecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);
							nccohreco_blip_counter_radius++; 
							NCCOHRecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v),weight);
							NCCOHRecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v),weight);	
							NCCOHRecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v),weight);							
						
						}

					}
					
					NCCOHReconBlips_radiusPlot->Fill(nccohreco_blip_counter_radius,weight);

					//------------------------------//

			// 							NCCOHReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
			// 		NCCOHReconumu_scorePlot->Fill(wc_numu_score,weight);
			// 					Recokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
			// 		NCCOHRecoPi0CosThetaPlot2D->Fill(True_CandidateMu_CosTheta->at(0),reco_Pmu_cos_theta);
			// 		NCCOHRecoSingleBinPlot2D->Fill(0.5,0.5);
			// 		NCCOHRecoPi0MomentumPlot2D->Fill(true_pmiss,pmiss);			
					
			// 					Recog1MomentumPlot->Fill(g1_p,weight);
			// Recog1CosThetaPlot->Fill(g1_costheta,weight);
			// Recog2MomentumPlot->Fill(g2_p,weight);
			// Recog2CosThetaPlot->Fill(g2_costheta,weight);
			// Recotwo_shower_anglePlot->Fill(two_shower_angle,weight);
			//Recokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
			// Recosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
			// Recosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
			// Recosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
			// Recosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					NCCOHReconBlips_savedPlot2D->Fill(nBlips_saved, nBlips_saved);
					NCCOHReconBlips_radiusPlot2D->Fill(nccohreco_blip_counter_radius,nccohreco_blip_counter_radius,weight);

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						NCCOHRecoBlip_xPlot2D->Fill(Blip_x->at(iblip),Blip_x->at(iblip));
						NCCOHRecoBlip_yPlot2D->Fill(Blip_y->at(iblip),Blip_y->at(iblip));
						NCCOHRecoBlip_zPlot2D->Fill(Blip_z->at(iblip),Blip_z->at(iblip));
						NCCOHRecoBlip_energyPlot2D->Fill(Blip_energy->at(iblip),Blip_energy->at(iblip));
						NCCOHRecoBlip_proxtrkdistPlot2D->Fill(Blip_proxtrkdist->at(iblip),Blip_proxtrkdist->at(iblip));

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						NCCOHRecoblip_vrtPlot2D->Fill(blip_vrt_dist,blip_vrt_dist);
						NCCOHRecoblip_cos_alphapi0Plot2D->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v),tools.CosAlpha(blip_3d, RecoVertex, pi0_v));
						NCCOHRecoblip_cos_alphag1Plot2D->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v),tools.CosAlpha(blip_3d, RecoVertex, g1_v));	
						NCCOHRecoblip_cos_alphag2Plot2D->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v),tools.CosAlpha(blip_3d, RecoVertex, g2_v));

					}						

					// Blips
					POTScaledNCCOHReconBlips_savedPlot2D->Fill(nBlips_saved, nBlips_saved,weight);
					POTScaledNCCOHReconBlips_radiusPlot2D->Fill(nccohreco_blip_counter_radius,nccohreco_blip_counter_radius,weight);

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						POTScaledNCCOHRecoBlip_xPlot2D->Fill(Blip_x->at(iblip),Blip_x->at(iblip),weight);
						POTScaledNCCOHRecoBlip_yPlot2D->Fill(Blip_y->at(iblip),Blip_y->at(iblip),weight);
						POTScaledNCCOHRecoBlip_zPlot2D->Fill(Blip_z->at(iblip),Blip_z->at(iblip),weight);
						POTScaledNCCOHRecoBlip_energyPlot2D->Fill(Blip_energy->at(iblip),Blip_energy->at(iblip),weight);
						POTScaledNCCOHRecoBlip_proxtrkdistPlot2D->Fill(Blip_proxtrkdist->at(iblip),Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						POTScaledNCCOHRecoblip_vrtPlot2D->Fill(blip_vrt_dist, blip_vrt_dist,weight);
						POTScaledNCCOHRecoblip_cos_alphapi0Plot2D->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v),tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);
						POTScaledNCCOHRecoblip_cos_alphag1Plot2D->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v),tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
						POTScaledNCCOHRecoblip_cos_alphag2Plot2D->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v),tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);

					}						

				} // End of the NCCOH signal

				//----------------------------------------//

				// Non-NCCOH beam related background or EXT BNB

				else {

					//----------------------------------------//	

					if (bkg_0pi0_X) { bkg_0pi0_X_counter++; }
					else if (bkg_Mpi0_X) { bkg_Mpi0_X_counter++; } 
					else if (bkg_bwds_1pi0_X) { bkg_bwds_1pi0_X_counter++; }
					else if (bkg_1n_0p_1pi0_X) { bkg_1n_0p_1pi0_X_counter++; }
					else if (bkg_Nn_0p_1pi0_X) { bkg_Nn_0p_1pi0_X_counter++; }
					else if (bkg_1p_0n_1pi0_X) { bkg_1p_0n_1pi0_X_counter++; }
					else if (bkg_Np_0n_1pi0_X) { bkg_Np_0n_1pi0_X_counter++; }
					else if (bkg_1pi0_Npipm_X) { bkg_1pi0_Npipm_X_counter++; }		
					else if (bkg_1pi0_Np_Nn_0pipm_X) { bkg_1pi0_Np_Nn_0pipm_X_counter++; }
					else if (bkg_1pi0_Np_Nn_Npipm_X) { bkg_1pi0_Np_Nn_Npipm_X_counter++; }
					else if (bkg_1pi0_Nmh_X) { bkg_1pi0_Nmh_X_counter++; }		
					else if (bkg_1pi0_Nl_X) { bkg_1pi0_Nl_X_counter++; }
					else { bkg_other_counter++; }

					//----------------------------------------//					

					NonNCCOHReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
					NonNCCOHReconumu_scorePlot->Fill(wc_numu_score,weight);
					NonNCCOHRecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
					NonNCCOHRecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					NonNCCOHRecoSingleBinPlot->Fill(0.5,weight);
					NonNCCOHRecoPi0MomentumPlot->Fill(pi0_p,weight);	
					NonNCCOHRecog1MomentumPlot->Fill(g1_p,weight);
					NonNCCOHRecog1CosThetaPlot->Fill(g1_costheta,weight);
					NonNCCOHRecog2MomentumPlot->Fill(g2_p,weight);
					NonNCCOHRecog2CosThetaPlot->Fill(g2_costheta,weight);
					NonNCCOHRecotwo_shower_anglePlot->Fill(two_shower_angle,weight);
					NonNCCOHRecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
					NonNCCOHRecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					NonNCCOHRecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					NonNCCOHRecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					NonNCCOHRecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					NonNCCOHReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int nonnccohreco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						NonNCCOHRecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						NonNCCOHRecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						NonNCCOHRecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						NonNCCOHRecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						NonNCCOHRecoblip_vrtPlot->Fill(blip_vrt_dist,weight);					

						if (blip_vrt_dist < radius) { 
							
							NonNCCOHRecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);
							nonnccohreco_blip_counter_radius++; 
							NonNCCOHRecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);
							NonNCCOHRecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
							NonNCCOHRecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);	
						
						}

					}		
					
					NonNCCOHReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);

					//------------------------------//

				} // End of the Non-NCCOH beam related background

				//----------------------------------------//

				// QE

				if (genie_mode == 0) {

					QERecoPi0MomentumPlot->Fill(pi0_p,weight);
					QERecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					QERecoSingleBinPlot->Fill(0.5,weight);	
					QEReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
					QEReconumu_scorePlot->Fill(wc_numu_score,weight);
					QERecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);	
					QERecog1MomentumPlot->Fill(g1_p,weight);
					QERecog1CosThetaPlot->Fill(g1_costheta,weight);
					QERecog2MomentumPlot->Fill(g2_p,weight);
					QERecog2CosThetaPlot->Fill(g2_costheta,weight);
					QERecotwo_shower_anglePlot->Fill(two_shower_angle,weight);	
					QERecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);	
					QERecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					QERecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					QERecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					QERecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);	

					// Blips
					QEReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int qereco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						QERecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						QERecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						QERecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						QERecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						QERecoblip_vrtPlot->Fill(blip_vrt_dist,weight);							

						if (blip_vrt_dist < radius) { 
							
							QERecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);
							qereco_blip_counter_radius++; 
							QERecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);
							QERecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
							QERecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);
													
						}

					}			
					
					QEReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);

				} // End of QE selection

				//----------------------------------------//

				// MEC

				if (genie_mode == 10) {

					MECRecoPi0MomentumPlot->Fill(pi0_p,weight);
					MECRecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					MECRecoSingleBinPlot->Fill(0.5,weight);	
					MECReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);	
					MECReconumu_scorePlot->Fill(wc_numu_score,weight);
					MECRecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
					MECRecog1MomentumPlot->Fill(g1_p,weight);
					MECRecog1CosThetaPlot->Fill(g1_costheta,weight);
					MECRecog2MomentumPlot->Fill(g2_p,weight);
					MECRecog2CosThetaPlot->Fill(g2_costheta,weight);
					MECRecotwo_shower_anglePlot->Fill(two_shower_angle,weight);		
					MECRecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
					MECRecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					MECRecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					MECRecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					MECRecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					MECReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int mecreco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						MECRecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						MECRecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						MECRecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						MECRecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						MECRecoblip_vrtPlot->Fill(blip_vrt_dist,weight);						

						if (blip_vrt_dist < radius) { 
							
							MECRecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);
							mecreco_blip_counter_radius++; 
							MECRecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);
							MECRecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
							MECRecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);
													
						}

					}		
					
					MECReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);
		
				}

				//----------------------------------------//

				// RES

				if (genie_mode == 1) {

					RESRecoPi0MomentumPlot->Fill(pi0_p,weight);
					RESRecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					RESRecoSingleBinPlot->Fill(0.5,weight);		
					RESReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);
					RESReconumu_scorePlot->Fill(wc_numu_score,weight);
					RESRecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
					RESRecog1MomentumPlot->Fill(g1_p,weight);
					RESRecog1CosThetaPlot->Fill(g1_costheta,weight);
					RESRecog2MomentumPlot->Fill(g2_p,weight);
					RESRecog2CosThetaPlot->Fill(g2_costheta,weight);
					RESRecotwo_shower_anglePlot->Fill(two_shower_angle,weight);
					RESRecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
					RESRecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					RESRecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					RESRecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					RESRecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					RESReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int resreco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						RESRecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						RESRecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						RESRecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						RESRecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						RESRecoblip_vrtPlot->Fill(blip_vrt_dist,weight);						
						
						if (blip_vrt_dist < radius) { 
							
							RESRecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);							
							resreco_blip_counter_radius++; 
							RESRecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);	
							RESRecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
							RESRecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);	
						
						}

					}		
					
					RESReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);
	
				}

				//----------------------------------------//

				// DIS

				if (genie_mode == 2) {

					DISRecoPi0MomentumPlot->Fill(pi0_p,weight);
					DISRecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					DISRecoSingleBinPlot->Fill(0.5,weight);		
					DISReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);	
					DISReconumu_scorePlot->Fill(wc_numu_score,weight);
					DISRecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
					DISRecog1MomentumPlot->Fill(g1_p,weight);
					DISRecog1CosThetaPlot->Fill(g1_costheta,weight);
					DISRecog2MomentumPlot->Fill(g2_p,weight);
					DISRecog2CosThetaPlot->Fill(g2_costheta,weight);
					DISRecotwo_shower_anglePlot->Fill(two_shower_angle,weight);		
					DISRecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
					DISRecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					DISRecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					DISRecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					DISRecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					DISReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int disreco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						DISRecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						DISRecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						DISRecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						DISRecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						DISRecoblip_vrtPlot->Fill(blip_vrt_dist,weight);

						if (blip_vrt_dist < radius) { 
							
							DISRecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);							
							disreco_blip_counter_radius++; 
							DISRecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);
							DISRecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
							DISRecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);						
						
						}

					}	
					
					DISReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);

				}

				//----------------------------------------//

				// COH

				if (genie_mode == 3) {

					COHRecoPi0MomentumPlot->Fill(pi0_p,weight);
					COHRecoPi0CosThetaPlot->Fill(pi0_costheta,weight);
					COHRecoSingleBinPlot->Fill(0.5,weight);		
					COHReconc_pio_scorePlot->Fill(wc_nc_pio_score,weight);	
					COHReconumu_scorePlot->Fill(wc_numu_score,weight);	
					COHRecokine_pio_flagPlot->Fill(wc_kine_pio_flag,weight);
					COHRecog1MomentumPlot->Fill(g1_p,weight);
					COHRecog1CosThetaPlot->Fill(g1_costheta,weight);
					COHRecog2MomentumPlot->Fill(g2_p,weight);
					COHRecog2CosThetaPlot->Fill(g2_costheta,weight);
					COHRecotwo_shower_anglePlot->Fill(two_shower_angle,weight);		
					COHRecokine_pio_vtx_disPlot->Fill(wc_kine_pio_vtx_dis,weight);
					COHRecosingle_photon_numu_scorePlot->Fill(wc_single_photon_numu_score,weight);
					COHRecosingle_photon_other_scorePlot->Fill(wc_single_photon_other_score,weight);
					COHRecosingle_photon_ncpi0_scorePlot->Fill(wc_single_photon_ncpi0_score,weight);
					COHRecosingle_photon_nue_scorePlot->Fill(wc_single_photon_nue_score,weight);

					// Blips
					COHReconBlips_savedPlot->Fill(nBlips_saved,weight);
					int cohreco_blip_counter_radius = 0;

					for (int iblip = 0; iblip < nBlips_saved; iblip++) {

						COHRecoBlip_xPlot->Fill(Blip_x->at(iblip),weight);
						COHRecoBlip_yPlot->Fill(Blip_y->at(iblip),weight);
						COHRecoBlip_zPlot->Fill(Blip_z->at(iblip),weight);
						COHRecoBlip_proxtrkdistPlot->Fill(Blip_proxtrkdist->at(iblip),weight);

						TVector3 blip_3d(Blip_x->at(iblip),Blip_y->at(iblip),Blip_z->at(iblip) );
						double blip_vrt_dist = (RecoVertex - blip_3d).Mag();
						COHRecoblip_vrtPlot->Fill(blip_vrt_dist,weight);

						if (blip_vrt_dist < radius) { 
							
							COHRecoBlip_energyPlot->Fill(Blip_energy->at(iblip),weight);							
							cohreco_blip_counter_radius++; 
							COHRecoblip_cos_alphapi0Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, pi0_v), weight);
							COHRecoblip_cos_alphag1Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g1_v), weight);	
							COHRecoblip_cos_alphag2Plot->Fill( tools.CosAlpha(blip_3d, RecoVertex, g2_v), weight);
						
						}

					}
					
					COHReconBlips_radiusPlot->Fill(reco_blip_counter_radius,weight);

				}


			} // End of the Overlay case and the breakdown into NCCOH/NonNCCOH & QE,MEC,RES,DIS

		} // End of the loop over the events

		std::cout << std::endl << "Created file: " << FileName << std::endl << std::endl;
		std::cout << "candidate events: " << pass_selection_counter << " [" << std::setprecision(2) << double(pass_selection_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "coh events: " << coh_counter << " [" << std::setprecision(2) << double(coh_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "signal events: " << signal_counter << " [" << std::setprecision(2) << double(signal_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl << std::endl;

		std::cout << "bkg_0pi0_X events: " << bkg_0pi0_X_counter << " [" << std::setprecision(2) << double(bkg_0pi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_Mpi0_X events: " << bkg_Mpi0_X_counter << " [" << std::setprecision(2) << double(bkg_Mpi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_bwds_1pi0_X events: " << bkg_bwds_1pi0_X_counter << " [" << std::setprecision(2) << double(bkg_bwds_1pi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1n_0p_1pi0_X events: " << bkg_1n_0p_1pi0_X_counter << " [" << std::setprecision(2) << double(bkg_1n_0p_1pi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_Nn_0p_1pi0_X events: " << bkg_Nn_0p_1pi0_X_counter << " [" << std::setprecision(2) << double(bkg_Nn_0p_1pi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1p_0n_1pi0_X events: " << bkg_1p_0n_1pi0_X_counter << " [" << std::setprecision(2) << double(bkg_1p_0n_1pi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_Np_0n_1pi0_X events: " << bkg_Np_0n_1pi0_X_counter << " [" << std::setprecision(2) << double(bkg_Np_0n_1pi0_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1pi0_Npipm_X events: " << bkg_1pi0_Npipm_X_counter << " [" << std::setprecision(2) << double(bkg_1pi0_Npipm_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1pi0_Np_Nn_0pipm_X events: " << bkg_1pi0_Np_Nn_0pipm_X_counter << " [" << std::setprecision(2) << double(bkg_1pi0_Np_Nn_0pipm_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1pi0_Np_Nn_Npipm_X events: " << bkg_1pi0_Np_Nn_Npipm_X_counter << " [" << std::setprecision(2) << double(bkg_1pi0_Np_Nn_Npipm_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1pi0_Nmh_X events: " << bkg_1pi0_Nmh_X_counter << " [" << std::setprecision(2) << double(bkg_1pi0_Nmh_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
		std::cout << "bkg_1pi0_Nl_X events: " << bkg_1pi0_Nl_X_counter << " [" << std::setprecision(2) << double(bkg_1pi0_Nl_X_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;								
		std::cout << "bkg_other events: " << bkg_other_counter << " [" << std::setprecision(2) << double(bkg_other_counter)/double(pass_selection_counter) *100. <<"%]" << std::endl;
	
		//----------------------------------------//	

		file->cd();
		file->Write();
		file->Close();
		fFile->Close();

		//----------------------------------------//	

//	} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

} // End of the program