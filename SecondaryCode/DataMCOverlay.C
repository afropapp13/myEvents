#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TTree.h>

// using namespace std;

//This defines our current settings for the fiducial volume (in cm)

double FVx = 256.35;
double FVy = 233;
double FVz = 1036.8;
double borderx = 10.;
double bordery = 20.;
double borderz = 10.;

//This function returns if a 3D point is within the fiducial volume
bool inFV(double x, double y, double z) {

	if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
	else return false;
};

bool inTPC(double x, double y, double z) {
  
	if(x<-50 || x>300) return false;
	if(y<-100 || y>100) return false;
	if(z<0.01 || z>1050) return false;
	else return true;
};

//This function returns the distance between a flash and a track (in one dimension, here used only for z direction)
double FlashTrackDist(double flash, double start, double end) {

	if(end >= start) {
	
		if(flash < end && flash > start) return 0;
		else return TMath::Min(fabs(flash-start), fabs(flash-end));
	} else {
	
		if(flash > end && flash < start) return 0;
		else return TMath::Min(fabs(flash-start), fabs(flash-end));
	};
};

int Topology(int nnumuCC, int nnueCC, int nNC, int cosmicflag, int OOFVflag) {

	//// This function return the true topology of the event, numu & anti-numu
	////1. numu CC
	////2. nue CC
	////3. NC
	////4. OOFV
	////5. cosmic
	////6. Other

	if (cosmicflag > 0 ) return 5;
	if (OOFVflag > 0 ) return 4;
	if (cosmicflag == 0 && OOFVflag == 0 && (nnumuCC + nnueCC) ==0 && nNC >0) return 3;
	if (cosmicflag== 0 && OOFVflag == 0 &&nNC == 0 && nnueCC == 0 && nnumuCC >0) return 1; 
	if (cosmicflag== 0 && OOFVflag == 0 &&nNC == 0 && nnueCC > 0 && nnumuCC == 0) return 2; 

	else return 6;
};

// Main function
int DataMCOverlay() {

	vector <string> samples;
	samples.push_back("Overlay");
	//samples.push_back("OverlayOld");
	//samples.push_back("MCC84");

	// Opening a text file to write purities & efficiencies
	ofstream outfile;
	outfile.open("/uboone/app/users/adi/cosmicMuons/cosmicMuonsOverlay/Eff_DataCosmicUnbiased_BNBMC.txt");

	// Initialize and fill track reco product names
	string TrackingName = "pandoraNu";
	string VertexingName = "pandoraNu";

 
	for (int isample = 0 ; isample < samples.size() ; isample++) {   
	
		TChain *treenc = new TChain("analysistree/anatree");

		//treenc -> Add("/pnfs/uboone/data/uboone/root-tuple/anatree_reco_extbnb_data_mcc8/anatree/prod_v06_26_01_06/00/00/63/40/ana_hist_ffad46ca-89a6-4521-b854-533cfa34ac1a.root");
		//TChain ch("MeritTuple")
		//TFileCollection fc("dum","","/uboone/app/users/adi/cosmicMuons/cosmicMuonsOverlay/newFiles.txt");
//		TFileCollection fc("dum","",Form("/uboone/app/users/adi/cosmicMuons/cosmicMuonsOverlay/filesList%s.txt",samples[isample].c_str()));

		TFileCollection fc("dum","",Form("/uboone/app/users/apapadop/LArSoft_v06_26_01_09/MyFirstList.txt"));

		treenc->AddFileInfoList(fc.GetList());

		//maximum array sizes
		const int maxentries = 50000;
		const int maxtracks = 10000;
		const int maxvtx = 5000;
		const int maxnu = 10;
		const int maxmc = 10;
		const int kMaxFlashes = 5000;
		const int maxshower = 1000;
		const int kMaxHits = 40000;

		//Define variables to read from Tree
		
		Int_t           event;
		Int_t           run;
		Int_t           subrun;
		Int_t           triggerbits; //this is only important for data. 2048 = BNB stream
		Double_t        potbnb; //this is only important for data
		Char_t          isdata;            //flag, 0=MC 1=data

		Short_t         ntracks_reco; //number of reco tracks
		Short_t         trkbestplane[maxtracks]; //plane that has most hits for a given track
		Float_t         trklen[maxtracks]; //track length
		Float_t         trkstartx[maxtracks];
		Float_t         trkstarty[maxtracks];
		Float_t         trkstartz[maxtracks];
		Float_t         trkendx[maxtracks];
		Float_t         trkendy[maxtracks];
		Float_t         trkendz[maxtracks];
		Float_t         trktheta[maxtracks];
		Float_t         trkphi[maxtracks];
		Float_t         trkmomrange[maxtracks]; //track momentum calculated from track range
		Short_t         trkId[maxtracks];
		Int_t           trkpdgtruth[maxtracks][3];////
		Short_t         trkorigin[maxtracks][3]; //for MC only: which true particle contributes most hits to the reco track: 2 = cosmic. 1 = neutrino
		Int_t           TrackIDTruth[maxtracks][3]; // MC id matched with reco track
		bool            vertexatstart[maxtracks]; //for analysis: is the vertex at start of the track?
		bool            vertexatend[maxtracks]; //for analysis: ist the vertex at end of track?

		Short_t         ntrkhits[maxtracks];
		Float_t         trkefftruth[maxtracks]; //completeness
		Float_t         trkpurtruth[maxtracks]; //purity of track
		Float_t         trkpurity[maxtracks];    // track purity based on hit information
		Float_t         trkcompleteness[maxtracks]; //track completeness based on hit information
		Float_t         trkpidpida[maxtracks][3];
		Float_t         trkdedx[maxtracks][3];
		Float_t         trkdqdx[maxtracks][3];

		//// hit info///

		Int_t           no_hits;                  //number of hits
		Int_t           hit_multiplicity;                  //number of hits
		Int_t           no_hits_stored;           //number of hits actually stored in the tree    
		Short_t         hit_plane[kMaxHits];      //plane number
		Short_t         hit_wire[kMaxHits];       //wire number
		Short_t         hit_channel[kMaxHits];    //channel ID
		Float_t         hit_peakT[kMaxHits];      //peak time
		Float_t         hit_charge[kMaxHits];     //charge (area)
		Float_t         hit_ph[kMaxHits];         //amplitude
		Float_t         hit_startT[kMaxHits];     //hit start time
		Float_t         hit_endT[kMaxHits];       //hit end time

		////
		
		Short_t         nshowers;
		Short_t         showerID;
		Float_t         shwr_startx[maxshower];
		Float_t         shwr_starty[maxshower];
		Float_t         shwr_startz[maxshower];
    
		Short_t         nfls_simpleFlashBeam; //number of flashes in the event
		Short_t         nfls_simpleFlashCosmic; //number of flashes in the event
		Float_t         flsTime_simpleFlashBeam[kMaxFlashes]; //flash time (in microseconds)
		Float_t         flsTime_simpleFlashCosmic[kMaxFlashes]; //flash time (in microseconds)
		Float_t         flsPe_simpleFlashBeam[kMaxFlashes]; //total number of photoelectrons corresponding to the flash
		Float_t         flsPe_simpleFlashCosmic[kMaxFlashes]; //total number of photoelectrons corresponding to the flash
		Float_t         flsZcenter_simpleFlashBeam[kMaxFlashes]; //z center of flash (in cm)
		Float_t         flsYcenter_simpleFlashBeam[kMaxFlashes]; //y center of flash (in cm)
		Float_t         flsXcenter_simpleFlashBeam[kMaxFlashes]; //x center of flash (in cm)
		Float_t         flsZcenter_simpleFlashCosmic[kMaxFlashes]; //z center of flash (in cm)
		Float_t         flsYcenter_simpleFlashCosmic[kMaxFlashes]; //y center of flash (in cm)
		Float_t         flsXcenter_simpleFlashCosmic[kMaxFlashes]; //x center of flash (in cm)
		Double_t        flsPePerOpDet_simpleFlashBeam[kMaxFlashes][32]; //PE per Optical detector
		Double_t        flsPePerOpDet_simpleFlashCosmic[kMaxFlashes][32]; //PE per Optical detector
		Float_t         flsZwidth_simpleFlashBeam[kMaxFlashes]; //z width of flash (in cm)
		Float_t         flsYwidth_simpleFlashBeam[kMaxFlashes]; //y width of flash (in cm)
		//Float_t         flash_xwidth[kMaxFlashes]; //x width of flash (in cm)
		Float_t         flsTwidth_simpleFlashBeam[kMaxFlashes]; //t width of flash (in ticks? how many us?)
		Float_t         flsZwidth_simpleFlashCosmic[kMaxFlashes]; //z width of flash (in cm)
		Float_t         flsYwidth_simpleFlashCosmic[kMaxFlashes]; //y width of flash (in cm)
		//Float_t         flash_xwidth[kMaxFlashes]; //x width of flash (in cm)
		Float_t         flsTwidth_simpleFlashCosmic[kMaxFlashes]; //t width of flash (in ticks? how many us?)

		Short_t         nvtx;
		Float_t         vtxx[maxvtx];
		Float_t         vtxy[maxvtx];
		Float_t         vtxz[maxvtx];

		//finding candidate nu interaction vertex in event
	
		bool            candvertex[maxvtx];
		Short_t         candtrack[maxvtx];
		Float_t         candlength[maxvtx];
		Short_t         numuvertex = -1;
		Short_t         mutrack = -1;
		Float_t         mutracklength = 0;

		//MC truth
		
		Int_t           mcevts_truth; //neutrino interactions per event
		Float_t         nuvtxx_truth[maxmc]; //true vertex x (in cm)
		Float_t         nuvtxy_truth[maxmc];
		Float_t         nuvtxz_truth[maxmc];
		Int_t           ccnc_truth[maxmc]; //CC = 0, NC = 1
		Int_t           nuPDG_truth[maxmc]; //true neutrino pdg code. numu = 14
		Int_t           mode_truth[maxmc]; //QE = 0, RES = 1, DIS = 2
		Int_t	        PDG_truth[maxtracks];
		Float_t         enu_truth[maxmc];
		Float_t         lep_mom_truth[maxmc];
		Float_t         lep_dcosz_truth[maxmc];

		Int_t           NumberOfMCTracks;

		//geant_list_size
		
		Float_t         XMCTrackStart[maxtracks];
		Float_t         YMCTrackStart[maxtracks];
		Float_t         ZMCTrackStart[maxtracks];

		Float_t         XMCTrackEnd[maxtracks];
		Float_t         YMCTrackEnd[maxtracks];
		Float_t         ZMCTrackEnd[maxtracks];
		Float_t         MCPx[maxtracks];
		Float_t         MCPy[maxtracks];
		Float_t         MCPz[maxtracks];
		Float_t         MCtheta[maxtracks];
		Float_t         MCphi[maxtracks];
		Float_t         MCmom[maxtracks];
		Int_t           MCstatus[maxtracks];//// MC tracks that are actually tracked (otherwise they are too low energetic, how much?) 
		Int_t           MCorigin[maxtracks];
		Int_t           MCprocess_primary[maxtracks];     
		Int_t           MCTrackID[maxtracks];
		Float_t         MCLenght[maxtracks];
		Float_t         MCpathlen[maxtracks];

		//define cut variables
		
		double flashwidth = 80; //cm. Distance flash-track
		double distcut = 5; //cm. Distance track start/end to vertex/// *** check it!!!
		double lengthcut = 75; //cm. Length of longest track/// ***check it!!
		double beammin = 3.2/*-0.36*/; //us. Beam window start///***only BNB?
		double beammax = 4.8/*-0.36*/; //us. Beam window end///***only BNB?
		double PEthresh = 50; //PE
		double MCTrackToMCVtxDist = 1; //cm. distance between mc track start and mc vertex///***what does it mean?
		double TrackToMCDist = 5; //cm. Distance track start/end to mcvertex///***what does it mean?

		treenc -> SetBranchAddress("event", &event);
		treenc -> SetBranchAddress("potbnb", &potbnb);
		
		////MCC8
		////using simpleFlash algorithm only
		
		treenc -> SetBranchAddress("nfls_simpleFlashBeam", &nfls_simpleFlashBeam);
		treenc -> SetBranchAddress("flsTime_simpleFlashBeam", flsTime_simpleFlashBeam);
		treenc -> SetBranchAddress("flsPe_simpleFlashBeam", flsPe_simpleFlashBeam);
		treenc -> SetBranchAddress("flsPePerOpDet_simpleFlashBeam", flsPePerOpDet_simpleFlashBeam);
		//treenc -> SetBranchAddress("flsXcenter_simpleFlashBeam", flsXcenter_simpleFlashBeam);
		treenc -> SetBranchAddress("flsZcenter_simpleFlashBeam", flsZcenter_simpleFlashBeam);
		treenc -> SetBranchAddress("flsYcenter_simpleFlashBeam", flsYcenter_simpleFlashBeam);
		//treenc -> SetBranchAddress("flsXwidth_simpleFlashBeam", flash_xwidth);
		treenc -> SetBranchAddress("flsYwidth_simpleFlashBeam", flsYwidth_simpleFlashBeam);
		treenc -> SetBranchAddress("flsZwidth_simpleFlashBeam", flsZwidth_simpleFlashBeam);
		treenc -> SetBranchAddress("flsTwidth_simpleFlashBeam", flsTwidth_simpleFlashBeam);

		treenc -> SetBranchAddress("nfls_simpleFlashCosmic", &nfls_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsTime_simpleFlashCosmic", flsTime_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsPe_simpleFlashCosmic", flsPe_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsPePerOpDet_simpleFlashCosmic", flsPePerOpDet_simpleFlashCosmic);
		//treenc -> SetBranchAddress("flsXcenter_simpleFlashCosmic", flsXcenter_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsZcenter_simpleFlashCosmic", flsZcenter_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsYcenter_simpleFlashCosmic", flsYcenter_simpleFlashCosmic);
		//treenc -> SetBranchAddress("flsXwidth_simpleFlashCosmic", flash_xwidth);
		treenc -> SetBranchAddress("flsYwidth_simpleFlashCosmic", flsYwidth_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsZwidth_simpleFlashCosmic", flsZwidth_simpleFlashCosmic);
		treenc -> SetBranchAddress("flsTwidth_simpleFlashCosmic", flsTwidth_simpleFlashCosmic);
		
		//////
    
		//// MCC7
		/*
		treenc -> SetBranchAddress("no_flashes", &no_flashes);
		treenc -> SetBranchAddress("flash_time", flash_time);
		treenc -> SetBranchAddress("flash_pe", flash_pe);
		treenc -> SetBranchAddress("flash_zcenter", flash_zcenter);
		treenc -> SetBranchAddress("flash_ycenter", flash_ycenter);
		*/

		///////
    
		treenc -> SetBranchAddress("isdata", &isdata);
		treenc -> SetBranchAddress("mcevts_truth", &mcevts_truth);
		treenc -> SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
		treenc -> SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
		treenc -> SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
		treenc -> SetBranchAddress("ccnc_truth", ccnc_truth);
		treenc -> SetBranchAddress("nuPDG_truth", nuPDG_truth);
		treenc -> SetBranchAddress("enu_truth", enu_truth);
		treenc -> SetBranchAddress("pdg", PDG_truth);
		treenc -> SetBranchAddress("mode_truth", mode_truth);
		treenc -> SetBranchAddress("lep_mom_truth", lep_mom_truth);
		treenc -> SetBranchAddress("lep_dcosz_truth", lep_dcosz_truth);
		treenc -> SetBranchAddress("geant_list_size", &NumberOfMCTracks);
		treenc -> SetBranchAddress("StartPointx", XMCTrackStart);
		treenc -> SetBranchAddress("StartPointy", YMCTrackStart);
		treenc -> SetBranchAddress("StartPointz", ZMCTrackStart);
		treenc -> SetBranchAddress("EndPointx", XMCTrackEnd);
		treenc -> SetBranchAddress("EndPointy", YMCTrackEnd);
		treenc -> SetBranchAddress("EndPointz", ZMCTrackEnd);
		treenc -> SetBranchAddress("theta", MCtheta);
		treenc -> SetBranchAddress("phi", MCphi); 
		treenc -> SetBranchAddress("P", MCmom);
		treenc -> SetBranchAddress("Px",MCPx);
		treenc -> SetBranchAddress("Py",MCPy);
		treenc -> SetBranchAddress("Pz",MCPz);
		treenc -> SetBranchAddress("status", MCstatus);      
		treenc -> SetBranchAddress("origin", MCorigin); 
		treenc -> SetBranchAddress("process_primary", MCprocess_primary);           
		treenc -> SetBranchAddress("TrackId", MCTrackID);
		treenc -> SetBranchAddress("mctrk_len_drifted", MCLenght);
		treenc -> SetBranchAddress("pathlen_drifted", MCpathlen);

		// Product specific stuff
 
		treenc -> SetBranchAddress("nshowers_showerrecopandora", &nshowers);
		treenc -> SetBranchAddress("showerID_showerrecopandora", &showerID);
		treenc -> SetBranchAddress("shwr_startx_showerrecopandora", shwr_startx);
		treenc -> SetBranchAddress("shwr_starty_showerrecopandora", shwr_starty);
		treenc -> SetBranchAddress("shwr_startz_showerrecopandora", shwr_startz);    
		treenc -> SetBranchAddress(("ntracks_"+TrackingName).c_str(),&ntracks_reco);
		treenc -> SetBranchAddress(("trklen_"+TrackingName).c_str(), trklen);
		treenc -> SetBranchAddress(("trkstartx_"+TrackingName).c_str(),trkstartx);
		treenc -> SetBranchAddress(("trkstarty_"+TrackingName).c_str(),trkstarty);
		treenc -> SetBranchAddress(("trkstartz_"+TrackingName).c_str(),trkstartz);
		treenc -> SetBranchAddress(("trkendx_"+TrackingName).c_str(),trkendx);
		treenc -> SetBranchAddress(("trkendy_"+TrackingName).c_str(),trkendy);
		treenc -> SetBranchAddress(("trkendz_"+TrackingName).c_str(),trkendz);
		treenc -> SetBranchAddress(("trktheta_"+TrackingName).c_str(),trktheta);
		treenc -> SetBranchAddress(("trkphi_"+TrackingName).c_str(),trkphi);
		treenc -> SetBranchAddress(("trkorigin_"+TrackingName).c_str(),trkorigin);
		treenc -> SetBranchAddress(("trkidtruth_"+TrackingName).c_str(),TrackIDTruth);
		treenc -> SetBranchAddress(("trkpdgtruth_"+TrackingName).c_str(),trkpdgtruth);
		treenc -> SetBranchAddress(("trkpidbestplane_"+TrackingName).c_str(), trkbestplane);
		treenc -> SetBranchAddress(("trkmomrange_"+TrackingName).c_str(), trkmomrange);///*** mom by range
		treenc -> SetBranchAddress(("trkpidpida_"+TrackingName).c_str(), trkpidpida);

		treenc -> SetBranchAddress(("nvtx_"+VertexingName).c_str(), &nvtx);
		treenc -> SetBranchAddress(("vtxx_"+VertexingName).c_str(), vtxx);
		treenc -> SetBranchAddress(("vtxy_"+VertexingName).c_str(), vtxy);
		treenc -> SetBranchAddress(("vtxz_"+VertexingName).c_str(), vtxz);

		/// hit info
		
		treenc -> SetBranchAddress(("ntrkhits_"+TrackingName).c_str(), ntrkhits);
		treenc -> SetBranchAddress(("trkefftruth_"+TrackingName).c_str(), trkefftruth);
		treenc -> SetBranchAddress(("trkpurtruth_"+TrackingName).c_str(), trkpurtruth);
		treenc -> SetBranchAddress(("trkpurity_"+TrackingName).c_str(), trkpurity);
		treenc -> SetBranchAddress(("trkcompleteness_"+TrackingName).c_str(), trkcompleteness);
		treenc -> SetBranchAddress(("trkdedx_"+TrackingName).c_str(), trkdedx);
		treenc -> SetBranchAddress(("trkdqdx_"+TrackingName).c_str(), trkdqdx);

		treenc -> SetBranchAddress("no_hits", &no_hits);
		treenc -> SetBranchAddress("no_hits_stored", &no_hits_stored);
		treenc -> SetBranchAddress("hit_plane", hit_plane);
		treenc -> SetBranchAddress("hit_wire", hit_wire);
		treenc -> SetBranchAddress("hit_channel", hit_channel);
		treenc -> SetBranchAddress("hit_peakT", hit_peakT);
		treenc -> SetBranchAddress("hit_charge", hit_charge);
		treenc -> SetBranchAddress("hit_ph", hit_ph);
		treenc -> SetBranchAddress("hit_startT", hit_startT);
		treenc -> SetBranchAddress("hit_endT", hit_endT);
    
		//////
		//unsigned long int Size = treenc -> GetEntries();
		unsigned long int Size = 1500;

		int theflash = -1;

		double diststart = 0;
		double distend = 0;
		double TrackRange = 0;

		int numbertrk = 0;

		double TotalPOT = 0.0;
		bool signal = false;
		bool intpc = false;


		TH1D *no_flashes_evt_cosmicData = new TH1D(Form("no_flashes_evt_cosmicData%s",samples[isample].c_str()), ";number of flashes;events",5,0,5);
		TH1D *no_flashes_evt_BNBMC = new TH1D(Form("no_flashes_evt_BNBMC%s",samples[isample].c_str()), ";number of flashes;events",5,0,5);
		TH1D *no_flashes_evt_timeCut_BNBMC = new TH1D(Form("no_flashes_evt_timeCut_BNBMC%s",samples[isample].c_str()), ";number of flashes;events",5,0,5);

		TH1D *flashPE_cosmicData = new TH1D(Form("flashPE_cosmicData%s",samples[isample].c_str()), ";flash PE;events",100,0,15000);
		TH1D *flashPE_BNBMC = new TH1D(Form("flashPE_BNBMC%s",samples[isample].c_str()), ";flash PE;events",100,0,15000);
		TH1D *flashPE_timeCut_BNBMC = new TH1D(Form("flashPE_timeCut_BNBMC%s",samples[isample].c_str()), ";flash PE;events",100,0,15000);

		TH1D *flash_time_cosmicData = new TH1D(Form("flash_time_cosmicData%s",samples[isample].c_str()), "flash time",100,-6000,6000);
		TH1D *flash_time_BNBMC = new TH1D(Form("flash_time_BNBMC%s",samples[isample].c_str()), "flash time",200,0,30);

		TH2D *flashPE_vs_time_BNBMC = new TH2D(Form("flashPE_vs_time_BNBMC%s",samples[isample].c_str()), ";flash PE;flash time",100,0,1000,100,0,30);
		TH2D *flash_PE_vs_OptDet_BNBMC = new TH2D(Form("flash_PE_vs_OptDet_BNBMC%s",samples[isample].c_str()), ";flash PE;#PMT",32,0,32,100,0,100);
		TH2D *flash_time_vs_OptDet_BNBMC = new TH2D(Form("flash_time_vs_OptDet_BNBMC%s",samples[isample].c_str()), ";flash time;#PMT",32,0,32,32,0,32);
		TH2D *flashPE_vs_time_cosmicData = new TH2D(Form("flashPE_vs_time_cosmicData%s",samples[isample].c_str()), ";flash PE; flash time",100,0,1000,100,0,30);
		TH2D *flash_PE_vs_OptDet_cosmicData = new TH2D(Form("flash_PE_vs_OptDet_cosmicData%s",samples[isample].c_str()), ";flash PE;#PMT",32,0,32,100,0,100);
		TH2D *flash_time_vs_OptDet_cosmicData = new TH2D(Form("flash_time_vs_OptDet_cosmicData%s",samples[isample].c_str()), ";flash time;#PMT",32,0,32,32,0,32);

		TH1D *flashZcenter_cosmicData = new TH1D(Form("flashZcenter_cosmicData%s",samples[isample].c_str()), "; flash Z center;events",100,0,1000);
		TH1D *flashYcenter_cosmicData = new TH1D(Form("flashYcenter_cosmicData%s",samples[isample].c_str()), "; flash Y center;events",100,-60,60);
		TH1D *flashZwidth_cosmicData = new TH1D(Form("flashZwidth_cosmicData%s",samples[isample].c_str()), "; flash Z width;events",100,0,500);

		/////Checking the flashes with weird timing in MCC84
		
		TH1D *flashZwidth_ZcenterCut_cosmicData = new TH1D(Form("flashZwidth_ZCenterCut_cosmicData%s",samples[isample].c_str()), "; flash Z width;events",100,0,500);
		TH2D *flashtime_flashZwidth_cosmicData = new TH2D(Form("flashtime_flashZwidth_cosmicData%s",samples[isample].c_str()),"; flash time ; flash Z width ", 100, -2500,2500,100,0,500);
		TH2D *flashtime_flashZcenter_cosmicData = new TH2D(Form("flashtime_flashZcenter_cosmicData%s",samples[isample].c_str()),"; flash time ; flash Z center ", 100, -2500,2500,100,0,1000);
		TH1D *flashYwidth_cosmicData = new TH1D(Form("flashYwidth_cosmicData%s",samples[isample].c_str()), "; flash Y width;events",100,0,100);
		TH1D *flashZwidth_BNBMC = new TH1D(Form("flashZwidth_BNBMC%s",samples[isample].c_str()), "; flash Z width;events",100,0,500);
		TH1D *flashYwidth_BNBMC = new TH1D(Form("flashYwidth_BNBMC%s",samples[isample].c_str()), "; flash Y width;events",100,0,100);
		TH1D *flashZcenter_BNBMC = new TH1D(Form("flashZcenter_BNBMC%s",samples[isample].c_str()), "; flash Z center;events",100,0,1000);
		TH1D *flashYcenter_BNBMC = new TH1D(Form("flashYcenter_BNBMC%s",samples[isample].c_str()), "; flash Y center;events",100,-60,60);
		TH1D *flashTwidth_BNBMC = new TH1D(Form("flashTwidth_BNBMC%s",samples[isample].c_str()), "; flash T width;events",100,3.98,4.02);
		TH1D *flashZwidth_timeCut_BNBMC = new TH1D(Form("flashZwidth_timeCut_BNBMC%s",samples[isample].c_str()), "; flash Z width;events",100,0,500);
		TH1D *flashYwidth_timeCut_BNBMC = new TH1D(Form("flashYwidth_timeCut_BNBMC%s",samples[isample].c_str()), "; flash Y width;events",100,0,100);
		TH1D *flashZcenter_timeCut_BNBMC = new TH1D(Form("flashZcenter_timeCut_BNBMC%s",samples[isample].c_str()), "; flash Z center;events",100,0,1000);
		TH1D *flashYcenter_timeCut_BNBMC = new TH1D(Form("flashYcenter_timeCut_BNBMC%s",samples[isample].c_str()), "; flash Y center;events",100,-60,60);
		TH1D *flashTwidth_timeCut_BNBMC = new TH1D(Form("flashTwidth_timeCut_BNBMC%s",samples[isample].c_str()), "; flash T width;events",100,3.98,4.02);
		TH1D *flashTwidth_cosmicData = new TH1D(Form("flashTwidth_cosmicData%s",samples[isample].c_str()), "; flash T width;events",100,0.27,0.29);
		TH1D *hitPeakT = new TH1D(Form("hit_peakT%s",samples[isample].c_str()), "hit peakT",100,0,7000);
		TH1D *hitPeakTLow = new TH1D(Form("hit_peakTLow%s",samples[isample].c_str()), "hit peakT",100,0,100);
		TH1D *hitMultiplicity = new TH1D(Form("hit_multiplicity%s",samples[isample].c_str()), "hit multiplicity",100,0,100000);
		TH1D *hitCharge = new TH1D(Form("hit_charge%s",samples[isample].c_str()), "hit charge",100,-100,1400);

		TH1D *hitAmplitudPlane0 = new TH1D(Form("hitAmplitud_plane0%s",samples[isample].c_str()), "hit amplitud",60,0,120);
		TH1D *hitAmplitudPlane1 = new TH1D(Form("hitAmplitud_plane1%s",samples[isample].c_str()), "hit amplitud",60,0,120);
		//TH1D *hitAmplitudPlane2 = new TH1D(Form("hitAmplitud_plane2%s",samples[isample].c_str()), "hit amplitud",60,0,120);
		
		TH1D *mcevts_truth_Plot = new TH1D(Form("mcevts_truth%s",samples[isample].c_str()), "neutrino interactions per event",100,0,100);		
		TH1D *lep_mom_truth_Plot = new TH1D(Form("lep_mom_truth%s",samples[isample].c_str()), "lepton momentum",100,0,10);

		//////////

		Int_t NCC1piG4=0;
		Int_t NCC0pi2protonsG4=0;
		Int_t NCCG4=0;
		Int_t NCC1pG4=0;
		Int_t NCCNpG4=0;
		Int_t CCGENIE=0;

		Int_t NCC0Pion0Proton=0;
		Int_t NCC0Pion1Proton=0;
		Int_t NCC0Pion2Proton=0;
		Int_t NCC0PionNProton=0;
		Int_t NCC1PionNProton=0;
		Int_t NCCNPionNProton=0;
		Int_t NCCNue=0;
		Int_t NNC=0;
		Int_t NOOFV=0;
		Int_t NEXT=0;
		Int_t NOther=0;
    
		TFile *f1 = new TFile(Form("/uboone/data/users/apapadop/MCC84_bnb_nu_uB_overlay/myPlots/histos%s.root",samples[isample].c_str()), "RECREATE");
		TList *list = new TList();

		list->Add(no_flashes_evt_cosmicData);
		list->Add(no_flashes_evt_BNBMC);
		list->Add(no_flashes_evt_timeCut_BNBMC);
		list->Add(flashPE_cosmicData);
		list->Add(flashPE_BNBMC);
		list->Add(flashPE_timeCut_BNBMC);
		list->Add(flash_time_cosmicData);
		list->Add(flash_time_BNBMC);
		list->Add(flashPE_vs_time_cosmicData);
		list->Add(flashPE_vs_time_BNBMC);
		list->Add(flash_PE_vs_OptDet_cosmicData);
		list->Add(flash_time_vs_OptDet_cosmicData);
		list->Add(flash_PE_vs_OptDet_BNBMC);
		list->Add(flash_time_vs_OptDet_BNBMC);
		list->Add(flashZcenter_cosmicData);
		list->Add(flashYcenter_cosmicData);
		list->Add(flashZwidth_cosmicData);
		list->Add(flashZwidth_ZcenterCut_cosmicData);
		list->Add(flashtime_flashZwidth_cosmicData);
		list->Add(flashtime_flashZcenter_cosmicData);
		list->Add(flashYwidth_cosmicData);
		list->Add(flashTwidth_cosmicData);
		list->Add(flashZcenter_BNBMC);
		list->Add(flashYcenter_BNBMC);
		list->Add(flashZwidth_BNBMC);
		list->Add(flashYwidth_BNBMC);
		list->Add(flashTwidth_BNBMC);
		list->Add(flashZcenter_timeCut_BNBMC);
		list->Add(flashYcenter_timeCut_BNBMC);
		list->Add(flashZwidth_timeCut_BNBMC);
		list->Add(flashYwidth_timeCut_BNBMC);
		list->Add(flashTwidth_timeCut_BNBMC);
		list->Add(hitPeakT);
		list->Add(hitPeakTLow);
		list->Add(hitMultiplicity);
		list->Add(hitCharge);
		list->Add(hitAmplitudPlane0);
		list->Add(hitAmplitudPlane1);
		//list->Add(hitAmplitudPlane2);
		list->Add(mcevts_truth_Plot);
		list->Add(lep_mom_truth_Plot);

		bool properFile = false;

		//Event Loop
		TFile* currentFile = nullptr;

		for(int i = 0; i < Size; i++) {
		
			if(i%1000 == 0) std::cout << "\t... " << i << std::endl;
			if(i<0) continue;
			// Get tree entries
			properFile = false;
			TFile* thisFile = treenc->GetFile();
    
			if (thisFile != currentFile) {
			
				currentFile = thisFile;
				std::cout << "New file '" << currentFile->GetPath() << "' at entry #" << i << std::endl;
				properFile = true;
			};
			
			treenc -> GetEntry(i);
			if(!treenc) continue;

			bool flashtag = false;
			float flashmax = 0;

			if(isdata == 0) {          //flag, 0=MC 1=data

				if(inTPC(nuvtxx_truth[0],nuvtxy_truth[0],nuvtxz_truth[0])) intpc = true;
				else intpc = false;
			};

			////////Topology nnumu vs nue vs NC vs cosmic vs OOFV
     
			bool CCevt = false;
			bool NCevt = false;
			bool numuevt = false;
			bool nueevt = false;
			bool OOFVevt = false;
			bool cosmicevt = false;

			if(isdata == 1) cosmicevt = true;/// cosmic data unbiased

			else if(isdata == 0) {

				if (mcevts_truth == 0)  std::cout<<"cosmic event in BNB MC"<<std::endl;
				if (!intpc) OOFVevt = true;
				
				if( intpc && ccnc_truth[0] == 0) {

					CCevt = true;
					if(abs(nuPDG_truth[0]) == 14) numuevt = true;
					if(abs(nuPDG_truth[0]) == 12) nueevt = true;
				};
				
				if(intpc && ccnc_truth[0] == 1) NCevt = true;
			};
     
			//cosmicevt = true;
			int topo_truth = -1;
 
			topo_truth = Topology(numuevt, nueevt, NCevt, cosmicevt, OOFVevt);

			/////// *******flash info for data & MC and overlay dataMC******** //////

			no_flashes_evt_BNBMC->Fill(nfls_simpleFlashBeam);
			no_flashes_evt_cosmicData->Fill(nfls_simpleFlashCosmic);
			//std::cout<<" number of flashes "<< no_flashes<<std::endl;
			int goodTimingFlash =0;
			
			for (int f = 0; f < nfls_simpleFlashBeam; f++) {
			
				flash_time_BNBMC->Fill(flsTime_simpleFlashBeam[f]);
				flashPE_BNBMC->Fill(flsPe_simpleFlashBeam[f]);
				flashZcenter_BNBMC->Fill(flsZcenter_simpleFlashBeam[f]);
				flashYcenter_BNBMC->Fill(flsYcenter_simpleFlashBeam[f]);
				flashZwidth_BNBMC->Fill(flsZwidth_simpleFlashBeam[f]);
				flashYwidth_BNBMC->Fill(flsYwidth_simpleFlashBeam[f]);
				flashTwidth_BNBMC->Fill(flsTwidth_simpleFlashBeam[f]);
				
				if (flsTime_simpleFlashBeam[f] < 23) { 
				
					flashPE_timeCut_BNBMC->Fill(flsPe_simpleFlashBeam[f]);
				 	flashZcenter_timeCut_BNBMC->Fill(flsZcenter_simpleFlashBeam[f]);
				 	flashYcenter_timeCut_BNBMC->Fill(flsYcenter_simpleFlashBeam[f]);
				 	flashZwidth_timeCut_BNBMC->Fill(flsZwidth_simpleFlashBeam[f]);
				 	flashYwidth_timeCut_BNBMC->Fill(flsYwidth_simpleFlashBeam[f]);
				 	flashTwidth_timeCut_BNBMC->Fill(flsTwidth_simpleFlashBeam[f]);
					goodTimingFlash++;
				};
				
				flashPE_vs_time_BNBMC->Fill(flsPe_simpleFlashBeam[f],flsTime_simpleFlashBeam[f]);
	
				for(int opt =0; opt< 32; opt++) {
				
					flash_PE_vs_OptDet_BNBMC->Fill(flsPe_simpleFlashBeam[f],flsPePerOpDet_simpleFlashBeam[f][opt]);
					flash_time_vs_OptDet_BNBMC->Fill(flsTime_simpleFlashBeam[f],flsPePerOpDet_simpleFlashBeam[f][opt]);
				};
			};
			
			no_flashes_evt_timeCut_BNBMC->Fill(goodTimingFlash);
			
			for(int f = 0; f < nfls_simpleFlashCosmic; f++) {
			
				flash_time_cosmicData->Fill(flsTime_simpleFlashCosmic[f]);
				flashPE_cosmicData->Fill(flsPe_simpleFlashCosmic[f]);
				flashPE_vs_time_cosmicData->Fill(flsPe_simpleFlashCosmic[f],flsTime_simpleFlashCosmic[f]);
				flashZcenter_cosmicData->Fill(flsZcenter_simpleFlashCosmic[f]);
				flashYcenter_cosmicData->Fill(flsYcenter_simpleFlashCosmic[f]);
				flashZwidth_cosmicData->Fill(flsZwidth_simpleFlashCosmic[f]);
				flashtime_flashZwidth_cosmicData->Fill(flsTime_simpleFlashCosmic[f],flsZwidth_simpleFlashCosmic[f]);
				flashtime_flashZcenter_cosmicData->Fill(flsTime_simpleFlashCosmic[f],flsZcenter_simpleFlashCosmic[f]);
				if (flsZcenter_simpleFlashCosmic[f] > 200) flashZwidth_ZcenterCut_cosmicData->Fill(flsZwidth_simpleFlashCosmic[f]);
				flashYwidth_cosmicData->Fill(flsYwidth_simpleFlashCosmic[f]);
				flashTwidth_cosmicData->Fill(flsTwidth_simpleFlashCosmic[f]);
	
				for(int opt = 0; opt < 32; opt++) {
		
					flash_PE_vs_OptDet_cosmicData->Fill(flsPe_simpleFlashBeam[f],flsPePerOpDet_simpleFlashCosmic[f][opt]);
					flash_time_vs_OptDet_cosmicData->Fill(flsTime_simpleFlashBeam[f],flsPePerOpDet_simpleFlashCosmic[f][opt]);
				};
			};
			
// ______________________________________________________________________________________________________________________________________________________________________________________________________

			for (int WhichMC = 0; WhichMC < mcevts_truth; WhichMC ++) {
			
				lep_mom_truth_Plot->Fill(lep_mom_truth[WhichMC]);
			};
			
			mcevts_truth_Plot->Fill(mcevts_truth);			

			////////////////////////////////////////////////////////////
			//////////////////////// Hit info///////////////////////////
			////////////////////////////////////////////////////////////


			// If the flash tag is true and we have POT/// event has a flash in the time window
			// for the moment, all events will pass flashtag (let's see later if we want to add a cut (beam window, min. PE,..
			flashtag = true;

			if(flashtag) {
			 
				///// loop over the hits
				for(int ihit = 0; ihit < no_hits; ihit++) {
				
					hitPeakT->Fill(hit_peakT[ihit]);
					hitPeakTLow->Fill(hit_peakT[ihit]);
					hitCharge->Fill(hit_charge[ihit]);

					/////per plane

					if (hit_plane[ihit] == 0) { hitAmplitudPlane0->Fill(hit_ph[ihit],hit_wire[ihit]); }
					else if (hit_plane[ihit] == 1) { hitAmplitudPlane1->Fill(hit_ph[ihit],hit_wire[ihit]); }
					//else if (hit_plane[ihit] == 2) { hitAmplitudPlane2->Fill(hit_ph[ihit],hit_wire[ihit]); }

				}; // end loop over hits
 
				hitMultiplicity->Fill(no_hits);

				////plane0_  /// per track!! trkdqdx[maxtracks][3];
			}; // end flash tag
		}; // end loop on events

		outfile << "--------------------------------------------------------------------------------------------" << std::endl;
		outfile << std::endl;
		outfile << "Track Reco Product Name : " << TrackingName << "  Vertex Reco Product Name : " << VertexingName << std::endl;
		outfile << "Number of events we are running over: " << Size << std::endl;

		//////

		//// write histos!!!
		f1->cd();
		list->Write();
		f1->Write();
		f1->Close();

		treenc -> ResetBranchAddresses();

	}; //end loop on sample
	return 0;
}; // end main function
