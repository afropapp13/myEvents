{
	const int NFiles = 5;
	TString File[NFiles]  = {"MCC84","Overlay","Adam","BNBPlusEXTBNB","EXTBNB"};
	TString TrackingName = "pandoraNu";
	TString VertexingName = "pandoraNu";

	for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {

		// Open the corresponding file
		TFile *_file0 = TFile::Open("./SamplesToUse/Full"+File[WhichFile]+"1500.root");

		cout << endl << File[WhichFile]+"1500.root is being created" << endl << endl;

		// Open the corresponding tree
		TTree* oldtree = (TTree*)(_file0->Get("analysistree/anatree"));

		// Deactivate all the branches
		oldtree->SetBranchStatus("*",0);

		// Activate only those that you want to copy to the new tree
//		oldtree -> SetBranchStatus("event",1);
//		oldtree -> SetBranchStatus("potbnb",1);
//		oldtree -> SetBranchStatus("nfls_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsTime_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsPe_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsPePerOpDet_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsZcenter_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsYcenter_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsYwidth_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsZwidth_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("flsTwidth_simpleFlashBeam",1);
//		oldtree -> SetBranchStatus("nfls_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsTime_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsPe_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsPePerOpDet_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsZcenter_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsYcenter_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsYwidth_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsZwidth_simpleFlashCosmic",1);
//		oldtree -> SetBranchStatus("flsTwidth_simpleFlashCosmic",1);
		oldtree -> SetBranchStatus("isdata",1);
		oldtree -> SetBranchStatus("mcevts_truth",1);
		oldtree -> SetBranchStatus("mctrk_TrackId",1);
		oldtree -> SetBranchStatus("mctrk_pdg",1);
		oldtree -> SetBranchStatus("mctrk_origin",1);
		oldtree -> SetBranchStatus("mctrk_Motherpdg",1);
		oldtree -> SetBranchStatus("nuvtxx_truth",1);
		oldtree -> SetBranchStatus("nuvtxy_truth",1);
		oldtree -> SetBranchStatus("nuvtxz_truth",1);
		oldtree -> SetBranchStatus("ccnc_truth",1);
		oldtree -> SetBranchStatus("nuPDG_truth",1);
		oldtree -> SetBranchStatus("enu_truth",1);
		oldtree -> SetBranchStatus("pdg",1);
		oldtree -> SetBranchStatus("no_mctracks",1);
		oldtree -> SetBranchStatus("mode_truth",1);
		oldtree -> SetBranchStatus("lep_mom_truth",1);
		oldtree -> SetBranchStatus("lep_dcosz_truth",1);
		oldtree -> SetBranchStatus("geant_list_size",1);
		oldtree -> SetBranchStatus("no_primaries",1);
//		oldtree -> SetBranchStatus("StartPointx",1);
//		oldtree -> SetBranchStatus("StartPointy",1);
//		oldtree -> SetBranchStatus("StartPointz",1);
//		oldtree -> SetBranchStatus("EndPointx",1);
//		oldtree -> SetBranchStatus("EndPointy",1);
//		oldtree -> SetBranchStatus("EndPointz",1);
		oldtree -> SetBranchStatus("theta",1);
		oldtree -> SetBranchStatus("phi",1); 
		oldtree -> SetBranchStatus("P",1);
		oldtree -> SetBranchStatus("Px",1);
		oldtree -> SetBranchStatus("Py",1);
		oldtree -> SetBranchStatus("Pz",1);
		oldtree -> SetBranchStatus("status",1);      
		oldtree -> SetBranchStatus("origin",1); 
//		oldtree -> SetBranchStatus("process_primary",1);           
		oldtree -> SetBranchStatus("TrackId",1);
//		oldtree -> SetBranchStatus("mctrk_len_drifted",1);
//		oldtree -> SetBranchStatus("pathlen_drifted",1);
//		oldtree -> SetBranchStatus("nshowers_showerrecopandora",1);
//		oldtree -> SetBranchStatus("showerID_showerrecopandora",1);
//		oldtree -> SetBranchStatus("shwr_startx_showerrecopandora",1);
//		oldtree -> SetBranchStatus("shwr_starty_showerrecopandora",1);
//		oldtree -> SetBranchStatus("shwr_startz_showerrecopandora",1);    
		oldtree -> SetBranchStatus(("ntracks_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkg4id_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trklen_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkId_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkke_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkpidpdg_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkstartx_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkstarty_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkstartz_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkendx_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkendy_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkendz_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trktheta_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkphi_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkorigin_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkidtruth_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkpdgtruth_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkpidbestplane_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkmomrange_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkmom_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkpidpida_"+TrackingName),1);
		oldtree -> SetBranchStatus(("nvtx_"+VertexingName),1);
		oldtree -> SetBranchStatus(("vtxx_"+VertexingName),1);
		oldtree -> SetBranchStatus(("vtxy_"+VertexingName),1);
		oldtree -> SetBranchStatus(("vtxz_"+VertexingName),1);
		oldtree -> SetBranchStatus(("ntrkhits_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkefftruth_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkpurtruth_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkpurity_"+TrackingName),1);
//		oldtree -> SetBranchStatus(("trkcompleteness_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkdedx_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkdqdx_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkresrg_"+TrackingName),1);

		oldtree -> SetBranchStatus(("trkstartx_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkstarty_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkstartz_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkstartd_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkendx_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkendy_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkendz_"+TrackingName),1);
		oldtree -> SetBranchStatus(("trkendd_"+TrackingName),1);

		oldtree -> SetBranchStatus("no_hits",1);
		oldtree -> SetBranchStatus("no_hits_stored",1);
		oldtree -> SetBranchStatus("hit_plane",1);
		oldtree -> SetBranchStatus("hit_wire",1);
//		oldtree -> SetBranchStatus("hit_channel",1);
//		oldtree -> SetBranchStatus("hit_peakT",1);
		oldtree -> SetBranchStatus("hit_charge",1);
		oldtree -> SetBranchStatus("hit_ph",1);
//		oldtree -> SetBranchStatus("hit_startT",1);
//		oldtree -> SetBranchStatus("hit_endT",1);

		oldtree -> SetBranchStatus(("nnuvtx_"+TrackingName),1);		
		oldtree -> SetBranchStatus(("nuvtxId_"+TrackingName),1);
		oldtree -> SetBranchStatus(("nuvtxx_"+TrackingName),1);
		oldtree -> SetBranchStatus(("nuvtxy_"+TrackingName),1);	
		oldtree -> SetBranchStatus(("nuvtxz_"+TrackingName),1);	
		oldtree -> SetBranchStatus(("nuvtxpdg_"+TrackingName),1);	
		oldtree -> SetBranchStatus(("nuvtxhasPFParticle_"+TrackingName),1);	
		oldtree -> SetBranchStatus(("nuvtxPFParticleID_"+TrackingName),1);	

		oldtree -> SetBranchStatus("nPFParticles",1);
		oldtree -> SetBranchStatus("pfp_selfID",1);
		oldtree -> SetBranchStatus("pfp_isPrimary",1);
		oldtree -> SetBranchStatus("pfp_numDaughters",1);
		oldtree -> SetBranchStatus("pfp_daughterIDs",1);
		oldtree -> SetBranchStatus("pfp_parentID",1);
		oldtree -> SetBranchStatus("pfp_vertexID",1);
		oldtree -> SetBranchStatus("pfp_isShower",1);
		oldtree -> SetBranchStatus("pfp_isTrack",1);
		//oldtree -> SetBranchStatus("pfp_isVertex",1);
		oldtree -> SetBranchStatus("pfp_trackID",1);
		oldtree -> SetBranchStatus("pfp_showerID",1);
		oldtree -> SetBranchStatus("pfp_pdgCode",1);
		oldtree -> SetBranchStatus("pfp_numClusters",1);
		oldtree -> SetBranchStatus("pfp_clusterIDs",1);
		oldtree -> SetBranchStatus("pfp_isNeutrino",1);
		oldtree -> SetBranchStatus("pfp_numNeutrinos",1);
		oldtree -> SetBranchStatus("pfp_neutrinoIDs",1);

		oldtree -> SetBranchStatus("trkcosmicscore_tagger_pandoraNu",1);
		oldtree -> SetBranchStatus("trkcosmicscore_flashmatch_pandoraNu",1);


		// Create a new file to save the new tree
		TFile* myFile = new TFile("./SamplesToUse/"+File[WhichFile]+"1500.root","RECREATE");
		// Clone the tree
		TTree* newtree = (TTree*)oldtree->CloneTree(15000); // Only 15000 events, set to 0 to copy all the events
		// Write the new tree in the new file
		newtree->Print();
		// Write the file
		myFile->Write();
		// Clean Up

		cout << endl << File[WhichFile]+"1500.root has been created" << endl << endl;

		delete _file0;
		delete myFile;
	};
};
