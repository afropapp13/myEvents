{

	gROOT->ProcessLine(".L ../../../generators/Tools.cxx++");
	gROOT->ProcessLine(".L mcc9_10_reco_selection.cxx++");

	//--------------------//

	//Flux uncertainties/Flux

	vector<TString> FluxWhichSampleArray;
	vector<TString> FluxEventWeightLabels;
	vector<int> FluxUniverses;
	int NFluxUniverses = 100;

	FluxWhichSampleArray.push_back("mcc9_10_Overlay9_Run4b_unified");
	FluxWhichSampleArray.push_back("mcc9_10_OverlayDirt9_Run4b_unified");	

	FluxEventWeightLabels.push_back("fluxes"); 
	FluxUniverses.push_back(NFluxUniverses);

	for (int i = 0;i < (int)(FluxWhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(FluxEventWeightLabels.size()); j++) {

			for (int k = 0; k < FluxUniverses[j]; k++) {	

				gROOT->ProcessLine("mcc9_10_reco_selection(\""+FluxWhichSampleArray[i]+"\",\"\",\""+FluxEventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	//--------------------//

	// G4 uncertainties

	vector<TString> G4WhichSampleArray;
	vector<TString> G4EventWeightLabels;
	vector<int> G4Universes;
	int NG4Universes = 100;

	G4WhichSampleArray.push_back("mcc9_10_Overlay9_Run4b_unified");
	G4WhichSampleArray.push_back("mcc9_10_OverlayDirt9_Run4b_unified");	

	G4EventWeightLabels.push_back("reinteractions"); 
	G4Universes.push_back(NG4Universes);

	for (int i = 0;i < (int)(G4WhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(G4EventWeightLabels.size()); j++) {

			for (int k = 0; k < G4Universes[j]; k++) {	

				gROOT->ProcessLine("mcc9_10_reco_selection(\""+G4WhichSampleArray[i]+"\",\"\",\""+G4EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	//--------------------//
	
	vector<TString> XSecWhichSampleArray;
	vector<TString> XSecEventWeightLabels;
	vector<int> XSecUniverses;
	int NXSecUniverses = 100;

	XSecWhichSampleArray.push_back("mcc9_10_Overlay9_Run4b_unified");
	XSecWhichSampleArray.push_back("mcc9_10_OverlayDirt9_Run4b_unified");	

	XSecEventWeightLabels.push_back("AxFFCCQEshape_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("DecayAngMEC_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("NormCCCOH_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("NormNCCOH_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("RPA_CCQE_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("ThetaDelta2NRad_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("Theta_Delta2Npi_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("VecFFCCQEshape_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("XSecShape_CCMEC_UBGenie"); XSecUniverses.push_back(2);
	XSecEventWeightLabels.push_back("All_UBGenie"); XSecUniverses.push_back(NXSecUniverses);

	for (int i = 0;i < (int)(XSecWhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(XSecEventWeightLabels.size()); j++) {

			for (int k = 0; k < XSecUniverses[j]; k++) {	

				gROOT->ProcessLine("mcc9_10_reco_selection(\""+XSecWhichSampleArray[i]+"\",\"\",\""+XSecEventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	}
	
	//--------------------//
	
	vector<TString> MCStatWhichSampleArray;
	vector<TString> MCStatEventWeightLabels;
	vector<int> MCStatUniverses;
	int NMCStatUniverses = 100;

	MCStatWhichSampleArray.push_back("mcc9_10_Overlay9_Run4b_unified");
	MCStatWhichSampleArray.push_back("mcc9_10_OverlayDirt9_Run4b_unified");	

	MCStatEventWeightLabels.push_back("MC_Stat"); 
	MCStatUniverses.push_back(NMCStatUniverses);

	for (int i = 0;i < (int)(MCStatWhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(MCStatEventWeightLabels.size()); j++) {

			for (int k = 0; k < MCStatUniverses[j]; k++) {	

				gROOT->ProcessLine("mcc9_10_reco_selection(\""+MCStatWhichSampleArray[i]+"\",\"\",\""+MCStatEventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	//--------------------//

	vector<TString> FDSWhichSampleArray;
			
	FDSWhichSampleArray.push_back("mcc9_10_Overlay9_Run4b_unified");	
	FDSWhichSampleArray.push_back("mcc9_10_OverlayDirt9_Run4b_unified");		

	for (int i = 0;i < (int)(FDSWhichSampleArray.size()); i++) {

		gROOT->ProcessLine("mcc9_10_reco_selection(\""+FDSWhichSampleArray[i]+"\",\"RS\").Loop()");

	}	

	//--------------------//
	
	gROOT->ProcessLine(".q");

	//--------------------//
	
}