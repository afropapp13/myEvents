{

	vector<TString> WhichSampleArray;
	vector<TString> EventWeightLabels;
	vector<int> Universes;
	int NGENIEUniverses = 100;

	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("Overlay9_Run1");
//	WhichSampleArray.push_back("Overlay9_Run2");
	WhichSampleArray.push_back("Overlay9_Run3");
//	WhichSampleArray.push_back("Overlay9_Run4");
//	WhichSampleArray.push_back("Overlay9_Run5");

	// -----------------------------------------------------------------------------------------

	EventWeightLabels.push_back("AxFFCCQEshape_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("DecayAngMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NormCCCOH_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NormNCCOH_UBGenie"); Universes.push_back(2);
//	EventWeightLabels.push_back("RPA_CCQE_Reduced_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("RPA_CCQE_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("ThetaDelta2NRad_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("Theta_Delta2Npi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("VecFFCCQEshape_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("XSecShape_CCMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("All_UBGenie"); Universes.push_back(NGENIEUniverses);

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx++");
	//gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L myRecoAnalysis.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(EventWeightLabels.size()); j++) {

			for (int k = 0; k < Universes[j]; k++) {	

				gROOT->ProcessLine("myRecoAnalysis(\""+WhichSampleArray[i]+"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

				if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
				  { gROOT->ProcessLine("myTrueAnalysis(\""+WhichSampleArray[i]+"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()"); } 

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
