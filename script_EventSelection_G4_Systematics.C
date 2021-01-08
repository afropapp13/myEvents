{

	vector<TString> WhichSampleArray;
	vector<TString> EventWeightLabels;
	vector<int> Universes;
//	int NG4Universes = 1000;	
	int NG4Universes = 100;

	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("Overlay9_Run1");
//	WhichSampleArray.push_back("Overlay9_Run2");
	WhichSampleArray.push_back("Overlay9_Run3");
//	WhichSampleArray.push_back("Overlay9_Run4");
//	WhichSampleArray.push_back("Overlay9_Run5");				
	
	// -----------------------------------------------------------------------------------------
	
	EventWeightLabels.push_back("reinteractions_piminus_Geant4"); Universes.push_back(NG4Universes);
	EventWeightLabels.push_back("reinteractions_piplus_Geant4"); Universes.push_back(NG4Universes);
	EventWeightLabels.push_back("reinteractions_proton_Geant4"); Universes.push_back(NG4Universes);
//	EventWeightLabels.push_back("xsr_scc_Fa3_SCC"); Universes.push_back(10);
//	EventWeightLabels.push_back("xsr_scc_Fv3_SCC"); Universes.push_back(10);		

	// -----------------------------------------------------------------------------------------

//	On the gpvm's'
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L myRecoAnalysis.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");

	for (int i = 0; i < (int)(WhichSampleArray.size()); i++) {
	
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
