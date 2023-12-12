{

	vector<TString> WhichSampleArray;
	vector<TString> EventWeightLabels;
	vector<int> Universes;
//	int MCStatUniverses = 1000;
	int MCStatUniverses = 100;

	// -----------------------------------------------------------------------------------------

//	WhichSampleArray.push_back("Overlay9_Run1");
//	WhichSampleArray.push_back("Overlay9_Run2");
//	WhichSampleArray.push_back("Overlay9_Run3");
//	WhichSampleArray.push_back("Overlay9_Run4");
//	WhichSampleArray.push_back("Overlay9_Run5");
	WhichSampleArray.push_back("Overlay9_Combined");

	// -----------------------------------------------------------------------------------------

	EventWeightLabels.push_back("MC_Stat"); Universes.push_back(MCStatUniverses);

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L reco_selection.cxx+");
	gROOT->ProcessLine(".L true_selection.cxx+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(EventWeightLabels.size()); j++) {

			for (int k = 0; k < Universes[j]; k++) {	

				gROOT->ProcessLine("reco_selection(\""+WhichSampleArray[i]+"\",\"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

				if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) { 
					
					gROOT->ProcessLine("true_selection(\""+WhichSampleArray[i]+"\",\"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");
					
				} 

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
