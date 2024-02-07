{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// NuWro

	WhichSampleArray.push_back("Overlay9NuWro_Run1");
	WhichSampleArray.push_back("Overlay9NuWro_Run2");
	WhichSampleArray.push_back("Overlay9NuWro_Run3");
	WhichSampleArray.push_back("Overlay9NuWro_Run4b");
	WhichSampleArray.push_back("Overlay9NuWro_Run4c");
	WhichSampleArray.push_back("Overlay9NuWro_Run4d");
	WhichSampleArray.push_back("Overlay9NuWro_Combined");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L reco_selection.cxx++");
	gROOT->ProcessLine(".L true_selection.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("reco_selection(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("true_selection(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
