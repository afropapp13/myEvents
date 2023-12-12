{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 3 Detector Variations

	WhichSampleArray.push_back("Overlay9_Run3_CV");
	WhichSampleArray.push_back("Overlay9_Run3_CVextra");
        WhichSampleArray.push_back("Overlay9_Run3_LYDown");
	WhichSampleArray.push_back("Overlay9_Run3_LYRayleigh");
	WhichSampleArray.push_back("Overlay9_Run3_LYAttenuation");

	WhichSampleArray.push_back("Overlay9_Run3_X");
	WhichSampleArray.push_back("Overlay9_Run3_YZ");
	WhichSampleArray.push_back("Overlay9_Run3_ThetaYZ");
	WhichSampleArray.push_back("Overlay9_Run3_ThetaXZ");
	WhichSampleArray.push_back("Overlay9_Run3_Recombination2");
	WhichSampleArray.push_back("Overlay9_Run3_SCE");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx++");		

	gROOT->ProcessLine(".L reco_selection.cxx+");
	gROOT->ProcessLine(".L true_selection.cxx+");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("reco_selection(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("true_selection(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
