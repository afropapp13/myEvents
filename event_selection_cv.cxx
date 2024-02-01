{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 3 // detector variations
	
        WhichSampleArray.push_back("BeamOn9_Run3");
        WhichSampleArray.push_back("ExtBNB9_Run3");
        WhichSampleArray.push_back("OverlayDirt9_Run3");
        WhichSampleArray.push_back("Overlay9_Run3");

	// Combined
	
        WhichSampleArray.push_back("BeamOn9_Combined");
        WhichSampleArray.push_back("ExtBNB9_Combined");
        WhichSampleArray.push_back("OverlayDirt9_Combined");
        WhichSampleArray.push_back("Overlay9_Combined");

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
