{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Combined

	WhichSampleArray.push_back("OverlayDirt9_Combined");
	WhichSampleArray.push_back("Overlay9_Combined");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx++");

	gROOT->ProcessLine(".L reco_selection.cxx+");
	gROOT->ProcessLine(".L true_selection.cxx+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		// GENIE w/o T2K tune
	
		gROOT->ProcessLine("reco_selection(\""+WhichSampleArray[i]+"\",\"NoTune\").Loop()");
		gROOT->ProcessLine("reco_selection(\""+WhichSampleArray[i]+"\",\"TwiceMEC\").Loop()");		

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) { 

			gROOT->ProcessLine("true_selection(\""+WhichSampleArray[i]+"\",\"NoTune\").Loop()"); 
		  	gROOT->ProcessLine("true_selection(\""+WhichSampleArray[i]+"\",\"TwiceMEC\").Loop()"); 

		}		   

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
