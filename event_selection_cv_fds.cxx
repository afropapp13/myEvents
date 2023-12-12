{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1

	WhichSampleArray.push_back("OverlayDirt9_Run1");
	WhichSampleArray.push_back("Overlay9_Run1");

	// Run 2

	WhichSampleArray.push_back("OverlayDirt9_Run2");
	WhichSampleArray.push_back("Overlay9_Run2");
	
	// Run 3

	WhichSampleArray.push_back("OverlayDirt9_Run3");
	WhichSampleArray.push_back("Overlay9_Run3");

	// Run 4a

	WhichSampleArray.push_back("OverlayDirt9_Run4a");
	WhichSampleArray.push_back("Overlay9_Run4a");	

	// Run 4b

	WhichSampleArray.push_back("OverlayDirt9_Run4b");
	WhichSampleArray.push_back("Overlay9_Run4b");	

	// Run 4c

	WhichSampleArray.push_back("OverlayDirt9_Run4c");
	WhichSampleArray.push_back("Overlay9_Run4c");	

	// Run 4d

	WhichSampleArray.push_back("OverlayDirt9_Run4d");
	WhichSampleArray.push_back("Overlay9_Run4d");	

	// Run 5

	WhichSampleArray.push_back("OverlayDirt9_Run5");
	WhichSampleArray.push_back("Overlay9_Run5");

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
