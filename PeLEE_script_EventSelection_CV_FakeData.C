{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1

//	WhichSampleArray.push_back("OverlayDirt9_Run1");
//	WhichSampleArray.push_back("Overlay9_Run1");

	// Run 2

//	WhichSampleArray.push_back("OverlayDirt9_Run2");
//	WhichSampleArray.push_back("Overlay9_Run2");
	
	// Run 3

	WhichSampleArray.push_back("OverlayDirt9_Run3");
	WhichSampleArray.push_back("Overlay9_Run3");

	// Run 4

//	WhichSampleArray.push_back("OverlayDirt9_Run4");
//	WhichSampleArray.push_back("Overlay9_Run4");	

	// Run 5

//	WhichSampleArray.push_back("OverlayDirt9_Run5");
//	WhichSampleArray.push_back("Overlay9_Run5");

	// Combined

	WhichSampleArray.push_back("OverlayDirt9_Combined");
	WhichSampleArray.push_back("Overlay9_Combined");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");

	gROOT->ProcessLine(".L PeLEE_myRecoAnalysis.C+");
	gROOT->ProcessLine(".L PeLEE_myTrueAnalysis.C+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		// GENIE w/o T2K tune
	
		gROOT->ProcessLine("PeLEE_myRecoAnalysis(\""+WhichSampleArray[i]+"\",\"NoTune\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("PeLEE_myTrueAnalysis(\""+WhichSampleArray[i]+"\",\"NoTune\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
