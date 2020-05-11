{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1 Systematics

	WhichSampleArray.push_back("Overlay9_Run1_CV");
	WhichSampleArray.push_back("Overlay9_Run1_X");
	WhichSampleArray.push_back("Overlay9_Run1_YZ");
	WhichSampleArray.push_back("Overlay9_Run1_LY");
	WhichSampleArray.push_back("Overlay9_Run1_LYRayleigh");
//	WhichSampleArray.push_back("Overlay9_Run1_LowE");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClass/Tools.cxx+");
	gROOT->ProcessLine(".L t.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("t(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("myTrueAnalysis(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
