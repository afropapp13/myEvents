{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// NuWro

	WhichSampleArray.push_back("Overlay9NuWro_Run1");
	WhichSampleArray.push_back("Overlay9NuWro_Run2");
	WhichSampleArray.push_back("Overlay9NuWro_Run3");
	WhichSampleArray.push_back("Overlay9NuWro_Combined");

	//WhichSampleArray.push_back("GENIEv2Overlay9_Combined");	

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L PeLEE_myCCQERecoAnalysis.C++");
	gROOT->ProcessLine(".L PeLEE_myCCQETrueAnalysis.C++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("PeLEE_myCCQERecoAnalysis(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("PeLEE_myCCQETrueAnalysis(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
