{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("Overlay9_Combined");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");

	gROOT->ProcessLine(".L PeLEE_myRecoAnalysis.C+");
	gROOT->ProcessLine(".L PeLEE_myTrueAnalysis.C+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		// GENIE w/o T2K tune
	
		gROOT->ProcessLine("PeLEE_myRecoAnalysis(\""+WhichSampleArray[i]+"\",\"GENIEv2\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) { 

		  	gROOT->ProcessLine("PeLEE_myTrueAnalysis(\""+WhichSampleArray[i]+"\",\"GENIEv2\").Loop()"); 

		}		   

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
