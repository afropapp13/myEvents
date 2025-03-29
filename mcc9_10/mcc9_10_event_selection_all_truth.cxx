{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

//	WhichSampleArray.push_back("mcc9_10_Overlay9_Run1");
	WhichSampleArray.push_back("mcc9_10_Overlay9_Run4b");
	WhichSampleArray.push_back("mcc9_10_Overlay9_Run4b_noweights");	

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L mcc9_10_true_selection.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("mcc9_10_true_selection(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}
	//--------------------//
	
	gROOT->ProcessLine(".q");

	//--------------------//
	
}
