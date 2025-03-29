{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------
	
	WhichSampleArray.push_back("mcc9_10_Overlay9_Run1");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L mcc9_10_reco_selection.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("mcc9_10_reco_selection(\""+WhichSampleArray[i]+"\").Loop()");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
