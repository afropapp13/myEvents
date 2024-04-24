{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Combined

	WhichSampleArray.push_back("BeamOn9_Combined");
	WhichSampleArray.push_back("ExtBNB9_Combined");
	WhichSampleArray.push_back("OverlayDirt9_Combined");
	WhichSampleArray.push_back("Overlay9_Combined");

	// BNB to Honda reweight
	//WhichSampleArray.push_back("Overlay9BNBToHonda_Combined");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L reco_selection.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("reco_selection(\""+WhichSampleArray[i]+"\").Loop()");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
