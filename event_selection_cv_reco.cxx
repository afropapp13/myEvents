{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Combined
	//WhichSampleArray.push_back("BeamOn9_Run1");
	
	WhichSampleArray.push_back("BeamOn9_Combined");
	//WhichSampleArray.push_back("ExtBNB9_Combined");
	//WhichSampleArray.push_back("OverlayDirt9_Combined");
	WhichSampleArray.push_back("Overlay9_Combined");

	// BNB to Honda reweight using true Enu
	//WhichSampleArray.push_back("BeamOn9BNBToHonda_Combined");
	//WhichSampleArray.push_back("ExtBNB9BNBToHonda_Combined");
	//WhichSampleArray.push_back("Overlay9BNBToHonda_Combined");
	//WhichSampleArray.push_back("OverlayDirt9BNBToHonda_Combined");

	// BNB to Honda reweight using ECal
	//WhichSampleArray.push_back("BeamOn9BNBToHondaECal_Combined");
	//WhichSampleArray.push_back("ExtBNB9BNBToHondaECal_Combined");
	//WhichSampleArray.push_back("Overlay9BNBToHondaECal_Combined");
	//WhichSampleArray.push_back("OverlayDirt9BNBToHondaECal_Combined");
	
	//WhichSampleArray.push_back("BeamOn9_Run1A_open_trigger");
	//WhichSampleArray.push_back("ExtBNB9_Run1A_open_trigger");
	//WhichSampleArray.push_back("OverlayDirt9_Run1A_open_trigger");
	//WhichSampleArray.push_back("Overlay9_Run1A_open_trigger");

	//WhichSampleArray.push_back("BeamOn9_Run1B_open_trigger");
	//WhichSampleArray.push_back("ExtBNB9_Run1B_open_trigger");
	//WhichSampleArray.push_back("OverlayDirt9_Run1B_open_trigger");
	//WhichSampleArray.push_back("Overlay9_Run1B_open_trigger");


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
