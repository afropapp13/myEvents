{

	vector<TString> WhichSampleArray;
	WhichSampleArray.push_back("Run1Data9");
	WhichSampleArray.push_back("ExtBNB9");
	WhichSampleArray.push_back("OverlayDirt9");

	WhichSampleArray.push_back("Overlay9");
//	WhichSampleArray.push_back("Overlay9_SCE");
//	WhichSampleArray.push_back("Overlay9_DLdown");

	gROOT->ProcessLine(".L t.C+");
	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("t(WhichSampleArray[i]).Loop()");

	}

	gROOT->ProcessLine(".q");
};
