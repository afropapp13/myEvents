{

	vector<TString> WhichSampleArray;

	WhichSampleArray.push_back("Overlay9");
	WhichSampleArray.push_back("Overlay9_SCE");
	WhichSampleArray.push_back("Overlay9_DLdown");

	gROOT->ProcessLine(".L ../../myClass/Tools.cxx+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");
	//gROOT->ProcessLine("myTrueAnalysis().Loop()");

	gROOT->ProcessLine(".L myTrueAnalysis.C+");
	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("myTrueAnalysis(WhichSampleArray[i]).Loop()");

	}

	gROOT->ProcessLine(".q");
};
