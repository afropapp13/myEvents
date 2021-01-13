{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1

	WhichSampleArray.push_back("BeamOn9_Run1");
	WhichSampleArray.push_back("ExtBNB9_Run1");
	WhichSampleArray.push_back("OverlayDirt9_Run1");
	WhichSampleArray.push_back("Overlay9_Run1");

	// Run 2

//	WhichSampleArray.push_back("BeamOn9_Run2");
//	WhichSampleArray.push_back("ExtBNB9_Run2");
//	WhichSampleArray.push_back("OverlayDirt9_Run2");
//	WhichSampleArray.push_back("Overlay9_Run2");
	
	// Run 3

	WhichSampleArray.push_back("BeamOn9_Run3");
	WhichSampleArray.push_back("ExtBNB9_Run3");
	WhichSampleArray.push_back("OverlayDirt9_Run3");
	WhichSampleArray.push_back("Overlay9_Run3");

	// Run 4

//	WhichSampleArray.push_back("BeamOn9_Run4");
//	WhichSampleArray.push_back("ExtBNB9_Run4");
//	WhichSampleArray.push_back("OverlayDirt9_Run4");
//	WhichSampleArray.push_back("Overlay9_Run4");	

	// Run 5

//	WhichSampleArray.push_back("BeamOn9_Run5");
//	WhichSampleArray.push_back("ExtBNB9_Run5");
//	WhichSampleArray.push_back("OverlayDirt9_Run5");
//	WhichSampleArray.push_back("Overlay9_Run5");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L myRecoAnalysis.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");


	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("myRecoAnalysis(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("myTrueAnalysis(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
