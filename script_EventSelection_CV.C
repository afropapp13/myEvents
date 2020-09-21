{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1

	WhichSampleArray.push_back("BeamOn9_Run1");
	WhichSampleArray.push_back("ExtBNB9_Run1");
	WhichSampleArray.push_back("OverlayDirt9_Run1");
	WhichSampleArray.push_back("Overlay9_Run1");
	
	// Run 3

	WhichSampleArray.push_back("BeamOn9_Run3");
	WhichSampleArray.push_back("ExtBNB9_Run3");
	WhichSampleArray.push_back("OverlayDirt9_Run3");
	WhichSampleArray.push_back("Overlay9_Run3");	

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");	

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
