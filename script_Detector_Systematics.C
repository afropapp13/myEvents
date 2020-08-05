{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1 Detector Variations
	
	WhichSampleArray.push_back("Overlay9_Run1_CV");
        WhichSampleArray.push_back("Overlay9_Run1_LYDown");
	WhichSampleArray.push_back("Overlay9_Run1_LYRayleigh");
	WhichSampleArray.push_back("Overlay9_Run1_LYAttenuation");

	// Run 3 Detector Variations

	WhichSampleArray.push_back("Overlay9_Run3_CV");
        WhichSampleArray.push_back("Overlay9_Run3_LYDown");
	WhichSampleArray.push_back("Overlay9_Run3_LYRayleigh");
	WhichSampleArray.push_back("Overlay9_Run3_LYAttenuation");

	WhichSampleArray.push_back("Overlay9_Run3_WireModX");
	WhichSampleArray.push_back("Overlay9_Run3_WireModYZ");
	WhichSampleArray.push_back("Overlay9_Run3_WireModThetaYZ");
	WhichSampleArray.push_back("Overlay9_Run3_WireModThetaXZ");
	WhichSampleArray.push_back("Overlay9_Run3_dEdx");
	WhichSampleArray.push_back("Overlay9_Run3_Recombination2");
	WhichSampleArray.push_back("Overlay9_Run3_SCE");	
	

	// -----------------------------------------------------------------------------------------

//	Locally
//	gROOT->ProcessLine(".L ../../myClass/Tools.cxx++");
//	On the gpvm's'
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");		

	gROOT->ProcessLine(".L t.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("t(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("myTrueAnalysis(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
