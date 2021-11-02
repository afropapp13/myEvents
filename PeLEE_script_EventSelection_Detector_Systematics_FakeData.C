{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1 Detector Variations
	
//	WhichSampleArray.push_back("Overlay9_Run1_CV");
//        WhichSampleArray.push_back("Overlay9_Run1_LYDown");
//	WhichSampleArray.push_back("Overlay9_Run1_LYRayleigh");
//	WhichSampleArray.push_back("Overlay9_Run1_LYAttenuation");

//	WhichSampleArray.push_back("Overlay9_Run1_WireModX");
//	WhichSampleArray.push_back("Overlay9_Run1_WireModYZ");
//	WhichSampleArray.push_back("Overlay9_Run1_WireModThetaYZ");
//	WhichSampleArray.push_back("Overlay9_Run1_WireModThetaXZ");
//	WhichSampleArray.push_back("Overlay9_Run1_dEdx");
//	WhichSampleArray.push_back("Overlay9_Run1_Recombination2");
//	WhichSampleArray.push_back("Overlay9_Run1_SCE");

//	WhichSampleArray.push_back("Overlay9_Run1_lowE");
//	WhichSampleArray.push_back("Overlay9_Run1_altBE");

	// Run 3 Detector Variations

	WhichSampleArray.push_back("Overlay9_Run3_CV");
	WhichSampleArray.push_back("Overlay9_Run3_CVextra");
    WhichSampleArray.push_back("Overlay9_Run3_LYDown");
	WhichSampleArray.push_back("Overlay9_Run3_LYRayleigh");
	WhichSampleArray.push_back("Overlay9_Run3_LYAttenuation");

	WhichSampleArray.push_back("Overlay9_Run3_X");
	WhichSampleArray.push_back("Overlay9_Run3_YZ");
	WhichSampleArray.push_back("Overlay9_Run3_ThetaYZ");
	WhichSampleArray.push_back("Overlay9_Run3_ThetaXZ");
	WhichSampleArray.push_back("Overlay9_Run3_dEdx");
	WhichSampleArray.push_back("Overlay9_Run3_Recombination2");
	WhichSampleArray.push_back("Overlay9_Run3_SCE");	

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");		

	gROOT->ProcessLine(".L PeLEE_myRecoAnalysis.C+");
	gROOT->ProcessLine(".L PeLEE_myTrueAnalysis.C+");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("PeLEE_myRecoAnalysis(\""+WhichSampleArray[i]+"\",\"NoTune\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("PeLEE_myTrueAnalysis(\""+WhichSampleArray[i]+"\",\"NoTune\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
