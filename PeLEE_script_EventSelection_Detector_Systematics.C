{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1 Detector Variations
	
	//WhichSampleArray.push_back("Overlay9_Run1_CV");
        //WhichSampleArray.push_back("Overlay9_Run1_LYDown");
	//WhichSampleArray.push_back("Overlay9_Run1_LYRayleigh");
	//WhichSampleArray.push_back("Overlay9_Run1_LYAttenuation");

	//WhichSampleArray.push_back("Overlay9_Run1_X");
	//WhichSampleArray.push_back("Overlay9_Run1_YZ");
	//WhichSampleArray.push_back("Overlay9_Run1_ThetaYZ");
	//WhichSampleArray.push_back("Overlay9_Run1_ThetaXZ");
//	WhichSampleArray.push_back("Overlay9_Run1_dEdx"); // Not be used, covered by the recomb2 variation
	//WhichSampleArray.push_back("Overlay9_Run1_Recombination2");
	//WhichSampleArray.push_back("Overlay9_Run1_SCE");

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
//	WhichSampleArray.push_back("Overlay9_Run3_dEdx"); // Not be used, covered by the recomb2 variation
	WhichSampleArray.push_back("Overlay9_Run3_Recombination2");
	WhichSampleArray.push_back("Overlay9_Run3_SCE");

//	WhichSampleArray.push_back("Overlay9_Run3_lowE");

//	WhichSampleArray.push_back("Overlay9NuWro_Run1");		

	// -----------------------------------------------------------------------------------------

//	On the gpvm's
	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");		

	gROOT->ProcessLine(".L PeLEE_myCCQERecoAnalysis.C+");
	gROOT->ProcessLine(".L PeLEE_myCCQETrueAnalysis.C+");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("PeLEE_myCCQERecoAnalysis(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("PeLEE_myCCQETrueAnalysis(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
