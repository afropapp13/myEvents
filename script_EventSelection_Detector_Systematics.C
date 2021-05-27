{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1 Detector Variations
	
	WhichSampleArray.push_back("Overlay9_Run1_CV");
        WhichSampleArray.push_back("Overlay9_Run1_LYDown");
	WhichSampleArray.push_back("Overlay9_Run1_LYRayleigh");
	WhichSampleArray.push_back("Overlay9_Run1_LYAttenuation");

	WhichSampleArray.push_back("Overlay9_Run1_WireModX");
	WhichSampleArray.push_back("Overlay9_Run1_WireModYZ");
	WhichSampleArray.push_back("Overlay9_Run1_WireModThetaYZ");
	WhichSampleArray.push_back("Overlay9_Run1_WireModThetaXZ");
	WhichSampleArray.push_back("Overlay9_Run1_dEdx");
	WhichSampleArray.push_back("Overlay9_Run1_Recombination2");
	WhichSampleArray.push_back("Overlay9_Run1_SCE");

//	WhichSampleArray.push_back("Overlay9_Run1_lowE");
//	WhichSampleArray.push_back("Overlay9_Run1_altBE");

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

//	WhichSampleArray.push_back("Overlay9_Run3_lowE");	
	

	// -----------------------------------------------------------------------------------------

//	On the gpvm's'
	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/STV_Tools.cxx++");		

	gROOT->ProcessLine(".L myRecoAnalysis.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("myRecoAnalysis(\""+WhichSampleArray[i]+"\").Loop()");

		if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
		  { gROOT->ProcessLine("myTrueAnalysis(\""+WhichSampleArray[i]+"\").Loop()"); } 

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
