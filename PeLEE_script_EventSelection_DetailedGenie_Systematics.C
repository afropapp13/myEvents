{

	vector<TString> WhichSampleArray;
	vector<TString> EventWeightLabels;
	vector<int> Universes;

	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("Overlay9_Run1_DecompXSecUnc");

	// -----------------------------------------------------------------------------------------

	EventWeightLabels.push_back("AGKYpT1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("AGKYxF1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("AhtBY_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("BhtBY_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("CV1uBY_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("CV2uBY_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("EtaNCEL_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrAbs_N_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrAbs_pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrCEx_N_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrCEx_pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrInel_N_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrInel_pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrPiProd_N_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FrPiProd_pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FracDelta_CCMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("FracPN_CCMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MFP_N_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MFP_pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MaCCQE_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MaCCRES_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MaNCEL_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MaNCRES_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MvCCRES_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("MvNCRES_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarnCC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarnCC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarnNC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarnNC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarpCC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarpCC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarpNC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvbarpNC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvnCC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvnCC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvnNC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvnNC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvpCC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvpCC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvpNC1pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NonRESBGvpNC2pi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NormCCMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("NormNCMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("RDecBR1eta_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("RDecBR1gamma_UBGenie"); Universes.push_back(2);	

	// Unisims										

	EventWeightLabels.push_back("UnShortAxFFCCQEshape_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("UnShortDecayAngMEC_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("UnShortRPA_CCQE_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("UnShortTheta_Delta2Npi_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("UnShortVecFFCCQEshape_UBGenie"); Universes.push_back(2);
	EventWeightLabels.push_back("UnShortXSecShape_CCMEC_UBGenie"); Universes.push_back(2);						

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L PeLEE_myRecoAnalysis.C+");
	gROOT->ProcessLine(".L PeLEE_myTrueAnalysis.C+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(EventWeightLabels.size()); j++) {

			for (int k = 0; k < Universes[j]; k++) {	

				gROOT->ProcessLine("PeLEE_myRecoAnalysis(\""+WhichSampleArray[i]+"\",\"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

				if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
				  { gROOT->ProcessLine("PeLEE_myTrueAnalysis(\""+WhichSampleArray[i]+"\",\"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()"); } 

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
