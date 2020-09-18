{

	vector<TString> WhichSampleArray;
	vector<TString> EventWeightLabels;
	vector<int> Universes;
//	int NFluxUniverses = 1000;
	int NFluxUniverses = 100;

	// -----------------------------------------------------------------------------------------

	// Run 1 Systematics

	WhichSampleArray.push_back("Overlay9_Run1");

	EventWeightLabels.push_back("horncurrent_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("kminus_PrimaryHadronNormalization"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("kplus_PrimaryHadronFeynmanScaling"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("kzero_PrimaryHadronSanfordWang"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("nucleoninexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("nucleonqexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("nucleontotxsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("piminus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("pioninexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("pionqexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("piontotxsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("piplus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(NFluxUniverses);
	EventWeightLabels.push_back("expskin_FluxUnisim"); Universes.push_back(10); // Watch out, different models

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/STV_Tools.cxx++");	

	gROOT->ProcessLine(".L myRecoAnalysis.C+");
	gROOT->ProcessLine(".L myTrueAnalysis.C+");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		for (int j = 0; j < (int)(EventWeightLabels.size()); j++) {

			for (int k = 0; k < Universes[j]; k++) {	

				gROOT->ProcessLine("myRecoAnalysis(\""+WhichSampleArray[i]+"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()");

				if (string(WhichSampleArray[i]).find("Overlay9") != std::string::npos) 
				  { gROOT->ProcessLine("myTrueAnalysis(\""+WhichSampleArray[i]+"\",\""+EventWeightLabels[j]+"\","+TString(std::to_string(k))+").Loop()"); } 

			} // End of the loop over the universes
			  
		} // End of the loop over the reinteraction labels	  

	} // End of the loop over the samples

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
