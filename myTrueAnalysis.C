#define myTrueAnalysis_cxx
#include "myTrueAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TProfile.h>

#include <iostream>
#include <iomanip>
#include <vector>

#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Tools.h"
#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/STV_Tools.h"

using namespace std;

TString TrueToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}


void myTrueAnalysis::Loop() {

	if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
	TH1D::SetDefaultSumw2();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Determine if CCQELike or STV analysis

	TString Extention = "";
	if (CCQElike) { Extention = "CCQE_"; }

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	TString Extension = "";

	// For overlays only for genie, flux and reinteraction uncertainties

	if (fUniverseIndex != -1) {

		Extension = "_"+fEventWeightLabel+"_"+TrueToStringInt(fUniverseIndex); 

	}

	// Output Files

	TString FileName = PathToFiles+Extention+"TruthSTVAnalysis_"+fWhichSample+Extension+"_"+UBCodeVersion+".root";	
	TFile* OutputFile = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------

	// Txt file to keep track of the event reduction at each stage

	TString TxtName = "/uboone/data/users/apapadop/myEvents/myTxtFiles/"+UBCodeVersion+"/"+Extention+"TxtmyTrueEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);

	// Txt file to keep track of the run/subrun/event of the candidate events

	TString RunTxtName = "/uboone/data/users/apapadop/myEvents/myTxtFiles/"+UBCodeVersion+"/"+Extention+"TxtmyTrueRunSubRunEvents_"+fWhichSample+"_"+UBCodeVersion+".txt";
	ofstream myRunTxtFile;
	myRunTxtFile.open(RunTxtName);
	myRunTxtFile << std::fixed << std::setprecision(2);
	myRunTxtFile << fWhichSample;

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// 1D True Variables

	TH1D* TruePi0Plot = new TH1D("TruePi0Plot",";# #pi^{0}",4,-0.5,3.5);
	TH1D* TrueNeutronPlot = new TH1D("TrueNeutronPlot",";# neutrons",6,-0.5,5.5);
	
	TH1D* TruekMissPlot = new TH1D("TruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
	TH1D* TruePMissMinusPlot = new TH1D("TruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
	TH1D* TruePMissPlot = new TH1D("TruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);	

	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);

	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);

	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueEQEPlot = new TH1D("TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	TH1D* TrueVertexXPlot = new TH1D("TrueVertexXPlot",RecoLabelXAxisVertexX,NBinsVertexX,MinVertexX,MaxVertexX);
	TH1D* TrueVertexYPlot = new TH1D("TrueVertexYPlot",RecoLabelXAxisVertexY,NBinsVertexY,MinVertexY,MaxVertexY);
	TH1D* TrueVertexZPlot = new TH1D("TrueVertexZPlot",RecoLabelXAxisVertexZ,NBinsVertexZ,MinVertexZ,MaxVertexZ);

	TH1D* TrueEvPlot = new TH1D("TrueEvPlot",RecoLabelXAxisEv,NBinsEv,MinEv,MaxEv);
	TH1D* TrueNuPlot = new TH1D("TrueNuPlot",RecoLabelXAxisNu,NBinsNu,MinNu,MaxNu);
	
	// 2D Analysis

	// For now and until box opening
	int NBins2DAnalysis = 4;
		
	TH2D* TrueCosThetaMuPmuPlot = new TH2D("TrueCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
	TH2D* TrueCosThetaPPpPlot = new TH2D("TrueCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);


	// --------------------------------------------------------------------------------------------------------------------------------

	TH1D* TrueMuonTrueMomentumLongitudinalRatio = new TH1D("TrueMuonTrueMomentumLongitudinalRatio",";#frac{P^{true}_{#mu,||}}{P^{true}_{#mu}}",25,0.,1.);		
	TH1D* TrueProtonTrueMomentumLongitudinalRatio = new TH1D("TrueProtonTrueMomentumLongitudinalRatio",";#frac{P^{true}_{p,||}}{P^{true}_{p}}",25,0.,1.);

	TH1D* TrueMuonTrueMomentumTransverseRatio = new TH1D("TrueMuonTrueMomentumTransverseRatio",";#frac{P^{true}_{#mu,T}}{P^{true}_{#mu}}",25,0.,1.);		
	TH1D* TrueProtonTrueMomentumTransverseRatio = new TH1D("TrueProtonTrueMomentumTransverseRatio",";#frac{P^{true}_{p,T}}{P^{true}_{p}}",25,0.,1.);	

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	Tools tools;

	// -------------------------------------------------------------------------------------------------------------------------------------------

	int TrueCC1pCounter = 0;
	int TrueCCQElikeCounter = 0;
	
	double SumWeights = 0.;
	
	// --------------------------------------------------------------------------------------------------------------------------------

	TH1D* POTScalePlot = new TH1D("POTScalePlot","",1,0,1);
	TH1D* NEventsPlot = new TH1D("NEventsPlot","",1,0,1);
	TH1D* NSelectedPlot = new TH1D("NSelectedPlot","",1,0,1);	
	
	// --------------------------------------------------------------------------------------------------------------------------------

	// POT Counting

	double POTCount = -99.;

	if (string(fWhichSample).find("Overlay") != std::string::npos) {

		TString PathToPOTFile = "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root";
		TFile* POTFile = TFile::Open(PathToPOTFile,"readonly");
		TH1D* POTCountHist = (TH1D*)(POTFile->Get("POTCountHist"));
		POTCount = POTCountHist->GetBinContent(1);
		POTFile->Close();

	}	
	
	// ------------------------------------------------------------------------------------------------------------------

	// POT Scaling

	double POTScalingFactor = 1.;

	double tor860_wcut = 1;
	double E1DCNT_wcut = 1.;
	double EXT = 1.;

	if (string(fWhichSample).find("Run1") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run1;
		E1DCNT_wcut = E1DCNT_wcut_Run1;
		EXT = EXT_Run1;

	}
	
	if (string(fWhichSample).find("Run2") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run2;
		E1DCNT_wcut = E1DCNT_wcut_Run2;
		EXT = EXT_Run2;

	}	
	
	if (string(fWhichSample).find("Run3") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run3;
		E1DCNT_wcut = E1DCNT_wcut_Run3;
		EXT = EXT_Run3;

	}	
	
	if (string(fWhichSample).find("Run4") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run4;
		E1DCNT_wcut = E1DCNT_wcut_Run4;
		EXT = EXT_Run4;

	}		

	if (string(fWhichSample).find("Run5") != std::string::npos) {

		tor860_wcut = tor860_wcut_Run5;
		E1DCNT_wcut = E1DCNT_wcut_Run5;
		EXT = EXT_Run5;

	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;	

		// ------------------------------------------------------------------------------------------------------------------------------

//		double weight = 1.;
//		double T2Kweight = 1.;
		if (Weight <= 0 || Weight > 30) { continue; }
		if (T2KWeight <= 0 || T2KWeight > 30) { continue; }		
		// Weight from v3.0.4 to v.3.0.6 * weight from application of T2K tune
		POTScalingFactor = tor860_wcut / POTCount;
		double weight = POTScalingFactor * Weight * T2KWeight * ROOTinoWeight;

		// -------------------------------------------------------------------------------------------------------------------------------------

		// Genie, flux & reinteraction weights for systematics

		if ( fUniverseIndex != -1 && (fWhichSample == "Overlay9_Run1" || fWhichSample == "Overlay9_Run2" || fWhichSample == "Overlay9_Run3" 
		|| fWhichSample == "Overlay9_Run4" || fWhichSample == "Overlay9_Run5") ) {

			// Genie weights

			if (fEventWeightLabel == "All_UBGenie") { weight = weight*All_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "AxFFCCQEshape_UBGenie") { weight = weight*AxFFCCQEshape_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "DecayAngMEC_UBGenie") { weight = weight*DecayAngMEC_UBGenie->at(fUniverseIndex) / T2KWeight; }
			if (fEventWeightLabel == "NormCCCOH_UBGenie") { weight = weight*NormCCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "NormNCCOH_UBGenie") { weight = weight*NormNCCOH_UBGenie->at(fUniverseIndex)/ T2KWeight; }
//			if (fEventWeightLabel == "RPA_CCQE_Reduced_UBGenie") { weight = weight*RPA_CCQE_Reduced_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "RPA_CCQE_UBGenie") { weight = weight*RPA_CCQE_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "ThetaDelta2NRad_UBGenie") { weight = weight*ThetaDelta2NRad_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "Theta_Delta2Npi_UBGenie") { weight = weight*Theta_Delta2Npi_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "VecFFCCQEshape_UBGenie") { weight = weight*VecFFCCQEshape_UBGenie->at(fUniverseIndex)/ T2KWeight; }
			if (fEventWeightLabel == "XSecShape_CCMEC_UBGenie") { weight = weight*XSecShape_CCMEC_UBGenie->at(fUniverseIndex)/ T2KWeight; }

			// Flux weights

			if (fEventWeightLabel == "expskin_FluxUnisim") { weight = weight*expskin_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "horncurrent_FluxUnisim") { weight = weight*horncurrent_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "kminus_PrimaryHadronNormalization") 
				{ weight = weight*kminus_PrimaryHadronNormalization->at(fUniverseIndex); }
			if (fEventWeightLabel == "kplus_PrimaryHadronFeynmanScaling") 
				{ weight = weight*kplus_PrimaryHadronFeynmanScaling->at(fUniverseIndex); }
			if (fEventWeightLabel == "kzero_PrimaryHadronSanfordWang") 
				{ weight = weight*kzero_PrimaryHadronSanfordWang->at(fUniverseIndex); }
			if (fEventWeightLabel == "nucleoninexsec_FluxUnisim") { weight = weight*nucleoninexsec_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "nucleonqexsec_FluxUnisim") { weight = weight*nucleonqexsec_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "nucleontotxsec_FluxUnisim") { weight = weight*nucleontotxsec_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "piminus_PrimaryHadronSWCentralSplineVariation") 
				{ weight = weight*piminus_PrimaryHadronSWCentralSplineVariation->at(fUniverseIndex); }
			if (fEventWeightLabel == "pioninexsec_FluxUnisim") { weight = weight*pioninexsec_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "pionqexsec_FluxUnisim") { weight = weight*pionqexsec_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "piontotxsec_FluxUnisim") { weight = weight*piontotxsec_FluxUnisim->at(fUniverseIndex); }
			if (fEventWeightLabel == "piplus_PrimaryHadronSWCentralSplineVariation") 
				{ weight = weight*piplus_PrimaryHadronSWCentralSplineVariation->at(fUniverseIndex); }

			// Reinteraction weights

			if (fEventWeightLabel == "reinteractions_piminus_Geant4") 
				{ weight = weight*reinteractions_piminus_Geant4->at(fUniverseIndex); }
			if (fEventWeightLabel == "reinteractions_piplus_Geant4") 
				{ weight = weight*reinteractions_piplus_Geant4->at(fUniverseIndex); }
			if (fEventWeightLabel == "reinteractions_proton_Geant4") 
				{ weight = weight*reinteractions_proton_Geant4->at(fUniverseIndex); }
//			if (fEventWeightLabel == "xsr_scc_Fa3_SCC") { weight = weight*xsr_scc_Fa3_SCC->at(fUniverseIndex); }
//			if (fEventWeightLabel == "xsr_scc_Fv3_SCC") { weight = weight*xsr_scc_Fv3_SCC->at(fUniverseIndex); }		

		}

		// ---------------------------------------------------------------------------------------------------------------------------------

		SumWeights += weight / POTScalingFactor;					

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Analysis over the simb::MCParticles

		std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
		std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();

		int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0;
		bool TrueCC1pEvent = false;
		bool TrueCCQElikeEvent = false;

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & 0 pi+/- (Ppi > 70 MeV / c)

		if (CC1p == 1 && NumberPi0 == 0) {
		
			// --------------------------------------------------------------------------------------------------------------------		
			
			// Containment of the vertex has already been demanded in PreTruthSelection
			// Soft fiducial volume for true vertex

			//if (Muon_MCParticle_StartContainment->at(0) == 0) { continue; }

			// --------------------------------------------------------------------------------------------------------------------		

			// True muon

			double TrueMuonCosTheta = Muon_MCParticle_CosTheta->at(0);
			double TrueMuonTheta = TMath::ACos(TrueMuonCosTheta);
			double TrueMuonPhi_Deg = Muon_MCParticle_Phi->at(0);
			double TrueMuonPhi = TrueMuonPhi_Deg * TMath::Pi() / 180.;
			double TrueMuonMomentum_GeV = Muon_MCParticle_Mom->at(0); // GeV
//			double TrueMuonMomentum_MeV = 1000. * TrueMuonMomentum_GeV; // MeV
//			double TrueMuon_KE_MeV = tools.PToKE(MuonPdg,TrueMuonMomentum_MeV); // MeV
//			double TrueMuon_KE_GeV = TrueMuon_KE_MeV / 1000.; // GeV
//			double TrueMuon_E_GeV = TrueMuon_KE_GeV + MuonMass_GeV; // GeV
			double TrueMuon_E_GeV = TMath::Sqrt( TMath::Power(TrueMuonMomentum_GeV,2.) + TMath::Power(MuonMass_GeV,2.) ); // GeV
//			int TrueMuonStartContainment = Muon_MCParticle_StartContainment->at(0);
//			int TrueMuonEndContainment = Muon_MCParticle_EndContainment->at(0);

//			TVector3 TVector3TrueMuon;
//			TVector3TrueMuon.SetMagThetaPhi(TrueMuonMomentum_GeV,TrueMuonTheta,TrueMuonPhi);
//			TLorentzVector TrueMuon4V(TVector3TrueMuon,TrueMuon_E_GeV);

			// True proton

			double TrueProtonCosTheta = Proton_MCParticle_CosTheta->at(0);
			double TrueProtonTheta = TMath::ACos(TrueProtonCosTheta);
			double TrueProtonPhi_Deg = Proton_MCParticle_Phi->at(0);
			double TrueProtonPhi = TrueProtonPhi_Deg * TMath::Pi() / 180.;
			double TrueProtonMomentum_GeV = Proton_MCParticle_Mom->at(0); // GeV
//			double TrueProtonMomentum_MeV = 1000. * TrueProtonMomentum_GeV; // MeV
//			double TrueProton_KE_MeV = tools.PToKE(ProtonPdg,TrueProtonMomentum_MeV); // MeV
//			double TrueProton_KE_GeV = TrueProton_KE_MeV / 1000.; // GeV
//			double TrueProton_E_GeV = TrueProton_KE_GeV + ProtonMass_GeV; // GeV
			double TrueProton_E_GeV = TMath::Sqrt( TMath::Power(TrueProtonMomentum_GeV,2.) + TMath::Power(ProtonMass_GeV,2.) ); // GeV
//			int TrueProtonStartContainment = Proton_MCParticle_StartContainment->at(0);
//			int TrueProtonEndContainment = Proton_MCParticle_EndContainment->at(0);

//			TVector3 TVector3TrueProton;
//			TVector3TrueProton.SetMagThetaPhi(TrueProtonMomentum_GeV,TrueProtonTheta,TrueProtonPhi);
//			TLorentzVector TrueProton4V(TVector3TrueProton,TrueProton_E_GeV);

			double TrueDeltaPhiProtonMuon_Deg = True_DeltaPhi->at(0);
			double TrueDeltaThetaProtonMuon_Deg = True_DeltaTheta->at(0);

			// Reconstructed calorimetric energy using true level info / MCParticles

			double TrueRecoECal = True_ECal->at(0); // GeV

			// Reconstructed QE energy using true level info / MCParticles

			double TrueRecoEQE = True_EQE->at(0);

			// Definition of Transverse Variables

			double TruekMiss = True_kMiss->at(0);
			double TruePMissMinus = True_PMissMinus->at(0);
			double TrueMissMomentum = True_PMiss->at(0);

			double TrueTransMissMomentum = True_Pt->at(0);

			double TrueDeltaAlphaT = True_DeltaAlphaT->at(0);
			double TrueDeltaPhiT = True_DeltaPhiT->at(0);

			// True Vertex

//			TVector3 TrueVertex(Muon_MCParticle_StartX->at(0),Muon_MCParticle_StartY->at(0),Muon_MCParticle_StartZ->at(0));
			TVector3 TrueVertex(True_Vx,True_Vy,True_Vz);
			
			// -----------------------------------------------------------------------------------------------------------------
			
			// CCQElike analysis
			// Use the DeltaTheta, DeltaPhi & Pt cuts
				
			if (CCQElike) {
			
				if ( !(TMath::Abs(TrueDeltaPhiProtonMuon_Deg - 180.) < 35.) ) { continue; }
				if ( !(TMath::Abs(TrueDeltaThetaProtonMuon_Deg - 90.) < 55.) ) { continue; }	
				if ( NumberPi0 != 0 ) { continue; }					
				if ( !(TrueTransMissMomentum < 0.35) ) { continue; }							
			
			}	
				
			// -----------------------------------------------------------------------------------------------------------------	

			// Reconstructed q vector using true level info / MCParticles

			double RecoTrueQ2 = True_Q2->at(0);

			// Demand that the true muon / proton start points and the true proton end point are contained
			// Angle selection cuts: collinearity to reject broken tracks
			// Momentum threshold
			// Same events fill all the plots

			if (
			    /*TrueMuonStartContainment == true 
			    && TrueProtonStartContainment == true 
			    &&*/ TrueMuonMomentum_GeV > ArrayNBinsMuonMomentum[0]
			    && TrueProtonMomentum_GeV > ArrayNBinsProtonMomentum[0]
			) {

				// -----------------------------------------------------------------------------------------------------------------

				// STV analysis

				if (
				    TrueTransMissMomentum > ArrayNBinsDeltaPT[0] 
				    && TrueTransMissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
				    && TrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] 
				    && TrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
				    && TrueDeltaPhiT > ArrayNBinsDeltaPhiT[0] 
				    && TrueDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
				    
				    && TrueMuonMomentum_GeV < ArrayNBinsMuonMomentum[NBinsMuonMomentum]
				    && TrueProtonMomentum_GeV < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
				    
				    && TrueMuonCosTheta > ArrayNBinsMuonCosTheta[0]
				    && TrueMuonCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
				    && TrueProtonCosTheta > ArrayNBinsProtonCosTheta[0]
				    && TrueProtonCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
				    
//				    && TrueRecoECal > ArrayNBinsECal[0]
//				    && TrueRecoEQE > ArrayNBinsEQE[0]
//				    && RecoTrueQ2 > ArrayNBinsQ2[0]
//				    && TrueRecoECal < ArrayNBinsECal[NBinsECal]
//				    && TrueRecoEQE < ArrayNBinsEQE[NBinsEQE]
//				    && RecoTrueQ2 < ArrayNBinsQ2[NBinsQ2]
				) {

					// No weight to be applied in the multiplicity plots

					TruePi0Plot->Fill(NumberPi0);
					TrueNeutronPlot->Fill(NumberNeutrons);

					// True CC1p event

					TrueCC1pEvent = true;
					TrueCC1pCounter++;

					// STV

					TruekMissPlot->Fill(TruekMiss,weight);
					TruePMissMinusPlot->Fill(TruePMissMinus,weight);
					TruePMissPlot->Fill(TrueMissMomentum,weight);

					TrueDeltaPTPlot->Fill(TrueTransMissMomentum,weight);
					TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
					TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);

					// Reconstructed energy & Q2

					TrueECalPlot->Fill(TrueRecoECal,weight);
					TrueEQEPlot->Fill(TrueRecoEQE,weight);
					TrueQ2Plot->Fill(RecoTrueQ2,weight);

					// Kinematic variables

					TrueMuonMomentumPlot->Fill(TrueMuonMomentum_GeV,weight);
					TrueMuonPhiPlot->Fill(TrueMuonPhi_Deg,weight);
					TrueMuonCosThetaPlot->Fill(TrueMuonCosTheta,weight);

					TrueProtonMomentumPlot->Fill(TrueProtonMomentum_GeV,weight);
					TrueProtonPhiPlot->Fill(TrueProtonPhi_Deg,weight);
					TrueProtonCosThetaPlot->Fill(TrueProtonCosTheta,weight);

					// True Vertex

					TrueVertexXPlot->Fill(TrueVertex.X(),weight);
					TrueVertexYPlot->Fill(TrueVertex.Y(),weight);
					TrueVertexZPlot->Fill(TrueVertex.Z(),weight);

					// True Energy

					TrueEvPlot->Fill(True_Ev,weight);

					double true_MuonEnergy = TMath::Sqrt( TMath::Power(MuonMass_GeV,2.) + TMath::Power(TrueMuonMomentum_GeV,2.) );
					double true_Nu = True_Ev - true_MuonEnergy;
					TrueNuPlot->Fill(true_Nu,weight);
					
					// 2D Analysis
		
					TrueCosThetaMuPmuPlot->Fill(TrueMuonCosTheta,TrueMuonMomentum_GeV,weight);
					TrueCosThetaPPpPlot->Fill(TrueProtonCosTheta,TrueProtonMomentum_GeV,weight);

					// Playground for CC1p true momenta (longitudinal & perpendicular) ratios

					TVector3 TrueCandidateMuon(1,1,1);
					TrueCandidateMuon.SetMag(TrueMuonMomentum_GeV);
					TrueCandidateMuon.SetPhi(TrueMuonPhi);
					TrueCandidateMuon.SetTheta(TMath::ACos(TrueMuonCosTheta));

					TVector3 TrueCandidateProton(1,1,1);
					TrueCandidateProton.SetMag(TrueProtonMomentum_GeV);
					TrueCandidateProton.SetPhi(TrueProtonPhi);
					TrueCandidateProton.SetTheta(TMath::ACos(TrueProtonCosTheta));	

					TrueMuonTrueMomentumLongitudinalRatio->Fill(TrueCandidateMuon.Z() / TrueCandidateMuon.Mag(),weight);
					TrueMuonTrueMomentumTransverseRatio->Fill(TrueCandidateMuon.Pt() / TrueCandidateMuon.Mag(),weight);

					TrueProtonTrueMomentumLongitudinalRatio->Fill(TrueCandidateProton.Z() / TrueCandidateProton.Mag(),weight);
					TrueProtonTrueMomentumTransverseRatio->Fill(TrueCandidateProton.Pt() / TrueCandidateProton.Mag(),weight);				
					
				}

				// -----------------------------------------------------------------------------------------------------------------

			} // End of the angle selection cuts && the demand that we fill the plots with the same events

		} // End of signal definition: 1 mu (Pmu > 100 MeV / c), 1p (Pp > 200 MeV / c) & pi+/- (Ppi > 70 MeV / c)

		// -------------------------------------------------------------------------------------------------------------------------

		// Storing the run/subrun/event of the candidate events 

		myRunTxtFile << "Run = " << Run << ", SubRun = " << SubRun << ", Event = " << Event ;

	} // End of the loop over the events

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// STV-CC1p analysis summary

	std::cout << std::endl;

	// -------------------------------------------------------------------------------------------------------------------------

	// Samdef events

//	TFile* fFile = new TFile(fPathToFile);
//	TH1D* SamdefEventsPlot = (TH1D*)(fFile->Get("SamdefEventPlot"));
//	double SamdefEvents = SamdefEventsPlot->GetBinContent(1);

	// -------------------------------------------------------------------------------------------------------------------------

	if ( string(fWhichSample).find("Overlay9") != std::string::npos ) {

		//std::cout << std::endl << "Samdef events = " << SamdefEvents << std::endl;
		std::cout << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

		//myTxtFile << "Samdef events = " << SamdefEvents << std::endl;
		myTxtFile << std::endl << "True CC1p events = " << TrueCC1pCounter << std::endl;

	}
	
	POTScalePlot->SetBinContent(1,POTScalingFactor);
	POTScalePlot->SetBinError(1,0);

	NEventsPlot->SetBinContent(1,nentries);
	NEventsPlot->SetBinError(1,TMath::Sqrt(nentries));

	NSelectedPlot->SetBinContent(1,TrueCC1pCounter);
	NSelectedPlot->SetBinError(1,sqrt(TrueCC1pCounter));	

	// --------------------------------------------------------------------------------------------------------------------------------------------

//	double ScaleDueToWeights = double(fChain->GetEntries()) / double(SumWeights);

	// --------------------------------------------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "File " << FileName << " has been created"<< std::endl << std::endl;
	OutputFile->cd();
	OutputFile->Write();
	OutputFile->Close();
	myTxtFile.close();

} // End of the program
