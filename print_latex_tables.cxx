#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>

#include <iomanip>
#include <iostream>
#include <vector>

//#include "../Secondary_Code/myFunctions.cpp"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

void PeLEE_PrintEvents(vector<TFile*> FileSample,vector<double> FilePOT,vector<double> FileSamdef,vector<double> FilePreSel, TString SelectionStage, bool PreSelection = true) {

	// 0: BeamOn
	// 1: MC Overlay
	// 2: ExtBNB
	// 3: Dirt

	TH1D* BeamOnPlot = (TH1D*)(FileSample[0]->Get(SelectionStage));
	TH1D* MCPlot = (TH1D*)(FileSample[1]->Get(SelectionStage)); 
	TH1D* ExtBNBPlot = (TH1D*)(FileSample[2]->Get(SelectionStage)); 
	TH1D* DirtPlot = (TH1D*)(FileSample[3]->Get(SelectionStage)); 

	// All the counter plots have 1 bin
	// Grab the bin entry to get the number of events that survive a given preselection / selection cut

	double BeamOnEvents = BeamOnPlot->GetBinContent(1);
	double MCEvents = MCPlot->GetBinContent(1);
	double ExtBNBEvents = ExtBNBPlot->GetBinContent(1);
	double DirtEvents = DirtPlot->GetBinContent(1);

	// Fraction of events with respect to initial number of events in samdef

	double BeamOnFracSamdef = BeamOnEvents / FileSamdef[0] * 100.;
	double MCFracSamdef = MCEvents / FileSamdef[1] * 100.;
	double ExtBNBFracSamdef = ExtBNBEvents / FileSamdef[2] * 100.;
	double DirtFracSamdef = DirtEvents / FileSamdef[3] * 100.;

	// Fraction of events with respect to number of events passing preselection

	double BeamOnFracPreSel = BeamOnEvents / FilePreSel[0] * 100.;
	double MCFracPreSel = MCEvents / FilePreSel[1] * 100.;
	double ExtBNBFracPreSel = ExtBNBEvents / FilePreSel[2] * 100.;
	double DirtFracPreSel = DirtEvents / FilePreSel[3] * 100.;

	if (PreSelection) {

		cout << BeamOnEvents << " [" << BeamOnEvents*FilePOT[0] << "] (" << BeamOnFracSamdef << "\\%)&";
		cout << ExtBNBEvents << " [" << ExtBNBEvents*FilePOT[2] << "] (" << ExtBNBFracSamdef << "\\%)&";
		cout << MCEvents << " [" << MCEvents*FilePOT[1] << "] (" << MCFracSamdef << "\\%)&";
		cout << DirtEvents << " [" << DirtEvents*FilePOT[3] << "] (" << DirtFracSamdef << "\\%)";

	} else {

		cout << BeamOnEvents << " [" << BeamOnEvents*FilePOT[0] << "] (" << BeamOnFracSamdef << "\\%/" << BeamOnFracPreSel << "\\%)&";
		cout << ExtBNBEvents << " [" << ExtBNBEvents*FilePOT[2] << "] (" << ExtBNBFracSamdef << "\\%/" << ExtBNBFracPreSel << "\\%)&";
		cout << MCEvents << " [" << MCEvents*FilePOT[1] << "] (" << MCFracSamdef << "\\%/" << MCFracPreSel << "\\%)&";
		cout << DirtEvents << " [" << DirtEvents*FilePOT[3] << "] (" << DirtFracSamdef << "\\%/" << DirtFracPreSel << "\\%)";

	}

	cout << endl << "\\tabularnewline \\hline"<< endl;

}

void print_latex_tables(TString BaseMC = "", bool PrintStats = false, bool PrintPurEff = false, bool PrintCosmicCont = false, bool PrintIntBreakDown = false) {

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::cout << std::fixed << std::setprecision(2);
	gStyle->SetOptStat(0);

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	VectorCuts.push_back("_PID_NuScore_CRT");

	int NCuts = (int)(VectorCuts.size());	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
	Runs.push_back("Run2");
	Runs.push_back("Run3");
	//Runs.push_back("Run4a");
	Runs.push_back("Run4b");
	Runs.push_back("Run4c");
	Runs.push_back("Run4d");
	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
//	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {												

		// -----------------------------------------------------------------------------------------------------------------------------------------
			
		//cout << "// -------------------------------------------------" << endl << endl;
		//cout << Runs[WhichRun] << endl << endl;

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];		

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles+"/"+Cuts+"/";

			TH1D::SetDefaultSumw2();

			// ---------------------------------------------------------------------------------------------------------------------

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;

			// 0: BeamOn
			// 1: Overlay
			// 2: ExtBNB
			// 3: Dirt
		
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");

			if (BaseMC == "") { NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }
			else if (BaseMC == "Overlay9NuWro") { NameOfSamples.push_back("STVStudies_Overlay9NuWro_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("NuWro MC"); }

			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");
			NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt");

			// ---------------------------------------------------------------------------------------------------------------------

			// Truth level plot for efficiency & purity 

			TFile* TruthCC1pFile = TFile::Open(PathToFiles+"TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root");
			TH1D* hTruthCC1pOverlay = (TH1D*)(TruthCC1pFile->Get("TrueMuonCosThetaSingleBinPlot"));

			// ---------------------------------------------------------------------------------------------------------------------

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();
			vector<double> FilePOT; FilePOT.clear();
			vector<double> FileSamdef; FileSamdef.clear();
			vector<double> FilePreSel; FilePreSel.clear();
			vector<double> FileFinal; FileFinal.clear();
			vector<double> FileFinalError; FileFinalError.clear();
			vector<double> FileCC1p; FileCC1p.clear();
			vector<double> FileCC1pError; FileCC1pError.clear();

			vector<TH1D*> Plots; Plots.resize(NSamples);
			vector<TH1D*> CC1pPlots; CC1pPlots.resize(NSamples);
			vector<TH1D*> NonCC1pPlots; NonCC1pPlots.resize(NSamples);
			vector<TH1D*> NonCC1pPlotsEvents; NonCC1pPlotsEvents.resize(NSamples);

			vector<TH1D*> CCQEPlots; CCQEPlots.resize(NSamples);
			vector<TH1D*> CCMECPlots; CCMECPlots.resize(NSamples);
			vector<TH1D*> CCRESPlots; CCRESPlots.resize(NSamples);
			vector<TH1D*> CCDISPlots; CCDISPlots.resize(NSamples);

			// -----------------------------------------------------------------------------------------------------------------------

			// Grab the reference plots

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// --------------------------------------------------------------------------------------------

				FileSample.push_back(TFile::Open(PathToFilesCut+NameOfSamples[WhichSample]));

				// --------------------------------------------------------------------------------------------
	
				Plots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("RecoMuonCosThetaSingleBinPlot"));
				CC1pPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CC1pRecoMuonCosThetaSingleBinPlot"));
				NonCC1pPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pRecoMuonCosThetaSingleBinPlot"));
				NonCC1pPlotsEvents[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pEventPlot"));

				CCQEPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCQERecoMuonCosThetaSingleBinPlot"));
				CCMECPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCMECRecoMuonCosThetaSingleBinPlot"));
				CCRESPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCRESRecoMuonCosThetaSingleBinPlot"));
				CCDISPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCDISRecoMuonCosThetaSingleBinPlot"));

				// --------------------------------------------------------------------------------------------

				TH1D* FinalPlot = (TH1D*)(FileSample[WhichSample]->Get("RecoMuonCosThetaSingleBinPlot"));
				double FinalCount = FinalPlot->GetBinContent(1);
				FileFinal.push_back(FinalCount);

				double FinalCountError = FinalPlot->GetBinError(1);
				FileFinalError.push_back(FinalCountError);

				TH1D* CC1pPlot = (TH1D*)(FileSample[WhichSample]->Get("CC1pRecoMuonCosThetaSingleBinPlot"));
				double CC1pCount = CC1pPlot->GetBinContent(1);
				FileCC1p.push_back(CC1pCount);

				double CC1pCountError = CC1pPlot->GetBinError(1);
				FileCC1pError.push_back(CC1pCountError);

				//-------------------//

			} // End of the loop over the samples

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			if (PrintStats) {

				// Print summary of final number of selected events for all samples & POT equivalent

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {		

					cout << LabelsOfSamples[WhichSample] << " & " << FileFinal[WhichSample] << " $\\pm$ " << FileFinalError[WhichSample];
					cout << " \\tabularnewline \\hline" << endl;

					// Special case for MC, print also the CC1p event count

					if (string(LabelsOfSamples[WhichSample]).find("MC") != std::string::npos) {

						cout << "\\Signal " << LabelsOfSamples[WhichSample] << " & " << FileCC1p[WhichSample] << " $\\pm$ " << FileCC1pError[WhichSample];
						cout << " \\tabularnewline \\hline" << endl;

					}

				}

				cout << endl << endl;

			}

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Purity & Efficiency study

			// -----------------------------------------------------------------------------------------------------------------------

			// Efficiency

			// 1: Overlay

			double CC1p = CC1pPlots[1]->GetBinContent(1);
			double CC1pError = CC1pPlots[1]->GetBinError(1);

			double TrueCC1p = hTruthCC1pOverlay->GetBinContent(1);
			double TrueCC1pError = hTruthCC1pOverlay->GetBinError(1);

			double Efficiency = CC1p / TrueCC1p * 100.;
			double EfficiencyError = Efficiency * TMath::Sqrt( TMath::Power(CC1pError/CC1p,2.) + TMath::Power(TrueCC1pError/TrueCC1p,2.) );

			// -----------------------------------------------------------------------------------------------------------------------
			// Purity

			TH1D* SumNonBeamOn = nullptr;

			double SumErrors = 0.;

			// Skip BeamOn, thus start from index = 1

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample ++) {		

				if (WhichSample == 1) { SumNonBeamOn = (TH1D*)(Plots[WhichSample]->Clone()); }
				else { SumNonBeamOn->Add(Plots[WhichSample]); }

				SumErrors += TMath::Power(Plots[WhichSample]->GetBinError(1),2.);

			}

			SumErrors = TMath::Sqrt(SumErrors);
			double Sum = SumNonBeamOn->Integral();

			double Purity = CC1p / Sum * 100.;
			double PurityError = Purity * TMath::Sqrt( TMath::Power(CC1pError/CC1p,2.) + TMath::Power(SumErrors/Sum,2.) );

			// -----------------------------------------------------------------------------------------------------------------------
	
			if (PrintPurEff) {

				cout << Runs[WhichRun] << " & " << Purity << " $\\pm$ " << PurityError << " & " << Efficiency << " $\\pm$ " << EfficiencyError  << " \\tabularnewline \\hline" << endl;
		
			}

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Cosmic & dirt contamination

			double Cosmics = Plots[2]->GetBinContent(1);
			double CosmicsError = Plots[2]->GetBinError(1) ;
			double CosmicsFrac = Cosmics / Sum * 100.;
			double CosmicsFracError = CosmicsFrac * TMath::Sqrt( TMath::Power(CosmicsError/Cosmics,2.) + TMath::Power(SumErrors/Sum,2.) );

			double Dirt = Plots[3]->GetBinContent(1);
			double DirtError = Plots[3]->GetBinError(1) ;
			double DirtFrac = Dirt / Sum * 100.;
			double DirtFracError = DirtFrac * TMath::Sqrt( TMath::Power(DirtError/Dirt,2.) + TMath::Power(SumErrors/Sum,2.) );

			if (PrintCosmicCont) {

				cout << Runs[WhichRun] << " & " << CosmicsFrac << " $\\pm$ " << CosmicsFracError;
				cout << " & " << DirtFrac << " $\\pm$ " << DirtFracError << " \\tabularnewline \\hline" << endl;
				//cout << endl << endl;

			}

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Interaction break down

			double QE = CCQEPlots[1]->Integral() / Plots[1]->Integral() * 100.;
			double MEC = CCMECPlots[1]->Integral() / Plots[1]->Integral() * 100.;
			double RES = CCRESPlots[1]->Integral() / Plots[1]->Integral() * 100.;
			double DIS = CCDISPlots[1]->Integral() / Plots[1]->Integral() * 100.;

			double sum = QE + MEC + RES + DIS;
			if (sum > 105 || sum < 95) { cout << "UNITARITY CHECKED FAILED! sum = " << sum << endl; }

			if (PrintIntBreakDown) {

	 			cout << Runs[WhichRun] << " & " << QE << " & " << MEC << " & " << RES << " & " << DIS << " \\tabularnewline \\hline" << endl;
				//cout << endl << endl;

			}

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

		if (PrintStats) { cout << endl << endl; }
	
	} // End of the loop over the runs

} // End of the program 
