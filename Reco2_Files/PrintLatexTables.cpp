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

#include "../Secondary_Code/myFunctions.cpp"
#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void PrintMCBkgEvents(TFile* FileSample,double FilePOT,double FileFinal,TString MCBkg, TString MCBkgLabel, double NonCC1pFinal) {

	TH1D* MCPlot = (TH1D*)(FileSample->Get(MCBkg));

	double MCEvents = MCPlot->GetBinContent(1);
	double MCEventsError = MCPlot->GetBinError(1);

	if (MCBkg == "NonCC1pEventPlot") { cout << "\\hline"<< endl; }

	cout << MCBkgLabel << " & " << MCEvents << " $\\pm$ " << MCEventsError << " & " << MCEvents*FilePOT << " $\\pm$ " << MCEventsError*FilePOT << " & ";
	cout << MCEvents / FileFinal * 100. << " & " << MCEvents / NonCC1pFinal * 100.;


	cout << " \\tabularnewline \\hline"<< endl;

}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------


void PrintEvents(vector<TFile*> FileSample,vector<double> FilePOT,vector<double> FileSamdef,vector<double> FilePreSel, TString SelectionStage, bool PreSelection = true) {

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

void PrintLatexTables(TString BaseMC = "") {

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::cout << std::fixed << std::setprecision(2);
	gStyle->SetOptStat(0);

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::vector<TString> PlotNames; PlotNames.clear();
	std::vector<TString> PlotLabels; PlotLabels.clear();

	PlotNames.push_back("SamdefEventPlot"); PlotLabels.push_back("Initial");
	PlotNames.push_back("SWTriggerEventPlot"); PlotLabels.push_back("SW Trigger");
	PlotNames.push_back("OneNuMuPFParticleEventPlot"); PlotLabels.push_back("1 $\\nu_{\\mu}$ PFParticle");
	PlotNames.push_back("MatchedTrackLikeDaughterEventPlot"); PlotLabels.push_back("2 track-like daughters");
	PlotNames.push_back("OneBeamFlashEventPlot"); PlotLabels.push_back("1 beam flash");
	PlotNames.push_back("MomentumThresholdEventPlot"); PlotLabels.push_back("P threshold");
	PlotNames.push_back("StartPointContainmentEventPlot"); PlotLabels.push_back("Start point in FV");
	PlotNames.push_back("ProtonEndPointContainmentEventPlot"); PlotLabels.push_back("End point p in FV");
	PlotNames.push_back("VertexContainmentEventPlot"); PlotLabels.push_back("Vertex in FV");
	PlotNames.push_back("MuonQualityEventPlot"); PlotLabels.push_back("Contained $\\mu$ quality cut");
	PlotNames.push_back("NoFlippedTrackEventPlot"); PlotLabels.push_back("No flipped tracks");
	PlotNames.push_back("SumTruncdEdxTrackEventPlot"); PlotLabels.push_back("Trunc dEdx sum");
	PlotNames.push_back("LowMCSQualityEventPlot"); PlotLabels.push_back("Exiting $\\mu$ length threshold");
	PlotNames.push_back("PidEventPlot"); PlotLabels.push_back("Calorimetry");
	PlotNames.push_back("NuScoreEventPlot"); PlotLabels.push_back("$\\nu$ score");
	PlotNames.push_back("CommonEventPlot"); PlotLabels.push_back("Common events");

	const int N1DPlots = PlotNames.size();
//	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::vector<TString> MCBkgPlotNames; MCBkgPlotNames.clear();
	std::vector<TString> MCBkgPlotLabels; MCBkgPlotLabels.clear();

	MCBkgPlotNames.push_back("CC2pEventPlot"); MCBkgPlotLabels.push_back("CC2p");
	MCBkgPlotNames.push_back("OutCommonRangeEventPlot"); MCBkgPlotLabels.push_back("Truth out of common range");
	MCBkgPlotNames.push_back("PiPEventPlot"); MCBkgPlotLabels.push_back("$\\pi$-p");
	MCBkgPlotNames.push_back("NeutralPiEventPlot"); MCBkgPlotLabels.push_back("$\\pi^{0}$ production");
	MCBkgPlotNames.push_back("CC1p1piEventPlot"); MCBkgPlotLabels.push_back("CC1p1$\\pi$");
	MCBkgPlotNames.push_back("CC3pEventPlot"); MCBkgPlotLabels.push_back("CC3p");
	MCBkgPlotNames.push_back("TrueVertexOutFVEventPlot"); MCBkgPlotLabels.push_back("True vertex outside FV");
	MCBkgPlotNames.push_back("PPEventPlot"); MCBkgPlotLabels.push_back("p-p");
	MCBkgPlotNames.push_back("BrokenMuEventPlot"); MCBkgPlotLabels.push_back("Broken $\\mu$ track");
	MCBkgPlotNames.push_back("CC4pEventPlot"); MCBkgPlotLabels.push_back("CC4p");
	MCBkgPlotNames.push_back("NCEventPlot"); MCBkgPlotLabels.push_back("NC");
	MCBkgPlotNames.push_back("CC2p1piEventPlot"); MCBkgPlotLabels.push_back("CC2p1$\\pi$");
	MCBkgPlotNames.push_back("InTimeCosmicsEventPlot"); MCBkgPlotLabels.push_back("In-Time Cosmics");
	MCBkgPlotNames.push_back("AntiMuPEventPlot"); MCBkgPlotLabels.push_back("$\\mu^{+}$-p");
	MCBkgPlotNames.push_back("MuPiEventPlot"); MCBkgPlotLabels.push_back("$\\mu$-$\\pi$");
	MCBkgPlotNames.push_back("MuEEventPlot"); MCBkgPlotLabels.push_back("$\\mu$-e");
	MCBkgPlotNames.push_back("BrokenPEventPlot"); MCBkgPlotLabels.push_back("Broken p track");
	MCBkgPlotNames.push_back("MultipleVerticesEventPlot"); MCBkgPlotLabels.push_back("Multiple vertices");
	MCBkgPlotNames.push_back("OtherMCBkgEventPlot"); MCBkgPlotLabels.push_back("Other");
	MCBkgPlotNames.push_back("NonCC1pEventPlot"); MCBkgPlotLabels.push_back("Total Non-\\Signal");

	const int MCBkgPlots = MCBkgPlotNames.size();
//	cout << "Number of MCBkg 1D Plots = " << MCBkgPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	//VectorCuts.push_back("");
	//VectorCuts.push_back("_PID");
	VectorCuts.push_back("_PID_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
//	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {												

		// -----------------------------------------------------------------------------------------------------------------------------------------
			
		cout << "// -------------------------------------------------" << endl << endl;
		cout << Runs[WhichRun] << endl << endl;

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
			TH1D* hTruthCC1pOverlay = (TH1D*)(TruthCC1pFile->Get("TrueMuonCosThetaPlot"));

			// ---------------------------------------------------------------------------------------------------------------------

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();
			vector<double> FilePOT; FilePOT.clear();
			vector<double> FileSamdef; FileSamdef.clear();
			vector<double> FilePreSel; FilePreSel.clear();
			vector<double> FileFinal; FileFinal.clear();
			vector<double> FileCC1p; FileCC1p.clear();

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

				// Numbers to keep track of

				TH1D* POTPlot = (TH1D*)(FileSample[WhichSample]->Get("POTScalePlot"));
				double POTCount = POTPlot->GetBinContent(1);
				FilePOT.push_back(POTCount);

				TH1D* SamdefPlot = (TH1D*)(FileSample[WhichSample]->Get("SamdefEventPlot"));
				double SamdefCount = SamdefPlot->GetBinContent(1);
				FileSamdef.push_back(SamdefCount);

				TH1D* PreSelPlot = (TH1D*)(FileSample[WhichSample]->Get("VertexContainmentEventPlot"));
				double PreSelCount = PreSelPlot->GetBinContent(1);
				FilePreSel.push_back(PreSelCount);

				TH1D* FinalPlot = (TH1D*)(FileSample[WhichSample]->Get("CommonEventPlot"));
				double FinalCount = FinalPlot->GetBinContent(1);
				FileFinal.push_back(FinalCount);

				TH1D* CC1pPlot = (TH1D*)(FileSample[WhichSample]->Get("NCC1pPlot"));
				double CC1pCount = CC1pPlot->GetBinContent(1);
				FileCC1p.push_back(CC1pCount);

				// --------------------------------------------------------------------------------------------
	
				Plots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("RecoMuonCosThetaPlot"));
				CC1pPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CC1pRecoMuonCosThetaPlot"));
				NonCC1pPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pRecoMuonCosThetaPlot"));
				NonCC1pPlotsEvents[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pEventPlot"));

				CCQEPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCQERecoMuonCosThetaPlot"));
				CCMECPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCMECRecoMuonCosThetaPlot"));
				CCRESPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCRESRecoMuonCosThetaPlot"));
				CCDISPlots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("CCDISRecoMuonCosThetaPlot"));

				// --------------------------------------------------------------------------------------------

			} // End of the loop over the samples

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Loop over the plots for event loss at each step

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				bool PreSelection = true;
				//if (PlotLabels[WhichPlot] == "Muon quality cut" || PlotLabels[WhichPlot] == "Calorimetry" 
				// || PlotLabels[WhichPlot] == "\\nu score" || PlotLabels[WhichPlot] == "Common events") { PreSelection = false; }


				cout << PlotLabels[WhichPlot] << " & ";
				PrintEvents(FileSample,FilePOT,FileSamdef,FilePreSel,PlotNames[WhichPlot],PreSelection);

				if (PlotLabels[WhichPlot] == "Initial" || PlotLabels[WhichPlot] == "Vertex in FV"  || PlotLabels[WhichPlot] == "Common ranges") 
					{ cout << "\\hline" << endl; }

				if (WhichPlot == N1DPlots-1) {

					cout << "Final & ";
					PrintEvents(FileSample,FilePOT,FileSamdef,FilePreSel,PlotNames[WhichPlot],PreSelection);

				}
		
			} // End of the loop over the plots for event loss at each step	

			cout << endl << endl;

			// --------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Print summary of final number of selected events for all samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {		

				cout << LabelsOfSamples[WhichSample] << " & " << FileFinal[WhichSample] << " $\\pm$ " << TMath::Sqrt(FileFinal[WhichSample]);
				cout << " & " << FileFinal[WhichSample] * FilePOT[WhichSample] << " $\\pm$ " << TMath::Sqrt(FileFinal[WhichSample] * FilePOT[WhichSample]);
				cout << " \\tabularnewline \\hline" << endl;

				// Special case for MC, print also the CC1p event count

				if (string(LabelsOfSamples[WhichSample]).find("MC") != std::string::npos) {

					cout << "\\Signal " << LabelsOfSamples[WhichSample] << " & " << FileCC1p[WhichSample] << " $\\pm$ " << TMath::Sqrt(FileCC1p[WhichSample]);
					cout << " & " << FileCC1p[WhichSample] * FilePOT[WhichSample] << " $\\pm$ " << TMath::Sqrt(FileCC1p[WhichSample] * FilePOT[WhichSample]);
					cout << " \\tabularnewline \\hline" << endl;

				}

			}

			cout << endl << endl;

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Purity & Efficiency study

			// -----------------------------------------------------------------------------------------------------------------------

			// Efficiency

			// 1: Overlay

			double CC1p = CC1pPlots[1]->Integral();
			double CC1pError = TMath::Sqrt(CC1pPlots[1]->Integral() * FilePOT[1]);

			double TrueCC1p = hTruthCC1pOverlay->Integral();
			double TrueCC1pError = TMath::Sqrt(TrueCC1p * FilePOT[1]);

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

				SumErrors += Plots[WhichSample]->Integral() * FilePOT[WhichSample];

			}

			SumErrors = TMath::Sqrt(SumErrors);
			double Sum = SumNonBeamOn->Integral();

			double Purity = CC1p / Sum * 100.;
			double PurityError = Purity * TMath::Sqrt( TMath::Power(CC1pError/CC1p,2.) + TMath::Power(SumErrors/Sum,2.) );

			// -----------------------------------------------------------------------------------------------------------------------
	
			cout << "\\Signal & " << Purity << " $\\pm$ " << PurityError << " & " << Efficiency << " $\\pm$ " << EfficiencyError;
			cout << " \\tabularnewline \\hline" << endl;
			cout << endl << endl;
		
			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Cosmic & dirt contamination

			double Cosmics = Plots[2]->Integral();
			double CosmicsError = TMath::Sqrt(Plots[2]->Integral() * FilePOT[2]) ;
			double CosmicsFrac = Cosmics / Sum * 100.;
			double CosmicsFracError = CosmicsFrac * TMath::Sqrt( TMath::Power(CosmicsError/Cosmics,2.) + TMath::Power(SumErrors/Sum,2.) );

			double Dirt = Plots[3]->Integral();
			double DirtError = TMath::Sqrt(Plots[3]->Integral() * FilePOT[3]) ;
			double DirtFrac = Dirt / Sum * 100.;
			double DirtFracError = DirtFrac * TMath::Sqrt( TMath::Power(DirtError/Dirt,2.) + TMath::Power(SumErrors/Sum,2.) );

			cout << "Cosmics & " << CosmicsFrac << " $\\pm$ " << CosmicsFracError << " \\tabularnewline \\hline" << endl;
			cout << "Dirt & " << DirtFrac << " $\\pm$ " << DirtFracError << " \\tabularnewline \\hline" << endl;
			cout << endl << endl;

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// Interaction break down

			double QE = CCQEPlots[1]->Integral() / Plots[1]->Integral() * 100.;
			double MEC = CCMECPlots[1]->Integral() / Plots[1]->Integral() * 100.;
			double RES = CCRESPlots[1]->Integral() / Plots[1]->Integral() * 100.;
			double DIS = CCDISPlots[1]->Integral() / Plots[1]->Integral() * 100.;

			double sum = QE + MEC + RES + DIS;
			if (sum > 105 || sum < 95) { cout << "UNITARITY CHECKED FAILED! sum = " << sum << endl; }

 			cout << "QE & " << QE << " \\tabularnewline \\hline" << endl;
 			cout << "MEC & " << MEC << " \\tabularnewline \\hline" << endl;
 			cout << "RES & " << RES << " \\tabularnewline \\hline" << endl;
 			cout << "DIS & " << DIS << " \\tabularnewline \\hline" << endl;

			cout << endl << endl;

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

			// MC Background Breakdown

			for (int WhichPlot = 0; WhichPlot < MCBkgPlots; WhichPlot ++) {

				PrintMCBkgEvents(FileSample[1],FilePOT[1],FileFinal[1],MCBkgPlotNames[WhichPlot],MCBkgPlotLabels[WhichPlot],NonCC1pPlotsEvents[1]->GetBinContent(1));
		
			}

			// -----------------------------------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------------

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

		cout << endl << endl;
	
	} // End of the loop over the runs

} // End of the program 
