#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>

#include <iostream>
#include <vector>

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void CC1pPurityEfficiency() {

	std::vector<TString> PlotNames; PlotNames.clear();

	// -----------------------------------------------------------------------------------------------------------------------------------------

	PlotNames.push_back("MuonMomentumPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToTrueFiles = "OutputFiles/"+UBCodeVersion+"/";
			TString PathToFiles = "OutputFiles/"+UBCodeVersion+"/"+Cuts+"/";

			TH1D::SetDefaultSumw2();

			// ---------------------------------------------------------------------------------------------------------------------

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > CC1pPlots; CC1pPlots.clear();
			vector<vector<TH1D*> > NonCC1pPlots; NonCC1pPlots.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;
			vector<double> POTScaleOfSamples; POTScaleOfSamples.clear();

			// 0: BeamOn
			// 1: Overlay
			// 2: ExtBNB
			// 3: Dirt
		
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");
			NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay");
			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");
			NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt");
			NameOfSamples.push_back("TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"); LabelsOfSamples.push_back("Overlay");

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				if (string(NameOfSamples[WhichSample]).find("Truth") != std::string::npos) { FileSample.push_back(TFile::Open(PathToTrueFiles+NameOfSamples[WhichSample])); }
				else { FileSample.push_back(TFile::Open(PathToFiles+NameOfSamples[WhichSample])); }

				vector<TH1D*> CurrentPlots; CurrentPlots.clear();
				vector<TH1D*> CC1pCurrentPlots; CC1pCurrentPlots.clear();
				vector<TH1D*> NonCC1pCurrentPlots; NonCC1pCurrentPlots.clear();

				TH1D* POTScalePlot = (TH1D*)(FileSample[WhichSample]->Get("POTScalePlot"));
				double POTScale = POTScalePlot->GetBinContent(1);
				POTScaleOfSamples.push_back(POTScale);
				//cout << NameOfSamples[WhichSample] << " POTScale Factor = " << POTScale << endl;

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					TH1D* CC1phist = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
					if (string(NameOfSamples[WhichSample]).find("Truth") != std::string::npos) { CC1phist = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot])); }
					TH1D* NonCC1phist = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pReco"+PlotNames[WhichPlot]));

					CurrentPlots.push_back(hist);
					CC1pCurrentPlots.push_back(CC1phist);
					NonCC1pCurrentPlots.push_back(NonCC1phist);
			
				}

				Plots.push_back(CurrentPlots);
				CC1pPlots.push_back(CC1pCurrentPlots);
				NonCC1pPlots.push_back(NonCC1pCurrentPlots);

			}

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				// --------------------------------------------------------------------------------------

				// True CC1p

				double TrueCC1pEventsPOTScaled = CC1pPlots[4][WhichPlot]->Integral();
				double TrueCC1pEvents = CC1pPlots[4][WhichPlot]->Integral() / POTScaleOfSamples[4];
				double TrueCC1pError = TMath::Sqrt(TrueCC1pEvents);
				double TrueCC1pErrorPOTScaled = TrueCC1pError * POTScaleOfSamples[4];

				// Reco CC1p

				double CC1pEventsPOTScaled = CC1pPlots[1][WhichPlot]->Integral();
				double CC1pEvents = CC1pEventsPOTScaled / POTScaleOfSamples[1];
				double CC1pError = TMath::Sqrt(CC1pEvents);
				double CC1pErrorPOTScaled = CC1pError * POTScaleOfSamples[1];
				//cout << "CC1p events = " << CC1pEventsPOTScaled << " CC1p Error = " << CC1pErrorPOTScaled << endl;

				// Reco Overlay

				double OverlayEventsPOTScaled = Plots[1][WhichPlot]->Integral();
				double OverlayEvents = OverlayEventsPOTScaled / POTScaleOfSamples[1];
				double OverlayError = TMath::Sqrt(OverlayEvents);
				double OverlayErrorPOTScaled = OverlayError * POTScaleOfSamples[1];

				// Reco Dirt

				double DirtEventsPOTScaled = Plots[3][WhichPlot]->Integral();
				double DirtEvents = Plots[3][WhichPlot]->Integral() / POTScaleOfSamples[3];
				double DirtError = TMath::Sqrt(DirtEvents);
				double DirtErrorPOTScaled = DirtError * POTScaleOfSamples[3];

				// Reco ExtBNB

				double ExtBNBEventsPOTScaled = Plots[2][WhichPlot]->Integral();
				double ExtBNBEvents = ExtBNBEventsPOTScaled / POTScaleOfSamples[2];
				double ExtBNBError = TMath::Sqrt(ExtBNBEvents);
				double ExtBNBErrorPOTScaled = ExtBNBError * POTScaleOfSamples[2];

				// Sum of NonBeamOn Samples

				TH1D* SumNonBeamOn = (TH1D*)Plots[1][WhichPlot]->Clone(); // Overlay
				SumNonBeamOn->Add(Plots[2][WhichPlot]); // ExtBNB
				SumNonBeamOn->Add(Plots[3][WhichPlot]); // Dirt
				double SumNonBeamOnEventsPOTScaled = SumNonBeamOn->Integral();
				double SumNonBeamOnEventsErrorPOTScaled = TMath::Sqrt(\
										TMath::Power(OverlayErrorPOTScaled,2.) + \
										TMath::Power(DirtErrorPOTScaled,2.) + \
										TMath::Power(ExtBNBErrorPOTScaled,2.) );

				// --------------------------------------------------------------------------------------

				// CC1p Efficiency 

				int CC1pEfficiency = CC1pEventsPOTScaled / TrueCC1pEventsPOTScaled * 1000.;
				int CC1pEfficiencyError = CC1pEfficiency * TMath::Sqrt( TMath::Power(CC1pErrorPOTScaled/CC1pEventsPOTScaled,2.) + \
										TMath::Power(TrueCC1pErrorPOTScaled/TrueCC1pEventsPOTScaled,2.) );

				TString LabelEfficiency = "CC1p Efficiency = " + ToString(CC1pEfficiency/10.) + " +/- "+ ToString(CC1pEfficiencyError/10.) +" %";

				// --------------------------------------------------------------------------------------

				// CC1p Purity 

				int CC1pPurity = CC1pEventsPOTScaled / SumNonBeamOnEventsPOTScaled * 1000.;
				int CC1pPurityError = CC1pPurity * \
					TMath::Sqrt( TMath::Power(CC1pErrorPOTScaled/CC1pEventsPOTScaled,2.) + TMath::Power(SumNonBeamOnEventsErrorPOTScaled/SumNonBeamOnEventsPOTScaled,2.) );

				TString LabelPurity = "CC1p Purity = " + ToString(CC1pPurity/10.) + " +/- "+ ToString(CC1pPurityError/10.) +" %";

				// --------------------------------------------------------------------------------------

				// Cosmic Contamination

				int CosmicContamination = ExtBNBEventsPOTScaled / SumNonBeamOnEventsPOTScaled * 1000.;
				int CosmicContaminationError = CosmicContamination * \
					TMath::Sqrt( TMath::Power(ExtBNBErrorPOTScaled/ExtBNBEventsPOTScaled,2.) + TMath::Power(SumNonBeamOnEventsErrorPOTScaled/SumNonBeamOnEventsPOTScaled,2.) );

				TString LabelCosmic = "Cosmics = " + ToString(CosmicContamination/10.) + " +/- "+ ToString(CosmicContaminationError/10.) +" %";

				// --------------------------------------------------------------------------------------

				std::cout << "Applied cuts = " << Cuts << std::endl;
				std::cout << LabelEfficiency << std::endl;
				std::cout << LabelPurity << std::endl;
				std::cout << LabelCosmic << std::endl << std::endl;

				// --------------------------------------------------------------------------------------

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs

} // End of the program 
