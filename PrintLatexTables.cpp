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

void PrintEvents(vector<TFile*> FileSample,vector<double> FilePOT,vector<double> FileSamdef,vector<double> FilePreSel, TString SelectionStage, bool PreSelection = true) {

	std::cout << std::fixed << std::setprecision(1);

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

		cout << BeamOnEvents << " [" << BeamOnEvents*FilePOT[0] << "] (" << BeamOnFracSamdef << "\\%) & ";
		cout << ExtBNBEvents << " [" << ExtBNBEvents*FilePOT[2] << "] (" << ExtBNBFracSamdef << "\\%) & ";
		cout << MCEvents << " [" << MCEvents*FilePOT[1] << "] (" << MCFracSamdef << "\\%) & ";
		cout << DirtEvents << " [" << DirtEvents*FilePOT[3] << "] (" << DirtFracSamdef << "\\%)";

	} else {

		cout << BeamOnEvents << " [" << BeamOnEvents*FilePOT[0] << "] (" << BeamOnFracSamdef << "\\%/" << BeamOnFracPreSel << "\\%) & ";
		cout << ExtBNBEvents << " [" << ExtBNBEvents*FilePOT[2] << "] (" << ExtBNBFracSamdef << "\\%/" << ExtBNBFracPreSel << "\\%) & ";
		cout << MCEvents << " [" << MCEvents*FilePOT[1] << "] (" << MCFracSamdef << "\\%/" << MCFracPreSel << "\\%) & ";
		cout << DirtEvents << " [" << DirtEvents*FilePOT[3] << "] (" << DirtFracSamdef << "\\%/" << DirtFracPreSel << "\\%)";

	}

	cout << " \\tabularnewline \\hline"<< endl;

}

void PrintLatexTables(TString BaseMC = "") {

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::vector<TString> PlotNames; PlotNames.clear();
	std::vector<TString> PlotLabels; PlotLabels.clear();

	PlotNames.push_back("SamdefEventPlot"); PlotLabels.push_back("Initial");
	PlotNames.push_back("SWTriggerEventPlot"); PlotLabels.push_back("SW Trigger");
	PlotNames.push_back("OneNuMuPFParticleEventPlot"); PlotLabels.push_back("1 $\\nu_{\\mu}$ PFParticle");
	PlotNames.push_back("MatchedTrackLikeDaughterEventPlot"); PlotLabels.push_back("2 daughters");
	PlotNames.push_back("OneBeamFlashEventPlot"); PlotLabels.push_back("1 beam flash");
	PlotNames.push_back("MomentumThresholdEventPlot"); PlotLabels.push_back("P threshold");
	PlotNames.push_back("StartPointContainmentEventPlot"); PlotLabels.push_back("Start point in FV");
	PlotNames.push_back("ProtonEndPointContainmentEventPlot"); PlotLabels.push_back("End point p in FV");
	PlotNames.push_back("VertexContainmentEventPlot"); PlotLabels.push_back("Vertex in FV");
	PlotNames.push_back("MuonQualityEventPlot"); PlotLabels.push_back("Muon quality cut");
	PlotNames.push_back("PidEventPlot"); PlotLabels.push_back("Calorimetry");
	PlotNames.push_back("NuScoreEventPlot"); PlotLabels.push_back("$\\nu$ score");
	PlotNames.push_back("CommonEventPlot"); PlotLabels.push_back("Common events");

	const int N1DPlots = PlotNames.size();
//	cout << "Number of 1D Plots = " << N1DPlots << endl;

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
			
		cout << endl << endl;

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];		

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles+"/"+Cuts+"/";

			TH1D::SetDefaultSumw2();

			// ---------------------------------------------------------------------------------------------------------------------

//			vector<vector<TH1D*> > Plots; Plots.clear();
//			vector<vector<TH1D*> > CC1pPlots; CC1pPlots.clear();
//			vector<vector<TH1D*> > NonCC1pPlots; NonCC1pPlots.clear();

//			vector<vector<TH1D*> > hratio;  hratio.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;

			// 0: BeamOn
			// 1: Overlay
			// 2: ExtBNB
			// 3: Dirt
		
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");

			if (BaseMC == "") { NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }
			else if (BaseMC == "Overlay9NuWro") { NameOfSamples.push_back("STVStudies_Overlay9NuWro_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("NuWro"); }

			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");
			NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt");

			// ---------------------------------------------------------------------------------------------------------------------

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();
			vector<double> FilePOT; FilePOT.clear();
			vector<double> FileSamdef; FileSamdef.clear();
			vector<double> FilePreSel; FilePreSel.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				FileSample.push_back(TFile::Open(PathToFilesCut+NameOfSamples[WhichSample]));

				TH1D* POTPlot = (TH1D*)(FileSample[WhichSample]->Get("POTScalePlot"));
				double POTCount = POTPlot->GetBinContent(1);
				FilePOT.push_back(POTCount);

				TH1D* SamdefPlot = (TH1D*)(FileSample[WhichSample]->Get("SamdefEventPlot"));
				double SamdefCount = SamdefPlot->GetBinContent(1);
				FileSamdef.push_back(SamdefCount);

				TH1D* PreSelPlot = (TH1D*)(FileSample[WhichSample]->Get("VertexContainmentEventPlot"));
				double PreSelCount = PreSelPlot->GetBinContent(1);
				FilePreSel.push_back(PreSelCount);

//				vector<TH1D*> CurrentPlots; CurrentPlots.clear();
//				vector<TH1D*> CC1pCurrentPlots; CC1pCurrentPlots.clear();
//				vector<TH1D*> NonCC1pCurrentPlots; NonCC1pCurrentPlots.clear();

//				vector<TH1D*> Currenthratio;  Currenthratio.clear();

//				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

//					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
//					TH1D* CC1phist = (TH1D*)(FileSample[WhichSample]->Get("CC1p"+PlotNames[WhichPlot]));
//					TH1D* NonCC1phist = (TH1D*)(FileSample[WhichSample]->Get("NonCC1p"+PlotNames[WhichPlot]));

//					CurrentPlots.push_back(hist);
//					CC1pCurrentPlots.push_back(CC1phist);
//					NonCC1pCurrentPlots.push_back(NonCC1phist);

//					Currenthratio.push_back((TH1D*)hist->Clone());
			
//				}

//				Plots.push_back(CurrentPlots);
//				CC1pPlots.push_back(CC1pCurrentPlots);
//				NonCC1pPlots.push_back(NonCC1pCurrentPlots);

//				hratio.push_back(Currenthratio);

			} // End of the loop over the samples

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				bool PreSelection = true;
				if (PlotLabels[WhichPlot] == "Muon quality cut" || PlotLabels[WhichPlot] == "Calorimetry" 
				 || PlotLabels[WhichPlot] == "\\nu score" || PlotLabels[WhichPlot] == "Common events") { PreSelection = false; }


				cout << PlotLabels[WhichPlot] << " & ";
				PrintEvents(FileSample,FilePOT,FileSamdef,FilePreSel,PlotNames[WhichPlot],PreSelection);

				if (PlotLabels[WhichPlot] == "Initial" || PlotLabels[WhichPlot] == "Vertex in FV"  || PlotLabels[WhichPlot] == "Common events") 
					{ cout << "\\hline" << endl; }

				if (WhichPlot == N1DPlots-1) {

					cout << "Final & ";
					PrintEvents(FileSample,FilePOT,FileSamdef,FilePreSel,PlotNames[WhichPlot],PreSelection);

				}
		

/*
				// -------------------------------------------------------------------------------------------------------------------

				hratio[1][WhichPlot]->Add(hratio[2][WhichPlot]);
				hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
				hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

				hratio[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{BeamOn}{MC+ExtBNB}");
				hratio[0][WhichPlot]->GetXaxis()->SetTitle(Plots[0][WhichPlot]->GetXaxis()->GetTitle());
				hratio[0][WhichPlot]->GetXaxis()->SetTitleSize(0.13);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelSize(0.12);
				hratio[0][WhichPlot]->GetXaxis()->SetTitleOffset(0.88);
				hratio[0][WhichPlot]->GetXaxis()->SetNdivisions(8);

				hratio[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.9*hratio[0][WhichPlot]->GetMinimum(),1.1*hratio[0][WhichPlot]->GetMaximum());
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(6);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.35);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleSize(0.1);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				hratio[0][WhichPlot]->Draw("e1 hist same");

				double RatioMin = hratio[0][WhichPlot]->GetXaxis()->GetXmin();
				double RatioMax = hratio[0][WhichPlot]->GetXaxis()->GetXmax();
				double YRatioCoord = 1.2;
				TLine* RatioLine = new TLine(RatioMin,YRatioCoord,RatioMax,YRatioCoord);
				RatioLine->SetLineWidth(4);
				RatioLine->SetLineColor(kPink+8);
				RatioLine->SetLineStyle(4);
				//RatioLine->Draw("same");
			
				topPad->cd();
				leg[WhichPlot]->SetTextSize(0.5);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->Draw();

				// --------------------------------------------------------------------------------------

				// Sum of NonBeamOn Samples

				TH1D* SumNonBeamOn = (TH1D*)Plots[1][WhichPlot]->Clone(); // ExtBNB
				SumNonBeamOn->Add(Plots[2][WhichPlot]); // Overlay
				SumNonBeamOn->Add(Plots[3][WhichPlot]); // Dirt

				// CC1p Purity 

				int CC1pPurity = CC1pPlots[1][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexPurity;
				latexPurity.SetTextFont(FontStyle);
				latexPurity.SetTextSize(0.09);
				TString LabelPurity = "CC1p = " + ToString(CC1pPurity/10.) + " %";
				latexPurity.DrawLatexNDC(0.59,0.9, LabelPurity);

				// --------------------------------------------------------------------------------------

				// Cosmic Contamination

				int CosmicContamination = Plots[2][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexCosmic;
				latexCosmic.SetTextFont(FontStyle);
				latexCosmic.SetTextSize(0.09);
				TString LabelCosmic = "Cosmics = " + ToString(CosmicContamination/10.) + " %";
				latexCosmic.DrawLatexNDC(0.59,0.8, LabelCosmic);

*/
			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

		cout << endl << endl;
	
	} // End of the loop over the runs

} // End of the program 
