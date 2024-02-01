#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void print_2d() {

	// -------------------------------------------------------------------------------------

	TH2D::SetDefaultSumw2();
	
	double TextSize = 0.07;
	const Int_t NCont = 999;	
	gStyle->SetPalette(55); 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_CRT");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		CutExtension = CutExtension + VectorCuts[i];

	}

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;

	PlotNames.push_back("RecoThetaZRecoECalPlot"); 
	PlotNames.push_back("RecoECalTrueECalPlot"); 
	PlotNames.push_back("RecoECalTrueEnuPlot"); 

	PlotNames.push_back("TrueThetaZTrueECalPlot"); 
	PlotNames.push_back("TrueThetaZTrueEnuPlot"); 
	
	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	vector<TString> NameOfSamples;
	NameOfSamples.push_back("Overlay9");
	const int NSamples = NameOfSamples.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");				

	const int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;	

	// -------------------------------------------------------------------------------------

	vector< vector<TFile*>> FileSample;
	FileSample.resize(NSamples, vector<TFile*>(NRuns));
	vector< vector <TH2D*> > Plots;
	Plots.resize(NSamples, vector<TH2D*>(N2DPlots));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -------------------------------------------------------------------------------------

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			TString ExactFileLocation = PathToFiles+CutExtension;

			FileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+\
							  Runs[WhichRun]+CutExtension+".root");

		}

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {


				Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample][WhichRun]->Get("POTScaledCC1p"+PlotNames[WhichPlot]+"2D"));
				
				// ---------------------------------------------------------------------------------------				
	
				TString PlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+NameOfSamples[WhichSample];
				TCanvas* PlotCanvas = new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768);
				PlotCanvas->cd();
				PlotCanvas->SetBottomMargin(0.16);
				PlotCanvas->SetLeftMargin(0.15);
				PlotCanvas->SetRightMargin(0.15);				
					
				gStyle->SetMarkerSize(1.5);
				gStyle->SetPaintTextFormat("4.2f");				
					
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);				
				Plots[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(8);
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelOffset(0.02);				
					
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);				
				Plots[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(8);
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.);				
									
				Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(TextSize);
				Plots[WhichSample][WhichPlot]->GetZaxis()->SetNdivisions(8);				
//				Plots[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(-0.1,1.);

				Plots[WhichSample][WhichPlot]->SetMarkerColor(kWhite);				
				Plots[WhichSample][WhichPlot]->SetMarkerSize(0.9);
				//Plots[WhichSample][WhichPlot]->SetTitle(Runs[WhichRun]);					
				Plots[WhichSample][WhichPlot]->Draw("colz"); 
					
				PlotCanvas->SaveAs(PlotPath+NameOfSamples[0]+"/Atmospherics_2D_"+PlotNames[WhichPlot]
						+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
					
				delete PlotCanvas;				
				 
			} // End of the loop over the plots
			
			FileSample[WhichSample][WhichRun]->Close();

		} // End of the loop over the samples

	} // End of the loop over the runs	

} // End of the program
