#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

#include "../myCCQEAnalysis/Constants.h"

using namespace std;
using namespace Constants;

void Create2DPlots() {

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
//	VectorCuts.push_back("_NuScore");
//	VectorCuts.push_back("_ThreePlaneLogChi2");
//	VectorCuts.push_back("_MatchedFlash");
//	VectorCuts.push_back("_Collinearity");

//	VectorCuts.push_back("_Chi2");
//	VectorCuts.push_back("_Collinearity");
//	VectorCuts.push_back("_MatchedFlash");
//	VectorCuts.push_back("_Distance");
//	VectorCuts.push_back("_Coplanarity");
//	VectorCuts.push_back("_TransImb");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		Cuts = Cuts + VectorCuts[i];

	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	TString PathToFiles = "OutputFiles/"+UBCodeVersion+"/";

	TH2D::SetDefaultSumw2();

	vector<TString> PlotNames;

	PlotNames.push_back("RecoChi2vsThreePlaneChi2TPlot"); 

	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); gStyle->SetTitleFont(FontStyle,"t"); SetOffsetAndSize();

	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9");

	const int NSamples = NameOfSamples.size();
	TCanvas* PlotCanvas[NSamples][N2DPlots] = {}; TFile* FileSample[NSamples] = {};
	TH2D* Plots[NSamples][N2DPlots] = {};

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		FileSample[WhichSample] = TFile::Open(PathToFiles+"CCQEStudies_"+NameOfSamples[WhichSample]+Cuts+".root");

	}

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {

			TString LogScale = "";

			PlotCanvas[WhichSample][WhichPlot] = new TCanvas(PlotNames[WhichPlot]+NameOfSamples[WhichSample],
					    PlotNames[WhichPlot]+NameOfSamples[WhichSample],205,34,1024,768);
			PlotCanvas[WhichSample][WhichPlot]->cd();
			Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(0.03);
			double ScalingFactor = double(Plots[0][WhichPlot]->GetEntries()) / double(Plots[WhichSample][WhichPlot]->GetEntries());
			Plots[WhichSample][WhichPlot]->Scale(ScalingFactor);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(0,Plots[0][WhichPlot]->GetMaximum());

			CenterAxisTitle(Plots[WhichSample][WhichPlot]);
			PlotCanvas[WhichSample][WhichPlot]->SetLogz();
			Plots[WhichSample][WhichPlot]->Draw("colz");


			PlotCanvas[WhichSample][WhichPlot]->SaveAs("./myPlots/pdf/2D/"+UBCodeVersion+"/"+
								PlotNames[WhichPlot]+NameOfSamples[WhichSample]+LogScale+"_"+WhichRun+"_"+UBCodeVersion+".pdf");
			//delete PlotCanvas[WhichSample][WhichPlot];

		} // End of the loop over the plots

	} // End of the loop over the samples


} // End of the program
