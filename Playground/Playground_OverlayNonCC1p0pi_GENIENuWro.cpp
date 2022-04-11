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

#include "../../Secondary_Code/myFunctions.cpp"
#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void Playground_OverlayNonCC1p0pi_GENIENuWro(TString BaseMC = "") {

	//----------------------------------------//

	double TextSize = 0.07;
	int TextFont = 132;	

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetTitleSize(TextSize,"t");	
	gStyle->SetTitleFont(TextFont,"t");	

	//----------------------------------------//

	vector<TString> PlotNames; PlotNames.clear();
	PlotNames.push_back("RecoDeltaAlphaTPlot");

	const int N1DPlots = PlotNames.size();

	//----------------------------------------//

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	//----------------------------------------//

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//----------------------------------------//

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);	

		TFile* GENIEFile = new TFile("/uboone/data/users/apapadop/myEvents/OutputFiles/v08_00_00_52/_NoCuts_PID_NuScore/STVStudies_Overlay9_Combined_NoCuts_PID_NuScore.root","readonly");
		TFile* NuWroFile = new TFile("/uboone/data/users/apapadop/myEvents/OutputFiles/v08_00_00_52/_NoCuts_PID_NuScore/STVStudies_Overlay9NuWro_Combined_NoCuts_PID_NuScore.root","readonly");					

		//----------------------------------------//

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			//----------------------------------------//			

			TH1D* GENIEHisto = (TH1D*)(GENIEFile->Get("NonCC1p"+PlotNames[WhichPlot]));
			TH1D* NuWroHisto = (TH1D*)(NuWroFile->Get("NonCC1p"+PlotNames[WhichPlot]));

			//----------------------------------------//

			TString PlotCanvasName = "Canvas_"+PlotNames[WhichPlot];
			TCanvas* PlotCanvas = new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetLeftMargin(0.15);

			//----------------------------------------//

			TLegend* leg = new TLegend(0.25,0.75,0.35,0.85);
			leg->SetBorderSize(0);
			leg->SetNColumns(1);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(TextFont);															

			//----------------------------------------//	

			GENIEHisto->SetTitle("NonCC1p0#pi Background");
			GENIEHisto->SetLineWidth(3);			

			GENIEHisto->GetXaxis()->CenterTitle();
			GENIEHisto->GetXaxis()->SetTitleSize(TextSize);
			GENIEHisto->GetXaxis()->SetLabelSize(TextSize);
			GENIEHisto->GetXaxis()->SetTitleFont(TextFont);
			GENIEHisto->GetXaxis()->SetLabelFont(TextFont);
			GENIEHisto->GetXaxis()->SetNdivisions(8);

			GENIEHisto->GetYaxis()->CenterTitle();
			GENIEHisto->GetYaxis()->SetTitleSize(TextSize);
			GENIEHisto->GetYaxis()->SetLabelSize(TextSize);
			GENIEHisto->GetYaxis()->SetTitleFont(TextFont);
			GENIEHisto->GetYaxis()->SetLabelFont(TextFont);
			GENIEHisto->GetYaxis()->SetTitleOffset(0.95);			
			GENIEHisto->GetYaxis()->SetNdivisions(8);
			GENIEHisto->GetYaxis()->SetRangeUser(0.,1.15*GENIEHisto->GetMaximum());
			GENIEHisto->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT));													

			GENIEHisto->SetLineColor(OverlayColor);		
			GENIEHisto->Draw("hist same");
			leg->AddEntry(GENIEHisto,"GENIE","l");

			NuWroHisto->SetLineWidth(3);
			NuWroHisto->SetLineColor(kOrange+7);		
			NuWroHisto->Draw("hist same");
			leg->AddEntry(NuWroHisto,"NuWro","l");			

			//----------------------------------------//

			leg->Draw("same");

			//----------------------------------------//			

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
