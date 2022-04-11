#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TLatex.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

#include "../../Secondary_Code/myFunctions.cpp"
#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void GENIENuWro_Efficiencies() {

	//----------------------------------------//

	double TextSize = 0.07;
	int TextFont = 132;

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	gStyle->SetTitleSize(TextSize,"t");	
	gStyle->SetTitleFont(TextFont,"t");	

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaAlphaTPlot"); 

	//----------------------------------------//

	const int N1DPlots = PlotNames.size();

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

		TFile* TrueGENIEFile = new TFile("/uboone/data/users/apapadop/myEvents/OutputFiles/v08_00_00_52/TruthSTVAnalysis_Overlay9_Combined_v08_00_00_52.root","readonly");
		TFile* TrueNuWroFile = new TFile("/uboone/data/users/apapadop/myEvents/OutputFiles/v08_00_00_52/TruthSTVAnalysis_Overlay9NuWro_Combined_v08_00_00_52.root","readonly");							

		//----------------------------------------//

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			//----------------------------------------//

			TString PlotCanvasName = "Canvas_"+PlotNames[WhichPlot];
			TCanvas* PlotCanvas = new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetLeftMargin(0.15);

			//----------------------------------------//

			TLegend* leg = new TLegend(0.25,0.25,0.35,0.35);
			leg->SetBorderSize(0);
			leg->SetNColumns(1);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(TextFont);																		

			//----------------------------------------//			

			TH1D* GENIEHisto = (TH1D*)(GENIEFile->Get("CC1pTrue"+PlotNames[WhichPlot]));
			TH1D* NuWroHisto = (TH1D*)(NuWroFile->Get("CC1pTrue"+PlotNames[WhichPlot]));

			TH1D* TrueGENIEHisto = (TH1D*)(TrueGENIEFile->Get("True"+PlotNames[WhichPlot]));
			TH1D* TrueNuWroHisto = (TH1D*)(TrueNuWroFile->Get("True"+PlotNames[WhichPlot]));			

			GENIEHisto->Divide(TrueGENIEHisto);
			NuWroHisto->Divide(TrueNuWroHisto);	

			//----------------------------------------//	

			GENIEHisto->SetTitle("CC1p0#pi Efficiency");
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
			GENIEHisto->GetYaxis()->SetTitleOffset(1.15);			
			GENIEHisto->GetYaxis()->SetNdivisions(8);
			GENIEHisto->GetYaxis()->SetRangeUser(0.,1.25*GENIEHisto->GetMaximum());
			GENIEHisto->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT));													

			GENIEHisto->SetMarkerColor(OverlayColor);
			GENIEHisto->SetMarkerSize(2.);
			GENIEHisto->SetMarkerStyle(20);						
			GENIEHisto->SetLineColor(OverlayColor);		
			GENIEHisto->Draw("e same");
			leg->AddEntry(GENIEHisto,"GENIE","lep");

			NuWroHisto->SetLineWidth(3);
			NuWroHisto->SetLineColor(kOrange+7);
			NuWroHisto->SetMarkerColor(kOrange+7);			
			NuWroHisto->SetMarkerSize(2.);
			NuWroHisto->SetMarkerStyle(20);								
			NuWroHisto->Draw("e same");
			leg->AddEntry(NuWroHisto,"NuWro","lep");			

			//----------------------------------------//

			leg->Draw("same");

			//----------------------------------------//								

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
