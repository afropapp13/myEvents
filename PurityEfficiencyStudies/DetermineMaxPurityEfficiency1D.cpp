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
#include <algorithm>
#include <cmath>

#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

std::vector<double> GetMax(double Array[], int SizeArray){

	double max = -99999.;
	double index = -99999.;	

	for (int i = 0; i < SizeArray; i++) {
	
		if (Array[i] > max) { max = Array[i]; index = i;}
	
	}
	
	std::vector<double> MaxQuantities{max,index};
	
	return MaxQuantities;
	
}

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

void DetermineMaxPurityEfficiency1D() {

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	

	int TextFont = 132;
	double TextSize = 0.07;
	int MarkerColor = 610;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> RunNumber;
	RunNumber.push_back("Run1");
	RunNumber.push_back("Run3");

	vector<TString> CutName; vector<TString> Cuts;
	CutName.push_back("LLP"); Cuts.push_back("_NoCuts");
	CutName.push_back("NuScore"); Cuts.push_back("_NoCuts_PID");

	// ----------------------------------------------------------------------------------------------------------------------------------------

	int NRuns = RunNumber.size();
	int NCuts = CutName.size();

	// ----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		for (int WhichCut = 0; WhichCut < NCuts; WhichCut++) {

			// --------------------------------------------------------------------------------------------------------------------------------

			int NBins = -99; double Min = -99., Max = 99.;

			if (CutName[WhichCut] == "LLP") {

				NBins = 50;
				Min = -5., Max = 5.;

			}

			if (CutName[WhichCut] == "LLMu") {

				NBins = 20;
				Min = -5., Max = 5.;

			}

			if (CutName[WhichCut] == "NuScore") {

				NBins = 20;
				Min = 0., Max = 1.;

			}

			if (CutName[WhichCut] == "Length") {

				NBins = 20;
				Min = -150., Max = 500.;

			}

			if (CutName[WhichCut] == "PMissMinus") {

				NBins = 20;
				Min = 0., Max = 1.50;

			}

			if (CutName[WhichCut] == "kMiss") {

				NBins = 20;
				Min = 0., Max = 1.05;

			}

			if (CutName[WhichCut] == "DeltaTheta") {

				NBins = 19;
				Min = 5., Max = 90.;

			}

			double Step = (Max - Min) / double(NBins);

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TFile* TruthCC1pFile = TFile::Open("OutputFiles/TruthPurityEfficiencyStudies_Overlay9_"+RunNumber[WhichRun]+Cuts[0]+".root");
			TH1D* hTruthCC1pOverlay = (TH1D*)(TruthCC1pFile->Get("TrueMuonCosThetaPlot"));
			double NTruthCC1pOverlay = hTruthCC1pOverlay->Integral();

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TFile* RecoOverlayFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_Overlay9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");
			TFile* RecoOverlayDirtFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_OverlayDirt9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");
			TFile* RecoBeamOnFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_BeamOn9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");
			TFile* RecoExtBNBFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_ExtBNB9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");

			// ----------------------------------------------------------------------------------------------------------------------------------------

			double XArray[NBins];

			double PurityArray[NBins];
			double EfficiencyArray[NBins];
			double ProductArray[NBins];

			// ----------------------------------------------------------------------------------------------------------------------------------------

			for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {

				double LocalThres = Min + Step * WhichBin;
				XArray[WhichBin] = LocalThres;

				// TString to grab the plot of interest (1D grid)

				TString PlotThres = "RecoMuonCosThetaPlot_"+CutName[WhichCut]+"Thres_"+TString(std::to_string(WhichBin));

				TH1D* hRecoCC1pOverlay = (TH1D*)(RecoOverlayFile->Get("CC1p"+PlotThres));
				TH1D* hRecoOverlay = (TH1D*)(RecoOverlayFile->Get(PlotThres));
				TH1D* hRecoDirt = (TH1D*)(RecoOverlayDirtFile->Get(PlotThres));
				TH1D* hRecoBeamOn = (TH1D*)(RecoBeamOnFile->Get(PlotThres));
				TH1D* hRecoExtBNB = (TH1D*)(RecoExtBNBFile->Get(PlotThres));

				double NRecoCC1pOverlay = hRecoCC1pOverlay->Integral();	
				double NRecoOverlay = hRecoOverlay->Integral();
				double NRecoDirt = hRecoDirt->Integral();
				double NRecoBeamOn = hRecoBeamOn->Integral();
				double NRecoExtBNB = hRecoExtBNB->Integral();

				double purity = NRecoCC1pOverlay / (NRecoOverlay + NRecoDirt + NRecoExtBNB) ;
				double efficiency = NRecoCC1pOverlay / NTruthCC1pOverlay;
				double product = purity * efficiency;

				PurityArray[WhichBin] = purity;
				EfficiencyArray[WhichBin] = efficiency;
				ProductArray[WhichBin] = product;

			}

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TString TStringProductCanvas = RunNumber[WhichRun]+"_"+CutName[WhichCut]+"ProductCanvas";
			TCanvas* ProductCanvas = new TCanvas(TStringProductCanvas,TStringProductCanvas,205,34,1024,768);
			ProductCanvas->SetBottomMargin(0.15);
			ProductCanvas->SetLeftMargin(0.15);
				
			TGraph* ProductGraph = new TGraph(NBins,XArray,ProductArray);		
				
			std::vector<double> LocalMax = GetMax(ProductArray,NBins);
				
			double GlobalMax = LocalMax.at(0);
			int GlobalThresBin = LocalMax.at(1);

			ProductGraph->SetLineColor(MarkerColor);
			ProductGraph->SetMarkerStyle(8);
			ProductGraph->SetMarkerColor(MarkerColor);
			ProductGraph->SetMarkerSize(2.);

			ProductGraph->GetXaxis()->SetNdivisions(6);
			ProductGraph->GetXaxis()->SetRangeUser(Min,Max);
			ProductGraph->GetXaxis()->SetTitle(CutName[WhichCut]+" Threshold");
			ProductGraph->GetXaxis()->SetTitleFont(TextFont);
			ProductGraph->GetXaxis()->SetTitleSize(TextSize);
			ProductGraph->GetXaxis()->SetLabelFont(TextFont);
			ProductGraph->GetXaxis()->SetLabelSize(TextSize);

			ProductGraph->GetYaxis()->SetRangeUser(0,0.25);
			ProductGraph->GetYaxis()->SetNdivisions(6);
			ProductGraph->GetYaxis()->SetTitle("Purity x Efficiency");
			ProductGraph->GetYaxis()->SetTitleOffset(1.05);
			ProductGraph->GetYaxis()->SetTitleFont(TextFont);
			ProductGraph->GetYaxis()->SetTitleSize(TextSize);
			ProductGraph->GetYaxis()->SetLabelFont(TextFont);
			ProductGraph->GetYaxis()->SetLabelSize(TextSize);

			ProductGraph->SetTitle();
			ProductGraph->Draw("ap");

			// ----------------------------------------------------------------------------------------------------------------------------------------	
			
			double SelectedThres = Min + Step * GlobalThresBin;

			TLatex* latProduct = new TLatex();
			latProduct->SetTextFont(TextFont);
			latProduct->SetTextSize(TextSize);
			latProduct->DrawLatexNDC(0.2,0.82,RunNumber[WhichRun] + " Max at "+CutName[WhichCut]+" > "+ToString(SelectedThres));

			// ----------------------------------------------------------------------------------------------------------------------------------------
			
			// Beam On events that we would have with the candidate set of cuts
			
			TString PlotThres = "RecoMuonCosThetaPlot_"+CutName[WhichCut]+"Thres_"+TString(std::to_string(GlobalThresBin));
			TH1D* hRecoBeamOn = (TH1D*)(RecoBeamOnFile->Get(PlotThres));
			
			TLatex* latBeamOn = new TLatex();
			latBeamOn->SetTextFont(TextFont);
			latBeamOn->SetTextSize(TextSize);
			latBeamOn->DrawLatexNDC(0.2,0.75,"Beam On Events = "+ToString(hRecoBeamOn->Integral()));

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TH1D* hRecoExtBNB = (TH1D*)(RecoExtBNBFile->Get(PlotThres));
			TH1D* hRecoOverlay = (TH1D*)(RecoOverlayFile->Get(PlotThres));
			TH1D* hRecoDirt = (TH1D*)(RecoOverlayDirtFile->Get(PlotThres));

			// ----------------------------------------------------------------------------------------------------------------------------------------
			
			TLatex* latExtBNB = new TLatex();
			latExtBNB->SetTextFont(TextFont);
			latExtBNB->SetTextSize(TextSize);
			latExtBNB->DrawLatexNDC(0.2,0.68,"Cosmics = "+ToString(int(hRecoExtBNB->Integral() / ( hRecoExtBNB->Integral() + hRecoOverlay->Integral() + hRecoDirt->Integral() )*100.)) + "%");

			ProductCanvas->SaveAs(PlotPath+CutName[WhichCut]+"_TwoDScanProduct_"+RunNumber[WhichRun]+".pdf");
			delete ProductCanvas;

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TString TStringPurityCanvas = RunNumber[WhichRun]+"_"+CutName[WhichCut]+"PurityCanvas";
			TCanvas* PurityCanvas = new TCanvas(TStringPurityCanvas,TStringPurityCanvas,205,34,1024,768);
			PurityCanvas->SetBottomMargin(0.15);
			PurityCanvas->SetLeftMargin(0.15);

			TGraph* PurityGraph = new TGraph(NBins,XArray,PurityArray);

			PurityGraph->SetLineColor(MarkerColor);
			PurityGraph->SetMarkerStyle(8);
			PurityGraph->SetMarkerColor(MarkerColor);
			PurityGraph->SetMarkerSize(2.);

			PurityGraph->GetXaxis()->SetNdivisions(6);
			PurityGraph->GetXaxis()->SetRangeUser(Min,Max);
			PurityGraph->GetXaxis()->SetTitle(CutName[WhichCut]+" Threshold");
			PurityGraph->GetXaxis()->SetTitleFont(TextFont);
			PurityGraph->GetXaxis()->SetTitleSize(TextSize);
			PurityGraph->GetXaxis()->SetLabelFont(TextFont);
			PurityGraph->GetXaxis()->SetLabelSize(TextSize);

			PurityGraph->GetYaxis()->SetRangeUser(0,1.);
			PurityGraph->GetYaxis()->SetNdivisions(6);
			PurityGraph->GetYaxis()->SetTitle("Purity");
			PurityGraph->GetYaxis()->SetTitleOffset(0.9);
			PurityGraph->GetYaxis()->SetTitleFont(TextFont);
			PurityGraph->GetYaxis()->SetTitleSize(TextSize);
			PurityGraph->GetYaxis()->SetLabelFont(TextFont);
			PurityGraph->GetYaxis()->SetLabelSize(TextSize);

			PurityGraph->SetTitle();
			PurityGraph->Draw("ap");

			TLatex* latPurity = new TLatex();
			latPurity->SetTextFont(TextFont);
			latPurity->SetTextSize(TextSize);
			latPurity->DrawLatexNDC(0.2,0.82,RunNumber[WhichRun] + " Purity at max = "+ToString(int(PurityArray[GlobalThresBin]*100.)) + "%");

			PurityCanvas->SaveAs(PlotPath+CutName[WhichCut]+"_OneDScanPurity_"+RunNumber[WhichRun]+".pdf");
			delete PurityCanvas;

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TString TStringEfficiencyCanvas = RunNumber[WhichRun]+"_"+CutName[WhichCut]+"EfficiencyCanvas";
			TCanvas* EfficiencyCanvas = new TCanvas(TStringEfficiencyCanvas,TStringEfficiencyCanvas,205,34,1024,768);
			EfficiencyCanvas->SetBottomMargin(0.15);
			EfficiencyCanvas->SetLeftMargin(0.15);

			TGraph* EfficiencyGraph = new TGraph(NBins,XArray,EfficiencyArray);

			EfficiencyGraph->SetLineColor(MarkerColor);
			EfficiencyGraph->SetMarkerStyle(8);
			EfficiencyGraph->SetMarkerColor(MarkerColor);
			EfficiencyGraph->SetMarkerSize(2.);

			EfficiencyGraph->GetXaxis()->SetNdivisions(6);
			EfficiencyGraph->GetXaxis()->SetRangeUser(Min,Max);
			EfficiencyGraph->GetXaxis()->SetTitle(CutName[WhichCut]+" Threshold");
			EfficiencyGraph->GetXaxis()->SetTitleFont(TextFont);
			EfficiencyGraph->GetXaxis()->SetTitleSize(TextSize);
			EfficiencyGraph->GetXaxis()->SetLabelFont(TextFont);
			EfficiencyGraph->GetXaxis()->SetLabelSize(TextSize);

			EfficiencyGraph->GetYaxis()->SetRangeUser(0,0.5);
			EfficiencyGraph->GetYaxis()->SetNdivisions(6);
			EfficiencyGraph->GetYaxis()->SetTitle("Efficiency");
			EfficiencyGraph->GetYaxis()->SetTitleOffset(0.9);
			EfficiencyGraph->GetYaxis()->SetTitleFont(TextFont);
			EfficiencyGraph->GetYaxis()->SetTitleSize(TextSize);
			EfficiencyGraph->GetYaxis()->SetLabelFont(TextFont);
			EfficiencyGraph->GetYaxis()->SetLabelSize(TextSize);

			EfficiencyGraph->SetTitle();
			EfficiencyGraph->Draw("ap"); 

			TLatex* latEfficiency = new TLatex();
			latEfficiency->SetTextFont(TextFont);
			latEfficiency->SetTextSize(TextSize);
			latEfficiency->DrawLatexNDC(0.2,0.82,RunNumber[WhichRun] + " Efficiency at max = "+ToString(int(EfficiencyArray[GlobalThresBin]*100.)) + "%");

			EfficiencyCanvas->SaveAs(PlotPath+CutName[WhichCut]+"_OneDScanEfficiency_"+RunNumber[WhichRun]+".pdf");	
			delete EfficiencyCanvas;

			// ----------------------------------------------------------------------------------------------------------------------------------------

		} // En of the loop over the cuts

	} // En of the loop over the runs

} // End of the program 
