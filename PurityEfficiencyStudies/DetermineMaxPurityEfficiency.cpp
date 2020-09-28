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

#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

void DetermineMaxPurityEfficiency() {

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	

	int TextFont = 132;
	double TextSize = 0.07;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString UserID = "apapadop";

	TString PlotsPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/Overlay9/"; 

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString RunNumber = "Run1";

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TFile* TruthCC1pFile = TFile::Open("OutputFiles/TruthPurityEfficiencyStudies_Overlay9_Run1_NoCuts.root");
	TH1D* hTruthCC1pOverlay = (TH1D*)(TruthCC1pFile->Get("TrueMuonCosThetaPlot"));
	double NTruthCC1pOverlay = hTruthCC1pOverlay->Integral();

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TFile* RecoOverlayFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_Overlay9_Run1_NoCuts.root");
	TFile* RecoOverlayDirtFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_OverlayDirt9_Run1_NoCuts.root");
	TFile* RecoBeamOnFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_BeamOn9_Run1_NoCuts.root");
	TFile* RecoExtBNBFile = TFile::Open("OutputFiles/PurityEfficiencyStudies_ExtBNB9_Run1_NoCuts.root");

	// --------------------------------------------------------------------------------------------------------------------------------

	int NBinsNuScore = 20;
	double MinNuScore = 0., MaxNuScore = 1.;
	double NuScoreStep = (MaxNuScore - MinNuScore) / double(NBinsNuScore);

	int NBinsLL = 20;
	double MinLL = -5., MaxLL = 5.;
	double LLStep = (MaxLL - MinLL) / double(NBinsLL);

	// ----------------------------------------------------------------------------------------------------------------------------------------

	double NuScoreArray[NBinsNuScore];
	double LLArray[NBinsLL];

	double PurityArray[NBinsNuScore][NBinsLL];
	double EfficiencyArray[NBinsNuScore][NBinsLL];
	double ProductArray[NBinsNuScore][NBinsLL];

	double InvertedPurityArray[NBinsLL][NBinsNuScore];
	double InvertedEfficiencyArray[NBinsLL][NBinsNuScore];
	double InvertedProductArray[NBinsLL][NBinsNuScore];

	// ----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichNuScore = 0; WhichNuScore < NBinsNuScore; WhichNuScore++) {

		double LocalNuScoreThres = MinNuScore + NuScoreStep * WhichNuScore;
		NuScoreArray[WhichNuScore] = LocalNuScoreThres;

		for (int WhichLL = 0; WhichLL < NBinsLL; WhichLL++) {

			double LocalLLThres = MinLL + LLStep * WhichLL;
			LLArray[WhichLL] = LocalLLThres;

			TString PlotNuScoreLLThres = "RecoMuonCosThetaPlot_NuScoreThres_"+TString(std::to_string(WhichNuScore))+"_LLThres_"+TString(std::to_string(WhichLL));

			TH1D* hRecoCC1pOverlay = (TH1D*)(RecoOverlayFile->Get("CC1p"+PlotNuScoreLLThres));
			TH1D* hRecoOverlay = (TH1D*)(RecoOverlayFile->Get(PlotNuScoreLLThres));
			TH1D* hRecoDirt = (TH1D*)(RecoOverlayDirtFile->Get(PlotNuScoreLLThres));
			TH1D* hRecoBeamOn = (TH1D*)(RecoBeamOnFile->Get(PlotNuScoreLLThres));
			TH1D* hRecoExtBNB = (TH1D*)(RecoExtBNBFile->Get(PlotNuScoreLLThres));

			double NRecoCC1pOverlay = hRecoCC1pOverlay->Integral();	
			double NRecoOverlay = hRecoOverlay->Integral();
			double NRecoDirt = hRecoDirt->Integral();
			double NRecoBeamOn = hRecoBeamOn->Integral();
			double NRecoExtBNB = hRecoExtBNB->Integral();

			double purity = NRecoCC1pOverlay / (NRecoOverlay + NRecoDirt + NRecoExtBNB) ;
			double efficiency = NRecoCC1pOverlay / NTruthCC1pOverlay;
			double product = purity * efficiency;
			
			PurityArray[WhichNuScore][WhichLL] = purity;
			EfficiencyArray[WhichNuScore][WhichLL] = efficiency;
			ProductArray[WhichNuScore][WhichLL] = product;

			InvertedPurityArray[WhichLL][WhichNuScore] = purity;
			InvertedEfficiencyArray[WhichLL][WhichNuScore] = efficiency;
			InvertedProductArray[WhichLL][WhichNuScore] = product;

		}

	}

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TCanvas* NuScorePurityCanvas = new TCanvas("NuScorePurityCanvas","NuScorePurityCanvas",205,34,1024,768);
	NuScorePurityCanvas->SetBottomMargin(0.15);
	NuScorePurityCanvas->SetLeftMargin(0.15);

	TLegend* legNuScorePurity = new TLegend(0.35,0.15,0.9,0.5);
	legNuScorePurity->SetNColumns(3);

	for (int WhichLL = 0; WhichLL < NBinsLL; WhichLL++) {

		double LocalLLThres = MinLL + LLStep * WhichLL;

		TGraph* NuScorePurityGraph = new TGraph(NBinsNuScore,NuScoreArray,InvertedPurityArray[WhichLL]);

		NuScorePurityGraph->SetLineColor(WhichLL+1);
		NuScorePurityGraph->SetMarkerStyle(8);
		NuScorePurityGraph->SetMarkerColor(5*WhichLL+1);
		NuScorePurityGraph->SetMarkerSize(2.);

		NuScorePurityGraph->GetXaxis()->SetNdivisions(6);
		NuScorePurityGraph->GetXaxis()->SetRangeUser(-0.05,1.);
		NuScorePurityGraph->GetXaxis()->SetTitle("#nu Score Threshold");
		NuScorePurityGraph->GetXaxis()->SetTitleFont(TextFont);
		NuScorePurityGraph->GetXaxis()->SetTitleSize(TextSize);
		NuScorePurityGraph->GetXaxis()->SetLabelFont(TextFont);
		NuScorePurityGraph->GetXaxis()->SetLabelSize(TextSize);

		NuScorePurityGraph->GetYaxis()->SetRangeUser(0,1.);
		NuScorePurityGraph->GetYaxis()->SetNdivisions(6);
		NuScorePurityGraph->GetYaxis()->SetTitle("Purity");
		NuScorePurityGraph->GetYaxis()->SetTitleOffset(0.9);
		NuScorePurityGraph->GetYaxis()->SetTitleFont(TextFont);
		NuScorePurityGraph->GetYaxis()->SetTitleSize(TextSize);
		NuScorePurityGraph->GetYaxis()->SetLabelFont(TextFont);
		NuScorePurityGraph->GetYaxis()->SetLabelSize(TextSize);

		NuScorePurityGraph->SetTitle();

		if (WhichLL == 0) { NuScorePurityGraph->Draw("ap"); }
		else { NuScorePurityGraph->Draw("p"); }
	
		TString LegendText = "LL > "+ToString(LocalLLThres);
		legNuScorePurity->AddEntry(NuScorePurityGraph,LegendText,"p");


	}

	legNuScorePurity->SetBorderSize(0);
	legNuScorePurity->SetTextSize(TextSize-0.03);
	legNuScorePurity->SetTextFont(TextFont);
	legNuScorePurity->Draw();

	NuScorePurityCanvas->SaveAs(PlotsPath+"TwoDScanPurity.pdf");

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TCanvas* NuScoreEfficiencyCanvas = new TCanvas("NuScoreEfficiencyCanvas","NuScoreEfficiencyCanvas",205,34,1024,768);
	NuScoreEfficiencyCanvas->SetBottomMargin(0.15);
	NuScoreEfficiencyCanvas->SetLeftMargin(0.15);

	TLegend* legNuScoreEfficiency = new TLegend(0.35,0.53,0.9,0.89);
	legNuScoreEfficiency->SetNColumns(3);

	for (int WhichLL = 0; WhichLL < NBinsLL; WhichLL++) {

		double LocalLLThres = MinLL + LLStep * WhichLL;

		TGraph* NuScoreEfficiencyGraph = new TGraph(NBinsNuScore,NuScoreArray,InvertedEfficiencyArray[WhichLL]);

		NuScoreEfficiencyGraph->SetLineColor(WhichLL+1);
		NuScoreEfficiencyGraph->SetMarkerStyle(8);
		NuScoreEfficiencyGraph->SetMarkerColor(5*WhichLL+1);
		NuScoreEfficiencyGraph->SetMarkerSize(2.);

		NuScoreEfficiencyGraph->GetXaxis()->SetNdivisions(6);
		NuScoreEfficiencyGraph->GetXaxis()->SetRangeUser(-0.05,1.);
		NuScoreEfficiencyGraph->GetXaxis()->SetTitle("#nu Score Threshold");
		NuScoreEfficiencyGraph->GetXaxis()->SetTitleFont(TextFont);
		NuScoreEfficiencyGraph->GetXaxis()->SetTitleSize(TextSize);
		NuScoreEfficiencyGraph->GetXaxis()->SetLabelFont(TextFont);
		NuScoreEfficiencyGraph->GetXaxis()->SetLabelSize(TextSize);

		NuScoreEfficiencyGraph->GetYaxis()->SetRangeUser(0,0.5);
		NuScoreEfficiencyGraph->GetYaxis()->SetNdivisions(6);
		NuScoreEfficiencyGraph->GetYaxis()->SetTitle("Efficiency");
		NuScoreEfficiencyGraph->GetYaxis()->SetTitleOffset(0.9);
		NuScoreEfficiencyGraph->GetYaxis()->SetTitleFont(TextFont);
		NuScoreEfficiencyGraph->GetYaxis()->SetTitleSize(TextSize);
		NuScoreEfficiencyGraph->GetYaxis()->SetLabelFont(TextFont);
		NuScoreEfficiencyGraph->GetYaxis()->SetLabelSize(TextSize);

		NuScoreEfficiencyGraph->SetTitle();

		if (WhichLL == 0) { NuScoreEfficiencyGraph->Draw("ap"); }
		else { NuScoreEfficiencyGraph->Draw("p"); }
	
		TString LegendText = "LL > "+ToString(LocalLLThres);
		legNuScoreEfficiency->AddEntry(NuScoreEfficiencyGraph,LegendText,"p");


	}

	legNuScoreEfficiency->SetBorderSize(0);
	legNuScoreEfficiency->SetTextSize(TextSize-0.03);
	legNuScoreEfficiency->SetTextFont(TextFont);
	legNuScoreEfficiency->Draw();

	NuScoreEfficiencyCanvas->SaveAs(PlotsPath+"TwoDScanEfficiency.pdf");

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TCanvas* NuScoreProductCanvas = new TCanvas("NuScoreProductCanvas","NuScoreProductCanvas",205,34,1024,768);
	NuScoreProductCanvas->SetBottomMargin(0.15);
	NuScoreProductCanvas->SetLeftMargin(0.15);

	TLegend* legNuScoreProduct = new TLegend(0.35,0.55,0.9,0.89);
	legNuScoreProduct->SetNColumns(3);

	for (int WhichLL = 0; WhichLL < NBinsLL; WhichLL++) {

		double LocalLLThres = MinLL + LLStep * WhichLL;

		TGraph* NuScoreProductGraph = new TGraph(NBinsNuScore,NuScoreArray,InvertedProductArray[WhichLL]);

		NuScoreProductGraph->SetLineColor(WhichLL+1);
		NuScoreProductGraph->SetMarkerStyle(8);
		NuScoreProductGraph->SetMarkerColor(5*WhichLL+1);
		NuScoreProductGraph->SetMarkerSize(2.);

		NuScoreProductGraph->GetXaxis()->SetNdivisions(6);
		NuScoreProductGraph->GetXaxis()->SetRangeUser(-0.05,1.);
		NuScoreProductGraph->GetXaxis()->SetTitle("#nu Score Threshold");
		NuScoreProductGraph->GetXaxis()->SetTitleFont(TextFont);
		NuScoreProductGraph->GetXaxis()->SetTitleSize(TextSize);
		NuScoreProductGraph->GetXaxis()->SetLabelFont(TextFont);
		NuScoreProductGraph->GetXaxis()->SetLabelSize(TextSize);

		NuScoreProductGraph->GetYaxis()->SetRangeUser(0,0.25);
		NuScoreProductGraph->GetYaxis()->SetNdivisions(6);
		NuScoreProductGraph->GetYaxis()->SetTitle("Purity*Efficiency");
		NuScoreProductGraph->GetYaxis()->SetTitleOffset(0.95);
		NuScoreProductGraph->GetYaxis()->SetTitleFont(TextFont);
		NuScoreProductGraph->GetYaxis()->SetTitleSize(TextSize);
		NuScoreProductGraph->GetYaxis()->SetLabelFont(TextFont);
		NuScoreProductGraph->GetYaxis()->SetLabelSize(TextSize);

		NuScoreProductGraph->SetTitle();

		if (WhichLL == 0) { NuScoreProductGraph->Draw("ap"); }
		else { NuScoreProductGraph->Draw("p"); }
	
		TString LegendText = "LL > "+ToString(LocalLLThres);
		legNuScoreProduct->AddEntry(NuScoreProductGraph,LegendText,"p");


	}

	legNuScoreProduct->SetBorderSize(0);
	legNuScoreProduct->SetTextSize(TextSize-0.03);
	legNuScoreProduct->SetTextFont(TextFont);
	legNuScoreProduct->Draw();

	NuScoreProductCanvas->SaveAs(PlotsPath+"TwoDScanProduct.pdf");

	// ----------------------------------------------------------------------------------------------------------------------------------------


} // End of the program 
