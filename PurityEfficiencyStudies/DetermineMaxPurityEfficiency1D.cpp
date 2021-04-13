#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "../../../myClasses/Constants.h"
#include "../../Secondary_Code/myFunctions.cpp"

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

				NBins = NBinsThreePlaneChi2LogLikelihood;
				Min = MinThreePlaneChi2LogLikelihood, Max = MaxThreePlaneChi2LogLikelihood;

			}

			if (CutName[WhichCut] == "NuScore") {

				NBins = NBinsNuScore;
				Min = MinNuScore, Max = MaxNuScore;

			}

			if (CutName[WhichCut] == "Length") {

				NBins = NBinsMuonLength;
				Min = MinMuonLength, Max = MaxMuonLength;

			}


			double Step = (Max - Min) / double(NBins);

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TFile* TruthCC1pFile = TFile::Open(PathToFiles+"TruthSTVAnalysis_Overlay9_"+RunNumber[WhichRun]+"_"+UBCodeVersion+".root");
			TH1D* hTruthCC1pOverlay = (TH1D*)(TruthCC1pFile->Get("TrueMuonCosThetaPlot"));
			double NTruthCC1pOverlay = hTruthCC1pOverlay->Integral();

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TFile* RecoOverlayFile = TFile::Open(PathToFiles+Cuts[WhichCut]+"/"+"STVStudies_Overlay9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");
			TFile* RecoOverlayDirtFile = TFile::Open(PathToFiles+Cuts[WhichCut]+"/"+"STVStudies_OverlayDirt9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");
			TFile* RecoBeamOnFile = TFile::Open(PathToFiles+Cuts[WhichCut]+"/"+"STVStudies_BeamOn9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");
			TFile* RecoExtBNBFile = TFile::Open(PathToFiles+Cuts[WhichCut]+"/"+"STVStudies_ExtBNB9_"+RunNumber[WhichRun]+Cuts[WhichCut]+".root");

			TH1D* hOverlayPOTScale = (TH1D*)(RecoOverlayFile->Get("POTScalePlot"));
			TH1D* hBeamOnPOTScale = (TH1D*)(RecoBeamOnFile->Get("POTScalePlot"));
			TH1D* hExtBNBPOTScale = (TH1D*)(RecoExtBNBFile->Get("POTScalePlot"));
			TH1D* hOverlayDirtPOTScale = (TH1D*)(RecoOverlayDirtFile->Get("POTScalePlot"));

			double OverlayPOTScale = hOverlayPOTScale->GetBinContent(1);
			double BeamOnPOTScale = hBeamOnPOTScale->GetBinContent(1);
			double ExtBNBPOTScale = hExtBNBPOTScale->GetBinContent(1);
			double OverlayDirtPOTScale = hOverlayDirtPOTScale->GetBinContent(1);

			double ErrorTrueCC1pOverlay = TMath::Sqrt(NTruthCC1pOverlay * OverlayPOTScale);

			// ----------------------------------------------------------------------------------------------------------------------------------------

			double XArray[NBins];

			double PurityArray[NBins];
			double EfficiencyArray[NBins];
			double ProductArray[NBins];

			double ErrorPurityArray[NBins];
			double ErrorEfficiencyArray[NBins];
			double ErrorProductArray[NBins];

			double BeamOn[NBins];
			double ExtBNB[NBins];
			double Dirt[NBins];
			double Overlay[NBins];
			double CC1p[NBins];
			double NonBeamOn[NBins];
			double CosmicFrac[NBins];

			double ErrorBeamOn[NBins];
			double ErrorExtBNB[NBins];
			double ErrorDirt[NBins];
			double ErrorOverlay[NBins];
			double ErrorCC1p[NBins];
			double ErrorNonBeamOn[NBins];
			double ErrorCosmicFrac[NBins];

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
				double BinNonBeamOn = NRecoOverlay + NRecoDirt + NRecoExtBNB;

				double BinErrorCC1pOverlay = TMath::Sqrt(NRecoCC1pOverlay * OverlayPOTScale);
				double BinErrorOverlay = TMath::Sqrt(NRecoOverlay * OverlayPOTScale);
				double BinErrorOverlayDirt = TMath::Sqrt(NRecoDirt * OverlayDirtPOTScale);
				double BinErrorBeamOn = TMath::Sqrt(NRecoBeamOn * BeamOnPOTScale);
				double BinErrorExtBNB = TMath::Sqrt(NRecoExtBNB * ExtBNBPOTScale);
				double BinErrorNonBeamOn = TMath::Sqrt( TMath::Power(BinErrorOverlay,2.) + TMath::Power(BinErrorOverlayDirt,2.) + TMath::Power(BinErrorExtBNB,2.) );

				double purity = NRecoCC1pOverlay / BinNonBeamOn ;
				double efficiency = NRecoCC1pOverlay / NTruthCC1pOverlay;
				double product = purity * efficiency;

				double PurityError = purity * TMath::Sqrt( TMath::Power(BinErrorCC1pOverlay/NRecoCC1pOverlay,2.) + TMath::Power(BinErrorNonBeamOn/BinNonBeamOn,2.) );
				double EfficiencyError = efficiency * TMath::Sqrt( TMath::Power(BinErrorCC1pOverlay/NRecoCC1pOverlay,2.) + TMath::Power(ErrorTrueCC1pOverlay/NTruthCC1pOverlay,2.) );
				double ProductError = TMath::Sqrt( TMath::Power(PurityError*efficiency,2.) + TMath::Power(EfficiencyError*purity,2.) ) ;

				PurityArray[WhichBin] = purity;
				EfficiencyArray[WhichBin] = efficiency;
				ProductArray[WhichBin] = product;

				ErrorPurityArray[WhichBin] = PurityError;
				ErrorEfficiencyArray[WhichBin] = EfficiencyError;
				ErrorProductArray[WhichBin] = ProductError;

				BeamOn[WhichBin] = NRecoBeamOn;
				ExtBNB[WhichBin] = NRecoExtBNB;
				Dirt[WhichBin] = NRecoDirt;
				Overlay[WhichBin] = NRecoOverlay;
				CC1p[WhichBin] = NRecoCC1pOverlay;
				NonBeamOn[WhichBin] = BinNonBeamOn;
				CosmicFrac[WhichBin] = NRecoExtBNB / BinNonBeamOn;

				ErrorBeamOn[WhichBin] = BinErrorBeamOn;
				ErrorExtBNB[WhichBin] = BinErrorExtBNB;
				ErrorDirt[WhichBin] = BinErrorOverlayDirt;
				ErrorOverlay[WhichBin] = BinErrorOverlay;
				ErrorCC1p[WhichBin] = BinErrorCC1pOverlay;
				ErrorNonBeamOn[WhichBin] = BinErrorNonBeamOn;
				ErrorCosmicFrac[WhichBin] = CosmicFrac[WhichBin] * TMath::Sqrt( TMath::Power( BinErrorExtBNB / NRecoExtBNB,2.) + TMath::Power( BinErrorNonBeamOn / BinNonBeamOn,2.) ) ;

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
			latProduct->SetTextSize(TextSize-0.01);
			latProduct->DrawLatexNDC(0.2,0.82,RunNumber[WhichRun] + " Max at "+CutName[WhichCut]+" > "+ToString(SelectedThres));

			// ----------------------------------------------------------------------------------------------------------------------------------------
			
			// Beam On events that we would have with the candidate set of cuts
			
			TString PlotThres = "RecoMuonCosThetaPlot_"+CutName[WhichCut]+"Thres_"+TString(std::to_string(GlobalThresBin));
			TH1D* hRecoBeamOn = (TH1D*)(RecoBeamOnFile->Get(PlotThres));
			
			TLatex* latBeamOn = new TLatex();
			latBeamOn->SetTextFont(TextFont);
			latBeamOn->SetTextSize(TextSize-0.01);
			latBeamOn->DrawLatexNDC(0.2,0.75,"Beam On Events = "+ToString(BeamOn[GlobalThresBin])+" #pm "+ ToString(round(ErrorBeamOn[GlobalThresBin],2)));

			// ----------------------------------------------------------------------------------------------------------------------------------------

			TH1D* hRecoExtBNB = (TH1D*)(RecoExtBNBFile->Get(PlotThres));
			TH1D* hRecoOverlay = (TH1D*)(RecoOverlayFile->Get(PlotThres));
			TH1D* hRecoDirt = (TH1D*)(RecoOverlayDirtFile->Get(PlotThres));

			// ----------------------------------------------------------------------------------------------------------------------------------------
			
			TLatex* latExtBNB = new TLatex();
			latExtBNB->SetTextFont(TextFont);
			latExtBNB->SetTextSize(TextSize-0.01);
			TString CosmicCont = ToString(round(CosmicFrac[GlobalThresBin]*100.,2.) );
			TString ErrorCosmicCont = ToString(round(ErrorCosmicFrac[GlobalThresBin]*100.,2.) );
			latExtBNB->DrawLatexNDC(0.2,0.68,"Cosmics = " + CosmicCont + " #pm " + ErrorCosmicCont + " %");

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
			latPurity->SetTextSize(TextSize-0.01);
			TString Purity = ToString(round(PurityArray[GlobalThresBin]*100.,2) );
			TString ErrorPurity = ToString(round(ErrorPurityArray[GlobalThresBin]*100.,2) );
			latPurity->DrawLatexNDC(0.2,0.82,RunNumber[WhichRun] + " Purity at max = "+ Purity + " #pm " + ErrorPurity + " %");

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
			latEfficiency->SetTextSize(TextSize-0.01);
			TString Efficiency = ToString(round(EfficiencyArray[GlobalThresBin]*100.,2) );
			TString ErrorEfficiency = ToString(round(ErrorEfficiencyArray[GlobalThresBin]*100.,2) );
			latEfficiency->DrawLatexNDC(0.2,0.82,RunNumber[WhichRun] + " Efficiency at max = " + Efficiency + " #pm " + ErrorEfficiency + " %");

			EfficiencyCanvas->SaveAs(PlotPath+CutName[WhichCut]+"_OneDScanEfficiency_"+RunNumber[WhichRun]+".pdf");	
			delete EfficiencyCanvas;

			// ----------------------------------------------------------------------------------------------------------------------------------------

		} // En of the loop over the cuts

	} // En of the loop over the runs

} // End of the program 
