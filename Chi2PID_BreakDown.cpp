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

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

#include "../myCCQEAnalysis/Constants.h"

using namespace std;
using namespace Constants;

void Chi2PID_BreakDown() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames; PlotNames.clear();

	TString StorePath = "myPlots/";
	TString PathToFiles = "OutputFiles/"+UBCodeVersion+"/";

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	PlotNames.push_back("RecoChi2Plot");
	PlotNames.push_back("RecoThreePlaneChi2Plot");
	PlotNames.push_back("RecoThreePlaneChi2LogLikelihoodPlot");

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
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

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
	vector<THStack*> THStacks; THStacks.clear();
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();
	vector<TLegend*> leg; leg.clear();

	vector<vector<TH1D*> > Plots; Plots.clear();
	vector<vector<TH1D*> > MuonPlots; MuonPlots.clear();
	vector<vector<TH1D*> > ProtonPlots; ProtonPlots.clear();
	vector<vector<TH1D*> > PionPlots; PionPlots.clear();
	vector<vector<TH1D*> > CosmicPlots; CosmicPlots.clear();

	vector<vector<TH1D*> > hratio;  hratio.clear();

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;
	
	NameOfSamples.push_back("CCQEStudies_Run1Data9"+Cuts+".root"); LabelsOfSamples.push_back("Data Run-1");
	NameOfSamples.push_back("CCQEStudies_Overlay9"+Cuts+".root"); LabelsOfSamples.push_back("Overlay");
	NameOfSamples.push_back("CCQEStudies_ExtBNB9"+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");
	NameOfSamples.push_back("CCQEStudies_OverlayDirt9"+Cuts+".root"); LabelsOfSamples.push_back("Dirt");

	vector<int> Colors; Colors.clear(); 
	Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kBlue); Colors.push_back(kMagenta);

	vector<int> ColorsOverlay; ColorsOverlay.clear(); 
	ColorsOverlay.push_back(kRed); ColorsOverlay.push_back(kCyan); ColorsOverlay.push_back(kGreen); ColorsOverlay.push_back(kMagenta); ColorsOverlay.push_back(kOrange+7);

	const int NSamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		FileSample.push_back(TFile::Open(PathToFiles+NameOfSamples[WhichSample]));

		vector<TH1D*> CurrentPlots; CurrentPlots.clear();
		vector<TH1D*> MuonCurrentPlots; MuonCurrentPlots.clear();
		vector<TH1D*> ProtonCurrentPlots; ProtonCurrentPlots.clear();
		vector<TH1D*> PionCurrentPlots; PionCurrentPlots.clear();
		vector<TH1D*> CosmicCurrentPlots; CosmicCurrentPlots.clear();

		vector<TH1D*> Currenthratio;  Currenthratio.clear();

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
			TH1D* Muonhist = (TH1D*)(FileSample[WhichSample]->Get("Muon"+PlotNames[WhichPlot]));
			TH1D* Protonhist = (TH1D*)(FileSample[WhichSample]->Get("Proton"+PlotNames[WhichPlot]));
			TH1D* Pionhist = (TH1D*)(FileSample[WhichSample]->Get("Pion"+PlotNames[WhichPlot]));
			TH1D* Cosmichist = (TH1D*)(FileSample[WhichSample]->Get("Cosmic"+PlotNames[WhichPlot]));

			CenterAxisTitle(hist);

			hist->SetLineColor(Colors[WhichSample]);
			if (LabelsOfSamples[WhichSample] == "Data Run-1") { hist->SetMarkerStyle(20); }

			CurrentPlots.push_back(hist);
			MuonCurrentPlots.push_back(Muonhist);
			ProtonCurrentPlots.push_back(Protonhist);
			PionCurrentPlots.push_back(Pionhist);
			CosmicCurrentPlots.push_back(Cosmichist);

			Currenthratio.push_back((TH1D*)hist->Clone());
		
		}

		Plots.push_back(CurrentPlots);
		MuonPlots.push_back(MuonCurrentPlots);
		ProtonPlots.push_back(ProtonCurrentPlots);
		PionPlots.push_back(PionCurrentPlots);
		CosmicPlots.push_back(CosmicCurrentPlots);

		hratio.push_back(Currenthratio);

	}

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
		PlotCanvas.push_back(new TCanvas(PlotNames[WhichPlot],PlotNames[WhichPlot],205,34,1024,768));
		PlotCanvas[WhichPlot]->cd();

		THStacks.push_back(new THStack(PlotNames[WhichPlot],""));

		TPad *topPad = new TPad("topPad", "", 0.005, 0.92, 0.995, 0.995);
		TPad *midPad = new TPad("midPad", "", 0.005, 0.3  , 0.995, 0.92);
		TPad *botPad = new TPad("botPad", "", 0.005, 0.005, 0.995, 0.3);
		topPad->SetTopMargin(0.3);
		topPad->SetBottomMargin(0.0);
		midPad->SetBottomMargin(0.03);
		midPad->SetTopMargin(0.0);
		botPad->SetTopMargin(0.);
		botPad->SetBottomMargin(0.25);
		botPad->SetGridy();
		topPad->Draw();
		midPad->Draw();
		botPad->Draw();

		leg.push_back(new TLegend(0.1,0.005,0.9,0.995));
		leg[WhichPlot]->SetBorderSize(0);
		leg[WhichPlot]->SetNColumns(3);

		double max = -99.;

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

			midPad->cd();
			Plots[WhichSample][WhichPlot]->SetTitle("");
			Plots[WhichSample][WhichPlot]->SetLineWidth(4);

			Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(5);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0);

			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(5);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# events");
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

			double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
			if (localmax > max) { max = localmax; }
			Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.2*max);

			if (LabelsOfSamples[WhichSample] == "Data Run-1") { 

				Plots[WhichSample][WhichPlot]->Draw("e1 same"); 
				leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample],"lep");

			}

			if (LabelsOfSamples[WhichSample] == "ExtBNB") {

					Plots[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
					Plots[WhichSample][WhichPlot]->SetFillColor(Colors[WhichSample]);
					THStacks[WhichPlot]->Add(Plots[WhichSample][WhichPlot],"hist");
					leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample],"f");
					THStacks[WhichPlot]->Draw("same");

			}

			if (LabelsOfSamples[WhichSample] == "Overlay") {

					MuonPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
					MuonPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
					MuonPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[0]);
					MuonPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[0]);
					THStacks[WhichPlot]->Add(MuonPlots[WhichSample][WhichPlot],"hist");
					THStacks[WhichPlot]->Add(MuonPlots[WhichSample+2][WhichPlot],"hist");

					ProtonPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
					ProtonPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
					ProtonPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[1]);
					ProtonPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[1]);
					THStacks[WhichPlot]->Add(ProtonPlots[WhichSample][WhichPlot],"hist");
					THStacks[WhichPlot]->Add(ProtonPlots[WhichSample+2][WhichPlot],"hist");

					PionPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
					PionPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
					PionPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[2]);
					PionPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[2]);
					THStacks[WhichPlot]->Add(PionPlots[WhichSample][WhichPlot],"hist");
					THStacks[WhichPlot]->Add(PionPlots[WhichSample+2][WhichPlot],"hist");

					CosmicPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
					CosmicPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
					CosmicPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[3]);
					CosmicPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[3]);
					THStacks[WhichPlot]->Add(CosmicPlots[WhichSample][WhichPlot],"hist");
					THStacks[WhichPlot]->Add(CosmicPlots[WhichSample+2][WhichPlot],"hist");

					leg[WhichPlot]->AddEntry(MuonPlots[WhichSample][WhichPlot],"#mu","f"); 
					leg[WhichPlot]->AddEntry(ProtonPlots[WhichSample][WhichPlot],"p","f"); 
					leg[WhichPlot]->AddEntry(PionPlots[WhichSample][WhichPlot],"#pi^{#pm}","f"); 
					leg[WhichPlot]->AddEntry(CosmicPlots[WhichSample][WhichPlot],"Cosmic","f"); 

					THStacks[WhichPlot]->Draw("same");
					
			}
				

		} // End of the loop over the samples

		Plots[0][WhichPlot]->Draw("e1 same");

//		hratio[1][WhichPlot]->Add(hratio[2][WhichPlot]);
//		hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
//		hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

		hratio[0][WhichPlot]->Add(hratio[2][WhichPlot],-1);
		hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
		hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

		hratio[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
		hratio[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
//		hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{Data}{Overlay + ExtBNB}");
		hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{Data - ExtBNB}{Overlay}");
		hratio[0][WhichPlot]->GetXaxis()->SetTitle(Plots[0][WhichPlot]->GetXaxis()->GetTitle());
		hratio[0][WhichPlot]->GetXaxis()->SetTitleSize(0.13);
		hratio[0][WhichPlot]->GetXaxis()->SetLabelSize(0.12);
		hratio[0][WhichPlot]->GetXaxis()->SetTitleOffset(0.88);
		hratio[0][WhichPlot]->GetXaxis()->SetNdivisions(5);

		hratio[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
		hratio[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
		hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.9*hratio[0][WhichPlot]->GetMinimum(),1.1*hratio[0][WhichPlot]->GetMaximum());
		hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(4);
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
		RatioLine->Draw("same");
		
		topPad->cd();
		leg[WhichPlot]->SetTextSize(0.5);
		leg[WhichPlot]->SetTextFont(FontStyle);
		leg[WhichPlot]->Draw();

		PlotCanvas[WhichPlot]->SaveAs("./myPlots/pdf/1D/"+UBCodeVersion+"/PID_BreakDown_"+PlotNames[WhichPlot]+Cuts+"_"+UBCodeVersion+".pdf");
		PlotCanvas[WhichPlot]->SaveAs("./myPlots/eps/1D/"+UBCodeVersion+"/PID_BreakDown_"+PlotNames[WhichPlot]+Cuts+"_"+UBCodeVersion+".eps");
		//delete PlotCanvas[WhichPlot];

	} // End of the loop over the plots

} // End of the program 
