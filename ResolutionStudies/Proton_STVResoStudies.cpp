#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

#include "../../Secondary_Code/myFunctions.cpp"
#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

   // Quadratic background function
   Double_t background(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
   }

   // Lorentzian Peak function
   Double_t lorentzianPeak(Double_t *x, Double_t *par) {
      return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,
      (x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
   }

   // Sum of background and peak function
   Double_t fitFunction(Double_t *x, Double_t *par) {
      return background(x,par) + lorentzianPeak(x,&par[3]);
   }

// ----------------------------------------------------------------------------------------------------------------

void Proton_STVResoStudies() {

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);
	TH1D::SetDefaultSumw2();

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// Selection cuts

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

	/*
	// up to v43
	VectorCuts.push_back("");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");
	*/

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		Cuts = Cuts + VectorCuts[i];

	}

	// ------------------------------------------------------------------------

	TString LocalPathToFiles = PathToFiles + Cuts + "/";
	TString NameOfSamples = "STVStudies_Overlay9_Run1"+Cuts+".root";
	TFile* FileSample = TFile::Open(LocalPathToFiles+NameOfSamples);

	// ------------------------------------------------------------------------

	std::vector<TString> PlotNames;
	std::vector<TString> PlotLabels;

	PlotNames.push_back("Playground_CC1pRecoDeltaPTPlot"); PlotLabels.push_back("All #delta P_{T}");
	PlotNames.push_back("Playground_CC1pRecoDeltaPTPlot_Slice_1"); PlotLabels.push_back("#delta P_{T} < 0.3 GeV/c");
	PlotNames.push_back("Playground_CC1pRecoDeltaPTPlot_Slice_2"); PlotLabels.push_back("0.3 < #delta P_{T} < 0.6 GeV/c");
	PlotNames.push_back("Playground_CC1pRecoDeltaPTPlot_Slice_3"); PlotLabels.push_back("#delta P_{T} > 0.6 GeV/c");

	PlotNames.push_back("Playground_CC1pRecoDeltaAlphaTPlot"); PlotLabels.push_back("All #delta#alpha_{T}");
	PlotNames.push_back("Playground_CC1pRecoDeltaAlphaTPlot_Slice_1"); PlotLabels.push_back("#delta#alpha_{T} < 60 deg");
	PlotNames.push_back("Playground_CC1pRecoDeltaAlphaTPlot_Slice_2"); PlotLabels.push_back("60 < #delta#alpha_{T} < 120 deg");
	PlotNames.push_back("Playground_CC1pRecoDeltaAlphaTPlot_Slice_3"); PlotLabels.push_back("#delta#alpha_{T} > 120 deg");

	PlotNames.push_back("Playground_CC1pRecoDeltaPhiTPlot"); PlotLabels.push_back("All #delta#phi_{T}");
	PlotNames.push_back("Playground_CC1pRecoDeltaPhiTPlot_Slice_1"); PlotLabels.push_back("#delta#phi_{T} < 30 deg");
	PlotNames.push_back("Playground_CC1pRecoDeltaPhiTPlot_Slice_2"); PlotLabels.push_back("30 < #delta#phi_{T} < 60 deg");
	PlotNames.push_back("Playground_CC1pRecoDeltaPhiTPlot_Slice_3"); PlotLabels.push_back("#delta#phi_{T} > 60 deg");

	// ------------------------------------------------------------------------

	const int NPlots = PlotNames.size();

	// ------------------------------------------------------------------------

	std::vector<TString> Discriminator;
	std::vector<TH1D*> Plots;
	std::vector<TString> LegendLabel;

	// ------------------------------------------------------------------------

	Discriminator.push_back("ShortProton"); LegendLabel.push_back("Short");
	Discriminator.push_back("MediumProton"); LegendLabel.push_back("Med");
	Discriminator.push_back("LongProton"); LegendLabel.push_back("Long");

	// ------------------------------------------------------------------------

	int NDiscriminators = Discriminator.size();

	// ------------------------------------------------------------------------

	for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

		// ------------------------------------------------------------------------

		TLegend* leg = new TLegend(0.12,0.91,0.9,0.99);
		leg->SetNColumns(2);
		leg->SetMargin(0.1);

		TCanvas* can = new TCanvas(PlotNames[WhichPlot],PlotNames[WhichPlot],205,34,1024,768);
		can->cd();

		Plots.clear();

		for (int WhichDiscriminator = 0; WhichDiscriminator < NDiscriminators; WhichDiscriminator++) {

			Plots.push_back((TH1D*)FileSample->Get(PlotNames[WhichPlot]+"_"+Discriminator[WhichDiscriminator]));

			for (int i = 0; i < 3; i++) { Plots[WhichDiscriminator]->Rebin(); }

			Plots[WhichDiscriminator]->GetXaxis()->CenterTitle();
			Plots[WhichDiscriminator]->GetXaxis()->SetNdivisions(10);
			Plots[WhichDiscriminator]->GetXaxis()->SetTitleFont(FontStyle);
			Plots[WhichDiscriminator]->GetXaxis()->SetLabelFont(FontStyle);

			Plots[WhichDiscriminator]->GetYaxis()->CenterTitle();
			Plots[WhichDiscriminator]->GetYaxis()->SetNdivisions(10);
			Plots[WhichDiscriminator]->GetYaxis()->SetTitleFont(FontStyle);
			Plots[WhichDiscriminator]->GetYaxis()->SetLabelFont(FontStyle);
			Plots[WhichDiscriminator]->GetYaxis()->SetTitle("Peak Normalized To 1");

			double SF = 1. / Plots[WhichDiscriminator]->GetMaximum();
			Plots[WhichDiscriminator]->Scale(SF);

			Plots[WhichDiscriminator]->GetYaxis()->SetRangeUser(0,1.1*Plots[0]->GetMaximum());

			Plots[WhichDiscriminator]->SetMarkerSize(2.);
			Plots[WhichDiscriminator]->SetMarkerStyle(20);
			Plots[WhichDiscriminator]->SetMarkerColor(WhichDiscriminator+1);
			Plots[WhichDiscriminator]->SetLineColor(WhichDiscriminator+1);
			Plots[WhichDiscriminator]->SetLineWidth(3);
			Plots[WhichDiscriminator]->Draw("p hist same");

			TF1* f = new TF1("f","gaus",-15,15);

			f->SetLineColor(WhichDiscriminator+1);
			Plots[WhichDiscriminator]->Fit(f,"R0Q");
			//f->Draw("same");

			leg->AddEntry(Plots[WhichDiscriminator],LegendLabel[WhichDiscriminator] + ", #mu = " + ToString(round(f->GetParameter(1),1)) + ", #sigma = " + ToString(round(f->GetParameter(2),1)),"p");
			

		} // End of the loop over the discriminators
	
		// ------------------------------------------------------------------------

		leg->SetTextSize(0.04);
		leg->SetTextFont(FontStyle);
		leg->SetBorderSize(0);
		leg->Draw();

		TLatex* lat = new TLatex(0.2,0.4,PlotLabels[WhichPlot]);
		lat->SetTextFont(FontStyle);
		lat->SetTextSize(0.04);
		lat->DrawLatexNDC(0.15,0.8,PlotLabels[WhichPlot]);

		can->SaveAs(PlotPath + "Proton_"+PlotNames[WhichPlot] + ".pdf");

		delete can;

	} // End of the loop over the plots

	// ------------------------------------------------------------------------

} // End of the program
