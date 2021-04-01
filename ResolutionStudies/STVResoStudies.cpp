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

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"

#include "../../myClasses/Constants.h"

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

void STVResoStudies() {

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

	TString PathToFiles = "../OutputFiles/" + UBCodeVersion + "/" + Cuts + "/";
	TString NameOfSamples = "STVStudies_Overlay9_Run1"+Cuts+".root";
	TFile* FileSample = TFile::Open(PathToFiles+NameOfSamples);

	// ------------------------------------------------------------------------

//	TString PlotName = "Playground_CC1pRecoMuonMomentumPlot"; 

	TString PlotName = "Playground_CC1pRecoDeltaPTPlot"; TString Label = "All #delta P_{T}";
//	TString PlotName = "Playground_CC1pRecoDeltaPTPlot_Slice_1"; TString Label = "#delta P_{T} < 0.3 GeV/c";
//	TString PlotName = "Playground_CC1pRecoDeltaPTPlot_Slice_2"; TString Label = "0.3 < #delta P_{T} < 0.6 GeV/c"; 
//	TString PlotName = "Playground_CC1pRecoDeltaPTPlot_Slice_3"; TString Label = "#delta P_{T} > 0.6 GeV/c"; 

//	TString PlotName = "Playground_CC1pRecoDeltaAlphaTPlot"; TString Label = "All #delta#alpha_{T}";
//	TString PlotName = "Playground_CC1pRecoDeltaAlphaTPlot_Slice_1"; TString Label = "#delta#alpha_{T} < 60 deg";
//	TString PlotName = "Playground_CC1pRecoDeltaAlphaTPlot_Slice_2"; TString Label = "60 < #delta#alpha_{T} < 120 deg";
//	TString PlotName = "Playground_CC1pRecoDeltaAlphaTPlot_Slice_3"; TString Label = "#delta#alpha_{T} > 120 deg"; 

//	TString PlotName = "Playground_CC1pRecoDeltaPhiTPlot"; TString Label = "All #delta#phi_{T}";
//	TString PlotName = "Playground_CC1pRecoDeltaPhiTPlot_Slice_1"; TString Label = "#delta#phi_{T} < 30 deg"; 
//	TString PlotName = "Playground_CC1pRecoDeltaPhiTPlot_Slice_2"; TString Label = "30 < #delta#phi_{T} < 60 deg";
//	TString PlotName = "Playground_CC1pRecoDeltaPhiTPlot_Slice_3"; TString Label = "#delta#phi_{T} > 60 deg"; 

	// ------------------------------------------------------------------------

	std::vector<TString> Discriminator;
	std::vector<TH1D*> Plots;
	std::vector<TString> LegendLabel;

	// ------------------------------------------------------------------------

	Discriminator.push_back("FullyContainedMuon"); LegendLabel.push_back("Contained");
	Discriminator.push_back("ExitingShortMuon"); LegendLabel.push_back("Exit Short");
	Discriminator.push_back("ExitingMediumMuon"); LegendLabel.push_back("Exit Med");
	Discriminator.push_back("ExitingLongMuon"); LegendLabel.push_back("Exit Long");

	// ------------------------------------------------------------------------

	int NDiscriminators = Discriminator.size();

	// ------------------------------------------------------------------------

	TCanvas* can = new TCanvas(PlotName,PlotName,205,34,1024,768);
	can->cd();

	// ------------------------------------------------------------------------

	TLegend* leg = new TLegend(0.12,0.91,0.9,0.99);
	leg->SetNColumns(2);
	leg->SetMargin(0.1);

	// ------------------------------------------------------------------------

	for (int WhichDiscriminator = 0; WhichDiscriminator < NDiscriminators; WhichDiscriminator++) {

		Plots.push_back((TH1D*)FileSample->Get(PlotName+"_"+Discriminator[WhichDiscriminator]));

		for (int i = 0; i < 3; i++) { Plots[WhichDiscriminator]->Rebin(); }

		Plots[WhichDiscriminator]->GetXaxis()->CenterTitle();
		Plots[WhichDiscriminator]->GetXaxis()->SetNdivisions(10);
		Plots[WhichDiscriminator]->GetXaxis()->SetTitleFont(FontStyle);
		Plots[WhichDiscriminator]->GetXaxis()->SetLabelFont(FontStyle);

		Plots[WhichDiscriminator]->GetYaxis()->CenterTitle();
		Plots[WhichDiscriminator]->GetYaxis()->SetNdivisions(10);
		Plots[WhichDiscriminator]->GetYaxis()->SetTitleFont(FontStyle);
		Plots[WhichDiscriminator]->GetYaxis()->SetLabelFont(FontStyle);
		Plots[WhichDiscriminator]->GetYaxis()->SetTitle("Arbitrary Units");

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
		

	}	
	// ------------------------------------------------------------------------

	leg->SetTextSize(0.04);
	leg->SetTextFont(FontStyle);
	leg->SetBorderSize(0);
	leg->Draw();

	TLatex* lat = new TLatex(0.2,0.4,Label);
	lat->SetTextFont(FontStyle);
	lat->SetTextSize(0.04);
	lat->DrawLatexNDC(0.15,0.8,Label);


	can->SaveAs("myResoPlots/"+UBCodeVersion+"/"+PlotName+".pdf");

	// ------------------------------------------------------------------------

} // End of the program
