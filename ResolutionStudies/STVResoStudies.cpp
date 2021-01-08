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

// ----------------------------------------------------------------------------------------------------------------

void STVResoStudies() {

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);
	TH1D::SetDefaultSumw2();

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// Selection cuts

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		Cuts = Cuts + VectorCuts[i];

	}

	// ------------------------------------------------------------------------

	TString PathToFiles = "../OutputFiles/" + UBCodeVersion + "/" + Cuts + "/";
	TString NameOfSamples = "STVStudies_Overlay9_Run1"+Cuts+".root";
	TFile* FileSample = TFile::Open(PathToFiles+NameOfSamples);

	// ------------------------------------------------------------------------

//	TString PlotName = "Playground_CC1pRecoDeltaPTPlot"; 
//	TString PlotName = "Playground_CC1pRecoDeltaAlphaTPlot"; 
	TString PlotName = "Playground_CC1pRecoDeltaPhiTPlot"; 

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

	TLegend* leg = new TLegend(0.5,0.7,0.8,0.85);

	// ------------------------------------------------------------------------

	for (int WhichDiscriminator = 0; WhichDiscriminator < NDiscriminators; WhichDiscriminator++) {

		Plots.push_back((TH1D*)FileSample->Get(PlotName+"_"+Discriminator[WhichDiscriminator]));

		Plots[WhichDiscriminator]->GetXaxis()->CenterTitle();
		Plots[WhichDiscriminator]->GetXaxis()->SetNdivisions(10);
		Plots[WhichDiscriminator]->GetXaxis()->SetTitleFont(FontStyle);
		Plots[WhichDiscriminator]->GetXaxis()->SetLabelFont(FontStyle);

		Plots[WhichDiscriminator]->GetYaxis()->CenterTitle();
		Plots[WhichDiscriminator]->GetYaxis()->SetNdivisions(10);
		Plots[WhichDiscriminator]->GetYaxis()->SetTitleFont(FontStyle);
		Plots[WhichDiscriminator]->GetYaxis()->SetLabelFont(FontStyle);
		Plots[WhichDiscriminator]->GetYaxis()->SetTitle("# Events");

		Plots[WhichDiscriminator]->SetLineColor(WhichDiscriminator+1);
		Plots[WhichDiscriminator]->SetLineWidth(3);
		Plots[WhichDiscriminator]->Draw("e same");

		TF1* f = new TF1("f","gaus",-100,100);
		f->SetLineColor(WhichDiscriminator+1);
		Plots[WhichDiscriminator]->Fit(f,"","",-15,15);
		f->Draw("same");

		leg->AddEntry(Plots[WhichDiscriminator],LegendLabel[WhichDiscriminator] + ", #mu = " + ToString(round(f->GetParameter(1),1)) + ", #sigma = " + ToString(round(f->GetParameter(2),1)),"l");
		

	}	
	// ------------------------------------------------------------------------

	leg->SetTextSize(0.04);
	leg->SetTextFont(FontStyle);
	leg->SetBorderSize(0);
	leg->Draw();

	// ------------------------------------------------------------------------

} // End of the program
