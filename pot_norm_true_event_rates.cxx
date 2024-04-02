#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TProfile.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include "../myClasses/Tools.h"
#include "../myClasses/Constants.h"
#include "../myClasses/myFunctions.cpp"

using namespace Constants;
using namespace std;

void myPrettyPlot(TH1D* h,int LineWidth = 2, int FontStyle = 132, int Ndivisions = 6, double TextSize = 0.06) {

	// ----------------------------------------------------------------------------------------------------------------

	h->SetLineWidth(LineWidth);

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	// ----------------------------------------------------------------------------------------------------------------

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetNdivisions(Ndivisions);
	h->GetYaxis()->SetTitle("True POT Normalized Events");

	return;	

}

void pot_norm_true_event_rates() {

	// ----------------------------------------------------------------------------------------------------------------------------------------

	int LineWidth = 3;
	int FontStyle = 132;
	int Ndivisions = 6; 
	double TextSize = 0.06; 

	gStyle->SetOptStat(0);	

	// ----------------------------------------------------------------------------------------------------------------------------------------

	const std::vector<int> Colors{kOrange+7,kAzure-4,kBlack,610,410,kRed+1,kGreen+3,kBlue};

	// -----------------------------------------------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();
	vector<TH1D*> Plots; Plots.clear();
	vector<double> pot; pot.clear();

	// -----------------------------------------------------------------------------------------------------------------------------

	TString PlotName = "TrueMuonCosThetaPlot";

	// -----------------------------------------------------------------------------------------------------------------------------

	FileNames.push_back("Overlay9_Run1"); Label.push_back("Run1"); pot.push_back(Fulltor860_wcut_Run1); 
	FileNames.push_back("Overlay9_Run2"); Label.push_back("Run2"); pot.push_back(Fulltor860_wcut_Run2);
	FileNames.push_back("Overlay9_Run3"); Label.push_back("Run3"); pot.push_back(Fulltor860_wcut_Run3);
	FileNames.push_back("Overlay9_Run4b"); Label.push_back("Run4b"); pot.push_back(Fulltor860_wcut_Run4b);
	FileNames.push_back("Overlay9_Run4c"); Label.push_back("Run4c"); pot.push_back(Fulltor860_wcut_Run4c);
	FileNames.push_back("Overlay9_Run4d"); Label.push_back("Run4d"); pot.push_back(Fulltor860_wcut_Run4d);
	FileNames.push_back("Overlay9_Run5"); Label.push_back("Run5"); pot.push_back(Fulltor860_wcut_Run5);

	const int NFiles = FileNames.size();

	// -----------------------------------------------------------------------------------------------------------------------------

	TString CanvasName = "OverlayCanvas";
	TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	can->SetBottomMargin(0.17);
	can->SetLeftMargin(0.15);

	// -----------------------------------------------------------------------------------------------------------------------------

	TLegend* leg = new TLegend(0.2,0.5,0.4,0.8);
	leg->SetNColumns(1);

	// -----------------------------------------------------------------------------------------------------------------------------

	for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {

		TFile* f = TFile::Open(PathToFiles+"/TruthSTVAnalysis_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root");
		TH1D* h = (TH1D*)(f->Get(PlotName));
		Plots.push_back(h);
			
		Plots[WhichFile]->SetLineColor(Colors[WhichFile]);
		myPrettyPlot(Plots[WhichFile]);
		Plots[WhichFile]->Scale(pot[0]/pot[WhichFile]);
		Plots[WhichFile]->Draw("e same");

		leg->AddEntry(Plots[WhichFile],Label[WhichFile] + " " + ToString(pot[WhichFile]),"l");

	}

	// -----------------------------------------------------------------------------------------------------------------------------

	leg->SetBorderSize(0);
	leg->SetTextSize(TextSize);
	leg->SetTextFont(FontStyle);
	leg->Draw();

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program
