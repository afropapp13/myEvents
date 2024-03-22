#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <THStack.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

using namespace std;

double TextSize = 0.08;
double FontStyle = 132;

//int DUNEColor = kBlue-9;
//int T2KColor = kOrange-3;
//int NOvAColor = kGreen+1;

int DUNEColor = kBlue-9;
int T2KColor = kOrange+1;
int NOvAColor = kOrange+5;

// -------------------------------------------------------------------------------------------------------------------------------------

void PrettyPlot(TH1D* h, int Color,bool shift) {

	h->SetLineColor(Color);
	h->SetLineWidth(3);
	h->Scale(1./h->GetMaximum());

	h->GetXaxis()->SetRangeUser(0.01,2.4);
	
	h->GetXaxis()->SetNdivisions(8);
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleOffset(0.9);
	h->GetXaxis()->SetTickSize(0.02);
	h->GetXaxis()->SetTitle("E_{#nu} [GeV]");

	h->GetYaxis()->SetRangeUser(0.,1.1);
	h->GetYaxis()->SetNdivisions(8);
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleSize(TextSize);
	h->GetYaxis()->SetLabelSize(0.);
	h->GetYaxis()->SetTickSize(0.);
	h->GetYaxis()->SetTitleOffset(0.5);
	h->GetYaxis()->SetTitle("#nu Flux [arb.]");

}

// -------------------------------------------------------------------------------------------------------------------------------------

void fluxes() {

	// -------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);
	TGaxis::SetMaxDigits(3);

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t"); gStyle->SetTitleFont(FontStyle,"t");

	// --------------------------------------------------------------------------------------------------------------------------------	

	TString CanvasName = "ThetaVis_FluxCanvas";
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->SetBottomMargin(0.15);

	// ---------------------------------------------------------------------------------------------------------------------------

	TFile* bnb_file = TFile::Open("../mySTVAnalysis/MCC9_FluxHist_volTPCActive.root");
	TH1D* bnb_flux = (TH1D*)( bnb_file->Get("hEnumu_cv") );

	TFile* honda_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/Fluxes/histogram.root");
	TH1D* honda_flux = (TH1D*)(honda_file->Get("h"));

	// --------------------------------------------------------------------------------------------------

	PrettyPlot(bnb_flux,DUNEColor,false);
	PrettyPlot(honda_flux,T2KColor,false);
							
	// --------------------------------------------------------------------------------------------------

	PlotCanvas->cd();

	bnb_flux->Draw("c hist same");
	honda_flux->Draw("c hist same");

	// --------------------------------------------------------------------------------------------------

	TLatex *bnb = new TLatex(); 
	bnb->SetTextFont(FontStyle); 
	bnb->SetTextColor(DUNEColor); 
	bnb->SetTextSize(TextSize);
	bnb->DrawLatexNDC(0.7,0.82,"BNB");

	TLatex *Honda = new TLatex(); 
	Honda->SetTextFont(FontStyle); 
	Honda->SetTextColor(T2KColor); 
	Honda->SetTextSize(TextSize);
	Honda->DrawLatexNDC(0.7,0.75,"Honda");

	PlotCanvas->SaveAs("/exp/uboone/data/users/apapadop/PeLEETuples_Atmospherics/FlatTTreePlots/"+CanvasName+".pdf");
	delete PlotCanvas;

	//---------------------//

} // End of the program
