#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TLatex.h>

#include <iostream>
#include <vector>

#include "../myClasses/Constants.h"
#include "../myClasses/myFunctions.cpp"

using namespace std;
using namespace Constants;

//----------------------------------------//

void print_1d() {

	//--------------------------------------//

	int font = 132;
	double size = 0.07;

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	gStyle->SetTitleSize(size,"xyzt");

	//--------------------------------------//

	vector<TString> plot_names; plot_names.clear();
	
	plot_names.push_back("CC1pRecoThetaBRTPlot");

	const int nplots = plot_names.size();

	//--------------------------------------//

	TString cut = "_NoCuts_PID_NuScore_CRT";
	TString run = "Combined";
	TString file_path = PathToFiles + "/"+cut+"/STVStudies_Overlay9_"+run+cut+".root";
	TFile* file = new TFile(file_path,"readonly");
	int colors = OverlayColor;

	//--------------------------------------//

	vector<TCanvas*> canvas; canvas.clear(); canvas.resize(nplots);
	vector<TH1D*>  plot; plot.clear(); plot.resize(nplots);

	//--------------------------------------//

	// Loop over the plots
	
	for (int iplot = 0; iplot < nplots; iplot ++) {

		TString canvas_name = "canvas_"+plot_names[iplot];
		canvas.at(iplot) = new TCanvas(canvas_name,canvas_name,205,34,1024,768);
		canvas.at(iplot)->SetBottomMargin(0.17);
		canvas.at(iplot)->SetTopMargin(0.18);
		canvas.at(iplot)->SetLeftMargin(0.15);
		
		plot.at(iplot) = (TH1D*)( file->Get( plot_names.at(iplot) ) );

		plot.at(iplot)->GetXaxis()->CenterTitle();
                plot.at(iplot)->GetXaxis()->SetTitleSize(size);
                plot.at(iplot)->GetXaxis()->SetLabelSize(size);
                plot.at(iplot)->GetXaxis()->SetTitleFont(font);
                plot.at(iplot)->GetXaxis()->SetLabelFont(font);
                plot.at(iplot)->GetXaxis()->SetTitleOffset(1.);
                plot.at(iplot)->GetXaxis()->SetNdivisions(8);

		plot.at(iplot)->GetYaxis()->CenterTitle();
                plot.at(iplot)->GetYaxis()->SetTitleSize(size);
                plot.at(iplot)->GetYaxis()->SetLabelSize(size);
                plot.at(iplot)->GetYaxis()->SetTitleFont(font);
                plot.at(iplot)->GetYaxis()->SetLabelFont(font);
                plot.at(iplot)->GetYaxis()->SetTitleOffset(1.1);
                plot.at(iplot)->GetYaxis()->SetNdivisions(6);

		plot.at(iplot)->SetLineColor(colors);
		plot.at(iplot)->SetLineWidth(3);
		plot.at(iplot)->Draw("hist");

		gPad->RedrawAxis();

		if (plot_names[iplot] == "CC1pRecoThetaBRTPlot") {

                	plot.at(iplot)->GetYaxis()->SetTitle("CC1p0#pi event count");
			
			double mean = plot.at(iplot)->GetMean();
			double median = GetMedian(plot.at(iplot));
			double sigma = plot.at(iplot)->GetRMS();
			double peak = FindOneDimHistoMaxValueBin(plot.at(iplot));
	
			TLatex latex;
			latex.SetTextFont(FontStyle);
			latex.SetTextSize(0.07);
			TString label = "#splitline{peak = " + to_string_with_precision(peak,2) + "^{o}, median = " + to_string_with_precision(median,2)  + "^{o}}{#mu = " + to_string_with_precision(mean,2) + "^{o}, #sigma = " + to_string_with_precision(sigma,2) + "^{o}}";
			latex.DrawLatexNDC(0.3,0.7, label);				

		}

	        TLatex *bnb = new TLatex();
        	bnb->SetTextFont(FontStyle);
        	bnb->SetTextSize(TextSize);
        	bnb->DrawLatexNDC(0.15,0.85,"MicroBooNE Simulation");

                TString canvas_path = PlotPath + cut + "/";
                TString canvas_export_name = "print_1d_"+plot_names.at(iplot)+"_"+run+"_"+UBCodeVersion+cut+".pdf";
                canvas.at(iplot)->SaveAs(canvas_path + canvas_export_name);
                delete canvas.at(iplot);

	} // End of the loop over the plots

} // End of the program 
