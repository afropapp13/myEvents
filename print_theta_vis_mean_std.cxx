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

void  print_theta_vis_mean_std() {

	//--------------------------------------//

	TH1D::SetDefaultSumw2();

	//--------------------------------------//

	vector<TString> plot_names; plot_names.clear();
	vector<TString> plot_label; plot_label.clear();
	
	plot_names.push_back("CC1pFineBinThetaVisPlot"); plot_label.push_back("All events");
	
	plot_names.push_back(""); plot_label.push_back("");
	
	plot_names.push_back("CC1pFineBinThetaVis_ECalSlices0_00To0_50Plot"); plot_label.push_back("$E_{\\mathrm{reco}} <$ 0.5\\,GeV");
	plot_names.push_back("CC1pFineBinThetaVis_ECalSlices0_50To0_80Plot"); plot_label.push_back("0.5 $ < E_{\\mathrm{reco}} <$ 0.8\\,GeV");
	plot_names.push_back("CC1pFineBinThetaVis_ECalSlices0_80To2_00Plot"); plot_label.push_back("$E_{\\mathrm{reco}} >$ 0.8\\,GeV");

	plot_names.push_back(""); plot_label.push_back("");
	
	plot_names.push_back("CC1pFineBinThetaVis_DeltaPnSlices0_00To0_20Plot"); plot_label.push_back("$p_{n} <$ 0.2\\,GeV/c");
	plot_names.push_back("CC1pFineBinThetaVis_DeltaPnSlices0_20To0_40Plot"); plot_label.push_back("0.2 $< p_{n} <$ 0.4\\,GeV/c");
	plot_names.push_back("CC1pFineBinThetaVis_DeltaPnSlices0_40To1_00Plot"); plot_label.push_back("$p_{n} >$ 0.4\\,GeV/c");
	
	plot_names.push_back(""); plot_label.push_back("");
	
	plot_names.push_back("CC1pThetaVis_PMissSlices0_00To0_12Plot"); plot_label.push_back("$|p_{\\mathrm{miss}}| <$ 0.12\\,GeV/c");
	plot_names.push_back("CC1pThetaVis_PMissSlices0_12To0_50Plot"); plot_label.push_back("$|p_{\\mathrm{miss}}| >$ 0.12\\,GeV/c");

	
	const int nplots = plot_names.size();

	//--------------------------------------//

	TString cut = "_NoCuts_PID_NuScore_CRT";
	TString run = "Combined";
	TString file_path = PathToFiles + "/"+cut+"/STVStudies_Overlay9_"+run+cut+".root";
	TFile* file = new TFile(file_path,"readonly");
	int colors = OverlayColor;

	//--------------------------------------//

	vector<TH1D*>  plot; plot.clear(); plot.resize(nplots);

	//--------------------------------------//

	// Loop over the plots
	
	for (int iplot = 0; iplot < nplots; iplot ++) {

		if (plot_names.at(iplot) == "") { 

			cout << " \\hline";

		} else {

			plot.at(iplot) = (TH1D*)( file->Get( plot_names.at(iplot) ) );

			double mean = plot.at(iplot)->GetMean();
			double sigma = plot.at(iplot)->GetRMS();
			double median = Median( plot.at(iplot) );
	
			//cout << plot_names.at(iplot) << " mean = " << mean << "  sigma = " << sigma << "  median = " << median << endl;
			cout << " " << plot_label.at(iplot) << " & " << to_string_with_precision(mean,1) << "  &  " << to_string_with_precision(sigma,1) << " & " << to_string_with_precision(median,1) << " \\tabularnewline \\hline" << endl;
		
		}

	} // End of the loop over the plots

} // End of the program 
