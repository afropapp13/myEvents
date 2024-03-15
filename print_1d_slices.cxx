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

using namespace std;
using namespace Constants;

//----------------------------------------//

TString to_string_with_precision(double value, const int n = 2) {

    std::ostringstream out;
    out.precision(n);
    out << std::fixed << value;
    return TString(out.str());

}

//----------------------------------------//

TString convert_to_string(double value) {

	TString StringValue = to_string_with_precision(value, 2);
	StringValue.ReplaceAll(".","_");
	StringValue.ReplaceAll("-","Minus");	

	return StringValue;

}

//----------------------------------------//

void print_1d_slices() {

	//--------------------------------------//

	int font = 132;
	double size = 0.07;

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	gStyle->SetTitleSize(size,"xyzt");

	//--------------------------------------//

	vector<TString> plot_names; plot_names.clear();
	vector< vector<double> > slices; slices.clear();
	vector< TString > variable; variable.clear();

	//--------------------------------------//

	// ECal
	
	plot_names.push_back("CC1pECalDiff_MuonMomentumSlices"); slices.push_back(TwoDArrayNBinsMuonMomentum); variable.push_back("P_{#mu} [GeV/c]");
	plot_names.push_back("CC1pECalReso_MuonMomentumSlices"); slices.push_back(TwoDArrayNBinsMuonMomentum); variable.push_back("P_{#mu} [GeV/c]");
	
	plot_names.push_back("CC1pECalDiff_MuonCosThetaSlices"); slices.push_back(TwoDArrayNBinsMuonCosTheta); variable.push_back("cos#theta_{#mu}");
	plot_names.push_back("CC1pECalReso_MuonCosThetaSlices"); slices.push_back(TwoDArrayNBinsMuonCosTheta); variable.push_back("cos#theta_{#mu}");

	plot_names.push_back("CC1pECalDiff_ProtonMomentumSlices"); slices.push_back(TwoDArrayNBinsProtonMomentum); variable.push_back("P_{p} [GeV/c]");
	plot_names.push_back("CC1pECalReso_ProtonMomentumSlices"); slices.push_back(TwoDArrayNBinsProtonMomentum); variable.push_back("P_{p} [GeV/c]");

	plot_names.push_back("CC1pECalDiff_ProtonCosThetaSlices"); slices.push_back(TwoDArrayNBinsProtonCosTheta); variable.push_back("cos#theta_{p}");
	plot_names.push_back("CC1pECalReso_ProtonCosThetaSlices"); slices.push_back(TwoDArrayNBinsProtonCosTheta); variable.push_back("cos#theta_{p}");

	plot_names.push_back("CC1pECalDiff_DeltaPTSlices"); slices.push_back(TwoDArrayNBinsDeltaPT); variable.push_back("#deltap_{T} [GeV/c]");
	plot_names.push_back("CC1pECalReso_DeltaPTSlices"); slices.push_back(TwoDArrayNBinsDeltaPT); variable.push_back("#deltap_{T} [GeV/c]");

	plot_names.push_back("CC1pECalDiff_DeltaPnSlices"); slices.push_back(TwoDArrayNBinsDeltaPn); variable.push_back("p_{n} [GeV/c]");
	plot_names.push_back("CC1pECalReso_DeltaPnSlices"); slices.push_back(TwoDArrayNBinsDeltaPn); variable.push_back("p_{n} [GeV/c]");

	plot_names.push_back("CC1pECalDiff_DeltaAlphaTSlices"); slices.push_back(TwoDArrayNBinsDeltaAlphaT); variable.push_back("#delta#alpha_{T} [deg]");
	plot_names.push_back("CC1pECalReso_DeltaAlphaTSlices"); slices.push_back(TwoDArrayNBinsDeltaAlphaT); variable.push_back("#delta#alpha_{T} [deg]");

	plot_names.push_back("CC1pECalDiff_DeltaAlpha3DSlices"); slices.push_back(TwoDArrayNBinsDeltaAlpha3D); variable.push_back("#alpha_{3D} [deg]");
	plot_names.push_back("CC1pECalReso_DeltaAlpha3DSlices"); slices.push_back(TwoDArrayNBinsDeltaAlpha3D); variable.push_back("#alpha_{3D} [deg]");
	
	//--------------------------------------//

	// ThetaVis

	plot_names.push_back("CC1pThetaVis_DeltaPnSlices"); slices.push_back(TwoDArrayNBinsDeltaPn); variable.push_back("p_{n}^{reco} [GeV/c]");	
	plot_names.push_back("CC1pThetaVis_DeltaPTSlices"); slices.push_back(TwoDArrayNBinsDeltaPT); variable.push_back("#deltap_{T}^{reco} [GeV/c]");	
	
	plot_names.push_back("CC1pThetaVis_ECalSlices"); slices.push_back(TwoDArrayNBinsECal); variable.push_back("E_{Cal}^{reco} [GeV]");	
	plot_names.push_back("CC1pTrueThetaVis_TrueECalSlices"); slices.push_back(TwoDArrayNBinsECal); variable.push_back("E_{Cal}^{true} [GeV]");	
	plot_names.push_back("CC1pTrueThetaVis_TrueEnuSlices"); slices.push_back(TwoDArrayNBinsECal); variable.push_back("E_{#nu}^{true} [GeV]");	
	plot_names.push_back("CC1pThetaVisDiff_ECalSlices"); slices.push_back(TwoDArrayNBinsECal); variable.push_back("E_{Cal} [GeV]");
	plot_names.push_back("CC1pThetaVisReso_ECalSlices"); slices.push_back(TwoDArrayNBinsECal); variable.push_back("E_{Cal} [GeV]");

	plot_names.push_back("CC1pThetaVisDiff_MuonMomentumSlices"); slices.push_back(TwoDArrayNBinsMuonMomentum); variable.push_back("p_{#mu} [GeV/c]");
	plot_names.push_back("CC1pThetaVisReso_MuonMomentumSlices"); slices.push_back(TwoDArrayNBinsMuonMomentum); variable.push_back("p_{#mu} [GeV/c]");

	plot_names.push_back("CC1pThetaVisDiff_ProtonMomentumSlices"); slices.push_back(TwoDArrayNBinsProtonMomentum); variable.push_back("p_{p} [GeV/c]");
	plot_names.push_back("CC1pThetaVisReso_ProtonMomentumSlices"); slices.push_back(TwoDArrayNBinsProtonMomentum); variable.push_back("p_{p} [GeV/c]");
	
	plot_names.push_back("CC1pThetaVisDiff_DeltaPTSlices"); slices.push_back(TwoDArrayNBinsDeltaPT); variable.push_back("#deltap_{T} [GeV/c]");
	plot_names.push_back("CC1pThetaVisReso_DeltaPTSlices"); slices.push_back(TwoDArrayNBinsDeltaPT); variable.push_back("#deltap_{T} [GeV/c]");
	
	plot_names.push_back("CC1pThetaVisDiff_DeltaAlphaTSlices"); slices.push_back(TwoDArrayNBinsDeltaAlphaT); variable.push_back("#delta#alpha_{T} [deg]");
	plot_names.push_back("CC1pThetaVisReso_DeltaAlphaTSlices"); slices.push_back(TwoDArrayNBinsDeltaAlphaT); variable.push_back("#delta#alpha_{T} [deg]");

	plot_names.push_back("CC1pThetaVisDiff_DeltaPnSlices"); slices.push_back(TwoDArrayNBinsDeltaPn); variable.push_back("p_{n} [GeV/c]");
	plot_names.push_back("CC1pThetaVisReso_DeltaPnSlices"); slices.push_back(TwoDArrayNBinsDeltaPn); variable.push_back("p_{n} [GeV/c]");
	
	plot_names.push_back("CC1pThetaVisDiff_DeltaAlpha3DSlices"); slices.push_back(TwoDArrayNBinsDeltaAlpha3D); variable.push_back("#alpha_{3D} [deg]");
	plot_names.push_back("CC1pThetaVisReso_DeltaAlpha3DSlices"); slices.push_back(TwoDArrayNBinsDeltaAlpha3D); variable.push_back("#alpha_{3D} [deg]");


	const int nplots = plot_names.size();
	const int nslices = slices.size();

	//--------------------------------------//

	TString cut = "_NoCuts_PID_NuScore_CRT";
	TString run = "Combined";
	TString file_path = PathToFiles + "/"+cut+"/STVStudies_Overlay9_"+run+cut+".root";
	TFile* file = new TFile(file_path,"readonly");
	vector<int> colors{OverlayColor,kOrange-3,kGreen+1,kRed+1,kBlue};
	vector<int> shapes{20,21,22,23,33};


	//--------------------------------------//

	vector<TLegend*> leg; leg.clear(); leg.resize(nplots);
	vector<TCanvas*> canvas; canvas.clear(); canvas.resize(nplots);

	// 1st index = plot name
	// 2nd index = slice
	vector<vector<TH1D*> > plot; plot.clear(); plot.resize(nplots);

	//--------------------------------------//

	// Loop over the plots
	
	for (int iplot = 0; iplot < nplots; iplot ++) {

		leg.at(iplot) = new TLegend(0.25,0.83,0.9,0.98);

		TString canvas_name = "canvas_"+plot_names[iplot];
		canvas.at(iplot) = new TCanvas(canvas_name,canvas_name,205,34,1024,768);
		canvas.at(iplot)->SetBottomMargin(0.16);
		canvas.at(iplot)->SetTopMargin(0.18);
		canvas.at(iplot)->SetLeftMargin(0.15);

		plot[iplot].resize( slices.at(iplot).size() );

		// Loop over the slices
		
		for (int islice = 0; islice < (int)slices.at(iplot).size()-1; islice++) {

			TString slice_name = convert_to_string(slices.at(iplot).at(islice))+"To"+convert_to_string(slices.at(iplot).at(islice+1));
			plot.at(iplot).at(islice) = (TH1D*)( file->Get( plot_names.at(iplot) + slice_name +"Plot") );

			plot.at(iplot).at(islice)->GetXaxis()->CenterTitle();
                        plot.at(iplot).at(islice)->GetXaxis()->SetTitleSize(size);
                        plot.at(iplot).at(islice)->GetXaxis()->SetLabelSize(size);
                        plot.at(iplot).at(islice)->GetXaxis()->SetTitleFont(font);
                        plot.at(iplot).at(islice)->GetXaxis()->SetLabelFont(font);
	                plot.at(iplot).at(islice)->GetXaxis()->SetTitleOffset(1.);
	                plot.at(iplot).at(islice)->GetXaxis()->SetNdivisions(8);

			plot.at(iplot).at(islice)->GetYaxis()->CenterTitle();
                        plot.at(iplot).at(islice)->GetYaxis()->SetTitleSize(size);
                        plot.at(iplot).at(islice)->GetYaxis()->SetLabelSize(size);
                        plot.at(iplot).at(islice)->GetYaxis()->SetTitleFont(font);
                        plot.at(iplot).at(islice)->GetYaxis()->SetLabelFont(font);
	                plot.at(iplot).at(islice)->GetYaxis()->SetTitleOffset(1.);
	                plot.at(iplot).at(islice)->GetYaxis()->SetNdivisions(8);

			// peak at 1
			plot.at(iplot).at(islice)->Scale(1./plot.at(iplot).at(islice)->GetMaximum());
	                plot.at(iplot).at(islice)->GetYaxis()->SetRangeUser(0.001,1.1);
	                plot.at(iplot).at(islice)->GetYaxis()->SetTitle("Peak normalized to 1");

			plot.at(iplot).at(islice)->SetLineColor(colors.at(islice));
			plot.at(iplot).at(islice)->SetMarkerColor(colors.at(islice));
			plot.at(iplot).at(islice)->SetMarkerStyle(shapes.at(islice));
			plot.at(iplot).at(islice)->SetMarkerSize(2.);

			plot.at(iplot).at(islice)->Draw("ep0 same");

			TString mean = to_string_with_precision( plot.at(iplot).at(islice)->GetMean() );
			TString std = to_string_with_precision( plot.at(iplot).at(islice)->GetStdDev() );


			TString leg_entry = to_string_with_precision(slices.at(iplot).at(islice)) +" < " + variable.at(iplot) + " < " + to_string_with_precision(slices.at(iplot).at(islice+1)) + ", #mu = " + mean + ", #sigma = " + std;
			TLegendEntry* lMC = leg.at(iplot)->AddEntry(plot.at(iplot).at(islice),leg_entry,"ep0");
			lMC->SetTextColor(colors.at(islice));	

		} // End of the loop over the slices

		gPad->RedrawAxis();

		leg.at(iplot)->SetBorderSize(0);
		leg.at(iplot)->SetTextSize(size-0.03);
		leg.at(iplot)->SetTextFont(font);
		leg.at(iplot)->SetNColumns(1);
		leg.at(iplot)->SetMargin(0.1);
		leg.at(iplot)->Draw();

                TString canvas_path = PlotPath + cut + "/";
                TString canvas_export_name = "Resolution_"+plot_names.at(iplot)+"_"+run+"_"+UBCodeVersion+cut+".pdf";
                canvas.at(iplot)->SaveAs(canvas_path + canvas_export_name);
                delete canvas.at(iplot);
	

	} // End of the loop over the plots

} // End of the program 
