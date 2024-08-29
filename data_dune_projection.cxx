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

#include "../myClasses/myFunctions.cpp"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void data_dune_projection(TString BaseMC = "") {

	//-------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(0.07,"t");

	//-------------------------------------//

	vector<TString> PlotNames; PlotNames.clear();

	PlotNames.push_back("RecoThetaVisPlot");

	//PlotNames.push_back("RecoThetaVis_DeltaPn_0_20To0_40Plot");	
	//PlotNames.push_back("RecoThetaVis_DeltaPn_0_40To1_00Plot");	
	//PlotNames.push_back("RecoSerialThetaVis_DeltaPnPlot");	
	
	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//-------------------------------------//

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("_PID_NuScore_CRT");

	int NCuts = (int)(VectorCuts.size());	

	//-------------------------------------//

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//-------------------------------------//

	vector<TString> samples;
	vector<TString> leg_samples;
	vector<TString> labels;
	vector<TString> reweight;
	vector<int> colors;

	samples.push_back("BeamOn9"); colors.push_back(kBlack); leg_samples.push_back("Data");
	samples.push_back("Overlay9"); colors.push_back(OverlayColor); leg_samples.push_back("MC");
	samples.push_back("OverlayDirt9"); colors.push_back(kGreen+1); leg_samples.push_back("Dirt");
	samples.push_back("ExtBNB9"); colors.push_back(kGray); leg_samples.push_back("Cosmics");

	reweight.push_back(""); labels.push_back("MicroBooNE");
	reweight.push_back("BNBToHondaECal"); labels.push_back("atmospheric DUNE projection");

	int nsamples = samples.size();
	int nreweight = reweight.size();	
	
	//-------------------------------------//

	// loop over the runs
	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//-------------------------------------//

		Cuts = "_NoCuts";

		// loop over the selection cuts
		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];
			TString PathToFilesCut = PathToFiles + "/"+Cuts+"/";

			//-------------------------------------//

			// loop over the reweighted samples
			for (int ireweight = 0; ireweight < nreweight; ireweight++) {

				//-------------------------------------//

				vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
				vector<TLegend*> leg; leg.clear();

				// loop over the plots
				for (int iplot = 0; iplot < N1DPlots; iplot ++) {

					TString PlotCanvasName = reweight[ireweight] + "_" + Runs[WhichRun]+"_"+PlotNames[iplot]+Cuts;
					PlotCanvas.push_back(new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768));
					PlotCanvas[iplot]->cd();

					TPad *topPad = new TPad("topPad", "", 0.005, 0.92, 0.995, 0.995);
					TPad *midPad = new TPad("midPad", "", 0.005, 0.01  , 0.995, 0.92);
					topPad->SetTopMargin(0.3);
					topPad->SetBottomMargin(0.0);
					midPad->SetBottomMargin(0.16);
					midPad->SetTopMargin(0.03);
					topPad->Draw();
					midPad->Draw();

					leg.push_back(new TLegend(0.1,0.005,0.93,0.995));
					leg[iplot]->SetBorderSize(0);
					leg[iplot]->SetNColumns(4);

					vector<TH1D*> Plots; Plots.clear();
					THStack* thstack = new THStack(PlotNames[iplot],"");	
	
					// loop over the samples
					for (int isample = 0; isample < nsamples; isample ++) {

						TString NameOfSamples = PathToFilesCut + "STVStudies_" + samples[isample] + reweight[ireweight] + "_"+Runs[WhichRun]+Cuts+".root"; 	
						TFile* FileSample = new TFile(NameOfSamples, "readonly");
						midPad->cd();

						Plots.push_back( (TH1D*)(FileSample->Get(PlotNames[iplot])) );
						Plots[isample]->SetLineColor( colors[isample] );
						Plots[isample]->SetFillColor( colors[isample] );
						Reweight(Plots[isample]);

						if ( string( samples[isample] ).find("BeamOn") != std::string::npos) {

							Plots[isample]->SetMarkerSize(1.);
							Plots[isample]->SetMarkerStyle(20);
							Plots[isample]->SetTitle("");
							Plots[isample]->SetLineWidth(1);
	
							Plots[isample]->GetXaxis()->SetTitleFont(FontStyle);
							Plots[isample]->GetXaxis()->SetTitleSize(TextSize);
							Plots[isample]->GetXaxis()->SetLabelFont(FontStyle);
							Plots[isample]->GetXaxis()->SetLabelSize(TextSize);
							Plots[isample]->GetXaxis()->SetNdivisions(8);
							Plots[isample]->GetXaxis()->CenterTitle();
	
							Plots[isample]->GetYaxis()->SetTitleFont(FontStyle);
							Plots[isample]->GetYaxis()->SetLabelFont(FontStyle);
							Plots[isample]->GetYaxis()->SetNdivisions(6);
							Plots[isample]->GetYaxis()->SetLabelSize(TextSize);
							Plots[isample]->GetYaxis()->SetTitle(Runs[WhichRun] + " events / bin");
							Plots[isample]->GetYaxis()->SetTitleSize(TextSize);
							Plots[isample]->GetYaxis()->SetTitleOffset(0.75);
							Plots[isample]->GetYaxis()->SetTickSize(0.01);
							Plots[isample]->GetYaxis()->CenterTitle();
							
							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[isample]->Draw("same e1");

							leg[iplot]->AddEntry(Plots[isample],leg_samples[isample],"pe"); 
	
						} else {

							thstack->Add(Plots[isample],"hist");
							//thstack->Draw("hist same");
						
							TH1D* stack = (TH1D*)(thstack->GetStack()->Last());
							stack->Draw("same hist");

							leg[iplot]->AddEntry(Plots[isample],leg_samples[isample],"f"); 
	
						}	
						
					} // end of the loop over the samples				


					thstack->Draw("same hist");	
					// Redraw the data
					Plots[0]->Draw("same e1");	
					gPad->RedrawAxis();
	
					//----------------------------------------//

					TH1D* stack = (TH1D*)(thstack->GetStack()->Last());
					
					double mean = Plots[0]->GetMean();
					double median = GetMedian(Plots[0]);
					double sigma = Plots[0]->GetRMS();
					double peak = FindOneDimHistoMaxValueBin(Plots[0]);
	
					TLatex latex;
					latex.SetTextFont(FontStyle);
					latex.SetTextSize(TextSize);
					TString label = "#splitline{peak = " + to_string_with_precision(peak,2) + "^{o}, median = " + to_string_with_precision(median,2)  + "^{o}}{#mu = " + to_string_with_precision(mean,2) + "^{o}, #sigma = " + to_string_with_precision(sigma,2) + "^{o}}";
					latex.DrawLatexNDC(0.3,0.7, label);				
					latex.DrawLatexNDC(0.3,0.8, "Data " + labels[ireweight]);				

					//----------------------------------------//

					topPad->cd();
					leg[iplot]->SetTextSize(0.6);
					leg[iplot]->SetTextFont(FontStyle);
					leg[iplot]->Draw();

					//----------------------------------------//

					TString CanvasPath = PlotPath + Cuts + "/InteractionBreakDown/";
					TString CanvasName = BaseMC + "data_dune_projection_"+reweight[ireweight] + "_"+PlotNames[iplot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
					PlotCanvas[iplot]->SaveAs(CanvasPath+CanvasName);
					delete PlotCanvas[iplot];

				} // End of the loop over the plots

			} // end of the loop over the reweighted samples

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs	

} // End of the program 
