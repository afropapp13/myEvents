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

void dune_projection(TString BaseMC = "") {

	//-------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(0.06,"t");

	double TextSize = 0.06;

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
	vector<TString> labels;
	
	samples.push_back("Overlay9"); labels.push_back("MicroBooNE");
	samples.push_back("Overlay9BNBToHonda"); labels.push_back("DUNE atmospheric projection");

	int nsamples = samples.size();
		
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

			// loop over the samples
			for (int isample = 0; isample < nsamples; isample++) {

				//-------------------------------------//

				TString NameOfSamples = PathToFilesCut + "STVStudies_" + samples[isample] + "_"+Runs[WhichRun]+Cuts+".root"; 	
				TFile* FileSample = new TFile(NameOfSamples, "readonly");

				vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
				vector<THStack*> THStacks; THStacks.clear();	
				vector<TLegend*> leg; leg.clear();
				vector<int> ColorsOverlay; ColorsOverlay.clear();

				vector<TH1D*> CCQEPlots; CCQEPlots.clear(); ColorsOverlay.push_back(OverlayColor);
				vector<TH1D*> CCMECPlots; CCMECPlots.clear(); ColorsOverlay.push_back(kOrange-3);
				vector<TH1D*> CCRESPlots; CCRESPlots.clear();  ColorsOverlay.push_back(kGreen+1); 
				vector<TH1D*> CCDISPlots; CCDISPlots.clear();  ColorsOverlay.push_back(kRed+1);

				for (int iplot = 0; iplot < N1DPlots; iplot ++) {

					CCQEPlots.push_back(  (TH1D*)(FileSample->Get("CCQE"+PlotNames[iplot])) );
					CCMECPlots.push_back( (TH1D*)(FileSample->Get("CCMEC"+PlotNames[iplot])) );
					CCRESPlots.push_back( (TH1D*)(FileSample->Get("CCRES"+PlotNames[iplot])) );
					CCDISPlots.push_back( (TH1D*)(FileSample->Get("CCDIS"+PlotNames[iplot])) );

					TString PlotCanvasName = samples[isample] + "_" + Runs[WhichRun]+"_"+PlotNames[iplot]+Cuts;
					PlotCanvas.push_back(new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768));
					PlotCanvas[iplot]->cd();

					THStacks.push_back(new THStack(PlotNames[iplot],""));

					TPad *topPad = new TPad("topPad", "", 0.005, 0.92, 0.995, 0.995);
					TPad *midPad = new TPad("midPad", "", 0.005, 0.01  , 0.995, 0.92);
					topPad->SetTopMargin(0.3);
					topPad->SetBottomMargin(0.0);
					midPad->SetBottomMargin(0.16);
					midPad->SetTopMargin(0.03);
					topPad->Draw();
					midPad->Draw();

					leg.push_back(new TLegend(0.09,0.005,0.9,0.995));
					leg[iplot]->SetBorderSize(0);
					leg[iplot]->SetNColumns(4);

					midPad->cd();

					// QE
					CCQEPlots[iplot]->SetLineColor(ColorsOverlay[0]);
					CCQEPlots[iplot]->SetFillColor(ColorsOverlay[0]);
					Reweight(CCQEPlots[iplot]);
					THStacks[iplot]->Add(CCQEPlots[iplot],"hist");

					// MEC
					CCMECPlots[iplot]->SetLineColor(ColorsOverlay[1]);
					CCMECPlots[iplot]->SetFillColor(ColorsOverlay[1]);
					Reweight(CCMECPlots[iplot]);
					THStacks[iplot]->Add(CCMECPlots[iplot],"hist");

					// RES
					CCRESPlots[iplot]->SetLineColor(ColorsOverlay[2]);
					CCRESPlots[iplot]->SetFillColor(ColorsOverlay[2]);
					Reweight(CCRESPlots[iplot]);
					THStacks[iplot]->Add(CCRESPlots[iplot],"hist");

					// DIS
					CCDISPlots[iplot]->SetLineColor(ColorsOverlay[3]);
					CCDISPlots[iplot]->SetFillColor(ColorsOverlay[3]);
					Reweight(CCDISPlots[iplot]);
					THStacks[iplot]->Add(CCDISPlots[iplot],"hist");

					// ---------------------------------------------//
					
					TH1D* stack = (TH1D*)(THStacks[iplot]->GetStack()->Last());
					stack->Draw("same hist");

					// ---------------------------------------------//

					TString qefrac = to_string_with_precision(CCQEPlots[iplot]->Integral()/stack->Integral()*100.,1);
					TLegendEntry* lqe = leg[iplot]->AddEntry(CCQEPlots[iplot],"QE (" + qefrac  + "%)","f");
					lqe->SetTextColor(ColorsOverlay[0]);

					TString mecfrac = to_string_with_precision(CCMECPlots[iplot]->Integral()/stack->Integral()*100.,1);
					TLegendEntry* lmec = leg[iplot]->AddEntry(CCMECPlots[iplot],"MEC (" + mecfrac  + "%)","f");
					lmec->SetTextColor(ColorsOverlay[1]);

					TString resfrac = to_string_with_precision(CCRESPlots[iplot]->Integral()/stack->Integral()*100.,1);
					TLegendEntry* lres = leg[iplot]->AddEntry(CCRESPlots[iplot],"RES (" + resfrac  + "%)","f");
					lres->SetTextColor(ColorsOverlay[2]);

					TString disfrac = to_string_with_precision(CCDISPlots[iplot]->Integral()/stack->Integral()*100.,1);
					TLegendEntry* ldis = leg[iplot]->AddEntry(CCDISPlots[iplot],"DIS (" + disfrac  + "%)","f");
					ldis->SetTextColor(ColorsOverlay[3]);
	
					// ---------------------------------------------//
					
					stack->SetTitle("");
					stack->SetLineWidth(1);
	
					stack->GetXaxis()->SetTitleFont(FontStyle);
					stack->GetXaxis()->SetTitleSize(TextSize);
					stack->GetXaxis()->SetLabelFont(FontStyle);
					stack->GetXaxis()->SetLabelSize(TextSize);
					stack->GetXaxis()->SetNdivisions(8);
					stack->GetXaxis()->CenterTitle();
	
					stack->GetYaxis()->SetTitleFont(FontStyle);
					stack->GetYaxis()->SetLabelFont(FontStyle);
					stack->GetYaxis()->SetNdivisions(6);
					stack->GetYaxis()->SetLabelSize(TextSize);

					if (Runs[WhichRun] == "Combined") {

						stack->GetYaxis()->SetTitle("Number of events / bin");

					} else {

						stack->GetYaxis()->SetTitle(Runs[WhichRun] + " events / bin");

					}

					stack->GetYaxis()->SetTitleSize(TextSize);
					stack->GetYaxis()->SetTitleOffset(0.85);
					stack->GetYaxis()->SetTickSize(0.01);
					stack->GetYaxis()->CenterTitle();
	
					THStacks[iplot]->Draw("same hist");		
					gPad->RedrawAxis();
	
					//----------------------------------------//

					double mean = stack->GetMean();
					double median = GetMedian(stack);
					double sigma = stack->GetRMS();
					double peak = FindOneDimHistoMaxValueBin(stack);
	
					TLatex latex;
					latex.SetTextFont(FontStyle);
					latex.SetTextSize(TextSize);
					TString label = "#splitline{p = " + to_string_with_precision(peak,2) + "^{o}, m = " + to_string_with_precision(median,2)  + "^{o}}{#mu = " + to_string_with_precision(mean,2) + "^{o}, #sigma = " + to_string_with_precision(sigma,2) + "^{o}}";
					latex.DrawLatexNDC(0.4,0.8, label);				
					latex.DrawLatexNDC(0.4,0.9, labels[isample]);				

					//----------------------------------------//

					topPad->cd();
					leg[iplot]->SetTextSize(0.6);
					leg[iplot]->SetTextFont(FontStyle);
					leg[iplot]->Draw();

					//----------------------------------------//

					TString CanvasPath = PlotPath + Cuts + "/InteractionBreakDown/";
					TString CanvasName = BaseMC + "dune_projection_"+samples[isample] + "_"+PlotNames[iplot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
					PlotCanvas[iplot]->SaveAs(CanvasPath+CanvasName);
					delete PlotCanvas[iplot];

				} // End of the loop over the plots

			} // end of the loop over the samples

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs	

} // End of the program 
