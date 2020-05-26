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

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void Chi2PID_BreakDown() {

	vector<TString> PlotNames; PlotNames.clear();

	TString StorePath = "myPlots/";

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Keep them for historical reasons but don't plot them
//	PlotNames.push_back("RecoChi2Plot");
//	PlotNames.push_back("RecoThreePlaneChi2Plot");
	PlotNames.push_back("RecoThreePlaneChi2LogLikelihoodPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// Selection cuts

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("");
//	VectorCuts.push_back("_NuScore");
//	VectorCuts.push_back("_ThreePlaneLogChi2");
//	VectorCuts.push_back("_Collinearity");

	int NCuts = (int)(VectorCuts.size());	

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// ----------------------------------------------------------------------------------------------------------------------------------------------------------------

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];

		//	} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFiles = "OutputFiles/" + UBCodeVersion + "/" + Cuts + "/";
			TH1D::SetDefaultSumw2();

			// ------------------------------------------------------------------------------------------------------------------------------------------------

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();
			vector<TLegend*> leg; leg.clear();

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > MuonPlots; MuonPlots.clear();
			vector<vector<TH1D*> > ProtonPlots; ProtonPlots.clear();
			vector<vector<TH1D*> > PionPlots; PionPlots.clear();
			vector<vector<TH1D*> > CosmicPlots; CosmicPlots.clear();

			vector<vector<TH1D*> > hratio;  hratio.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;
	
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");
			NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay");
			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");
			NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt");

			vector<int> Colors; Colors.clear(); 
			Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kGray+2); Colors.push_back(kMagenta);

			vector<int> ColorsOverlay; ColorsOverlay.clear(); 
			ColorsOverlay.push_back(kRed); ColorsOverlay.push_back(kBlue); ColorsOverlay.push_back(kCyan); ColorsOverlay.push_back(kMagenta); ColorsOverlay.push_back(kOrange+7);

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				FileSample.push_back(TFile::Open(PathToFiles+NameOfSamples[WhichSample]));

				vector<TH1D*> CurrentPlots; CurrentPlots.clear();
				vector<TH1D*> MuonCurrentPlots; MuonCurrentPlots.clear();
				vector<TH1D*> ProtonCurrentPlots; ProtonCurrentPlots.clear();
				vector<TH1D*> PionCurrentPlots; PionCurrentPlots.clear();
				vector<TH1D*> CosmicCurrentPlots; CosmicCurrentPlots.clear();

				vector<TH1D*> Currenthratio;  Currenthratio.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					TH1D* Muonhist = (TH1D*)(FileSample[WhichSample]->Get("Muon"+PlotNames[WhichPlot]));
					TH1D* Protonhist = (TH1D*)(FileSample[WhichSample]->Get("Proton"+PlotNames[WhichPlot]));
					TH1D* Pionhist = (TH1D*)(FileSample[WhichSample]->Get("Pion"+PlotNames[WhichPlot]));
					TH1D* Cosmichist = (TH1D*)(FileSample[WhichSample]->Get("Cosmic"+PlotNames[WhichPlot]));

					CenterAxisTitle(hist);

					hist->SetLineColor(Colors[WhichSample]);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						hist->SetMarkerSize(2.); 
					}

					CurrentPlots.push_back(hist);
					MuonCurrentPlots.push_back(Muonhist);
					ProtonCurrentPlots.push_back(Protonhist);
					PionCurrentPlots.push_back(Pionhist);
					CosmicCurrentPlots.push_back(Cosmichist);

					Currenthratio.push_back((TH1D*)hist->Clone());
		
				}

				Plots.push_back(CurrentPlots);
				MuonPlots.push_back(MuonCurrentPlots);
				ProtonPlots.push_back(ProtonCurrentPlots);
				PionPlots.push_back(PionCurrentPlots);
				CosmicPlots.push_back(CosmicCurrentPlots);

				hratio.push_back(Currenthratio);

			}

			// -------------------------------------------------------------------------------------------------------------------------------------------------------

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
				PlotCanvas.push_back(new TCanvas(PlotNames[WhichPlot]+Cuts,PlotNames[WhichPlot]+Cuts,205,34,1024,768));
				PlotCanvas[WhichPlot]->cd();

				THStacks.push_back(new THStack(PlotNames[WhichPlot],""));

				TPad *topPad = new TPad("topPad", "", 0.005, 0.92, 0.995, 0.995);
				TPad *midPad = new TPad("midPad", "", 0.005, 0.3  , 0.995, 0.92);
				TPad *botPad = new TPad("botPad", "", 0.005, 0.005, 0.995, 0.3);
				topPad->SetTopMargin(0.3);
				topPad->SetBottomMargin(0.0);
				midPad->SetBottomMargin(0.03);
				midPad->SetTopMargin(0.0);
				botPad->SetTopMargin(0.);
				botPad->SetBottomMargin(0.25);
				botPad->SetGridy();
				topPad->Draw();
				midPad->Draw();
				botPad->Draw();

				leg.push_back(new TLegend(0.1,0.005,0.9,0.995));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(3);

				double max = -99.;

				// Loop over the samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

					midPad->cd();
					Plots[WhichSample][WhichPlot]->SetTitle("");
					Plots[WhichSample][WhichPlot]->SetLineWidth(4);

					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(5);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0);

					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(5);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# events");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.2*max);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						Plots[WhichSample][WhichPlot]->Draw("e1 same"); 
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample],"ep");

					}

					if (LabelsOfSamples[WhichSample] == "ExtBNB") {

							Plots[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
							Plots[WhichSample][WhichPlot]->SetFillColor(Colors[WhichSample]);
							Plots[WhichSample][WhichPlot]->SetFillStyle(3004);
							Plots[WhichSample][WhichPlot]->SetLineWidth(1);

							THStacks[WhichPlot]->Add(Plots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Draw("same");

					}

					if (LabelsOfSamples[WhichSample] == "Overlay") {

							MuonPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							MuonPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							MuonPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							MuonPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							THStacks[WhichPlot]->Add(MuonPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(MuonPlots[WhichSample+2][WhichPlot],"hist");

							ProtonPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							ProtonPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							ProtonPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							ProtonPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							THStacks[WhichPlot]->Add(ProtonPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(ProtonPlots[WhichSample+2][WhichPlot],"hist");

							PionPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							PionPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							PionPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							PionPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							THStacks[WhichPlot]->Add(PionPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(PionPlots[WhichSample+2][WhichPlot],"hist");

							CosmicPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CosmicPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							CosmicPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CosmicPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[3]);

							THStacks[WhichPlot]->Add(CosmicPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CosmicPlots[WhichSample+2][WhichPlot],"hist");

							leg[WhichPlot]->AddEntry(MuonPlots[WhichSample][WhichPlot],"#mu","f"); 
							leg[WhichPlot]->AddEntry(ProtonPlots[WhichSample][WhichPlot],"p","f"); 

							// Tweak to get the ExtBNB label at the right place below BeamOn
							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],LabelsOfSamples[2],"f");

							leg[WhichPlot]->AddEntry(PionPlots[WhichSample][WhichPlot],"#pi^{#pm}","f"); 
							leg[WhichPlot]->AddEntry(CosmicPlots[WhichSample][WhichPlot],"Overlay Cosmic","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}
				

				} // End of the loop over the samples

				Plots[0][WhichPlot]->Draw("e1 same"); 

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.08);
				text->DrawTextNDC(0.74, 0.87, Runs[WhichRun]);

				// ---------------------------------------------------------------------------------------------------------------------------------------------------

				// Ratio plots

				hratio[0][WhichPlot]->Add(hratio[2][WhichPlot],-1);
				hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
				hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

				hratio[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{BeamOn - ExtBNB}{Overlay}");
				hratio[0][WhichPlot]->GetXaxis()->SetTitle(Plots[0][WhichPlot]->GetXaxis()->GetTitle());
				hratio[0][WhichPlot]->GetXaxis()->SetTitleSize(0.13);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelSize(0.12);
				hratio[0][WhichPlot]->GetXaxis()->SetTitleOffset(0.88);
				hratio[0][WhichPlot]->GetXaxis()->SetNdivisions(5);

				hratio[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.9*hratio[0][WhichPlot]->GetMinimum(),1.1*hratio[0][WhichPlot]->GetMaximum());
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(4);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.35);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleSize(0.1);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				hratio[0][WhichPlot]->Draw("e1 hist same");

				double RatioMin = hratio[0][WhichPlot]->GetXaxis()->GetXmin();
				double RatioMax = hratio[0][WhichPlot]->GetXaxis()->GetXmax();
				double YRatioCoord = 1.2;
				TLine* RatioLine = new TLine(RatioMin,YRatioCoord,RatioMax,YRatioCoord);
				RatioLine->SetLineWidth(4);
				RatioLine->SetLineColor(kPink+8);
				RatioLine->SetLineStyle(4);
				RatioLine->Draw("same");
		
				topPad->cd();
				leg[WhichPlot]->SetTextSize(0.5);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->SetMargin(0.6);
				leg[WhichPlot]->Draw();

				PlotCanvas[WhichPlot]->SaveAs("./myPlots/pdf/1D/"+UBCodeVersion+"/"+Cuts+"/PID_BreakDown_"+PlotNames[WhichPlot]+Cuts+"_"+UBCodeVersion+".pdf");
				//delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the run numbers

} // End of the program 
