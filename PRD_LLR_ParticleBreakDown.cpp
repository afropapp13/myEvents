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

#include "../Secondary_Code/myFunctions.cpp"
#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void PRD_LLR_ParticleBreakDown() {

	vector<TString> PlotNames; PlotNames.clear();

	TString StorePath = "myPlots/";
	gStyle->SetOptStat(0);	
	double TextSize = 0.06;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Keep them for historical reasons but don't plot them
//	PlotNames.push_back("RecoChi2Plot");
//	PlotNames.push_back("RecoThreePlaneChi2Plot");
//	PlotNames.push_back("RecoThreePlaneChi2LogLikelihoodPlot");
	PlotNames.push_back("RecoLLRPIDPlot");
	PlotNames.push_back("RecoMuonLLRPIDPlot");
	PlotNames.push_back("RecoProtonLLRPIDPlot");	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// Selection cuts

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("");
//	VectorCuts.push_back("_PID");
//	VectorCuts.push_back("_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	// -------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -------------------------------------------------------------------------------------------------------------------------------------

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);

		// -------------------------------------------------------------------------------------------------------------------------------------

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];

		//	} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles + "/" + Cuts + "/";
			TH1D::SetDefaultSumw2();

			// --------------------------------------------------------------------------------------------------------------------------

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
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

//			vector<int> ColorsOverlay{kBlue-5,kOrange+7,kYellow+1,kRed+1,kBlue};
			vector<int> ColorsOverlay{OverlayColor,kOrange-3,kGreen+1,kRed+1,kBlue};

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				FileSample.push_back(TFile::Open(PathToFilesCut+NameOfSamples[WhichSample]));

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

					//CenterAxisTitle(hist);

					hist->SetLineColor(Colors[WhichSample]);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						hist->SetMarkerSize(1.); 
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

			// ---------------------------------------------------------------------------------------------------------------------------

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
				TString CanvasNameString = PlotNames[WhichPlot]+Cuts + "_" + Runs[WhichRun];
				PlotCanvas.push_back(new TCanvas(CanvasNameString,CanvasNameString,205,34,1024,768));
				PlotCanvas[WhichPlot]->SetBottomMargin(0.12);
				PlotCanvas[WhichPlot]->SetLeftMargin(0.14);
				PlotCanvas[WhichPlot]->cd();

				THStacks.push_back(new THStack(PlotNames[WhichPlot],""));

				leg.push_back(new TLegend(0.12,0.91,0.9,0.995));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(5);

				double max = -99.;

				// Loop over the samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

					Plots[WhichSample][WhichPlot]->SetTitle("");
					Plots[WhichSample][WhichPlot]->SetLineWidth(1);

					Plots[WhichSample][WhichPlot]->GetXaxis()->SetRangeUser(-0.9,1.);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(8);					
					Plots[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);					
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);					

					Plots[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);					
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("Number of events / bin");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.2);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);					
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0.02);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.05*max);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						Plots[WhichSample][WhichPlot]->Draw("e same"); 
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],"BNB Data","ep");

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

							//THStacks[WhichPlot]->Add(CosmicPlots[WhichSample][WhichPlot],"hist");
							//THStacks[WhichPlot]->Add(CosmicPlots[WhichSample+2][WhichPlot],"hist");

							// Tweak to get the ExtBNB label at the right place next to BeamOn
							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic","f");

							leg[WhichPlot]->AddEntry(MuonPlots[WhichSample][WhichPlot],"#mu","f"); 
							leg[WhichPlot]->AddEntry(ProtonPlots[WhichSample][WhichPlot],"p","f"); 
							leg[WhichPlot]->AddEntry(PionPlots[WhichSample][WhichPlot],"#pi^{#pm}","f"); 
							//leg[WhichPlot]->AddEntry(CosmicPlots[WhichSample][WhichPlot],"Overlay Cosmic","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}
				

				} // End of the loop over the samples

				Plots[0][WhichPlot]->Draw("e same"); 

				// ---------------------------------------------------------------------------------------------------------------------

				leg[WhichPlot]->SetTextSize(TextSize);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->SetMargin(0.5);
				leg[WhichPlot]->Draw();

				//----------------------------------------//

				PlotCanvas[WhichPlot]->cd();
				TLatex *textPOT = new TLatex();
				textPOT->SetTextFont(FontStyle);
				textPOT->SetTextSize(TextSize);
				textPOT->DrawLatexNDC(0.2, 0.83,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");								

				//----------------------------------------//				

				gPad->RedrawAxis();

				TString CanvasPath = PlotPath + Cuts+"/";
				TString CanvasName = "Fig131.pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				//delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the run numbers

} // End of the program 
