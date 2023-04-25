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
	vector<TString> SaveNames; SaveNames.clear();	

	TString StorePath = "myPlots/";
	gStyle->SetOptStat(0);	
	double TextSize = 0.1;
	int MCUncColor = kGray;	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// Keep them for historical reasons but don't plot them
//	PlotNames.push_back("RecoChi2Plot");
//	PlotNames.push_back("RecoThreePlaneChi2Plot");
//	PlotNames.push_back("RecoThreePlaneChi2LogLikelihoodPlot");
	PlotNames.push_back("RecoLLRPIDPlot"); SaveNames.push_back("Fig131");
	PlotNames.push_back("RecoMuonLLRPIDPlot"); SaveNames.push_back("Fig139");
	PlotNames.push_back("RecoProtonLLRPIDPlot"); SaveNames.push_back("Fig140");	

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

		TString NameExtractedXSec = MigrationMatrixPath+"WienerSVD_Total_EventRateCovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* CovFile = new TFile(NameExtractedXSec,"readonly");

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
			vector<THStack*> THStacksMCUnc; THStacksMCUnc.clear();			
			gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();
			vector<TLegend*> legMC; legMC.clear();			

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

					if (PlotNames[WhichPlot] == "RecoLLRPIDPlot") { 
						
						hist->GetXaxis()->SetRangeUser(-0.9,1.); 
						
					}
					
					if (PlotNames[WhichPlot] == "RecoProtonLLRPIDPlot") { 
						
						hist->GetXaxis()->SetRangeUser(-0.9,0.95); 
						
					}

					if (PlotNames[WhichPlot] == "RecoMuonLLRPIDPlot") { 
						
						hist->GetXaxis()->SetRangeUser(0.35,1.); 
						
					}					

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
				
				TPad *topPad = new TPad("topPad", "", 0.005, 0.86, 0.995, 0.95);
				TPad *midPad = new TPad("midPad", "", 0.005, 0.3  , 0.995, 0.86);
				TPad *botPad = new TPad("botPad", "", 0.005, 0.01, 0.995, 0.3);
				topPad->SetTopMargin(0.);
				topPad->SetBottomMargin(0.0);
				midPad->SetBottomMargin(0.025);
				midPad->SetTopMargin(0.01);
				midPad->SetLeftMargin(0.13);				
				botPad->SetTopMargin(0.01);				
				botPad->SetBottomMargin(0.35);
				botPad->SetLeftMargin(0.13);				
				topPad->Draw();
				midPad->Draw();
				botPad->Draw();				

				THStacks.push_back(new THStack(PlotNames[WhichPlot],""));

				leg.push_back(new TLegend(0.1,0.005,0.9,0.95));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(5);

				double max = -99.;

				// Loop over the samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

					Plots[WhichSample][WhichPlot]->SetTitle("");
					Plots[WhichSample][WhichPlot]->SetLineWidth(1);	

					if (PlotNames[WhichPlot] == "RecoProtonLLRPIDPlot") { 
						
						Plots[WhichSample][WhichPlot]->GetXaxis()->SetRangeUser(-0.9,0.95); 
						
					}

					if (PlotNames[WhichPlot] == "RecoMuonLLRPIDPlot") { 
						
						Plots[WhichSample][WhichPlot]->GetXaxis()->SetRangeUser(0.35,1.); 
						
					}									

					midPad->cd();
					//Plots[WhichSample][WhichPlot]->GetXaxis()->SetRangeUser(-0.9,1.);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(8);					
					Plots[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);					
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0.);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);					

					Plots[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);					
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("Number of events / bin");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.7);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);					
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0.02);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.1,1.05*max);

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
							leg[WhichPlot]->AddEntry(PionPlots[WhichSample][WhichPlot],"#pi","f"); 
							//leg[WhichPlot]->AddEntry(CosmicPlots[WhichSample][WhichPlot],"Overlay Cosmic","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}
				

				} // End of the loop over the samples

				Plots[0][WhichPlot]->Draw("e same"); 

				// ---------------------------------------------------------------------------------------------------------------------

				topPad->cd();
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->SetMargin(0.5);
				leg[WhichPlot]->Draw();

				//----------------------------------------//

				midPad->cd();
				TLatex *textPOT = new TLatex();
				textPOT->SetTextFont(FontStyle);
				textPOT->SetTextSize(TextSize);
				textPOT->DrawLatexNDC(0.18, 0.89,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");								

				//----------------------------------------//				

				gPad->RedrawAxis();
				
				// -------------------------------------------------------------------- //				

				// Uncertainty band
				
				TString PlotNameDuplicate = PlotNames[WhichPlot];				
				TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
				//if (Runs[WhichRun] != "Combined") {  textSlice->DrawLatexNDC(0.115, 0.8, LatexLabel[ReducedPlotName]);	}				

				TH2D* CovMatrix = (TH2D*)(CovFile->Get("TotalCovariance_"+ReducedPlotName));
				TH2D* StatCovMatrix = (TH2D*)(CovFile->Get("StatCovariance_"+ReducedPlotName));		
				CovMatrix->Add(StatCovMatrix,-1);		
				// Sanity check, stat errors should be identical to the ones coming from the Stat covariances 
				//TH2D* CovMatrix = (TH2D*)(CovFile->Get("StatCovariance_"+ReducedPlotName));				
				//CovMatrix->Scale(TMath::Power( (IntegratedFlux*NTargets)/Units ,2.));

				int n = Plots[0][WhichPlot]->GetXaxis()->GetNbins();
				TH1D* MCUnc = (TH1D*)(Plots[0][WhichPlot]->Clone());				

				for (int i = 1; i <= n;i++ ) { 

					double MCCV = ( (TH1*) (THStacks[WhichPlot]->GetStack()->Last()) )->GetBinContent(i);
					double Unc = TMath::Sqrt( CovMatrix->GetBinContent(i,i) );

					MCUnc->SetBinContent(i,MCCV);					
					MCUnc->SetBinError(i, Unc);				

				}

				MCUnc->SetMarkerSize(0.);
				MCUnc->SetMarkerColor(MCUncColor);				
				MCUnc->SetLineColor(kWhite);
				MCUnc->SetLineWidth(1);				
				MCUnc->SetFillColor(MCUncColor);
				//MCUnc->SetFillStyle(3005);	
				//MCUnc->Draw("e2 same");										

				//gStyle->SetErrorX(0); // Removing the horizontal errors
				Plots[0][WhichPlot]->Draw("e same");				

				gPad->RedrawAxis();								

				//----------------------------------------//

				hratio[1][WhichPlot]->Add(hratio[2][WhichPlot]);
				hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
				hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

				hratio[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
				hratio[0][WhichPlot]->GetXaxis()->SetTitle(Plots[0][WhichPlot]->GetXaxis()->GetTitle());
				hratio[0][WhichPlot]->GetXaxis()->SetTitleSize(0.15);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelSize(0.12);
				hratio[0][WhichPlot]->GetXaxis()->CenterTitle();				

				hratio[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.);
				hratio[0][WhichPlot]->GetXaxis()->SetNdivisions(8);

				hratio[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.56,1.44);
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(6);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.3);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleSize(0.15);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				hratio[0][WhichPlot]->Draw("e same");

				double RatioMin = hratio[0][WhichPlot]->GetXaxis()->GetXmin();
				double RatioMax = hratio[0][WhichPlot]->GetXaxis()->GetXmax();

				if (PlotNames[WhichPlot] == "RecoLLRPIDPlot") { 
						
					RatioMin = -0.9; RatioMax = 1.; 
						
				}

				if (PlotNames[WhichPlot] == "RecoProtonLLRPIDPlot") { 
						
					RatioMin = -0.9; RatioMax = 0.95; 
						
				}

				if (PlotNames[WhichPlot] == "RecoMuonLLRPIDPlot") { 
						
					RatioMin = 0.35; RatioMax = 1.;
						
				}

				double YRatioCoord = 1.;
				TLine* RatioLine = new TLine(RatioMin,YRatioCoord,RatioMax,YRatioCoord);
				RatioLine->SetLineWidth(1);
				RatioLine->SetLineColor(kBlack);
				RatioLine->SetLineStyle(kDashed);

				//----------------------------------------//

				// Uncertainty band on ratio plot

				TH1D* MCUncDown = (TH1D*)MCUnc->Clone();
				TH1D* MCUncTwice = (TH1D*)MCUnc->Clone();				

				for (int WhichBin = 1; WhichBin <= n; WhichBin++) {

					double MCCV = MCUnc->GetBinContent(WhichBin);
					double Unc = MCUnc->GetBinError(WhichBin);
					double FracUnc = Unc / MCCV;

					MCUncDown->SetBinContent(WhichBin,1.-FracUnc);					
					MCUncTwice->SetBinContent(WhichBin,2*FracUnc);	

				}

				THStacksMCUnc.push_back(new THStack(PlotNames[WhichPlot] + "MCUnc",PlotNames[WhichPlot] + "MCUnc"));	

				MCUncDown->SetLineColor(MCUncColor);
				MCUncDown->SetFillColor(kWhite);
				MCUncDown->SetLineWidth(1);

				MCUncTwice->SetLineColor(MCUncColor);
				MCUncTwice->SetFillColor(MCUncColor);
				MCUncTwice->SetLineWidth(1);				

				botPad->cd();
				THStacksMCUnc[WhichPlot]->Add(MCUncDown,"hist");
				THStacksMCUnc[WhichPlot]->Draw("same");

				THStacksMCUnc[WhichPlot]->Add(MCUncTwice,"hist");
				THStacksMCUnc[WhichPlot]->Draw("same");					

				RatioLine->Draw("same");
				hratio[0][WhichPlot]->Draw("e same");	
				gPad->RedrawAxis();														

				//----------------------------------------//

				topPad->cd();
				leg[WhichPlot]->SetTextSize(0.6);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->Draw();

				//----------------------------------------//

				botPad->cd();

				if (PlotNames[WhichPlot] == "RecoMuonLLRPIDPlot") { legMC.push_back(new TLegend(0.15,0.88,0.35,0.98)); }
				if (PlotNames[WhichPlot] == "RecoProtonLLRPIDPlot") { legMC.push_back(new TLegend(0.55,0.88,0.75,0.98)); }
				if (PlotNames[WhichPlot] == "RecoLLRPIDPlot") { legMC.push_back(new TLegend(0.55,0.88,0.75,0.98)); }				

				legMC[WhichPlot]->SetBorderSize(0);
				legMC[WhichPlot]->SetNColumns(3);
				legMC[WhichPlot]->AddEntry(MCUnc,"Prediction Uncertainty","f");
				legMC[WhichPlot]->SetTextSize(0.1);
				legMC[WhichPlot]->SetMargin(.5);				
				legMC[WhichPlot]->SetTextFont(FontStyle);	
						
				legMC[WhichPlot]->Draw();
				
				//----------------------------------------//
				
				// Subfig label
				
				midPad->cd();
				TLatex *textsubfig = new TLatex();
				textsubfig->SetTextFont(FontStyle);
				textsubfig->SetTextSize(0.1);
				if (PlotNames[WhichPlot] == "RecoProtonLLRPIDPlot") { textsubfig->DrawLatexNDC(0.79, 0.89,"(a)"); }
				if (PlotNames[WhichPlot] == "RecoMuonLLRPIDPlot") { textsubfig->DrawLatexNDC(0.79, 0.89,"(b)"); }				
				
				//----------------------------------------//

				TString CanvasPath = PlotPath + Cuts+"/";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+SaveNames[WhichPlot]+".pdf");
				//delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the run numbers

} // End of the program 
