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

void PRD_InteractionBreakDown(TString BaseMC = "") {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	

	//int MCUncColor = kMagenta-9;
	int MCUncColor = kGray;	

	//----------------------------------------//

	vector<TString> PlotNames; PlotNames.clear();

	PlotNames.push_back("RecoDeltaPTPlot");
	PlotNames.push_back("RecoDeltaAlphaTPlot");
	PlotNames.push_back("RecoDeltaPhiTPlot");
	PlotNames.push_back("RecoDeltaPnPlot");
	PlotNames.push_back("RecoDeltaPtxPlot");
	PlotNames.push_back("RecoDeltaPtyPlot");
	PlotNames.push_back("RecoECalPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("_PID_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	//----------------------------------------//

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	TFile* FluxFile = TFile::Open("../mySTVAnalysis/MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));		

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//----------------------------------------//

		TString NameExtractedXSec = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* CovFile = new TFile(NameExtractedXSec,"readonly");		

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);	
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);			

		//----------------------------------------//

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];											

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles + "/"+Cuts+"/";

			//----------------------------------------//

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			vector<THStack*> THStacksMCUnc; THStacksMCUnc.clear();			
			gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();
			vector<TLegend*> legMC; legMC.clear();			

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > CCQEPlots; CCQEPlots.clear();
			vector<vector<TH1D*> > CCMECPlots; CCMECPlots.clear();
			vector<vector<TH1D*> > CCRESPlots; CCRESPlots.clear();
			vector<vector<TH1D*> > CCDISPlots; CCDISPlots.clear();
			vector<vector<TH1D*> > CCCCDISPlots; CCCCDISPlots.clear();

			vector<vector<TH1D*> > hratio;  hratio.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;
	
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");
			NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay");
			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");

			NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt");

			vector<int> Colors; Colors.clear(); 
			Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kGray+2); Colors.push_back(MCUncColor);

//			vector<int> ColorsOverlay{kBlue-5,kYellow+1,kOrange+7,kRed+1,kBlue};
			vector<int> ColorsOverlay{OverlayColor,kOrange-3,kGreen+1,kRed+1,kBlue};

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				FileSample.push_back(TFile::Open(PathToFilesCut+NameOfSamples[WhichSample]));

				vector<TH1D*> CurrentPlots; CurrentPlots.clear();
				vector<TH1D*> CCQECurrentPlots; CCQECurrentPlots.clear();
				vector<TH1D*> CCMECCurrentPlots; CCMECCurrentPlots.clear();
				vector<TH1D*> CCRESCurrentPlots; CCRESCurrentPlots.clear();
				vector<TH1D*> CCDISCurrentPlots; CCDISCurrentPlots.clear();

				vector<TH1D*> Currenthratio;  Currenthratio.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					TH1D* CCQEhist = (TH1D*)(FileSample[WhichSample]->Get("CCQE"+PlotNames[WhichPlot]));
					TH1D* CCMEChist = (TH1D*)(FileSample[WhichSample]->Get("CCMEC"+PlotNames[WhichPlot]));
					TH1D* CCREShist = (TH1D*)(FileSample[WhichSample]->Get("CCRES"+PlotNames[WhichPlot]));
					TH1D* CCDIShist = (TH1D*)(FileSample[WhichSample]->Get("CCDIS"+PlotNames[WhichPlot]));

					hist->GetXaxis()->CenterTitle();
					hist->GetYaxis()->CenterTitle();	

					TString TStringXaxis = 	hist->GetXaxis()->GetTitle();	
					hist->GetXaxis()->SetTitle("Reconstructed " + TStringXaxis);			

					//------------------------------//				

					hist->SetLineColor(Colors[WhichSample]);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						hist->SetMarkerSize(1.); 
					}

					CurrentPlots.push_back(hist);
					CCQECurrentPlots.push_back(CCQEhist);
					CCMECCurrentPlots.push_back(CCMEChist);
					CCRESCurrentPlots.push_back(CCREShist);
					CCDISCurrentPlots.push_back(CCDIShist);

					Currenthratio.push_back((TH1D*)hist->Clone());
		
				} // End of the loop over the plots

				Plots.push_back(CurrentPlots);
				CCQEPlots.push_back(CCQECurrentPlots);
				CCMECPlots.push_back(CCMECCurrentPlots);
				CCRESPlots.push_back(CCRESCurrentPlots);
				CCDISPlots.push_back(CCDISCurrentPlots);

				hratio.push_back(Currenthratio);

			} // End of the loop over the samples

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
				TString PlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+Cuts;
				PlotCanvas.push_back(new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768));
				PlotCanvas[WhichPlot]->cd();

				THStacks.push_back(new THStack(PlotNames[WhichPlot],""));

				TPad *topPad = new TPad("topPad", "", 0.005, 0.86, 0.995, 0.95);
				TPad *midPad = new TPad("midPad", "", 0.005, 0.3  , 0.995, 0.86);
				TPad *botPad = new TPad("botPad", "", 0.005, 0.01, 0.995, 0.3);
				topPad->SetTopMargin(0.);
				topPad->SetBottomMargin(0.0);
				midPad->SetBottomMargin(0.025);
				midPad->SetTopMargin(0.01);
				botPad->SetTopMargin(0.01);
				botPad->SetBottomMargin(0.35);
				topPad->Draw();
				midPad->Draw();
				botPad->Draw();

				leg.push_back(new TLegend(0.1,0.001,0.9,0.98));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(3);		

				double max = -99.;						

				// Loop over the samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

					midPad->cd();
					Plots[WhichSample][WhichPlot]->SetTitle("");
					Plots[WhichSample][WhichPlot]->SetLineWidth(1);

					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(0.);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(8);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0);

					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
//					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("Number of events / bin");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(1E-5,1.15*max);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						Plots[WhichSample][WhichPlot]->Draw("e same"); 
						TString NBeamOnEvents = ToString((int)(Plots[WhichSample][WhichPlot]->GetEntries()));
						//leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot], "MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT","");
						//leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot], "","");
						//leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot], "","");												
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot], "BNB Data","ep");		

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

							CCQEPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQEPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							CCQEPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQEPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							THStacks[WhichPlot]->Add(CCQEPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCQEPlots[WhichSample+2][WhichPlot],"hist");

							CCMECPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							CCMECPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							THStacks[WhichPlot]->Add(CCMECPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCMECPlots[WhichSample+2][WhichPlot],"hist");

							CCRESPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							CCRESPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							THStacks[WhichPlot]->Add(CCRESPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCRESPlots[WhichSample+2][WhichPlot],"hist");

							CCDISPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							CCDISPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							THStacks[WhichPlot]->Add(CCDISPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCDISPlots[WhichSample+2][WhichPlot],"hist");

							int IntNCCQEEvents = (int)(std::round(CCQEPlots[WhichSample][WhichPlot]->Integral()));
							int IntNCCMECEvents = (int)(std::round(CCMECPlots[WhichSample][WhichPlot]->Integral()));	
							int IntNCCRESEvents = (int)(std::round(CCRESPlots[WhichSample][WhichPlot]->Integral()));	
							int IntNCCDISEvents = (int)(std::round(CCDISPlots[WhichSample][WhichPlot]->Integral()));	
							int IntNExtBNBEvents = (int)(std::round(Plots[2][WhichPlot]->Integral()));
							int IntTotal = IntNCCQEEvents + IntNCCMECEvents + IntNCCRESEvents + IntNCCDISEvents + IntNExtBNBEvents;

							int PercNCCQEEvents	  = std::round((double)IntNCCQEEvents / (double)IntTotal * 100.);
							int PercNCCMECEvents  = std::round((double)IntNCCMECEvents / (double)IntTotal * 100.);
							int PercNCCRESEvents  = std::round((double)IntNCCRESEvents / (double)IntTotal * 100.);
							int PercNCCDISEvents  = std::round((double)IntNCCDISEvents / (double)IntTotal * 100.);
							int PercNExtBNBEvents = std::round((double)IntNExtBNBEvents / (double)IntTotal * 100.);																															

							TString NCCQEEvents = ToString(IntNCCQEEvents);
							TString NCCMECEvents = ToString(IntNCCMECEvents);	
							TString NCCRESEvents = ToString(IntNCCRESEvents);	
							TString NCCDISEvents = ToString(IntNCCDISEvents);	
							TString NExtBNBEvents = ToString(IntNExtBNBEvents);																						

							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic (" + ToString(PercNExtBNBEvents) + "%)","f"); // ExtBNB
							leg[WhichPlot]->AddEntry(CCQEPlots[WhichSample][WhichPlot],"MC QE (" + ToString(PercNCCQEEvents) + "%)","f"); 
							leg[WhichPlot]->AddEntry(CCMECPlots[WhichSample][WhichPlot],"MC MEC (" + ToString(PercNCCMECEvents) + "%)","f"); 
							leg[WhichPlot]->AddEntry(CCRESPlots[WhichSample][WhichPlot],"MC RES (" + ToString(PercNCCRESEvents) + "%)","f"); 
							leg[WhichPlot]->AddEntry(CCDISPlots[WhichSample][WhichPlot],"MC DIS (" + ToString(PercNCCDISEvents) + "%)","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}

				} // End of the loop over the samples								

				// -------------------------------------------------------------------- //				

				// Uncertainty band

				TString CopyPlotName = PlotNames[WhichPlot];
				TString ReducedPlotName = CopyPlotName.ReplaceAll("Reco","");
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
					double Unc = TMath::Sqrt( CovMatrix->GetBinContent(i,i) ) * (IntegratedFlux*NTargets)/Units;

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
				hratio[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.);
				hratio[0][WhichPlot]->GetXaxis()->SetNdivisions(8);

				hratio[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
//				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.81*hratio[0][WhichPlot]->GetMinimum(),1.19*hratio[0][WhichPlot]->GetMaximum());
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.61,1.39);
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(6);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.3);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleSize(0.15);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				hratio[0][WhichPlot]->Draw("e same");

				double RatioMin = hratio[0][WhichPlot]->GetXaxis()->GetXmin();
				double RatioMax = hratio[0][WhichPlot]->GetXaxis()->GetXmax();
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

				THStacksMCUnc[WhichPlot]->Add(MCUncDown,"hist");
				THStacksMCUnc[WhichPlot]->Draw("same");

				THStacksMCUnc[WhichPlot]->Add(MCUncTwice,"hist");
				THStacksMCUnc[WhichPlot]->Draw("same");					

				RatioLine->Draw("same");
				hratio[0][WhichPlot]->Draw("e same");	
				gPad->RedrawAxis();														

				//----------------------------------------//

				topPad->cd();
				leg[WhichPlot]->SetTextSize(0.5);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->Draw();

				//----------------------------------------//

				PlotCanvas[WhichPlot]->cd();
				TLatex *textPOT = new TLatex();
				textPOT->SetTextFont(FontStyle);
				textPOT->SetTextSize(0.045);

				if (PlotNames[WhichPlot] == "RecoDeltaPTPlot" || PlotNames[WhichPlot] == "RecoDeltaPhiTPlot"  || PlotNames[WhichPlot] == "RecoECalPlot") {

					textPOT->DrawLatexNDC(0.51, 0.81,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");				

				} else {

					textPOT->DrawLatexNDC(0.115, 0.81,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");

				}

				legMC.push_back(new TLegend(0.15,0.88,0.25,0.98));
				legMC[WhichPlot]->SetBorderSize(0);
				legMC[WhichPlot]->SetNColumns(3);
				legMC[WhichPlot]->AddEntry(MCUnc,"Prediction Uncertainty","f");
				legMC[WhichPlot]->SetTextSize(0.1);
				legMC[WhichPlot]->SetMargin(1.);				
				legMC[WhichPlot]->SetTextFont(FontStyle);	

				botPad->cd();							
				legMC[WhichPlot]->Draw();						

				//----------------------------------------//

				TString CanvasPath = PlotPath + Cuts + "/InteractionBreakDown/";
				TString CanvasName = BaseMC + "PRD_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

				//----------------------------------------//				

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs	

} // End of the program 
