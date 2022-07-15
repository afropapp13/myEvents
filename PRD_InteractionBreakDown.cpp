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

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//----------------------------------------//

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);		

		//----------------------------------------//

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];											

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles + "/"+Cuts+"/";

			//----------------------------------------//

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();

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
			Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kGray+2); Colors.push_back(kMagenta);

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

					//------------------------------//				

					hist->SetLineColor(Colors[WhichSample]);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						hist->SetMarkerSize(2.); 
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

				TPad *topPad = new TPad("topPad", "", 0.005, 0.92, 0.995, 0.995);
				TPad *midPad = new TPad("midPad", "", 0.005, 0.3  , 0.995, 0.92);
				TPad *botPad = new TPad("botPad", "", 0.005, 0.01, 0.995, 0.3);
				topPad->SetTopMargin(0.3);
				topPad->SetBottomMargin(0.0);
				midPad->SetBottomMargin(0.03);
				midPad->SetTopMargin(0.03);
				botPad->SetTopMargin(0.03);
				botPad->SetBottomMargin(0.26);
				botPad->SetGridx();
				botPad->SetGridy();
				topPad->Draw();
				midPad->Draw();
				botPad->Draw();

				leg.push_back(new TLegend(0.1,0.005,0.83,0.495));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(6);

				double max = -99.;

				// Loop over the samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

					midPad->cd();
					Plots[WhichSample][WhichPlot]->SetTitle("");
					Plots[WhichSample][WhichPlot]->SetLineWidth(4);

					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(8);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0);

					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT) + "  POT");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(1E-5,1.15*max);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						Plots[WhichSample][WhichPlot]->Draw("e1 same"); 
						TString NBeamOnEvents = ToString((int)(Plots[WhichSample][WhichPlot]->GetEntries()));
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot], "MicroBooNE Data","ep");

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

							TString NCCQEEvents = ToString((int)(CCQEPlots[WhichSample][WhichPlot]->Integral()));
							TString NCCMECEvents = ToString((int)(CCMECPlots[WhichSample][WhichPlot]->Integral()));	
							TString NCCRESEvents = ToString((int)(CCRESPlots[WhichSample][WhichPlot]->Integral()));	
							TString NCCDISEvents = ToString((int)(CCDISPlots[WhichSample][WhichPlot]->Integral()));	
							TString NExtBNBEvents = ToString((int)(Plots[2][WhichPlot]->Integral()));																						

							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic","f"); // ExtBNB
							leg[WhichPlot]->AddEntry(CCQEPlots[WhichSample][WhichPlot],"QE","f"); 
							leg[WhichPlot]->AddEntry(CCMECPlots[WhichSample][WhichPlot],"MEC","f"); 
							leg[WhichPlot]->AddEntry(CCRESPlots[WhichSample][WhichPlot],"RES","f"); 
							leg[WhichPlot]->AddEntry(CCDISPlots[WhichSample][WhichPlot],"DIS","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}
				

				} // End of the loop over the samples

				Plots[0][WhichPlot]->Draw("e1 same");

				gPad->RedrawAxis();				

				//----------------------------------------//

				hratio[1][WhichPlot]->Add(hratio[2][WhichPlot]);
				hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
				hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

				hratio[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
				hratio[0][WhichPlot]->GetXaxis()->SetTitle(Plots[0][WhichPlot]->GetXaxis()->GetTitle());
				hratio[0][WhichPlot]->GetXaxis()->SetTitleSize(0.13);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelSize(0.12);
				hratio[0][WhichPlot]->GetXaxis()->SetTitleOffset(0.88);
				hratio[0][WhichPlot]->GetXaxis()->SetNdivisions(8);

				hratio[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.9*hratio[0][WhichPlot]->GetMinimum(),1.1*hratio[0][WhichPlot]->GetMaximum());
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(6);
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
				//RatioLine->Draw("same");
		
				topPad->cd();
				leg[WhichPlot]->SetTextSize(0.5);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->Draw();

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
