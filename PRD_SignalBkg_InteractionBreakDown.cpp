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

void PRD_SignalBkg_InteractionBreakDown(TString BaseMC = "") {

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
	//PlotNames.push_back("RecoDeltaPnPlot");
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
			vector<vector<TH1D*> > CCQESignalPlots; CCQESignalPlots.clear();
			vector<vector<TH1D*> > CCQEBkgPlots; CCQEBkgPlots.clear();			
			vector<vector<TH1D*> > CCMECSignalPlots; CCMECSignalPlots.clear();
			vector<vector<TH1D*> > CCMECBkgPlots; CCMECBkgPlots.clear();			
			vector<vector<TH1D*> > CCRESSignalPlots; CCRESSignalPlots.clear();
			vector<vector<TH1D*> > CCRESBkgPlots; CCRESBkgPlots.clear();			
			vector<vector<TH1D*> > CCDISSignalPlots; CCDISSignalPlots.clear();
			vector<vector<TH1D*> > CCDISBkgPlots; CCDISBkgPlots.clear();			
			vector<vector<TH1D*> > CCCCDISSignalPlots; CCCCDISSignalPlots.clear();
			vector<vector<TH1D*> > CCCCDISBkgPlots; CCCCDISBkgPlots.clear();			

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
				vector<TH1D*> CCQESignalCurrentPlots; CCQESignalCurrentPlots.clear();
				vector<TH1D*> CCQEBkgCurrentPlots; CCQEBkgCurrentPlots.clear();				
				vector<TH1D*> CCMECSignalCurrentPlots; CCMECSignalCurrentPlots.clear();
				vector<TH1D*> CCMECBkgCurrentPlots; CCMECBkgCurrentPlots.clear();				
				vector<TH1D*> CCRESSignalCurrentPlots; CCRESSignalCurrentPlots.clear();
				vector<TH1D*> CCRESBkgCurrentPlots; CCRESBkgCurrentPlots.clear();				
				vector<TH1D*> CCDISSignalCurrentPlots; CCDISSignalCurrentPlots.clear();
				vector<TH1D*> CCDISBkgCurrentPlots; CCDISBkgCurrentPlots.clear();				

				vector<TH1D*> Currenthratio;  Currenthratio.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					TH1D* CCQESignalhist = (TH1D*)(FileSample[WhichSample]->Get("CCQESignal"+PlotNames[WhichPlot]));
					TH1D* CCQEBkghist = (TH1D*)(FileSample[WhichSample]->Get("CCQEBkg"+PlotNames[WhichPlot]));					
					TH1D* CCMECSignalhist = (TH1D*)(FileSample[WhichSample]->Get("CCMECSignal"+PlotNames[WhichPlot]));
					TH1D* CCMECBkghist = (TH1D*)(FileSample[WhichSample]->Get("CCMECBkg"+PlotNames[WhichPlot]));					
					TH1D* CCRESSignalhist = (TH1D*)(FileSample[WhichSample]->Get("CCRESSignal"+PlotNames[WhichPlot]));
					TH1D* CCRESBkghist = (TH1D*)(FileSample[WhichSample]->Get("CCRESBkg"+PlotNames[WhichPlot]));					
					TH1D* CCDISSignalhist = (TH1D*)(FileSample[WhichSample]->Get("CCDISSignal"+PlotNames[WhichPlot]));
					TH1D* CCDISBkghist = (TH1D*)(FileSample[WhichSample]->Get("CCDISBkg"+PlotNames[WhichPlot]));					

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
					CCQESignalCurrentPlots.push_back(CCQESignalhist);
					CCQEBkgCurrentPlots.push_back(CCQEBkghist);					
					CCMECSignalCurrentPlots.push_back(CCMECSignalhist);
					CCMECBkgCurrentPlots.push_back(CCMECBkghist);					
					CCRESSignalCurrentPlots.push_back(CCRESSignalhist);
					CCRESBkgCurrentPlots.push_back(CCRESBkghist);					
					CCDISSignalCurrentPlots.push_back(CCDISSignalhist);
					CCDISBkgCurrentPlots.push_back(CCDISBkghist);					

					Currenthratio.push_back((TH1D*)hist->Clone());
		
				} // End of the loop over the plots

				Plots.push_back(CurrentPlots);
				CCQESignalPlots.push_back(CCQESignalCurrentPlots);
				CCQEBkgPlots.push_back(CCQEBkgCurrentPlots);				
				CCMECSignalPlots.push_back(CCMECSignalCurrentPlots);
				CCMECBkgPlots.push_back(CCMECBkgCurrentPlots);				
				CCRESSignalPlots.push_back(CCRESSignalCurrentPlots);
				CCRESBkgPlots.push_back(CCRESBkgCurrentPlots);				
				CCDISSignalPlots.push_back(CCDISSignalCurrentPlots);
				CCDISBkgPlots.push_back(CCDISBkgCurrentPlots);				

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

				leg.push_back(new TLegend(0.08,0.001,0.91,0.98));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(5);		

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

							CCQESignalPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQESignalPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							CCQESignalPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQESignalPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							THStacks[WhichPlot]->Add(CCQESignalPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCQESignalPlots[WhichSample+2][WhichPlot],"hist");

							CCQEBkgPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQEBkgPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							CCQEBkgPlots[WhichSample][WhichPlot]->SetFillStyle(3001);	
							CCQEBkgPlots[WhichSample][WhichPlot]->SetLineWidth(1);													
							CCQEBkgPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQEBkgPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[0]);	
							CCQEBkgPlots[WhichSample+2][WhichPlot]->SetFillStyle(3001);	
							CCQEBkgPlots[WhichSample+2][WhichPlot]->SetLineWidth(1);														
							THStacks[WhichPlot]->Add(CCQEBkgPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCQEBkgPlots[WhichSample+2][WhichPlot],"hist");							

							CCMECSignalPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECSignalPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							CCMECSignalPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECSignalPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							THStacks[WhichPlot]->Add(CCMECSignalPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCMECSignalPlots[WhichSample+2][WhichPlot],"hist");

							CCMECBkgPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECBkgPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							CCMECBkgPlots[WhichSample][WhichPlot]->SetFillStyle(3001);	
							CCMECBkgPlots[WhichSample][WhichPlot]->SetLineWidth(1);							
							CCMECBkgPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECBkgPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							CCMECBkgPlots[WhichSample+2][WhichPlot]->SetFillStyle(3001);	
							CCMECBkgPlots[WhichSample+2][WhichPlot]->SetLineWidth(1);													
							THStacks[WhichPlot]->Add(CCMECBkgPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCMECBkgPlots[WhichSample+2][WhichPlot],"hist");							

							CCRESSignalPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESSignalPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);							
							CCRESSignalPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESSignalPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[2]);						
							THStacks[WhichPlot]->Add(CCRESSignalPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCRESSignalPlots[WhichSample+2][WhichPlot],"hist");

							CCRESBkgPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESBkgPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							CCRESBkgPlots[WhichSample][WhichPlot]->SetFillStyle(3001);	
							CCRESBkgPlots[WhichSample][WhichPlot]->SetLineWidth(1);							
							CCRESBkgPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESBkgPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							CCRESBkgPlots[WhichSample+2][WhichPlot]->SetFillStyle(3004);
							CCRESBkgPlots[WhichSample+2][WhichPlot]->SetFillStyle(3001);	
							CCRESBkgPlots[WhichSample+2][WhichPlot]->SetLineWidth(1);																				
							THStacks[WhichPlot]->Add(CCRESBkgPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCRESBkgPlots[WhichSample+2][WhichPlot],"hist");							

							CCDISSignalPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISSignalPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);								
							CCDISSignalPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISSignalPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							THStacks[WhichPlot]->Add(CCDISSignalPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCDISSignalPlots[WhichSample+2][WhichPlot],"hist");

							CCDISBkgPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISBkgPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							CCDISBkgPlots[WhichSample][WhichPlot]->SetFillStyle(3001);	
							CCDISBkgPlots[WhichSample][WhichPlot]->SetLineWidth(1);							
							CCDISBkgPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISBkgPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							CCDISBkgPlots[WhichSample+2][WhichPlot]->SetFillStyle(3004);
							CCDISBkgPlots[WhichSample+2][WhichPlot]->SetLineWidth(1);														
							THStacks[WhichPlot]->Add(CCDISBkgPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(CCDISBkgPlots[WhichSample+2][WhichPlot],"hist");							

							int IntNCCQESignalEvents = (int)(std::round(CCQESignalPlots[WhichSample][WhichPlot]->Integral()));
							int IntNCCQEBkgEvents = (int)(std::round(CCQEBkgPlots[WhichSample][WhichPlot]->Integral()));							
							int IntNCCMECSignalEvents = (int)(std::round(CCMECSignalPlots[WhichSample][WhichPlot]->Integral()));
							int IntNCCMECBkgEvents = (int)(std::round(CCMECBkgPlots[WhichSample][WhichPlot]->Integral()));								
							int IntNCCRESSignalEvents = (int)(std::round(CCRESSignalPlots[WhichSample][WhichPlot]->Integral()));
							int IntNCCRESBkgEvents = (int)(std::round(CCRESBkgPlots[WhichSample][WhichPlot]->Integral()));								
							int IntNCCDISSignalEvents = (int)(std::round(CCDISSignalPlots[WhichSample][WhichPlot]->Integral()));
							int IntNCCDISBkgEvents = (int)(std::round(CCDISBkgPlots[WhichSample][WhichPlot]->Integral()));								
							int IntNExtBNBEvents = (int)(std::round(Plots[2][WhichPlot]->Integral()));
							int IntTotal = IntNCCQESignalEvents + IntNCCMECSignalEvents + IntNCCRESSignalEvents + IntNCCDISSignalEvents + IntNCCQEBkgEvents + IntNCCMECBkgEvents + IntNCCRESBkgEvents + IntNCCDISBkgEvents + IntNExtBNBEvents;

							int PercNCCQESignalEvents	  = std::round((double)IntNCCQESignalEvents / (double)IntTotal * 100.);
							int PercNCCQEBkgEvents	  = std::round((double)IntNCCQEBkgEvents / (double)IntTotal * 100.);							
							int PercNCCMECSignalEvents  = std::round((double)IntNCCMECSignalEvents / (double)IntTotal * 100.);
							int PercNCCMECBkgEvents  = std::round((double)IntNCCMECBkgEvents / (double)IntTotal * 100.);							
							int PercNCCRESSignalEvents  = std::round((double)IntNCCRESSignalEvents / (double)IntTotal * 100.);
							int PercNCCRESBkgEvents  = std::round((double)IntNCCRESBkgEvents / (double)IntTotal * 100.);							
							int PercNCCDISSignalEvents  = std::round((double)IntNCCDISSignalEvents / (double)IntTotal * 100.);
							int PercNCCDISBkgEvents  = std::round((double)IntNCCDISBkgEvents / (double)IntTotal * 100.);							
							int PercNExtBNBEvents = std::round((double)IntNExtBNBEvents / (double)IntTotal * 100.);

							TString NCCQESignalEvents = ToString(IntNCCQESignalEvents);
							TString NCCQEBkgEvents = ToString(IntNCCQEBkgEvents);							
							TString NCCMECSignalEvents = ToString(IntNCCMECSignalEvents);
							TString NCCMECBkgEvents = ToString(IntNCCMECBkgEvents);								
							TString NCCRESSignalEvents = ToString(IntNCCRESSignalEvents);
							TString NCCRESBkgEvents = ToString(IntNCCRESBkgEvents);								
							TString NCCDISSignalEvents = ToString(IntNCCDISSignalEvents);
							TString NCCDISBkgEvents = ToString(IntNCCDISBkgEvents);								
							TString NExtBNBEvents = ToString(IntNExtBNBEvents);																						

							leg[WhichPlot]->AddEntry(CCQESignalPlots[WhichSample][WhichPlot],"S QE (" + ToString(PercNCCQESignalEvents) + "%)","f");
							leg[WhichPlot]->AddEntry(CCMECSignalPlots[WhichSample][WhichPlot],"S MEC (" + ToString(PercNCCMECSignalEvents) + "%)","f");	
							leg[WhichPlot]->AddEntry(CCRESSignalPlots[WhichSample][WhichPlot],"S RES (" + ToString(PercNCCRESSignalEvents) + "%)","f");	
							leg[WhichPlot]->AddEntry(CCDISSignalPlots[WhichSample][WhichPlot],"S DIS (" + ToString(PercNCCDISSignalEvents) + "%)","f");																			

							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic (" + ToString(PercNExtBNBEvents) + "%)","f"); // ExtBNB							
							leg[WhichPlot]->AddEntry(CCQEBkgPlots[WhichSample][WhichPlot],"B QE (" + ToString(PercNCCQEBkgEvents) + "%)","f");	
							leg[WhichPlot]->AddEntry(CCMECBkgPlots[WhichSample][WhichPlot],"B MEC (" + ToString(PercNCCMECBkgEvents) + "%)","f");
							leg[WhichPlot]->AddEntry(CCRESBkgPlots[WhichSample][WhichPlot],"B RES (" + ToString(PercNCCRESBkgEvents) + "%)","f");
							leg[WhichPlot]->AddEntry(CCDISBkgPlots[WhichSample][WhichPlot],"B DIS (" + ToString(PercNCCDISBkgEvents) + "%)","f");							 

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
				leg[WhichPlot]->SetTextSize(0.4);
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
				TString CanvasName = BaseMC + "PRD_SignalBkg_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

				//----------------------------------------//				

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs	

} // End of the program 
