#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>

#include <iostream>
#include <vector>

#include "../Secondary_Code/myFunctions.cpp"
#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"

void PeLEE_Create1DPlotsTHStack_SubSpectrum() {

	// -----------------------------------------------------------------------------------------------------------------------------------------

	GlobalSettings();
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(6);	

	// -----------------------------------------------------------------------------------------------------------------------------------------

//	std::vector<TString> PlotNames; PlotNames.clear();

//	PlotNames.push_back("MuonMomentumPlot");
//	PlotNames.push_back("ProtonMomentumPlot");
//	PlotNames.push_back("MuonCosThetaPlot");
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonPhiPlot");
//	PlotNames.push_back("DeltaPTPlot");
//	PlotNames.push_back("DeltaAlphaTPlot");
//	PlotNames.push_back("DeltaPhiTPlot");

//	PlotNames.push_back("PMissMinusPlot");
//	PlotNames.push_back("PMissPlot");
//	PlotNames.push_back("kMissPlot");

//	PlotNames.push_back("DeltaPLPlot");
//	PlotNames.push_back("DeltaPnPlot");
//	PlotNames.push_back("DeltaPtxPlot");
//	PlotNames.push_back("DeltaPtyPlot");
//	PlotNames.push_back("APlot");

//	PlotNames.push_back("CCQEMuonMomentumPlot");
//	PlotNames.push_back("CCQEProtonMomentumPlot");
//	PlotNames.push_back("CCQEMuonCosThetaPlot");
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	TFile* FluxFile = TFile::Open("../mySTVAnalysis/MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);		
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);		

		// -----------------------------------------------------------------------------------------------------------------------------------------

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];

		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles+"/"+Cuts+"/";

			// ---------------------------------------------------------------------------------------------------------------------

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			gStyle->SetPalette(55); const Int_t NCont = 999; 
			gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > CC1pPlots; CC1pPlots.clear();
			vector<vector<TH1D*> > NonCC1pPlots; NonCC1pPlots.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;

			// 0: BeamOn
			// 1: Overlay
			// 2: ExtBNB
			// 3: Dirt
		
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");
			NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay");
			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");
			NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt");

			vector<int> Colors; Colors.clear(); 
			Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kGray+2); Colors.push_back(kMagenta);

			vector<int> ColorsOverlay{kBlue-5,kYellow+1,kOrange+7,kRed+1,kBlue};

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				FileSample.push_back(TFile::Open(PathToFilesCut+NameOfSamples[WhichSample]));

				vector<TH1D*> CurrentPlots; CurrentPlots.clear();
				vector<TH1D*> CC1pCurrentPlots; CC1pCurrentPlots.clear();
				vector<TH1D*> NonCC1pCurrentPlots; NonCC1pCurrentPlots.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					TH1D* CC1phist = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
					TH1D* NonCC1phist = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pReco"+PlotNames[WhichPlot]));

					hist->GetXaxis()->CenterTitle();
					hist->GetYaxis()->CenterTitle();					

					hist->SetLineColor(Colors[WhichSample]);
				
					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						hist->SetMarkerSize(1.); 
					}

					CurrentPlots.push_back(hist);
					CC1pCurrentPlots.push_back(CC1phist);
					NonCC1pCurrentPlots.push_back(NonCC1phist);
			
				}

				Plots.push_back(CurrentPlots);
				CC1pPlots.push_back(CC1pCurrentPlots);
				NonCC1pPlots.push_back(NonCC1pCurrentPlots);

			}

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
				midPad->SetTopMargin(0.03);
				botPad->SetTopMargin(0.03);
				botPad->SetBottomMargin(0.26);
				botPad->SetGridx();
				botPad->SetGridy();
				topPad->Draw();
				midPad->Draw();
				botPad->Draw();

				leg.push_back(new TLegend(0.1,0.005,0.9,0.995));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(3);

				double max = -99.;

				// Loop over the samples

				TH1D* BeamOnClone = nullptr;

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

					midPad->cd();
					Plots[WhichSample][WhichPlot]->SetTitle("");
					Plots[WhichSample][WhichPlot]->SetLineWidth(1);

					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(8);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0);

					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT));
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					if (WhichSample == 0) { BeamOnClone = (TH1D*)Plots[0][WhichPlot]->Clone(); }

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						//gStyle->SetErrorX(0); // Removing the horizontal errors
						BeamOnClone->SetLineColor(kWhite);
						BeamOnClone->SetMarkerColor(kWhite);
						BeamOnClone->Draw("e1 same"); 
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample]+" (Stat Unc)","ep");

					}

					if (LabelsOfSamples[WhichSample] == "ExtBNB") {

							BeamOnClone->Add(Plots[WhichSample][WhichPlot],-1);
					}

					if (LabelsOfSamples[WhichSample] == "Dirt") {

							BeamOnClone->Add(Plots[WhichSample][WhichPlot],-1);

					}

					if (LabelsOfSamples[WhichSample] == "Overlay") {

							CC1pPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CC1pPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							THStacks[WhichPlot]->Add(CC1pPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(CC1pPlots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample] + " CC1p","f");
							THStacks[WhichPlot]->Draw("same");

							BeamOnClone->Add(NonCC1pPlots[WhichSample][WhichPlot],-1);

					}
					

				} // End of the loop over the samples

				BeamOnClone->GetYaxis()->SetRangeUser(0.,0.9*max);
				PlotCanvas[WhichPlot]->Update();
				BeamOnClone->SetLineColor(kBlack);
				BeamOnClone->SetMarkerColor(kBlack);
				BeamOnClone->Draw("e1 same"); 

				TString NameExtractedXSec = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TFile* CovFile = new TFile(NameExtractedXSec,"readonly");

				TString CopyPlotName = PlotNames[WhichPlot];
				TString ReducedPlotName = CopyPlotName.ReplaceAll("Reco","");
				TH2D* CovMatrix = (TH2D*)(CovFile->Get("TotalCovariance_"+ReducedPlotName));
				TH2D* StatCovMatrix = (TH2D*)(CovFile->Get("StatCovariance_"+ReducedPlotName));		
				CovMatrix->Add(StatCovMatrix,-1);		
				// Sanity check, stat errors should be identical to the ones coming from the Stat covariances 
				//TH2D* CovMatrix = (TH2D*)(CovFile->Get("StatCovariance_"+ReducedPlotName));				
				//CovMatrix->Scale(TMath::Power( (IntegratedFlux*NTargets)/Units ,2.));

				TH1D* DataClone = (TH1D*)(BeamOnClone->Clone());
				int n = DataClone->GetXaxis()->GetNbins();

				for (int i = 1; i <= n;i++ ) { 

					double MCCV = CC1pPlots[1][WhichPlot]->GetBinContent(i); // 1 = Overlay samples

					DataClone->SetBinContent(i,MCCV);
					DataClone->SetBinError(i, TMath::Sqrt( CovMatrix->GetBinContent(i,i) ) * (IntegratedFlux*NTargets)/Units);
				}

				BeamOnClone->Draw("e1 same"); // Data Stat Only			

				DataClone->SetMarkerSize(0.);
				DataClone->SetMarkerColor(ColorsOverlay[0]);				
				DataClone->SetLineColor(ColorsOverlay[0]);
				DataClone->SetFillColor(kOrange+7);
				DataClone->SetFillStyle(3004);								
				DataClone->Draw("e2 same");	// Full Syst	

				BeamOnClone->Draw("e1 same"); // Data Stat Only, draw again to be on top	

				leg[WhichPlot]->AddEntry(DataClone,"Syst Unc","f");												

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.09);
				text->DrawTextNDC(0.14, 0.9, Runs[WhichRun]);

				// -------------------------------------------------------------------------------------------------------------------

				TH1D* Ratio = (TH1D*)(BeamOnClone->Clone());

				Ratio->Divide(CC1pPlots[1][WhichPlot]);

				Ratio->GetXaxis()->SetTitleFont(FontStyle);
				Ratio->GetXaxis()->SetLabelFont(FontStyle);
				Ratio->GetYaxis()->SetTitle("#frac{BeamOn}{Overlay CC1p}");
				Ratio->GetXaxis()->SetTitle(Plots[0][WhichPlot]->GetXaxis()->GetTitle());
				Ratio->GetXaxis()->SetTitleSize(0.13);
				Ratio->GetXaxis()->SetLabelSize(0.12);
				Ratio->GetXaxis()->SetTitleOffset(0.88);
				Ratio->GetXaxis()->SetNdivisions(8);

				Ratio->GetYaxis()->SetTitleFont(FontStyle);
				Ratio->GetYaxis()->SetLabelFont(FontStyle);
				Ratio->GetYaxis()->SetRangeUser(0.51,1.49);
				Ratio->GetYaxis()->SetNdivisions(6);
				Ratio->GetYaxis()->SetTitleOffset(0.35);
				Ratio->GetYaxis()->SetTitleSize(0.1);
				Ratio->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				Ratio->SetLineColor(kBlack);
				Ratio->Draw("e1 hist same");

				double RatioMin = Ratio->GetXaxis()->GetXmin();
				double RatioMax = Ratio->GetXaxis()->GetXmax();
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

				// --------------------------------------------------------------------------------------

				TString CanvasPath = PlotPath + Cuts+"/TopologicalBreakDown/";
				TString CanvasName = "THStack_Syst_BreakDown_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		//} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs

} // End of the program 
