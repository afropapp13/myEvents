#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPad.h>

#include <iostream>
#include <vector>

#include "../Secondary_Code/myFunctions.cpp"
#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void PeLEE_Create1DPlotsTHStack_TopologicalBreakDown(TString BaseMC = "") {

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::vector<TString> PlotNames; PlotNames.clear();

	PlotNames.push_back("RecoMuonMomentumPlot");
	PlotNames.push_back("RecoProtonMomentumPlot");
	PlotNames.push_back("RecoMuonCosThetaPlot");
	PlotNames.push_back("RecoProtonCosThetaPlot");
	PlotNames.push_back("RecoMuonPhiPlot");
	PlotNames.push_back("RecoProtonPhiPlot");
	PlotNames.push_back("RecoDeltaPTPlot");
	PlotNames.push_back("RecoDeltaAlphaTPlot");
	PlotNames.push_back("RecoDeltaPhiTPlot");

	if (BaseMC == "") {

	PlotNames.push_back("RecoCCQEMuonMomentumPlot");
	PlotNames.push_back("RecoCCQEProtonMomentumPlot"); 
	PlotNames.push_back("RecoCCQEMuonCosThetaPlot");
	PlotNames.push_back("RecoCCQEProtonCosThetaPlot");
	PlotNames.push_back("RecoCCQEMuonPhiPlot");
	PlotNames.push_back("RecoCCQEProtonPhiPlot");

//	PlotNames.push_back("RecoPMissMinusPlot");
//	PlotNames.push_back("RecoPMissPlot");
//	PlotNames.push_back("RecokMissPlot");

	PlotNames.push_back("RecoNuScorePlot");
//	PlotNames.push_back("RecoFlashScorePlot");
//	PlotNames.push_back("RecoDistancePlot");
//	PlotNames.push_back("RecoLengthDifferencePlot");
//	PlotNames.push_back("RecodYZPlot");
//	PlotNames.push_back("RecoNPEPlot");

	PlotNames.push_back("RecoDeltaPhiPlot");
	PlotNames.push_back("RecoDeltaThetaPlot");
//	PlotNames.push_back("RecoDeltaForwardThetaPlot");
//	PlotNames.push_back("RecoDeltaBackwardThetaPlot");

//	PlotNames.push_back("RecoThreePlaneChi2LogLikelihoodCandidateMuonPlot");
//	PlotNames.push_back("RecoThreePlaneChi2LogLikelihoodCandidateProtonPlot");

	PlotNames.push_back("RecoMuonLLRPIDPlot");
	PlotNames.push_back("RecoProtonLLRPIDPlot");

	PlotNames.push_back("RecoECalPlot");
	PlotNames.push_back("RecoEQEPlot");
	PlotNames.push_back("RecoQ2Plot");

	PlotNames.push_back("RecoMuonLengthPlot");
	PlotNames.push_back("RecoProtonLengthPlot");
//	PlotNames.push_back("RecodMuonTracksScorePlot");
//	PlotNames.push_back("RecodProtonTracksScorePlot");
	PlotNames.push_back("RecodMuonVertexDistancePlot");
	PlotNames.push_back("RecodProtonVertexDistancePlot");
//	PlotNames.push_back("RecoVertexActivityPlot");
//	PlotNames.push_back("RecoNonZeroVertexActivityPlot");

	PlotNames.push_back("RecoVertexXPlot");
	PlotNames.push_back("RecoVertexYPlot");
	PlotNames.push_back("RecoVertexZPlot");

	PlotNames.push_back("RecoEvPlot");
	PlotNames.push_back("RecoNuPlot");

	//PlotNames.push_back("RecoContainedMuonMomentumPlot");
	//PlotNames.push_back("RecoUncontainedMuonMomentumPlot");
	// PlotNames.push_back("RecoContainedMuonLengthPlot");
	// PlotNames.push_back("RecoUncontainedMuonLengthPlot");

	}

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

	//vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// For the alternative models for now we have only Run1

		if (BaseMC != "" && Runs[WhichRun] == "Run3") { continue; }

		// -----------------------------------------------------------------------------------------------------------------------------------------

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);													

		// -----------------------------------------------------------------------------------------------------------------------------------------

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];

			// For the alternative MC, we want the figures after the application of all cuts
			if (BaseMC == "Overlay9NuWro" && i != NCuts-1) { continue; }	
			if (BaseMC == "NoTuneOverlay9" && i != NCuts-1) { continue; }		

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles+"/"+Cuts+"/";

			TH1D::SetDefaultSumw2();

			// ---------------------------------------------------------------------------------------------------------------------

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			gStyle->SetPalette(55); const Int_t NCont = 999; 
			gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > CC1pPlots; CC1pPlots.clear();
			vector<vector<TH1D*> > NonCC1pPlots; NonCC1pPlots.clear();

			vector<vector<TH1D*> > hratio;  hratio.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;

			// 0: BeamOn
			// 1: Overlay
			// 2: ExtBNB
			// 3: Dirt
		
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");

			if (BaseMC == "") { NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }
			else if (BaseMC == "Overlay9NuWro") { NameOfSamples.push_back("STVStudies_Overlay9NuWro_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }
			else if (BaseMC == "NoTuneOverlay9") { NameOfSamples.push_back("NoTuneSTVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }

			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");

			if (BaseMC != "NoTuneOverlay9") { NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
			else { NameOfSamples.push_back("NoTuneSTVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }

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

				vector<TH1D*> Currenthratio;  Currenthratio.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					TH1D* CC1phist = (TH1D*)(FileSample[WhichSample]->Get("CC1p"+PlotNames[WhichPlot]));
					TH1D* NonCC1phist = (TH1D*)(FileSample[WhichSample]->Get("NonCC1p"+PlotNames[WhichPlot]));

					hist->GetXaxis()->CenterTitle();
					hist->GetYaxis()->CenterTitle();					

					hist->SetLineColor(Colors[WhichSample]);
				
					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						hist->SetMarkerSize(2.); 
					}

					CurrentPlots.push_back(hist);
					CC1pCurrentPlots.push_back(CC1phist);
					NonCC1pCurrentPlots.push_back(NonCC1phist);

					Currenthratio.push_back((TH1D*)hist->Clone());
			
				}

				Plots.push_back(CurrentPlots);
				CC1pPlots.push_back(CC1pCurrentPlots);
				NonCC1pPlots.push_back(NonCC1pCurrentPlots);

				hratio.push_back(Currenthratio);

			}

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
		
				TString PlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+Cuts;
				PlotCanvas.push_back(new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768));
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
				botPad->SetBottomMargin(0.25);
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
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT));
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.3*max);

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
							leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample],"f");
							THStacks[WhichPlot]->Draw("same");

					}

					if (LabelsOfSamples[WhichSample] == "Dirt") {

							CC1pPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CC1pPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							THStacks[WhichPlot]->Add(CC1pPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(CC1pPlots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample] + " CC1p","f");
							THStacks[WhichPlot]->Draw("same");

							NonCC1pPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							NonCC1pPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							THStacks[WhichPlot]->Add(NonCC1pPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(NonCC1pPlots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample] + " NonCC1p","f");
							THStacks[WhichPlot]->Draw("same");

					}

					if (LabelsOfSamples[WhichSample] == "Overlay") {

							CC1pPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CC1pPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							THStacks[WhichPlot]->Add(CC1pPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(CC1pPlots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample] + " CC1p","f");
							THStacks[WhichPlot]->Draw("same");

							NonCC1pPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							NonCC1pPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							THStacks[WhichPlot]->Add(NonCC1pPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(NonCC1pPlots[WhichSample][WhichPlot],LabelsOfSamples[WhichSample] + " NonCC1p","f");
							THStacks[WhichPlot]->Draw("same");

					}
					

				} // End of the loop over the samples

				Plots[0][WhichPlot]->Draw("e1 same"); 

				gPad->RedrawAxis();

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.09);
				text->DrawTextNDC(0.14, 0.87, Runs[WhichRun]);

				// -------------------------------------------------------------------------------------------------------------------

				hratio[1][WhichPlot]->Add(hratio[2][WhichPlot]);
				hratio[1][WhichPlot]->Add(hratio[3][WhichPlot]);
				hratio[0][WhichPlot]->Divide(hratio[1][WhichPlot]);

				hratio[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				hratio[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				hratio[0][WhichPlot]->GetYaxis()->SetTitle("#frac{BeamOn}{MC+ExtBNB}");
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

				// --------------------------------------------------------------------------------------

				// Sum of NonBeamOn Samples

				TH1D* SumNonBeamOn = (TH1D*)Plots[1][WhichPlot]->Clone(); // ExtBNB
				SumNonBeamOn->Add(Plots[2][WhichPlot]); // Overlay
				SumNonBeamOn->Add(Plots[3][WhichPlot]); // Dirt

				// CC1p Purity 

				int CC1pPurity = CC1pPlots[1][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexPurity;
				latexPurity.SetTextFont(FontStyle);
				latexPurity.SetTextSize(0.09);
				TString LabelPurity = "CC1p = " + ToString(CC1pPurity/10.) + " %";
				latexPurity.DrawLatexNDC(0.59,0.9, LabelPurity);

				// --------------------------------------------------------------------------------------

				// Cosmic Contamination

				int CosmicContamination = Plots[2][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexCosmic;
				latexCosmic.SetTextFont(FontStyle);
				latexCosmic.SetTextSize(0.09);
				TString LabelCosmic = "Cosmics = " + ToString(CosmicContamination/10.) + " %";
				latexCosmic.DrawLatexNDC(0.59,0.8, LabelCosmic);

				// --------------------------------------------------------------------------------------

				TString CanvasPath = PlotPath + Cuts+"/TopologicalBreakDown/";
				TString CanvasName = BaseMC + "THStack_BreakDown_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs

} // End of the program 
