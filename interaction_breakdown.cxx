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

void interaction_breakdown(TString BaseMC = "") {

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	

	// ----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> PlotNames; PlotNames.clear();

	PlotNames.push_back("RecoMuonCosThetaPlot");
	PlotNames.push_back("RecoMuonCosThetaSingleBinPlot");
	PlotNames.push_back("RecoThetaVisPlot");
	PlotNames.push_back("RecoThetaVis_ECal_0_00To0_50Plot");	
	PlotNames.push_back("RecoThetaVis_ECal_0_50To0_80Plot");	
	PlotNames.push_back("RecoThetaVis_ECal_0_80To2_00Plot");	
	PlotNames.push_back("RecoSerialThetaVis_ECalPlot");	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	//VectorCuts.push_back("");
	VectorCuts.push_back("_PID_NuScore_CRT");

	int NCuts = (int)(VectorCuts.size());	

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	TFile* FluxFile = TFile::Open("../mySTVAnalysis/MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));		

	// -------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// We needs these for the uncertainty band

		TString NameExtractedXSec = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* CovFile = new TFile(NameExtractedXSec,"readonly");		

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);	

		// -------------------------------------------------------------------------------------------------------------------------------------

		bool plot_unc = false;
		if (BaseMC == "" && Runs[WhichRun] == "Combined") { plot_unc = true; }

		//-------------------------------------//

		Cuts = "_NoCuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];
	
			// For the alternative MC, we want the figures after the application of all cuts
			if (BaseMC == "Overlay9NuWro" && i != NCuts-1) { continue; }
			if (BaseMC == "NoTuneOverlay9" && i != NCuts-1) { continue; }
			if (BaseMC == "TwiceMECOverlay9" && i != NCuts-1) { continue; }		
			if (BaseMC == "GENIEv2Overlay9" && i != NCuts-1) { continue; }

			// GENIE v2 has only combined run
			if (BaseMC == "GENIEv2Overlay9" && Runs[WhichRun] != "Combined") { continue; }	

			// NuWro/Tweaked GENIE don't have Run 4a
			if (BaseMC == "Overlay9NuWro" && (Runs[WhichRun] == "Run5" || Runs[WhichRun] == "Run4a" || Runs[WhichRun] == "Run4b" || Runs[WhichRun] == "Run4aRutgers") ) { continue; }
			if (BaseMC == "NoTuneOverlay9" && (Runs[WhichRun] == "Run5" || Runs[WhichRun] == "Run4a" || Runs[WhichRun] == "Run4b" || Runs[WhichRun] == "Run4aRutgers") ) { continue; }
			if (BaseMC == "TwiceMECOverlay9" && (Runs[WhichRun] == "Run5" || Runs[WhichRun] == "Run4a" || Runs[WhichRun] == "Run4b" || Runs[WhichRun] == "Run4aRutgers") ) { continue; }												

//		} // If we want to run only on a specific cut combination, include this } and remove the one at the end of the program

			TString PathToFilesCut = PathToFiles + "/"+Cuts+"/";

			// ---------------------------------------------------------------------------------------------------------------------------

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			// Uncertainty band
			vector<THStack*> THStacksMCUnc; THStacksMCUnc.clear();			
			gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > CCQEPlots; CCQEPlots.clear();
			vector<vector<TH1D*> > CCMECPlots; CCMECPlots.clear();
			vector<vector<TH1D*> > CCRESPlots; CCRESPlots.clear();
			vector<vector<TH1D*> > CCDISPlots; CCDISPlots.clear();
			vector<vector<TH1D*> > CCCCDISPlots; CCCCDISPlots.clear();

			vector<vector<TH1D*> > bin_width_Plots; bin_width_Plots.clear();
			vector<vector<TH1D*> > bin_width_CCQEPlots; bin_width_CCQEPlots.clear();
			vector<vector<TH1D*> > bin_width_CCMECPlots; bin_width_CCMECPlots.clear();
			vector<vector<TH1D*> > bin_width_CCRESPlots; bin_width_CCRESPlots.clear();
			vector<vector<TH1D*> > bin_width_CCDISPlots; bin_width_CCDISPlots.clear();
			vector<vector<TH1D*> > bin_width_CCCCDISPlots; bin_width_CCCCDISPlots.clear();

			vector<vector<TH1D*> > hratio;  hratio.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;
	
			NameOfSamples.push_back("STVStudies_BeamOn9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");

			if (BaseMC == "") { NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }
			else if (BaseMC == "Overlay9NuWro") { NameOfSamples.push_back("STVStudies_Overlay9NuWro_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }			
			else if (BaseMC == "NoTuneOverlay9") { NameOfSamples.push_back("NoTuneSTVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }
			else if (BaseMC == "GENIEv2Overlay9") { NameOfSamples.push_back("GENIEv2STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }			
			else if (BaseMC == "TwiceMECOverlay9") { NameOfSamples.push_back("TwiceMECSTVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }			

			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");

			if (BaseMC == "NoTuneOverlay9") { NameOfSamples.push_back("NoTuneSTVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
			else if (BaseMC == "TwiceMECOverlay9") { NameOfSamples.push_back("TwiceMECSTVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
			else { NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
	
			vector<int> Colors; Colors.clear();
			// Unblind 
			Colors.push_back(kBlack); 
			// Blind
			//Colors.push_back(kWhite); 
			Colors.push_back(kRed); Colors.push_back(kGray+2); Colors.push_back(kMagenta);

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

					// The N-dimensional analysis has been developed based on the bin number, not the actual range

					if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

						TString XaxisTitle = hist->GetXaxis()->GetTitle();
						XaxisTitle.ReplaceAll("deg","bin #");
						XaxisTitle.ReplaceAll("GeV/c","bin #");
						XaxisTitle.ReplaceAll("GeV","bin #");				
						hist->GetXaxis()->SetTitle(XaxisTitle);

					}								

					//------------------------------//				

					hist->SetLineColor(Colors[WhichSample]);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 
				
						hist->SetMarkerStyle(20);
						// Unblind
						hist->SetMarkerSize(1.); 
						// Blind
						//hist->SetMarkerSize(0.); 

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

				bin_width_Plots.push_back(CurrentPlots);
				bin_width_CCQEPlots.push_back(CCQECurrentPlots);
				bin_width_CCMECPlots.push_back(CCMECCurrentPlots);
				bin_width_CCRESPlots.push_back(CCRESCurrentPlots);
				bin_width_CCDISPlots.push_back(CCDISCurrentPlots);

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
				TPad *botPad = new TPad("botPad", "", 0.005, 0.005, 0.995, 0.3);
				topPad->SetTopMargin(0.3);
				topPad->SetBottomMargin(0.0);
				midPad->SetBottomMargin(0.03);
				midPad->SetTopMargin(0.03);
				botPad->SetTopMargin(0.03);
				botPad->SetBottomMargin(0.3);
				//botPad->SetGridx();
				//botPad->SetGridy();
				topPad->Draw();
				midPad->Draw();
				botPad->Draw();

				leg.push_back(new TLegend(0.1,0.005,0.93,0.995));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(3);

				double max = -99.;

				// Loop over the samples

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
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle(Runs[WhichRun] + " events / bin");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.65);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					if (WhichSample == 0) { 

						bin_width_Plots[0][WhichPlot] = (TH1D*)(Plots[0][WhichPlot]->Clone()); 
						Reweight(bin_width_Plots[0][WhichPlot]); 
						max = FindOneDimHistoMaxValue(bin_width_Plots[0][WhichPlot]);
						bin_width_Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.3*max);

					}
					
					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						bin_width_Plots[WhichSample][WhichPlot]->Draw("e same"); 
						TString NBeamOnEvents = ToString((int)(Plots[WhichSample][WhichPlot]->GetEntries()));
						// Unblind
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],"BNB Data ("+NBeamOnEvents+")","ep");

					}

					if (LabelsOfSamples[WhichSample] == "ExtBNB") {

							Plots[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
							Plots[WhichSample][WhichPlot]->SetFillColor(Colors[WhichSample]);
							Plots[WhichSample][WhichPlot]->SetFillStyle(3004);
							Plots[WhichSample][WhichPlot]->SetLineWidth(1);

							bin_width_Plots[WhichSample][WhichPlot] = (TH1D*)(Plots[WhichSample][WhichPlot]->Clone());
							Reweight(bin_width_Plots[WhichSample][WhichPlot]);
	
							THStacks[WhichPlot]->Add(bin_width_Plots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Draw("same");

					}

					if (LabelsOfSamples[WhichSample] == "Overlay") {

							CCQEPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQEPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							bin_width_CCQEPlots[WhichSample][WhichPlot] = (TH1D*)(CCQEPlots[WhichSample][WhichPlot]->Clone());
							Reweight(bin_width_CCQEPlots[WhichSample][WhichPlot]);
							CCQEPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							CCQEPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							bin_width_CCQEPlots[WhichSample+2][WhichPlot] = (TH1D*)(CCQEPlots[WhichSample+2][WhichPlot]->Clone());
							Reweight(bin_width_CCQEPlots[WhichSample+2][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_CCQEPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(bin_width_CCQEPlots[WhichSample+2][WhichPlot],"hist");

							CCMECPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							bin_width_CCMECPlots[WhichSample][WhichPlot] = (TH1D*)(CCMECPlots[WhichSample][WhichPlot]->Clone());
							Reweight(bin_width_CCMECPlots[WhichSample][WhichPlot]);
							CCMECPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							CCMECPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							bin_width_CCMECPlots[WhichSample+2][WhichPlot] = (TH1D*)(CCMECPlots[WhichSample+2][WhichPlot]->Clone());
							Reweight(bin_width_CCMECPlots[WhichSample+2][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_CCMECPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(bin_width_CCMECPlots[WhichSample+2][WhichPlot],"hist");

							CCRESPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							bin_width_CCRESPlots[WhichSample][WhichPlot] = (TH1D*)(CCRESPlots[WhichSample][WhichPlot]->Clone());
							Reweight(bin_width_CCRESPlots[WhichSample][WhichPlot]);
							CCRESPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							CCRESPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							bin_width_CCRESPlots[WhichSample+2][WhichPlot] = (TH1D*)(CCRESPlots[WhichSample+2][WhichPlot]->Clone());
							Reweight(bin_width_CCRESPlots[WhichSample+2][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_CCRESPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(bin_width_CCRESPlots[WhichSample+2][WhichPlot],"hist");

							CCDISPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							bin_width_CCDISPlots[WhichSample][WhichPlot] = (TH1D*)(CCDISPlots[WhichSample][WhichPlot]->Clone());
							Reweight(bin_width_CCDISPlots[WhichSample][WhichPlot]);
							CCDISPlots[WhichSample+2][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							CCDISPlots[WhichSample+2][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							bin_width_CCDISPlots[WhichSample+2][WhichPlot] = (TH1D*)(CCDISPlots[WhichSample+2][WhichPlot]->Clone());
							Reweight(bin_width_CCDISPlots[WhichSample+2][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_CCDISPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Add(bin_width_CCDISPlots[WhichSample+2][WhichPlot],"hist");

							TString NCCQEEvents = ToString((int)(CCQEPlots[WhichSample][WhichPlot]->Integral()));
							TString NCCMECEvents = ToString((int)(CCMECPlots[WhichSample][WhichPlot]->Integral()));	
							TString NCCRESEvents = ToString((int)(CCRESPlots[WhichSample][WhichPlot]->Integral()));	
							TString NCCDISEvents = ToString((int)(CCDISPlots[WhichSample][WhichPlot]->Integral()));	
							TString NExtBNBEvents = ToString((int)(Plots[2][WhichPlot]->Integral()));																						

							leg[WhichPlot]->AddEntry(CCQEPlots[WhichSample][WhichPlot],"QE (" + NCCQEEvents + ")","f"); 
							leg[WhichPlot]->AddEntry(CCMECPlots[WhichSample][WhichPlot],"MEC (" + NCCMECEvents + ")","f"); 
							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic (" + NExtBNBEvents + ")","f"); // ExtBNB
							leg[WhichPlot]->AddEntry(CCRESPlots[WhichSample][WhichPlot],"RES (" + NCCRESEvents + ")","f"); 
							leg[WhichPlot]->AddEntry(CCDISPlots[WhichSample][WhichPlot],"DIS (" + NCCDISEvents + ")","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}
				

				} // End of the loop over the samples

				// Unblind
				bin_width_Plots[0][WhichPlot]->Draw("e same");

				gPad->RedrawAxis();

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.07);

				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(0.07);
				TString PlotNameDuplicate = PlotNames[WhichPlot];
				TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
				textSlice->DrawLatexNDC(0.115, 0.8, LatexLabel[ReducedPlotName]);	

				// --------------------------------------------------------------------------------------------------------

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
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.51,1.49);
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(6);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.35);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleSize(0.1);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				hratio[0][WhichPlot]->Draw("e same");

				double RatioMin = hratio[0][WhichPlot]->GetXaxis()->GetXmin();
				double RatioMax = hratio[0][WhichPlot]->GetXaxis()->GetXmax();
				double YRatioCoord = 1.;
				TLine* RatioLine = new TLine(RatioMin,YRatioCoord,RatioMax,YRatioCoord);
				RatioLine->SetLineWidth(2);
				RatioLine->SetLineColor(kBlack);
				RatioLine->SetLineStyle(kDashed);
		
				topPad->cd();
				leg[WhichPlot]->SetTextSize(0.5);
				leg[WhichPlot]->SetTextFont(FontStyle);
				leg[WhichPlot]->Draw();

				// --------------------------------------------------------------------------------------

				// Sum of NonBeamOn Samples

				TH1D* SumNonBeamOn = (TH1D*)Plots[1][WhichPlot]->Clone(); // ExtBNB
				SumNonBeamOn->Add(Plots[2][WhichPlot]); // Overlay
				SumNonBeamOn->Add(Plots[3][WhichPlot]); // Dirt

				// CCQE Purity 

				int CCQEPurity = CCQEPlots[1][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex* latexPurity = new TLatex();
				latexPurity->SetTextFont(FontStyle);
				latexPurity->SetTextSize(0.07);
				TString LabelPurity = "QE = " + ToString(CCQEPurity/10.) + " %";
				latexPurity->DrawLatexNDC(0.61,0.89, LabelPurity);

				// --------------------------------------------------------------------------------------

				// Cosmic Contamination 

				int CosmicContamination = Plots[2][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexCosmic;
				latexCosmic.SetTextFont(FontStyle);
				latexCosmic.SetTextSize(0.07);
				TString LabelCosmic = "Cosmics = " + ToString(CosmicContamination/10.) + " %";
				latexCosmic.DrawLatexNDC(0.61,0.8, LabelCosmic);
				
				//----------------------------------------//
				
				// POT label
				
				TLatex *textPOT = new TLatex();
				textPOT->SetTextFont(FontStyle);
				textPOT->SetTextSize(0.07);
				textPOT->DrawLatexNDC(0.115, 0.89,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");					

				// -------------------------------------------------------------------- //				

				if ( string(PlotNames[WhichPlot]).find("RecoThetaVis") != std::string::npos ) {

					TLatex latexDataStats;
					latexDataStats.SetTextFont(FontStyle);
					latexDataStats.SetTextSize(0.07);
					double data_peak = FindOneDimHistoMaxValueBin(Plots[0][WhichPlot]);
					double data_mean = bin_width_Plots[0][WhichPlot]->GetMean();
					double data_std = bin_width_Plots[0][WhichPlot]->GetRMS();
					TString LabelDataStats = "#splitline{Data peak = " + to_string_with_precision(data_peak,2) + "}{#mu = " + to_string_with_precision(data_mean,2) + ", #sigma = " + to_string_with_precision(data_std,2) + "}";
					latexDataStats.DrawLatexNDC(0.61,0.6, LabelDataStats);				

					TH1D* MC = (TH1D*) (THStacks[WhichPlot]->GetStack()->Last());
					TH1D* clone_MC = (TH1D*)(MC->Clone());
					rm_bin_width(clone_MC);
					TLatex latexMCStats;
					latexMCStats.SetTextFont(FontStyle);
					latexMCStats.SetTextSize(0.07);
					double mc_peak = FindOneDimHistoMaxValueBin(MC);
					double mc_mean = MC->GetMean();
					double mc_std = MC->GetRMS();
					TString LabelMCStats = "#splitline{MC peak = " + to_string_with_precision(mc_peak,2) + "}{#mu = " + to_string_with_precision(mc_mean,2) + ", #sigma = " + to_string_with_precision(mc_std,2) + "}";
					latexMCStats.DrawLatexNDC(0.61,0.4, LabelMCStats);				


				}

				//----------------------------------------//

				// Uncertainty band

				int n = Plots[0][WhichPlot]->GetXaxis()->GetNbins();
				TString CopyPlotName = PlotNames[WhichPlot];
				// Total covariance matrix
				TH2D* CovMatrix = (TH2D*)(CovFile->Get("TotalCovariance_"+ReducedPlotName));
				// Clone the covariance matrix, so that you can scale it to the correct units (events vs flux averaged)
				TH2D* CovMatrixEvents = (TH2D*)CovMatrix->Clone();

				for (int i = 1; i <= n;i++ ) { 

					for (int j = 1; j <= n;j++ ) { 

						double bin_entry = CovMatrix->GetBinContent(i,j);
						// Scale the covariances to events, not flux averaged events as they are right now
						double scaled_bin_entry = bin_entry * TMath::Power( (IntegratedFlux*NTargets)/Units, 2);

						CovMatrixEvents->SetBinContent(i,j,scaled_bin_entry);					

					}				

				}

				// Statistical covariance matrix
				TH2D* StatCovMatrix = (TH2D*)(CovFile->Get("StatCovariance_"+ReducedPlotName));		
				// Statistical covariance needs to be removed from total
				CovMatrix->Add(StatCovMatrix,-1);		
				// Sanity check, stat errors should be identical to the ones coming from the Stat covariances 
				//TH2D* CovMatrix = (TH2D*)(CovFile->Get("StatCovariance_"+ReducedPlotName));				
				//CovMatrix->Scale(TMath::Power( (IntegratedFlux*NTargets)/Units ,2.));

				TH1D* MCUnc = (TH1D*)(Plots[0][WhichPlot]->Clone());				
				TH1D* MCStack = (TH1D*) (THStacks[WhichPlot]->GetStack()->Last());
				TH1D* MCStackClone = (TH1D*)(MCStack->Clone());
				rm_bin_width(MCStackClone);

				for (int i = 1; i <= n;i++ ) { 

					double MCCV = MCStackClone->GetBinContent(i);
					// Scale the covariances to events, not flux averaged events as they are right now
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
				//Plots[0][WhichPlot]->Draw("e same");				

				//gPad->RedrawAxis();								

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
				if (plot_unc) { THStacksMCUnc[WhichPlot]->Draw("same"); }

				THStacksMCUnc[WhichPlot]->Add(MCUncTwice,"hist");
				if (plot_unc) { THStacksMCUnc[WhichPlot]->Draw(" same"); }					

				RatioLine->Draw("same");
				hratio[0][WhichPlot]->Draw("e same");	
				gPad->RedrawAxis();														

				//----------------------------------------//

				// Chi2, p-value, sigma

				double chi2, pval, sigma; int ndof;
				
				CalcChiSquared(Plots[0][WhichPlot],MCStackClone,CovMatrixEvents,chi2,ndof,pval,sigma);
				TString Chi2Ndof = "#chi^{2}/ndof = " + to_string_with_precision(chi2,1) + "/" + TString(std::to_string(ndof)) +", p = " + to_string_with_precision(pval,2) + ", " + to_string_with_precision(sigma,2) + "#sigma";

				TLatex latexChi2;
				latexChi2.SetTextFont(FontStyle);
				latexChi2.SetTextSize(0.1);
				if (plot_unc) { latexChi2.DrawLatexNDC(0.15,0.88,Chi2Ndof); }

				// --------------------------------------------------------------------------------------

				TString CanvasPath = PlotPath + Cuts + "/InteractionBreakDown/";
				TString CanvasName = BaseMC + "THStack_BreakDown_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs	

} // End of the program 
