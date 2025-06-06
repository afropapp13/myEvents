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

#include "../../../generators/constants.h"
#include "../../../generators/helper_functions.cxx"

using namespace std;
using namespace constants;

void mcc9_10_topological_breakdown(TString BaseMC = "") {

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);

	// -----------------------------------------------------------------------------------------------------------------------------------------

	std::vector<TString> PlotNames; PlotNames.clear();

	PlotNames.push_back("RecoSingleBinPlot");
	PlotNames.push_back("RecoPi0CosThetaPlot");
	PlotNames.push_back("RecoPi0MomentumPlot");
	PlotNames.push_back("Recog1CosThetaPlot");
	PlotNames.push_back("Recog1MomentumPlot");
	PlotNames.push_back("Recog2CosThetaPlot");
	PlotNames.push_back("Recog2MomentumPlot");
	PlotNames.push_back("Recotwo_shower_anglePlot");

	// Blips

	PlotNames.push_back("ReconBlips_radiusPlot");
	PlotNames.push_back("ReconBlips_savedPlot");
	PlotNames.push_back("RecoBlip_xPlot");
	PlotNames.push_back("RecoBlip_yPlot");
	PlotNames.push_back("RecoBlip_zPlot");
	PlotNames.push_back("RecoBlip_energyPlot");
	PlotNames.push_back("Recoblip_vrtPlot");
	PlotNames.push_back("Recoblip_cos_alphapi0Plot");
	PlotNames.push_back("Recoblip_cos_alphag1Plot");	
	PlotNames.push_back("Recoblip_cos_alphag2Plot");	
	
	PlotNames.push_back("Reconc_pio_scorePlot");
	PlotNames.push_back("Reconumu_scorePlot");
	PlotNames.push_back("Recokine_pio_flagPlot");
	PlotNames.push_back("Recokine_pio_vtx_disPlot");

	PlotNames.push_back("Recosingle_photon_numu_scorePlot");
	PlotNames.push_back("Recosingle_photon_other_scorePlot");
	PlotNames.push_back("Recosingle_photon_ncpi0_scorePlot");
	PlotNames.push_back("Recosingle_photon_nue_scorePlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_nocuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	//VectorCuts.push_back("_PID_NuScore_CRT");

	int NCuts = (int)(VectorCuts.size());	

	// -----------------------------------------------------------------------------------------------------------------------------------------

	int NRuns = (int)(xsec_Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));		

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// We needs these for the uncertainty band
		// event rate based
		/*TString NameExtractedXSec = MigrationMatrixPath+"ER_WienerSVD_Total_CovarianceMatrices_Overlay9_"+xsec_Runs[WhichRun]+".root";
		TFile* CovFile = new TFile(NameExtractedXSec,"readonly");		
*/
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(xsec_Runs[WhichRun]);
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);	
				
		// -----------------------------------------------------------------------------------------------------------------------------------------

		bool plot_unc = false;
		if (BaseMC == "") { plot_unc = true; }

		//-------------------------------------//

		Cuts = "_nocuts";

		for (int i = 0; i < NCuts; i++) {

			Cuts = Cuts + VectorCuts[i];

			// For the alternative MC, we want the figures after the application of all cuts
			if (BaseMC == "Overlay9NuWro" && i != NCuts-1) { continue; }	
			if (BaseMC == "GENIEv2Overlay9" && i != NCuts-1) { continue; }			
			if (BaseMC == "NoTuneOverlay9" && i != NCuts-1) { continue; }		
			if (BaseMC == "TwiceMECOverlay9" && i != NCuts-1) { continue; }		

			TString PathToFilesCut = event_selection_file_path+"/"+Cuts+"/";

			TH1D::SetDefaultSumw2();

			//--------------------------------------------------------------------------------

			vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
			vector<THStack*> THStacks; THStacks.clear();
			// Uncertainty band
			vector<THStack*> THStacksMCUnc; THStacksMCUnc.clear();			
			gStyle->SetPalette(55); const Int_t NCont = 999; 
			gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");
			vector<TLegend*> leg; leg.clear();

			vector<vector<TH1D*> > Plots; Plots.clear();
			vector<vector<TH1D*> > NCCOHPlots; NCCOHPlots.clear();
			vector<vector<TH1D*> > NonNCCOHPlots; NonNCCOHPlots.clear();

			vector<vector<TH1D*> > bin_width_Plots; bin_width_Plots.clear();
			vector<vector<TH1D*> > bin_width_NCCOHPlots; bin_width_NCCOHPlots.clear();
			vector<vector<TH1D*> > bin_width_NonNCCOHPlots; bin_width_NonNCCOHPlots.clear();

			vector<vector<TH1D*> > hratio;  hratio.clear();

			vector<TString> LabelsOfSamples;
			vector<TString> NameOfSamples;

			// 0: BeamOn
			// 1: Overlay
			// 2: ExtBNB
			// 3: Dirt
	
			// 0	
			NameOfSamples.push_back("ncpi0_mcc9_10_BeamOn9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("BeamOn");

			// 1
			if (BaseMC == "") { NameOfSamples.push_back("ncpi0_mcc9_10_Overlay9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }
			else if (BaseMC == "Overlay9NuWro") { NameOfSamples.push_back("ncpi0_mcc9_10_Overlay9NuWro_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }
			else if (BaseMC == "GENIEv2Overlay9") { NameOfSamples.push_back("GENIEv2ncpi0_mcc9_10_Overlay9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }		
			else if (BaseMC == "NoTuneOverlay9") { NameOfSamples.push_back("NoTunencpi0_mcc9_10_Overlay9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }
			else if (BaseMC == "TwiceMECOverlay9") { NameOfSamples.push_back("TwiceMECncpi0_mcc9_10_Overlay9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("MC"); }			

			// 2
			NameOfSamples.push_back("ncpi0_mcc9_10_ExtBNB9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");

			// 3
			if (BaseMC == "NoTuneOverlay9") { NameOfSamples.push_back("NoTunencpi0_mcc9_10_OverlayDirt9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
			else if (BaseMC == "TwiceMECOverlay9") { NameOfSamples.push_back("TwiceMECncpi0_mcc9_10_OverlayDirt9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
			else { NameOfSamples.push_back("ncpi0_mcc9_10_OverlayDirt9_"+xsec_Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
	
			vector<int> Colors; Colors.clear(); 
			// Unblind
			Colors.push_back(kBlack);
			// Blind 
			//Colors.push_back(kWhite); 
			Colors.push_back(kRed); Colors.push_back(kGray+2); Colors.push_back(kMagenta);

//			vector<int> ColorsOverlay{kBlue-5,kYellow+1,kOrange+7,kRed+1,kBlue};
			vector<int> ColorsOverlay{kMagenta,kAzure+7,kOrange-3,kGreen+1,kRed+1,kBlue};

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				FileSample.push_back(TFile::Open(PathToFilesCut+NameOfSamples[WhichSample]));

				vector<TH1D*> CurrentPlots; CurrentPlots.clear();
				vector<TH1D*> NCCOHCurrentPlots; NCCOHCurrentPlots.clear();
				vector<TH1D*> NonNCCOHCurrentPlots; NonNCCOHCurrentPlots.clear();

				vector<TH1D*> Currenthratio;  Currenthratio.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					TH1D* NCCOHhist = (TH1D*)(FileSample[WhichSample]->Get("NCCOH"+PlotNames[WhichPlot]));
					TH1D* NonNCCOHhist = (TH1D*)(FileSample[WhichSample]->Get("NonNCCOH"+PlotNames[WhichPlot]));

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
					NCCOHCurrentPlots.push_back(NCCOHhist);
					NonNCCOHCurrentPlots.push_back(NonNCCOHhist);

					Currenthratio.push_back((TH1D*)hist->Clone());
			
				}

				Plots.push_back(CurrentPlots);
				NCCOHPlots.push_back(NCCOHCurrentPlots);
				NonNCCOHPlots.push_back(NonNCCOHCurrentPlots);

				bin_width_Plots.push_back(CurrentPlots);
				bin_width_NCCOHPlots.push_back(NCCOHCurrentPlots);
				bin_width_NonNCCOHPlots.push_back(NonNCCOHCurrentPlots);

				hratio.push_back(Currenthratio);

			}
			
			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
		
				TString PlotCanvasName = xsec_Runs[WhichRun]+"_"+PlotNames[WhichPlot]+Cuts;
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

				leg.push_back(new TLegend(0.1,0.005,0.9,0.995));
				leg[WhichPlot]->SetBorderSize(0);
				leg[WhichPlot]->SetNColumns(3);
				leg[WhichPlot]->SetMargin(0.15);				

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

					if (xsec_Runs[WhichRun] == "Combined") {

						Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("Number of  events / bin");

					} else {

						Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle(xsec_Runs[WhichRun] + " events / bin");

					}

					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.65);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0.01);
		
					if (WhichSample == 0) { 

						bin_width_Plots[0][WhichPlot] = (TH1D*)(Plots[0][WhichPlot]->Clone()); 
						divide_bin_width(bin_width_Plots[0][WhichPlot]); 
						//max = find_bin_max_value(bin_width_Plots[0][WhichPlot]);
						//bin_width_Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.47*max);

					}
	
					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						bin_width_Plots[WhichSample][WhichPlot]->Draw("e same"); 
						TString NBeamOnEvents = ToString((int)(Plots[WhichSample][WhichPlot]->Integral()));
						// Unblind
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot], "BNB Data ("+NBeamOnEvents+")","ep");

					}

					if (LabelsOfSamples[WhichSample] == "ExtBNB") {

							Plots[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
							Plots[WhichSample][WhichPlot]->SetFillColor(Colors[WhichSample]);
							Plots[WhichSample][WhichPlot]->SetFillStyle(3004);
							Plots[WhichSample][WhichPlot]->SetLineWidth(1);

							bin_width_Plots[WhichSample][WhichPlot] = (TH1D*)(Plots[WhichSample][WhichPlot]->Clone());
							divide_bin_width(bin_width_Plots[WhichSample][WhichPlot]);
	
							TString NExtBNBEvents = ToString( (int)(Plots[WhichSample][WhichPlot]->Integral() ) );
							THStacks[WhichPlot]->Add(bin_width_Plots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Draw("same");

					}

					if (LabelsOfSamples[WhichSample] == "Dirt") {

							TString NNCCOHEvents = ToString( (int)(NCCOHPlots[WhichSample][WhichPlot]->Integral() ) );
							NCCOHPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[2]);
							NCCOHPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[2]);
							bin_width_NCCOHPlots[WhichSample][WhichPlot] = (TH1D*)(NCCOHPlots[WhichSample][WhichPlot]->Clone());
							divide_bin_width(bin_width_NCCOHPlots[WhichSample][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_NCCOHPlots[WhichSample][WhichPlot],"hist");
							THStacks[WhichPlot]->Draw("same");

							TString NNonNCCOHEvents = ToString( (int)(NonNCCOHPlots[WhichSample][WhichPlot]->Integral() ) );
							NonNCCOHPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[3]);
							NonNCCOHPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[3]);
							bin_width_NonNCCOHPlots[WhichSample][WhichPlot] = (TH1D*)(NonNCCOHPlots[WhichSample][WhichPlot]->Clone());
							divide_bin_width(bin_width_NonNCCOHPlots[WhichSample][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_NonNCCOHPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(NonNCCOHPlots[WhichSample][WhichPlot],"Out-of-cryo ("+NNonNCCOHEvents+")","f");
							THStacks[WhichPlot]->Draw("same");

					}

					if (LabelsOfSamples[WhichSample] == "MC") {

							TString NNCCOHEvents = ToString( (int)(NCCOHPlots[WhichSample][WhichPlot]->Integral() ) );
							NCCOHPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[0]);
							NCCOHPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[0]);
							bin_width_NCCOHPlots[WhichSample][WhichPlot] = (TH1D*)(NCCOHPlots[WhichSample][WhichPlot]->Clone());
							divide_bin_width(bin_width_NCCOHPlots[WhichSample][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_NCCOHPlots[WhichSample][WhichPlot],"hist");

							// add the cosmic label first
							TString NExtBNBEvents = ToString( (int)(Plots[2][WhichPlot]->Integral() ) );							
							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic ("+NExtBNBEvents+")","f");							
							// Unblind
							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"","");	 // blank space

							leg[WhichPlot]->AddEntry(NCCOHPlots[WhichSample][WhichPlot],"MC NCCOH-like ("+NNCCOHEvents+")","f");
							THStacks[WhichPlot]->Draw("same");

							TString NNonNCCOHEvents = ToString( (int)(NonNCCOHPlots[WhichSample][WhichPlot]->Integral() ) );
							NonNCCOHPlots[WhichSample][WhichPlot]->SetLineColor(ColorsOverlay[1]);
							NonNCCOHPlots[WhichSample][WhichPlot]->SetFillColor(ColorsOverlay[1]);
							bin_width_NonNCCOHPlots[WhichSample][WhichPlot] = (TH1D*)(NonNCCOHPlots[WhichSample][WhichPlot]->Clone());
							divide_bin_width(bin_width_NonNCCOHPlots[WhichSample][WhichPlot]);
							THStacks[WhichPlot]->Add(bin_width_NonNCCOHPlots[WhichSample][WhichPlot],"hist");
							leg[WhichPlot]->AddEntry(NonNCCOHPlots[WhichSample][WhichPlot],"MC nonNCCOH-like ("+NNonNCCOHEvents+")","f");
							THStacks[WhichPlot]->Draw("same");

					}
					

				} // End of the loop over the samples

				TH1D* stack_max = (TH1D*) (THStacks[WhichPlot]->GetStack()->Last());
				double m_stack = TMath::Max( find_bin_max_value(bin_width_Plots[0][WhichPlot]) , find_bin_max_value(stack_max) );
				bin_width_Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.35*m_stack);

				// Area normalize bc we don't have POT
				//bin_width_Plots[0][WhichPlot]->Scale( stack_max->Integral("width") / bin_width_Plots[0][WhichPlot]->Integral("width") );
				// Unblind and draw on top
				bin_width_Plots[0][WhichPlot]->Draw("e same"); 
				
				// -----------------------------------------------------------------------------------	

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

				// ------------------------------------------------------------------------------------

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

				// NCCOH Purity 

				int NCCOHPurity = NCCOHPlots[1][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				
				TLatex latexPurity;
				latexPurity.SetTextFont(FontStyle);
				latexPurity.SetTextSize(0.07);
				TString LabelPurity = "NCCOH-like = " + ToString(NCCOHPurity/10.) + " %";
				latexPurity.DrawLatexNDC(0.61,0.89, LabelPurity);
				
				//----------------------------------------//
				
				// POT label
				
				TLatex *textPOT = new TLatex();
				textPOT->SetTextFont(FontStyle);
				textPOT->SetTextSize(0.07);

				if (xsec_Runs[WhichRun] == "Combined") { 

					textPOT->DrawLatexNDC(0.115, 0.89,"MicroBooNE 1.30 #times 10^{21} POT");

				} else {
								
					textPOT->DrawLatexNDC(0.115, 0.89,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");								
				}

				//----------------------------------------//

				// Cosmic Contamination

				int CosmicContamination = Plots[2][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexCosmic;
				latexCosmic.SetTextFont(FontStyle);
				latexCosmic.SetTextSize(0.07);
				TString LabelCosmic = "Cosmics = " + ToString(CosmicContamination/10.) + " %";
				latexCosmic.DrawLatexNDC(0.61,0.8, LabelCosmic);				

				// -------------------------------------------------------------------- //				

				if ( string(PlotNames[WhichPlot]).find("RecoThetaVis") != std::string::npos ) {

					TLatex latexDataStats;
					latexDataStats.SetTextFont(FontStyle);
					latexDataStats.SetTextSize(0.07);
					double data_peak = find_bin_max_value(bin_width_Plots[0][WhichPlot]);
					double data_mean = bin_width_Plots[0][WhichPlot]->GetMean();
					double data_std = bin_width_Plots[0][WhichPlot]->GetRMS();
					TString LabelDataStats = "#splitline{Data peak = " + to_string_with_precision(data_peak,2) + "}{#mu = " + to_string_with_precision(data_mean,2) + ", #sigma' = " + to_string_with_precision(data_std,2) + "}";
					//latexDataStats.DrawLatexNDC(0.61,0.6, LabelDataStats);				

					TH1D* MC = (TH1D*) (THStacks[WhichPlot]->GetStack()->Last());
					TH1D* clone_MC = (TH1D*)(MC->Clone());
					rm_bin_width(clone_MC);
					TLatex latexMCStats;
					latexMCStats.SetTextFont(FontStyle);
					latexMCStats.SetTextSize(0.07);
					double mc_peak = find_bin_max_value(MC);
					double mc_mean = MC->GetMean();
					double mc_std = MC->GetRMS();
					TString LabelMCStats = "#splitline{MC peak = " + to_string_with_precision(mc_peak,2) + "}{#mu = " + to_string_with_precision(mc_mean,2) + ", #sigma' = " + to_string_with_precision(mc_std,2) + "}";
					//latexMCStats.DrawLatexNDC(0.61,0.4, LabelMCStats);				

				}

				//----------------------------------------//
/*
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
						//double scaled_bin_entry = bin_entry * TMath::Power( (IntegratedFlux*NTargets)/Units, 2);

						// no scaling needed for event rates
						double scaled_bin_entry = bin_entry;


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
				TH1D* MCStack = (TH1D*) (THStacks[WhichPlot]->GetStack())->Last();
				
				//TH1D* MCStackClone = (TH1D*)(MCStack->Clone());
				//rm_bin_width(MCStackClone);

				TH1D* MCStackClone = (TH1D*)(Plots[1][WhichPlot]->Clone()); // Overlay
				MCStackClone->Add(Plots[2][WhichPlot]); // ExtBNB
				MCStackClone->Add(Plots[3][WhichPlot]); // Dirt

				for (int i = 1; i <= n;i++ ) { 

					double MCCV = MCStackClone->GetBinContent(i);
					// Scale the covariances to events, not flux averaged events as they are right now
					//double Unc = TMath::Sqrt( CovMatrix->GetBinContent(i,i) ) * (IntegratedFlux*NTargets)/Units;

					//no scaling needed for event rates
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
				TString Chi2Ndof = "#chi^{2}/ndf = " + to_string_with_precision(chi2,1) + "/" + TString(std::to_string(ndof)) +", p = " + to_string_with_precision(pval,2) + ", " + to_string_with_precision(sigma,2) + "#sigma'";

				TLatex latexChi2;
				latexChi2.SetTextFont(FontStyle);
				latexChi2.SetTextSize(0.1);
				if (plot_unc) { latexChi2.DrawLatexNDC(0.15,0.88,Chi2Ndof); }				

				//----------------------------------------//

				// Plot vertical lines
				// Add latex label with phase space limits

				if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

					TString clone_name = PlotNames[WhichPlot];
					clone_name.ReplaceAll("Reco","");
					vector<int> bin_break_points = get_2d_bin_break_points( map_to_2d_bin.at(clone_name) );

					int nbreaks = bin_break_points.size() - 1;
					vector<TLine*> line; line.resize(nbreaks);

					for (int ipoint = 0; ipoint < nbreaks; ipoint ++) {

						line.at(ipoint) = new TLine( bin_break_points.at(ipoint) + 0.5,0., bin_break_points.at(ipoint) + 0.5, bin_width_Plots[0][WhichPlot]->GetMaximum() );
						midPad->cd();
						line.at(ipoint)->SetLineStyle(kDashed);
						line.at(ipoint)->Draw("same");

					}
	
					//----------------------------------------//

					vector<TLatex*> slice; slice.resize(nbreaks+1);

					for (int ipoint = 0; ipoint < nbreaks + 1; ipoint ++) {

			
						slice.at(ipoint) = new TLatex();
						slice.at(ipoint)->SetTextFont(FontStyle);
						slice.at(ipoint)->SetTextSize(0.04);
						TString phase_space = MapUncorCor[ clone_name + "_" + TString(std::to_string(ipoint) ) ];
						midPad->cd();
						if (ipoint == 0) { slice.at(ipoint)->DrawLatex( bin_break_points.at(ipoint) / 3. , 0.7 * bin_width_Plots[0][WhichPlot]->GetMaximum(), LatexLabel[phase_space ]); }
						else { slice.at(ipoint)->DrawLatex( bin_break_points.at(ipoint - 1) + ( bin_break_points.at(ipoint) - bin_break_points.at(ipoint-1) ) / 3. , 0.7 * bin_width_Plots[0][WhichPlot]->GetMaximum(), LatexLabel[phase_space ]); }


					}

				}
*/
				//----------------------------------------//

				TString CanvasPath = plot_path + Cuts+"/TopologicalBreakDown/";
				TString CanvasName = BaseMC + "mcc9_10_THStack_BreakDown_"+PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun]+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

				//----------------------------------------//

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs

} // End of the program 
