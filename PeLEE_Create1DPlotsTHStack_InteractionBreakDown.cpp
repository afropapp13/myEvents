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

void PeLEE_Create1DPlotsTHStack_InteractionBreakDown(TString BaseMC = "") {

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	

	// ----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> PlotNames; PlotNames.clear();

	PlotNames.push_back("RecoMuonMomentumPlot");
	PlotNames.push_back("RecoProtonMomentumPlot");
	PlotNames.push_back("RecoMuonCosThetaPlot");
	PlotNames.push_back("RecoProtonCosThetaPlot");
	PlotNames.push_back("RecoDeltaPTPlot");
	PlotNames.push_back("RecoDeltaAlphaTPlot");
	PlotNames.push_back("RecoDeltaAlpha3DqPlot");
	PlotNames.push_back("RecoDeltaPhiTPlot");
	PlotNames.push_back("RecoDeltaPhi3DPlot");	

	PlotNames.push_back("RecoDeltaPnPlot");
	PlotNames.push_back("RecoMuonCosThetaSingleBinPlot");	

	PlotNames.push_back("RecoECalPlot");

    PlotNames.push_back("RecoECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot");
    PlotNames.push_back("RecoECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot");

    PlotNames.push_back("RecoECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot");
    PlotNames.push_back("RecoECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot");		

    PlotNames.push_back("RecoECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot");
    PlotNames.push_back("RecoECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot");

    PlotNames.push_back("RecoSerialECal_DeltaPTDeltaAlphaTPlot");
    PlotNames.push_back("RecoSerialECal_MuonCosThetaMuonMomentumPlot");
    PlotNames.push_back("RecoSerialECal_ProtonCosThetaProtonMomentumPlot");					

	}

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	TString Cuts = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	//VectorCuts.push_back("");
	VectorCuts.push_back("_PID_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	// ------------------------------------------------------------------------------------------------------------------------------------------

	//vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
//	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -------------------------------------------------------------------------------------------------------------------------------------

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);		

		// -------------------------------------------------------------------------------------------------------------------------------------

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

			if (BaseMC == "") { NameOfSamples.push_back("STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }
			else if (BaseMC == "Overlay9NuWro") { NameOfSamples.push_back("STVStudies_Overlay9NuWro_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }			
			else if (BaseMC == "NoTuneOverlay9") { NameOfSamples.push_back("NoTuneSTVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }
			else if (BaseMC == "GENIEv2Overlay9") { NameOfSamples.push_back("GENIEv2STVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }			
			else if (BaseMC == "TwiceMECOverlay9") { NameOfSamples.push_back("TwiceMECSTVStudies_Overlay9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Overlay"); }			

			NameOfSamples.push_back("STVStudies_ExtBNB9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("ExtBNB");

			if (BaseMC != "NoTuneOverlay9") { NameOfSamples.push_back("STVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }
			else { NameOfSamples.push_back("NoTuneSTVStudies_OverlayDirt9_"+Runs[WhichRun]+Cuts+".root"); LabelsOfSamples.push_back("Dirt"); }

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
//				vector<TH1D*> CCCosmicCurrentPlots; CCCosmicCurrentPlots.clear();

				vector<TH1D*> Currenthratio;  Currenthratio.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* hist = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					TH1D* CCQEhist = (TH1D*)(FileSample[WhichSample]->Get("CCQE"+PlotNames[WhichPlot]));
					TH1D* CCMEChist = (TH1D*)(FileSample[WhichSample]->Get("CCMEC"+PlotNames[WhichPlot]));
					TH1D* CCREShist = (TH1D*)(FileSample[WhichSample]->Get("CCRES"+PlotNames[WhichPlot]));
					TH1D* CCDIShist = (TH1D*)(FileSample[WhichSample]->Get("CCDIS"+PlotNames[WhichPlot]));
//					TH1D* CCCosmichist = (TH1D*)(FileSample[WhichSample]->Get("Cosmic"+PlotNames[WhichPlot]));

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
						hist->SetMarkerSize(1.); 
					}

					CurrentPlots.push_back(hist);
					CCQECurrentPlots.push_back(CCQEhist);
					CCMECCurrentPlots.push_back(CCMEChist);
					CCRESCurrentPlots.push_back(CCREShist);
					CCDISCurrentPlots.push_back(CCDIShist);
//					CCCosmicCurrentPlots.push_back(CCCosmichist);

					Currenthratio.push_back((TH1D*)hist->Clone());
		
				} // End of the loop over the plots

				Plots.push_back(CurrentPlots);
				CCQEPlots.push_back(CCQECurrentPlots);
				CCMECPlots.push_back(CCMECCurrentPlots);
				CCRESPlots.push_back(CCRESCurrentPlots);
				CCDISPlots.push_back(CCDISCurrentPlots);
//				CCCCDISPlots.push_back(CCCosmicCurrentPlots);

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
				botPad->SetGridx();
				botPad->SetGridy();
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
//					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# Events / " + ToString(DataPOT));
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle("Number of " + Runs[WhichRun] + " events / bin");
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.08);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.6);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTickSize(0);

					double localmax = Plots[WhichSample][WhichPlot]->GetMaximum();
					if (localmax > max) { max = localmax; }
					Plots[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.3*max);

					if (LabelsOfSamples[WhichSample] == "BeamOn") { 

						gStyle->SetErrorX(0); // Removing the horizontal errors
						Plots[WhichSample][WhichPlot]->Draw("e same"); 
						TString NBeamOnEvents = ToString((int)(Plots[WhichSample][WhichPlot]->GetEntries()));
						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],"BNB Data ("+NBeamOnEvents+")","ep");
//						leg[WhichPlot]->AddEntry(Plots[WhichSample][WhichPlot],"BNB Data","ep");

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

							leg[WhichPlot]->AddEntry(CCQEPlots[WhichSample][WhichPlot],"QE (" + NCCQEEvents + ")","f"); 
							leg[WhichPlot]->AddEntry(CCMECPlots[WhichSample][WhichPlot],"MEC (" + NCCMECEvents + ")","f"); 
							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic (" + NExtBNBEvents + ")","f"); // ExtBNB
//							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic","f"); // ExtBNB
							leg[WhichPlot]->AddEntry(CCRESPlots[WhichSample][WhichPlot],"RES (" + NCCRESEvents + ")","f"); 
							leg[WhichPlot]->AddEntry(CCDISPlots[WhichSample][WhichPlot],"DIS (" + NCCDISEvents + ")","f"); 

// 							leg[WhichPlot]->AddEntry(CCQEPlots[WhichSample][WhichPlot],"QE","f"); 
// 							leg[WhichPlot]->AddEntry(CCMECPlots[WhichSample][WhichPlot],"MEC","f"); 
// //							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],LabelsOfSamples[2] + "","f"); // ExtBNB
// 							leg[WhichPlot]->AddEntry(Plots[2][WhichPlot],"Cosmic","f"); // ExtBNB
// 							leg[WhichPlot]->AddEntry(CCRESPlots[WhichSample][WhichPlot],"RES","f"); 
// 							leg[WhichPlot]->AddEntry(CCDISPlots[WhichSample][WhichPlot],"DIS","f"); 

							THStacks[WhichPlot]->Draw("same");
					
					}
				

				} // End of the loop over the samples

				Plots[0][WhichPlot]->Draw("e same");

				gPad->RedrawAxis();

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.09);
				//if (Runs[WhichRun] != "Combined") { text->DrawTextNDC(0.115, 0.9, Runs[WhichRun]); }

				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(0.09);
				TString PlotNameDuplicate = PlotNames[WhichPlot];
				TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
				if (Runs[WhichRun] != "Combined") {  textSlice->DrawLatexNDC(0.115, 0.8, LatexLabel[ReducedPlotName]);	}							

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
//				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.9*hratio[0][WhichPlot]->GetMinimum(),1.1*hratio[0][WhichPlot]->GetMaximum());
				hratio[0][WhichPlot]->GetYaxis()->SetRangeUser(0.1,1.9);
				hratio[0][WhichPlot]->GetYaxis()->SetNdivisions(6);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.35);
				hratio[0][WhichPlot]->GetYaxis()->SetTitleSize(0.1);
				hratio[0][WhichPlot]->GetYaxis()->SetLabelSize(0.11);

				botPad->cd();
				hratio[0][WhichPlot]->Draw("e same");

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

				// CCQE Purity 

				int CCQEPurity = CCQEPlots[1][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex* latexPurity = new TLatex();
				latexPurity->SetTextFont(FontStyle);
				latexPurity->SetTextSize(0.09);
				TString LabelPurity = "QE = " + ToString(CCQEPurity/10.) + " %";
				latexPurity->DrawLatexNDC(0.61,0.89, LabelPurity);

				// --------------------------------------------------------------------------------------

				// Cosmic Contamination 

				int CosmicContamination = Plots[2][WhichPlot]->Integral() / SumNonBeamOn->Integral() * 1000.;

				midPad->cd();
				TLatex latexCosmic;
				latexCosmic.SetTextFont(FontStyle);
				latexCosmic.SetTextSize(0.09);
				TString LabelCosmic = "Cosmics = " + ToString(CosmicContamination/10.) + " %";
				latexCosmic.DrawLatexNDC(0.61,0.8, LabelCosmic);
				
				//----------------------------------------//
				
				// POT label
				
				TLatex *textPOT = new TLatex();
				textPOT->SetTextFont(FontStyle);
				textPOT->SetTextSize(0.09);
				textPOT->DrawLatexNDC(0.115, 0.89,"MicroBooNE " + ToString(DataPOT).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT");					

				// --------------------------------------------------------------------------------------

				// Calculate the data / MC ratio

				if (PlotNames[WhichPlot] == "RecoProtonMomentumPlot") {

					cout << "Data / MC ratio = " << Plots[0][WhichPlot]->Integral() / SumNonBeamOn->Integral() << endl;

				}

				// --------------------------------------------------------------------------------------

				TString CanvasPath = PlotPath + Cuts + "/InteractionBreakDown/";
				TString CanvasName = BaseMC + "THStack_BreakDown_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Cuts+".pdf";
				PlotCanvas[WhichPlot]->SaveAs(CanvasPath+CanvasName);
				delete PlotCanvas[WhichPlot];

			} // End of the loop over the plots

		} // If we want to run on all cut combinations, include this } and remove the one at the beginning of the program

	} // End of the loop over the runs	

} // End of the program 
