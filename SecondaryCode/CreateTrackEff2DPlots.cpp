void CreateTrackEff2DPlots() {

	const int NSamples = 2, N2DPlots = 1; TString NameOfSamples[NSamples] = {"Overlay","MCC84"/*,"EXTUNB","BNBPlusEXTBNB"*/};
	TString LabelsOfSamples[NSamples] = {"BNB GENIE & EXT UNBIASED","BNB GENIE & COSMIC CORSIKA"/*,"EXT UNBIASED","BNB + EXTBNB"*/};
	TString NameOfPlots[N2DPlots] = {"NRecoTracks_vs_NPFTracks"};
	
	TCanvas* PlotCanvas[NSamples][N2DPlots] = {}; TH2D* Plots[NSamples][N2DPlots] = {}; TFile* FileSample[NSamples] = {};
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();
	
	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++){

		FileSample[WhichSample] = TFile::Open("./myFiles/histos"+NameOfSamples[WhichSample]+".root");

		for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++){

			PlotCanvas[WhichSample][WhichPlot] = new TCanvas(NameOfPlots[WhichPlot]+NameOfSamples[WhichSample],NameOfPlots[WhichPlot]+NameOfSamples[WhichSample],205,34,1024,768);
			PlotCanvas[WhichSample][WhichPlot]->cd();
			Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample]->Get(NameOfPlots[WhichPlot]+NameOfSamples[WhichSample]));
			Plots[WhichSample][WhichPlot]->SetTitle(LabelsOfSamples[WhichSample]);
			CenterAxisTitle(Plots[WhichSample][WhichPlot]);
			PlotCanvas[WhichSample][WhichPlot]->SetLogz();
			Plots[WhichSample][WhichPlot]->Draw("colz");
			PlotCanvas[WhichSample][WhichPlot]->SaveAs("./myPlots/pdf/2D/"+NameOfPlots[WhichPlot]+NameOfSamples[WhichSample]+".pdf");
			PlotCanvas[WhichSample][WhichPlot]->SaveAs("/home/afroditi/Dropbox/PhD/MIT_TAU/MITAU_GENIE_Adi/cosmicOverlay/note/plots/Afro/"
								   +NameOfPlots[WhichPlot]+NameOfSamples[WhichSample]+".pdf");
			PlotCanvas[WhichSample][WhichPlot]->SaveAs("/home/afroditi/Dropbox/PhD/LaTeX/Graduate_Studies/PhD_Dissertation/Pictures/MCC84VsOverlay/"
								   +NameOfPlots[WhichPlot]+NameOfSamples[WhichSample]+".pdf");
			PlotCanvas[WhichSample][WhichPlot]->SaveAs("./myPlots/eps/2D/"+NameOfPlots[WhichPlot]+NameOfSamples[WhichSample]+".eps");
			//delete PlotCanvas[WhichSample][WhichPlot];
		};
	};
};
