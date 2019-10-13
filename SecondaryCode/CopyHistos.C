void CopyHistos(){

	TFile* file_new = new TFile("SmallhistosMCC84.root","recreate");
	TFile* file_old = TFile::Open("histosMCC84.root");

	TH2D* Muon_dEdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Muon_dEdx_vs_ResidualRange"));
	TH2D* Proton_dEdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Proton_dEdx_vs_ResidualRange"));
	TH2D* Muon_dQdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Muon_dQdx_vs_ResidualRange"));
	TH2D* Proton_dQdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Proton_dqdx_vs_ResidualRange"));

	TH2D* Trunc_Muon_dEdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Trunc_Muon_dEdx_vs_ResidualRange"));
	TH2D* Trunc_Proton_dEdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Trunc_Proton_dEdx_vs_ResidualRange"));
	TH2D* Trunc_Muon_dQdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Trunc_Muon_dQdx_vs_ResidualRange"));
	TH2D* Trunc_Proton_dQdx_vs_ResidualRange = (TH2D*)(file_old->Get("AfroFinder/Trunc_Proton_dqdx_vs_ResidualRange"));

/*TCanvas* canvas = new TCanvas();
canvas->cd();
Muon_dEdx_vs_ResidualRange->Draw("coltz");*/

file_new->cd();

Muon_dEdx_vs_ResidualRange->Write();
Proton_dEdx_vs_ResidualRange->Write();
Muon_dQdx_vs_ResidualRange->Write();
Proton_dQdx_vs_ResidualRange->Write();

Trunc_Muon_dEdx_vs_ResidualRange->Write();
Trunc_Proton_dEdx_vs_ResidualRange->Write();
Trunc_Muon_dQdx_vs_ResidualRange->Write();
Trunc_Proton_dQdx_vs_ResidualRange->Write();


file_new->Close();
};
