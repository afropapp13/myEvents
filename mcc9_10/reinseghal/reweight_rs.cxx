#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

void reweight_rs() {

    TFile* f_bs = new TFile("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_6_0_G18_10a_02_11a/14_1000180400_NC_v3_6_0_G18_10a_02_11a.xml.root","readonly");
    TFile* f_rs = new TFile("/exp/uboone/data/users/apapadop/gLEE/NCPi0/generators/geniev3_reinsehgal/14_1000180400_NC_v3_6_0_G18_10a_02_11a.xml.root","readonly");     

    TGraph* g_bs = (TGraph*)(f_bs->Get("nu_mu_Ar40/coh_nc"));
    TGraph* g_rs = (TGraph*)(f_rs->Get("nu_mu_Ar40/coh_nc")); 

    //--------------------//

    TCanvas* cg = new TCanvas("","");
    cg->SetBottomMargin(0.15);
    cg->SetLeftMargin(0.15);    

    g_rs->GetXaxis()->SetRangeUser(0.,2.9);
    g_rs->GetXaxis()->CenterTitle();
    g_rs->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    g_rs->GetXaxis()->SetTitleSize(0.06);
    g_rs->GetXaxis()->SetTitleFont(132);  
    g_rs->GetXaxis()->SetLabelSize(0.06);
    g_rs->GetXaxis()->SetLabelFont(132);
    g_rs->GetXaxis()->SetNdivisions(8);  
    
    g_rs->GetYaxis()->CenterTitle();
    g_rs->GetYaxis()->SetTitleSize(0.06);
    g_rs->GetYaxis()->SetTitleFont(132);  
    g_rs->GetYaxis()->SetLabelSize(0.06);
    g_rs->GetYaxis()->SetLabelFont(132);
    g_rs->GetYaxis()->SetNdivisions(8);  
    g_rs->GetYaxis()->SetTitle("#sigma_{nccoh} (10^{-38} cm^{2})");        

    g_rs->SetTitle("");    
    g_rs->SetLineColor(kGreen+2);
    g_rs->SetMarkerColor(kGreen+2);   
    g_rs->SetMarkerStyle(20);
    g_rs->SetMarkerSize(2);            
    g_rs->Draw("A*");
 
    g_bs->SetLineColor(kOrange+7);
    g_bs->SetMarkerColor(kOrange+7);   
    g_bs->SetMarkerStyle(20);
    g_bs->SetMarkerSize(2);            
    g_bs->Draw("* same");    

    TLegend* leg = new TLegend(0.2,0.65,0.4,0.85);
	TLegendEntry* lrs = leg->AddEntry(g_rs,"RS","l");
	lrs->SetTextColor(kGreen+2);
	TLegendEntry* lbs = leg->AddEntry(g_bs,"BS","l");
	lbs->SetTextColor(kOrange+7);    	

    leg->SetTextFont(132);
    leg->SetTextSize(0.06);  
	leg->SetBorderSize(0);    
    leg->Draw();

    cg->SaveAs("xsec_graphs.pdf");

    //--------------------//    

    double min_e = 0.;
    double max_e = 3.;    
    int np_bs = g_bs->GetN();
    double step = (max_e - min_e) / np_bs;

    double x[np_bs];
    double y[np_bs];    

    for(int i=0; i < np_bs; ++i) {

        double e_bs = g_bs->Eval(min_e + i * step);
        double e_rs = g_rs->Eval(min_e + i * step);  
        double r = e_rs / e_bs;
        if ( fabs(r) != r ) { r = 0; }
    
        x[i] = min_e + i * step;
        y[i] = r; 

    } 
    
    TGraph* h_spline = new TGraph(np_bs,x,y);       

    TFile* f_out = new TFile("rs_spline.root","recreate");
    h_spline->Write("h_spline");

    f_out->Close();

    //--------------------//        
  
}