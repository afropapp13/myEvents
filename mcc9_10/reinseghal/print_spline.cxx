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

void print_spline() {

    TFile* f = new TFile("rs_spline.root","readonly");
    TGraph* g = (TGraph*)(f->Get("h_spline"));

    //--------------------//

    TCanvas* cg = new TCanvas("","");
    cg->SetBottomMargin(0.15);
    cg->SetLeftMargin(0.15);    

    g->GetXaxis()->SetRangeUser(0.,2.9);
    g->GetXaxis()->CenterTitle();
    g->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    g->GetXaxis()->SetTitleSize(0.06);
    g->GetXaxis()->SetTitleFont(132);  
    g->GetXaxis()->SetLabelSize(0.06);
    g->GetXaxis()->SetLabelFont(132);
    g->GetXaxis()->SetNdivisions(8);  
    
    g->GetYaxis()->CenterTitle();
    g->GetYaxis()->SetTitleSize(0.06);
    g->GetYaxis()->SetTitleFont(132);  
    g->GetYaxis()->SetLabelSize(0.06);
    g->GetYaxis()->SetLabelFont(132);
    g->GetYaxis()->SetNdivisions(8);  
    g->GetYaxis()->SetTitle("ratio RS/BS");        

    g->SetTitle("");    
    g->SetLineColor(kGreen+2);
    g->SetMarkerColor(kGreen+2);   
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);      
    g->Draw("AP");

    cg->SaveAs("xsec_graph_ratio.pdf");

    f->Close();
 
}