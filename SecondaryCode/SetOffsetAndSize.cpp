#include <TStyle.h>
#include <TROOT.h>

using namespace std;

void SetOffsetAndSize(){

	gStyle->SetTitleSize(0.07,"t");

	gStyle->SetTitleOffset(0.7,"x");
	gStyle->SetTitleSize(0.06,"x");

	gStyle->SetTitleOffset(0.7,"y");
	gStyle->SetTitleSize(0.06,"y");

	gStyle->SetStatX(0.9);                
	gStyle->SetStatY(0.9);  
	gStyle->SetStatH(0.2);

	gStyle->SetHistLineWidth(2);

	gStyle->SetOptStat(0);
	//gStyle->SetOptStat("m");

	gStyle->SetTitleAlign(23);

	gROOT->ForceStyle();

};
