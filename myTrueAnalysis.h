#ifndef myTrueAnalysis_h
#define myTrueAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <vector>
#include <vector>

#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Constants.h"

using namespace Constants;

class myTrueAnalysis {

private:
   TString fWhichSample;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   double          Weight;
   double          T2KWeight;
   int             CC1p;   
   Int_t           NumberMCParticles;

   vector<double>  *Muon_MCParticle_StartX;   
   vector<double>  *Muon_MCParticle_StartY;
   vector<double>  *Muon_MCParticle_StartZ;
   vector<double>  *Muon_MCParticle_EndX;   
   vector<double>  *Muon_MCParticle_EndY;
   vector<double>  *Muon_MCParticle_EndZ;         
   vector<double>  *Muon_MCParticle_Mom;
   vector<double>  *Muon_MCParticle_Phi;
   vector<double>  *Muon_MCParticle_CosTheta;
   vector<double>  *Muon_MCParticle_Length;
//   vector<int>     *Muon_MCParticle_StartContainment;
//   vector<int>     *Muon_MCParticle_EndContainment;
   vector<int>     *Muon_MCParticle_Pdg;

   vector<double>  *Proton_MCParticle_StartX;   
   vector<double>  *Proton_MCParticle_StartY;
   vector<double>  *Proton_MCParticle_StartZ;
   vector<double>  *Proton_MCParticle_EndX;   
   vector<double>  *Proton_MCParticle_EndY;
   vector<double>  *Proton_MCParticle_EndZ;   
   vector<double>  *Proton_MCParticle_Mom;
   vector<double>  *Proton_MCParticle_Phi;
   vector<double>  *Proton_MCParticle_CosTheta;
   vector<double>  *Proton_MCParticle_Length;
//   vector<int>     *Proton_MCParticle_StartContainment;
//   vector<int>     *Proton_MCParticle_EndContainment;
   vector<int>     *Proton_MCParticle_Pdg;   
   
   vector<double>  *True_Pt;
   vector<double>  *True_DeltaAlphaT;
   vector<double>  *True_DeltaPhiT;
   vector<double>  *True_ECal;
   vector<double>  *True_EQE;
   vector<double>  *True_Q2;         

   // List of branches
   TBranch        *b_Weight;   //!
   TBranch        *b_T2KWeight;   //!
   TBranch        *b_CC1p;   //!   
   TBranch        *b_NumberMCParticles;   //!

   TBranch        *b_Muon_MCParticle_StartX;   //!
   TBranch        *b_Muon_MCParticle_StartY;   //!
   TBranch        *b_Muon_MCParticle_StartZ;   //!
   TBranch        *b_Muon_MCParticle_EndX;   //!
   TBranch        *b_Muon_MCParticle_EndY;   //!
   TBranch        *b_Muon_MCParticle_EndZ;   //!            
   TBranch        *b_Muon_MCParticle_Mom;   //!
   TBranch        *b_Muon_MCParticle_Phi;   //!
   TBranch        *b_Muon_MCParticle_CosTheta;   //!
   TBranch        *b_Muon_MCParticle_Length;   //!
//   TBranch        *b_Muon_MCParticle_StartContainment;   //!
//   TBranch        *b_Muon_MCParticle_EndContainment;   //!
   TBranch        *b_Muon_MCParticle_Pdg;   //!

   TBranch        *b_Proton_MCParticle_StartX;   //!
   TBranch        *b_Proton_MCParticle_StartY;   //!
   TBranch        *b_Proton_MCParticle_StartZ;   //!
   TBranch        *b_Proton_MCParticle_EndX;   //!
   TBranch        *b_Proton_MCParticle_EndY;   //!
   TBranch        *b_Proton_MCParticle_EndZ;   //!   
   TBranch        *b_Proton_MCParticle_Mom;   //!
   TBranch        *b_Proton_MCParticle_Phi;   //!
   TBranch        *b_Proton_MCParticle_CosTheta;   //!
   TBranch        *b_Proton_MCParticle_Length;   //!
//   TBranch        *b_Proton_MCParticle_StartContainment;   //!
//   TBranch        *b_Proton_MCParticle_EndContainment;   //!
   TBranch        *b_Proton_MCParticle_Pdg;   //!   
   
   TBranch        *b_True_Pt;   //!
   TBranch        *b_True_DeltaAlphaT;   //!
   TBranch        *b_True_DeltaPhiT;   //!
   TBranch        *b_True_ECal;   //!
   TBranch        *b_True_EQE;   //!
   TBranch        *b_True_Q2;   //!                  

   myTrueAnalysis(TString WhichSample,TTree *tree=0);
   virtual ~myTrueAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif


#ifdef myTrueAnalysis_cxx
myTrueAnalysis::myTrueAnalysis(TString WhichSample, TTree *tree) : fChain(0) 
{

   fWhichSample = WhichSample;

// LOcally
//   TString PathToFile = "mySamples/"+UBCodeVersion+"/PreTruthSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

// On the gpvm's
   TString PathToFile = "mySamples/"+UBCodeVersion+"/PreTruthSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(PathToFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(PathToFile);
      }
      f->GetObject("myPreTruthSelection",tree);

   }
   Init(tree);
}

myTrueAnalysis::~myTrueAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myTrueAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myTrueAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myTrueAnalysis::Init(TTree *tree)
{

   // Set object pointer

   Muon_MCParticle_StartX = 0;
   Muon_MCParticle_StartY = 0;
   Muon_MCParticle_StartZ = 0;
   Muon_MCParticle_EndX = 0;
   Muon_MCParticle_EndY = 0;
   Muon_MCParticle_EndZ = 0;            
   Muon_MCParticle_Mom = 0;
   Muon_MCParticle_Phi = 0;
   Muon_MCParticle_CosTheta = 0;
   Muon_MCParticle_Length = 0;
//   Muon_MCParticle_StartContainment = 0;
//   Muon_MCParticle_EndContainment = 0;
   Muon_MCParticle_Pdg = 0;

   Proton_MCParticle_StartX = 0;
   Proton_MCParticle_StartY = 0;
   Proton_MCParticle_StartZ = 0;
   Proton_MCParticle_EndX = 0;
   Proton_MCParticle_EndY = 0;
   Proton_MCParticle_EndZ = 0;
   Proton_MCParticle_Mom = 0;
   Proton_MCParticle_Phi = 0;
   Proton_MCParticle_CosTheta = 0;
   Proton_MCParticle_Length = 0;
//   Proton_MCParticle_StartContainment = 0;
//   Proton_MCParticle_EndContainment = 0;
   Proton_MCParticle_Pdg = 0;
   
   True_Pt = 0;
   True_DeltaAlphaT = 0;
   True_DeltaPhiT = 0;
   True_ECal = 0;
   True_EQE = 0;
   True_Q2 = 0;                  

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("T2KWeight", &T2KWeight, &b_T2KWeight);
   fChain->SetBranchAddress("CC1p", &CC1p, &b_CC1p);   
   fChain->SetBranchAddress("NumberMCParticles", &NumberMCParticles, &b_NumberMCParticles);

   fChain->SetBranchAddress("Muon_MCParticle_StartX", &Muon_MCParticle_StartX, &b_Muon_MCParticle_StartX);
   fChain->SetBranchAddress("Muon_MCParticle_StartY", &Muon_MCParticle_StartY, &b_Muon_MCParticle_StartY);      
   fChain->SetBranchAddress("Muon_MCParticle_StartZ", &Muon_MCParticle_StartZ, &b_Muon_MCParticle_StartZ);
   fChain->SetBranchAddress("Muon_MCParticle_EndX", &Muon_MCParticle_EndX, &b_Muon_MCParticle_EndX);
   fChain->SetBranchAddress("Muon_MCParticle_EndY", &Muon_MCParticle_EndY, &b_Muon_MCParticle_EndY);      
   fChain->SetBranchAddress("Muon_MCParticle_EndZ", &Muon_MCParticle_EndZ, &b_Muon_MCParticle_EndZ);      
   fChain->SetBranchAddress("Muon_MCParticle_Mom", &Muon_MCParticle_Mom, &b_Muon_MCParticle_Mom);
   fChain->SetBranchAddress("Muon_MCParticle_Phi", &Muon_MCParticle_Phi, &b_Muon_MCParticle_Phi);
   fChain->SetBranchAddress("Muon_MCParticle_CosTheta", &Muon_MCParticle_CosTheta, &b_Muon_MCParticle_CosTheta);
   fChain->SetBranchAddress("Muon_MCParticle_Length", &Muon_MCParticle_Length, &b_Muon_MCParticle_Length);
//   fChain->SetBranchAddress("Muon_MCParticle_StartContainment", &Muon_MCParticle_StartContainment, &b_Muon_MCParticle_StartContainment);
//   fChain->SetBranchAddress("Muon_MCParticle_EndContainment", &Muon_MCParticle_EndContainment, &b_Muon_MCParticle_EndContainment);
   fChain->SetBranchAddress("Muon_MCParticle_Pdg", &Muon_MCParticle_Pdg, &b_Muon_MCParticle_Pdg);
   
   fChain->SetBranchAddress("Proton_MCParticle_StartX", &Proton_MCParticle_StartX, &b_Proton_MCParticle_StartX);
   fChain->SetBranchAddress("Proton_MCParticle_StartY", &Proton_MCParticle_StartY, &b_Proton_MCParticle_StartY);      
   fChain->SetBranchAddress("Proton_MCParticle_StartZ", &Proton_MCParticle_StartZ, &b_Proton_MCParticle_StartZ);
   fChain->SetBranchAddress("Proton_MCParticle_EndX", &Proton_MCParticle_EndX, &b_Proton_MCParticle_EndX);
   fChain->SetBranchAddress("Proton_MCParticle_EndY", &Proton_MCParticle_EndY, &b_Proton_MCParticle_EndY);      
   fChain->SetBranchAddress("Proton_MCParticle_EndZ", &Proton_MCParticle_EndZ, &b_Proton_MCParticle_EndZ);    
   fChain->SetBranchAddress("Proton_MCParticle_Mom", &Proton_MCParticle_Mom, &b_Proton_MCParticle_Mom);
   fChain->SetBranchAddress("Proton_MCParticle_Phi", &Proton_MCParticle_Phi, &b_Proton_MCParticle_Phi);
   fChain->SetBranchAddress("Proton_MCParticle_CosTheta", &Proton_MCParticle_CosTheta, &b_Proton_MCParticle_CosTheta);
   fChain->SetBranchAddress("Proton_MCParticle_Length", &Proton_MCParticle_Length, &b_Proton_MCParticle_Length);
//   fChain->SetBranchAddress("Proton_MCParticle_StartContainment", &Proton_MCParticle_StartContainment, &b_Proton_MCParticle_StartContainment);
//   fChain->SetBranchAddress("Proton_MCParticle_EndContainment", &Proton_MCParticle_EndContainment, &b_Proton_MCParticle_EndContainment);
   fChain->SetBranchAddress("Proton_MCParticle_Pdg", &Proton_MCParticle_Pdg, &b_Proton_MCParticle_Pdg);   
   
   fChain->SetBranchAddress("True_Pt", &True_Pt, &b_True_Pt);
   fChain->SetBranchAddress("True_DeltaAlphaT", &True_DeltaAlphaT, &b_True_DeltaAlphaT);
   fChain->SetBranchAddress("True_DeltaPhiT", &True_DeltaPhiT, &b_True_DeltaPhiT);
   fChain->SetBranchAddress("True_ECal", &True_ECal, &b_True_ECal);
   fChain->SetBranchAddress("True_EQE", &True_EQE, &b_True_EQE);      
   fChain->SetBranchAddress("True_Q2", &True_Q2, &b_True_Q2);   
   
   Notify();
}

Bool_t myTrueAnalysis::Notify()
{
   return kTRUE;
}

void myTrueAnalysis::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myTrueAnalysis::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef myTrueAnalysis_cxx
