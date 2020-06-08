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
   vector<double>  *MCParticle_Mom;
   vector<double>  *MCParticle_Phi;
   vector<double>  *MCParticle_CosTheta;
   vector<double>  *MCParticle_Length;
   vector<int>     *MCParticle_StartContainment;
   vector<int>     *MCParticle_EndContainment;
   vector<int>     *MCParticle_Pdg;

   // List of branches
   TBranch        *b_Weight;   //!
   TBranch        *b_T2KWeight;   //!
   TBranch        *b_CC1p;   //!   
   TBranch        *b_NumberMCParticles;   //!
   TBranch        *b_MCParticle_Mom;   //!
   TBranch        *b_MCParticle_Phi;   //!
   TBranch        *b_MCParticle_CosTheta;   //!
   TBranch        *b_MCParticle_Length;   //!
   TBranch        *b_MCParticle_StartContainment;   //!
   TBranch        *b_MCParticle_EndContainment;   //!
   TBranch        *b_MCParticle_Pdg;   //!

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
   MCParticle_Mom = 0;
   MCParticle_Phi = 0;
   MCParticle_CosTheta = 0;
   MCParticle_Length = 0;
   MCParticle_StartContainment = 0;
   MCParticle_EndContainment = 0;
   MCParticle_Pdg = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("T2KWeight", &T2KWeight, &b_T2KWeight);
   fChain->SetBranchAddress("CC1p", &CC1p, &b_CC1p);   
   fChain->SetBranchAddress("NumberMCParticles", &NumberMCParticles, &b_NumberMCParticles);
   fChain->SetBranchAddress("MCParticle_Mom", &MCParticle_Mom, &b_MCParticle_Mom);
   fChain->SetBranchAddress("MCParticle_Phi", &MCParticle_Phi, &b_MCParticle_Phi);
   fChain->SetBranchAddress("MCParticle_CosTheta", &MCParticle_CosTheta, &b_MCParticle_CosTheta);
   fChain->SetBranchAddress("MCParticle_Length", &MCParticle_Length, &b_MCParticle_Length);
   fChain->SetBranchAddress("MCParticle_StartContainment", &MCParticle_StartContainment, &b_MCParticle_StartContainment);
   fChain->SetBranchAddress("MCParticle_EndContainment", &MCParticle_EndContainment, &b_MCParticle_EndContainment);
   fChain->SetBranchAddress("MCParticle_Pdg", &MCParticle_Pdg, &b_MCParticle_Pdg);
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
