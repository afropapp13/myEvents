#ifndef myAnalysis_h
#define myAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <vector>
#include <vector>

#include "../myCCQEAnalysis/Constants.h"

TString WhichSample = "Run1Data9";
//TString WhichSample = "Overlay9";
//TString WhichSample = "ExtBNB9";
//TString WhichSample = "OverlayDirt9";

//TString WhichSample = "Overlay9_SCE";
//TString WhichSample = "Overlay9_DLdown";

using namespace Constants;


TString PathToFile = "mySamples/"+UBCodeVersion+"/PreSelection_"+WhichSample+"_"+UBCodeVersion+".root";


class myAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
//   Int_t           PassedSwTrigger;
   double           Weight;
   Int_t           NBeamFlashes;
   vector<double>  *BeamFlashes_YCenter;
   vector<double>  *BeamFlashes_ZCenter;
   vector<double>  *BeamFlashes_TotalPE;
   Int_t           NumberMCParticles;
   vector<double>  *CandidateMuP_Distance;
   vector<double>  *Vertex_X;
   vector<double>  *Vertex_Y;
   vector<double>  *Vertex_Z;
   vector<double>  *CandidateMu_P;
   vector<double>  *CandidateMu_Phi;
   vector<double>  *CandidateMu_CosTheta;
//   vector<double>  *CandidateMu_Length;
   vector<double>  *CandidateMu_Chi2_YPlane;
   vector<double>  *CandidateMu_ThreePlaneChi2;
   vector<int>     *CandidateMu_StartContainment;
   vector<int>     *CandidateMu_EndContainment;
   vector<int>     *CandidateMu_MCParticle_Pdg;
   vector<double>  *CandidateMu_MCParticle_Purity;
   vector<double>  *True_CandidateMu_P;
   vector<double>  *True_CandidateMu_Phi;
   vector<double>  *True_CandidateMu_CosTheta;
//   vector<double>  *True_CandidateMu_Length;
   vector<int>     *True_CandidateMu_StartContainment;
   vector<int>     *True_CandidateMu_EndContainment;
   vector<double>  *CandidateP_P;
   vector<double>  *CandidateP_Phi;
   vector<double>  *CandidateP_CosTheta;
//   vector<double>  *CandidateP_Length;
   vector<double>  *CandidateP_Chi2_YPlane;
   vector<double>  *CandidateP_ThreePlaneChi2;
   vector<int>     *CandidateP_StartContainment;
   vector<int>     *CandidateP_EndContainment;
   vector<int>     *CandidateP_MCParticle_Pdg;
   vector<double>  *CandidateP_MCParticle_Purity;
   vector<double>  *True_CandidateP_P;
   vector<double>  *True_CandidateP_Phi;
   vector<double>  *True_CandidateP_CosTheta;
//   vector<double>  *True_CandidateP_Length;
   vector<int>     *True_CandidateP_StartContainment;
   vector<int>     *True_CandidateP_EndContainment;

   // List of branches
//   TBranch        *b_PassedSwTrigger;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_NBeamFlashes;   //!
   TBranch        *b_BeamFlashes_YCenter;   //!
   TBranch        *b_BeamFlashes_ZCenter;   //!
   TBranch        *b_BeamFlashes_TotalPE;   //!
   TBranch        *b_NumberMCParticles;   //!
   TBranch        *b_CandidateMuP_Distance;   //!
   TBranch        *b_Vertex_X;   //!
   TBranch        *b_Vertex_Y;   //!
   TBranch        *b_Vertex_Z;   //!
   TBranch        *b_CandidateMu_P;   //!
   TBranch        *b_CandidateMu_Phi;   //!
   TBranch        *b_CandidateMu_CosTheta;   //!
//   TBranch        *b_CandidateMu_Length;   //!
   TBranch        *b_CandidateMu_Chi2_YPlane;   //!
   TBranch        *b_CandidateMu_ThreePlaneChi2;   //!
   TBranch        *b_CandidateMu_StartContainment;   //!
   TBranch        *b_CandidateMu_EndContainment;   //!
   TBranch        *b_CandidateMu_MCParticle_Pdg;   //!
   TBranch        *b_CandidateMu_MCParticle_Purity;   //!
   TBranch        *b_True_CandidateMu_P;   //!
   TBranch        *b_True_CandidateMu_Phi;   //!
   TBranch        *b_True_CandidateMu_CosTheta;   //!
//   TBranch        *b_True_CandidateMu_Length;   //!
   TBranch        *b_True_CandidateMu_StartContainment;   //!
   TBranch        *b_True_CandidateMu_EndContainment;   //!
   TBranch        *b_CandidateP_P;   //!
   TBranch        *b_CandidateP_Phi;   //!
   TBranch        *b_CandidateP_CosTheta;   //!
//   TBranch        *b_CandidateP_Length;   //!
   TBranch        *b_CandidateP_Chi2_YPlane;   //!
   TBranch        *b_CandidateP_ThreePlaneChi2;   //!
   TBranch        *b_CandidateP_StartContainment;   //!
   TBranch        *b_CandidateP_EndContainment;   //!
   TBranch        *b_CandidateP_MCParticle_Pdg;   //!
   TBranch        *b_CandidateP_MCParticle_Purity;   //!
   TBranch        *b_True_CandidateP_P;   //!
   TBranch        *b_True_CandidateP_Phi;   //!
   TBranch        *b_True_CandidateP_CosTheta;   //!
//   TBranch        *b_True_CandidateP_Length;   //!
   TBranch        *b_True_CandidateP_StartContainment;   //!
   TBranch        *b_True_CandidateP_EndContainment;   //!

   myAnalysis(TTree *tree=0);
   virtual ~myAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myAnalysis_cxx
myAnalysis::myAnalysis(TTree *tree) : fChain(0) 
{

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(PathToFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(PathToFile);
      }
      f->GetObject("myPreSelection",tree);

   }
   Init(tree);
}

myAnalysis::~myAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myAnalysis::LoadTree(Long64_t entry)
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

void myAnalysis::Init(TTree *tree)
{

   // Set object pointer
   BeamFlashes_YCenter = 0;
   BeamFlashes_ZCenter = 0;
   BeamFlashes_TotalPE = 0;
   CandidateMuP_Distance = 0;
   Vertex_X = 0;
   Vertex_Y = 0;
   Vertex_Z = 0;
   CandidateMu_P = 0;
   CandidateMu_Phi = 0;
   CandidateMu_CosTheta = 0;
//   CandidateMu_Length = 0;
   CandidateMu_Chi2_YPlane = 0;
   CandidateMu_ThreePlaneChi2 = 0;
   CandidateMu_StartContainment = 0;
   CandidateMu_EndContainment = 0;
   CandidateMu_MCParticle_Pdg = 0;
   CandidateMu_MCParticle_Purity = 0;
   True_CandidateMu_P = 0;
   True_CandidateMu_Phi = 0;
   True_CandidateMu_CosTheta = 0;
//   True_CandidateMu_Length = 0;
   True_CandidateMu_StartContainment = 0;
   True_CandidateMu_EndContainment = 0;
   CandidateP_P = 0;
   CandidateP_Phi = 0;
   CandidateP_CosTheta = 0;
//   CandidateP_Length = 0;
   CandidateP_Chi2_YPlane = 0;
   CandidateP_ThreePlaneChi2 = 0;
   CandidateP_StartContainment = 0;
   CandidateP_EndContainment = 0;
   CandidateP_MCParticle_Pdg = 0;
   CandidateP_MCParticle_Purity = 0;
   True_CandidateP_P = 0;
   True_CandidateP_Phi = 0;
   True_CandidateP_CosTheta = 0;
//   True_CandidateP_Length = 0;
   True_CandidateP_StartContainment = 0;
   True_CandidateP_EndContainment = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

//   fChain->SetBranchAddress("PassedSwTrigger", &PassedSwTrigger, &b_PassedSwTrigger);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("NBeamFlashes", &NBeamFlashes, &b_NBeamFlashes);
   fChain->SetBranchAddress("BeamFlashes_YCenter", &BeamFlashes_YCenter, &b_BeamFlashes_YCenter);
   fChain->SetBranchAddress("BeamFlashes_ZCenter", &BeamFlashes_ZCenter, &b_BeamFlashes_ZCenter);
   fChain->SetBranchAddress("BeamFlashes_TotalPE", &BeamFlashes_TotalPE, &b_BeamFlashes_TotalPE);
   fChain->SetBranchAddress("NumberMCParticles", &NumberMCParticles, &b_NumberMCParticles);
   fChain->SetBranchAddress("CandidateMuP_Distance", &CandidateMuP_Distance, &b_CandidateMuP_Distance);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("CandidateMu_P", &CandidateMu_P, &b_CandidateMu_P);
   fChain->SetBranchAddress("CandidateMu_Phi", &CandidateMu_Phi, &b_CandidateMu_Phi);
   fChain->SetBranchAddress("CandidateMu_CosTheta", &CandidateMu_CosTheta, &b_CandidateMu_CosTheta);
//   fChain->SetBranchAddress("CandidateMu_Length", &CandidateMu_Length, &b_CandidateMu_Length);
   fChain->SetBranchAddress("CandidateMu_Chi2_YPlane", &CandidateMu_Chi2_YPlane, &b_CandidateMu_Chi2_YPlane);
   fChain->SetBranchAddress("CandidateMu_ThreePlaneChi2", &CandidateMu_ThreePlaneChi2, &b_CandidateMu_ThreePlaneChi2);
   fChain->SetBranchAddress("CandidateMu_StartContainment", &CandidateMu_StartContainment, &b_CandidateMu_StartContainment);
   fChain->SetBranchAddress("CandidateMu_EndContainment", &CandidateMu_EndContainment, &b_CandidateMu_EndContainment);
   fChain->SetBranchAddress("CandidateMu_MCParticle_Pdg", &CandidateMu_MCParticle_Pdg, &b_CandidateMu_MCParticle_Pdg);
   fChain->SetBranchAddress("CandidateMu_MCParticle_Purity", &CandidateMu_MCParticle_Purity, &b_CandidateMu_MCParticle_Purity);
   fChain->SetBranchAddress("True_CandidateMu_P", &True_CandidateMu_P, &b_True_CandidateMu_P);
   fChain->SetBranchAddress("True_CandidateMu_Phi", &True_CandidateMu_Phi, &b_True_CandidateMu_Phi);
   fChain->SetBranchAddress("True_CandidateMu_CosTheta", &True_CandidateMu_CosTheta, &b_True_CandidateMu_CosTheta);
//   fChain->SetBranchAddress("True_CandidateMu_Length", &True_CandidateMu_Length, &b_True_CandidateMu_Length);
   fChain->SetBranchAddress("True_CandidateMu_StartContainment", &True_CandidateMu_StartContainment, &b_True_CandidateMu_StartContainment);
   fChain->SetBranchAddress("True_CandidateMu_EndContainment", &True_CandidateMu_EndContainment, &b_True_CandidateMu_EndContainment);
   fChain->SetBranchAddress("CandidateP_P", &CandidateP_P, &b_CandidateP_P);
   fChain->SetBranchAddress("CandidateP_Phi", &CandidateP_Phi, &b_CandidateP_Phi);
   fChain->SetBranchAddress("CandidateP_CosTheta", &CandidateP_CosTheta, &b_CandidateP_CosTheta);
//   fChain->SetBranchAddress("CandidateP_Length", &CandidateP_Length, &b_CandidateP_Length);
   fChain->SetBranchAddress("CandidateP_Chi2_YPlane", &CandidateP_Chi2_YPlane, &b_CandidateP_Chi2_YPlane);
   fChain->SetBranchAddress("CandidateP_ThreePlaneChi2", &CandidateP_ThreePlaneChi2, &b_CandidateP_ThreePlaneChi2);
   fChain->SetBranchAddress("CandidateP_StartContainment", &CandidateP_StartContainment, &b_CandidateP_StartContainment);
   fChain->SetBranchAddress("CandidateP_EndContainment", &CandidateP_EndContainment, &b_CandidateP_EndContainment);
   fChain->SetBranchAddress("CandidateP_MCParticle_Pdg", &CandidateP_MCParticle_Pdg, &b_CandidateP_MCParticle_Pdg);
   fChain->SetBranchAddress("CandidateP_MCParticle_Purity", &CandidateP_MCParticle_Purity, &b_CandidateP_MCParticle_Purity);
   fChain->SetBranchAddress("True_CandidateP_P", &True_CandidateP_P, &b_True_CandidateP_P);
   fChain->SetBranchAddress("True_CandidateP_Phi", &True_CandidateP_Phi, &b_True_CandidateP_Phi);
   fChain->SetBranchAddress("True_CandidateP_CosTheta", &True_CandidateP_CosTheta, &b_True_CandidateP_CosTheta);
//   fChain->SetBranchAddress("True_CandidateP_Length", &True_CandidateP_Length, &b_True_CandidateP_Length);
   fChain->SetBranchAddress("True_CandidateP_StartContainment", &True_CandidateP_StartContainment, &b_True_CandidateP_StartContainment);
   fChain->SetBranchAddress("True_CandidateP_EndContainment", &True_CandidateP_EndContainment, &b_True_CandidateP_EndContainment);
   Notify();
}

Bool_t myAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myAnalysis_cxx
