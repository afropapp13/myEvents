#ifndef t_h
#define t_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include "../myCCQEAnalysis/Constants.h"

#include <vector>
#include <vector>

using namespace Constants;

class t {

private:
	TString fPathToFile;
	TString fWhichSample;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
//   Int_t           PassedSwTrigger;
   double           Weight;
   int             CC1p;
   int             MCParticle_Mode;
   double          NuScore;
   double          FlashScore;
   Int_t           NBeamFlashes;
   vector<double>  *BeamFlashes_YCenter;
   vector<double>  *BeamFlashes_ZCenter;
   vector<double>  *BeamFlashes_TotalPE;
//   Int_t           NumberMCParticles;
//   vector<double>  *MCParticle_Mom;
//   vector<double>  *MCParticle_Mode;
//   vector<double>  *MCParticle_Phi;
//   vector<double>  *MCParticle_CosTheta;
////   vector<double>  *MCParticle_Length;
//   vector<int>     *MCParticle_StartContainment;
//   vector<int>     *MCParticle_EndContainment;
//   vector<int>     *MCParticle_Pdg;
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
   vector<double>  *CandidateMu_ThreePlaneLogLikelihood;
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
   vector<double>  *CandidateP_ThreePlaneLogLikelihood;
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

//   Int_t           NuMuPFParticles;
//   vector<int>     *PFParticle_NuMuDaughters;
//   vector<vector<int> > *PFParticle_NuMuDaughtersPdgCode;

   // List of branches
//   TBranch        *b_PassedSwTrigger;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_CC1p;   //!
   TBranch        *b_MCParticle_Mode;   //!
   TBranch        *b_NuScore;   //!
   TBranch        *b_FlashScore;   //!
   TBranch        *b_NBeamFlashes;   //!
   TBranch        *b_BeamFlashes_YCenter;   //!
   TBranch        *b_BeamFlashes_ZCenter;   //!
   TBranch        *b_BeamFlashes_TotalPE;   //!
//   TBranch        *b_NumberMCParticles;   //!
//   TBranch        *b_MCParticle_Mom;   //!
//   TBranch        *b_MCParticle_Mode;   //!
//   TBranch        *b_MCParticle_Phi;   //!
//   TBranch        *b_MCParticle_CosTheta;   //!
////   TBranch        *b_MCParticle_Length;   //!
//   TBranch        *b_MCParticle_StartContainment;   //!
//   TBranch        *b_MCParticle_EndContainment;   //!
//   TBranch        *b_MCParticle_Pdg;   //!
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
   TBranch        *b_CandidateMu_ThreePlaneLogLikelihood;   //!
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
   TBranch        *b_CandidateP_ThreePlaneLogLikelihood;   //!
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

//   TBranch        *b_NuMuPFParticles;   //!
//   TBranch        *b_PFParticle_NuMuDaughters;   //!
//   TBranch        *b_PFParticle_NuMuDaughtersPdgCode;   //!

   t(TString WhichSample, TTree *tree=0);
   virtual ~t();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef t_cxx
t::t(TString WhichSample, TTree *tree) : fChain(0) 
{

   fWhichSample = WhichSample;

   fPathToFile = "mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fPathToFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(fPathToFile);
      }
      f->GetObject("myPreSelection",tree);

   }
   Init(tree);
}

t::~t()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t t::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t t::LoadTree(Long64_t entry)
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

void t::Init(TTree *tree)
{

   // Set object pointer
   BeamFlashes_YCenter = 0;
   BeamFlashes_ZCenter = 0;
   BeamFlashes_TotalPE = 0;
//   MCParticle_Mom = 0;
//   MCParticle_Mode = 0;
//   MCParticle_Phi = 0;
//   MCParticle_CosTheta = 0;
////   MCParticle_Length = 0;
//   MCParticle_StartContainment = 0;
//   MCParticle_EndContainment = 0;
//   MCParticle_Pdg = 0;
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
   CandidateMu_ThreePlaneLogLikelihood = 0;
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
   CandidateP_ThreePlaneLogLikelihood = 0;
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

//   PFParticle_NuMuDaughters = 0;
//   PFParticle_NuMuDaughtersPdgCode = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

//   fChain->SetBranchAddress("PassedSwTrigger", &PassedSwTrigger, &b_PassedSwTrigger);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("CC1p", &CC1p, &b_CC1p);
   fChain->SetBranchAddress("MCParticle_Mode", &MCParticle_Mode, &b_MCParticle_Mode);
   fChain->SetBranchAddress("NuScore", &NuScore, &b_NuScore);
   fChain->SetBranchAddress("FlashScore", &FlashScore, &b_FlashScore);
   fChain->SetBranchAddress("NBeamFlashes", &NBeamFlashes, &b_NBeamFlashes);
   fChain->SetBranchAddress("BeamFlashes_YCenter", &BeamFlashes_YCenter, &b_BeamFlashes_YCenter);
   fChain->SetBranchAddress("BeamFlashes_ZCenter", &BeamFlashes_ZCenter, &b_BeamFlashes_ZCenter);
   fChain->SetBranchAddress("BeamFlashes_TotalPE", &BeamFlashes_TotalPE, &b_BeamFlashes_TotalPE);
//   fChain->SetBranchAddress("NumberMCParticles", &NumberMCParticles, &b_NumberMCParticles);
//   fChain->SetBranchAddress("MCParticle_Mom", &MCParticle_Mom, &b_MCParticle_Mom);
//   fChain->SetBranchAddress("MCParticle_Mode", &MCParticle_Mode, &b_MCParticle_Mode);
//   fChain->SetBranchAddress("MCParticle_Phi", &MCParticle_Phi, &b_MCParticle_Phi);
//   fChain->SetBranchAddress("MCParticle_CosTheta", &MCParticle_CosTheta, &b_MCParticle_CosTheta);
////   fChain->SetBranchAddress("MCParticle_Length", &MCParticle_Length, &b_MCParticle_Length);
//   fChain->SetBranchAddress("MCParticle_StartContainment", &MCParticle_StartContainment, &b_MCParticle_StartContainment);
//   fChain->SetBranchAddress("MCParticle_EndContainment", &MCParticle_EndContainment, &b_MCParticle_EndContainment);
//   fChain->SetBranchAddress("MCParticle_Pdg", &MCParticle_Pdg, &b_MCParticle_Pdg);
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
   fChain->SetBranchAddress("CandidateMu_ThreePlaneLogLikelihood", &CandidateMu_ThreePlaneLogLikelihood, &b_CandidateMu_ThreePlaneLogLikelihood);
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
   fChain->SetBranchAddress("CandidateP_ThreePlaneLogLikelihood", &CandidateP_ThreePlaneLogLikelihood, &b_CandidateP_ThreePlaneLogLikelihood);
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

//   fChain->SetBranchAddress("NuMuPFParticles", &NuMuPFParticles, &b_NuMuPFParticles);
//   fChain->SetBranchAddress("PFParticle_NuMuDaughters", &PFParticle_NuMuDaughters, &b_PFParticle_NuMuDaughters);
//   fChain->SetBranchAddress("PFParticle_NuMuDaughtersPdgCode", &PFParticle_NuMuDaughtersPdgCode, &b_PFParticle_NuMuDaughtersPdgCode);

   Notify();
}

Bool_t t::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void t::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t t::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef t_cxx
