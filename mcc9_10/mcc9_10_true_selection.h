#ifndef mcc9_10_true_selection_h
#define mcc9_10_true_selection_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <vector>
#include <vector>

#include "../../../generators/constants.h"

using namespace constants;

class mcc9_10_true_selection {

private:

   TString fTune;
   TString fPathToFile;
   TString fWhichSample;
   TString fEventWeightLabel;
   int     fUniverseIndex;  
   TFile* fFile; 

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        Weight;
   Double_t        T2KWeight;
   Double_t        ROOTinoWeight;
   Double_t        POTWeight;
   Int_t           Run;
   Int_t           SubRun;
   Int_t           Event;
   TString         *run_period;
   vector<unsigned short> *All_UBGenie;
   vector<double>  *AxFFCCQEshape_UBGenie;
   vector<double>  *DecayAngMEC_UBGenie;
   vector<double>  *NormCCCOH_UBGenie;
   vector<double>  *NormNCCOH_UBGenie;
   vector<double>  *RPA_CCQE_UBGenie;
   vector<double>  *ThetaDelta2NRad_UBGenie;
   vector<double>  *Theta_Delta2Npi_UBGenie;
   vector<double>  *VecFFCCQEshape_UBGenie;
   vector<double>  *XSecShape_CCMEC_UBGenie;
   vector<unsigned short> *fluxes;
   vector<unsigned short> *reinteractions;
   Double_t        True_Ev;
   Double_t        True_Vx;
   Double_t        True_Vy;
   Double_t        True_Vz;
   Int_t           signal;
   Int_t           NCcoh;
   Int_t           NCres;
   vector<int>     *pi0_MCParticle_Mode;
   vector<double>  *pi0_MCParticle_Mom;
   vector<double>  *pi0_MCParticle_Phi;
   vector<double>  *pi0_MCParticle_CosTheta;
   vector<double>  *pi0_MCParticle_StartX;
   vector<double>  *pi0_MCParticle_StartY;
   vector<double>  *pi0_MCParticle_StartZ;
   vector<int>     *pi0_MCParticle_StartContainment;
   vector<double>  *pi0_MCParticle_EndX;
   vector<double>  *pi0_MCParticle_EndY;
   vector<double>  *pi0_MCParticle_EndZ;
   vector<int>     *pi0_MCParticle_EndContainment;
   vector<int>     *pi0_MCParticle_Pdg;

   // List of branches
   TBranch        *b_Weight;   //!
   TBranch        *b_T2KWeight;   //!
   TBranch        *b_ROOTinoWeight;   //!
   TBranch        *b_POTWeight;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_run_period;   //!
   TBranch        *b_All_UBGenie;   //!
   TBranch        *b_AxFFCCQEshape_UBGenie;   //!
   TBranch        *b_DecayAngMEC_UBGenie;   //!
   TBranch        *b_NormCCCOH_UBGenie;   //!
   TBranch        *b_NormNCCOH_UBGenie;   //!
   TBranch        *b_RPA_CCQE_UBGenie;   //!
   TBranch        *b_ThetaDelta2NRad_UBGenie;   //!
   TBranch        *b_Theta_Delta2Npi_UBGenie;   //!
   TBranch        *b_VecFFCCQEshape_UBGenie;   //!
   TBranch        *b_XSecShape_CCMEC_UBGenie;   //!
   TBranch        *b_fluxes;   //!
   TBranch        *b_reinteractions;   //!
   TBranch        *b_True_Ev;   //!
   TBranch        *b_True_Vx;   //!
   TBranch        *b_True_Vy;   //!
   TBranch        *b_True_Vz;   //!
   TBranch        *b_signal;   //!
   TBranch        *b_NCcoh;   //!
   TBranch        *b_NCres;   //!
   TBranch        *b_pi0_MCParticle_Mode;   //!
   TBranch        *b_pi0_MCParticle_Mom;   //!
   TBranch        *b_pi0_MCParticle_Phi;   //!
   TBranch        *b_pi0_MCParticle_CosTheta;   //!
   TBranch        *b_pi0_MCParticle_StartX;   //!
   TBranch        *b_pi0_MCParticle_StartY;   //!
   TBranch        *b_pi0_MCParticle_StartZ;   //!
   TBranch        *b_pi0_MCParticle_StartContainment;   //!
   TBranch        *b_pi0_MCParticle_EndX;   //!
   TBranch        *b_pi0_MCParticle_EndY;   //!
   TBranch        *b_pi0_MCParticle_EndZ;   //!
   TBranch        *b_pi0_MCParticle_EndContainment;   //!
   TBranch        *b_pi0_MCParticle_Pdg;   //!

   mcc9_10_true_selection(TString WhichSample="",TString Tune="",TString WhichEventWeightLabel="", int UniverseIndex=-1,TTree *tree=0);
   virtual ~mcc9_10_true_selection();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif


#ifdef mcc9_10_true_selection_cxx
mcc9_10_true_selection::mcc9_10_true_selection(TString WhichSample, TString Tune, TString WhichEventWeightLabel, int UniverseIndex, TTree *tree) : fChain(0) 
{

   fTune = Tune;
   fWhichSample = WhichSample;
   fEventWeightLabel = WhichEventWeightLabel;
   fUniverseIndex = UniverseIndex;   

//	//pnfsToXRootD /pnfs/persistent/path/to/your/file
	fPathToFile = "/exp/uboone/data/users/"+UserID+"/ncpi0/PreTruthSelection_"+fWhichSample+".root";   
  
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fPathToFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(fPathToFile);
      }
      f->GetObject("myPreTruthSelection",tree);
      fFile = f;

   }
   Init(tree);
}

mcc9_10_true_selection::~mcc9_10_true_selection()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mcc9_10_true_selection::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mcc9_10_true_selection::LoadTree(Long64_t entry)
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

void mcc9_10_true_selection::Init(TTree *tree)
{

   // Set object pointer
   All_UBGenie = 0;
   AxFFCCQEshape_UBGenie = 0;
   DecayAngMEC_UBGenie = 0;
   NormCCCOH_UBGenie = 0;
   NormNCCOH_UBGenie = 0;
   RPA_CCQE_UBGenie = 0;
   ThetaDelta2NRad_UBGenie = 0;
   Theta_Delta2Npi_UBGenie = 0;
   VecFFCCQEshape_UBGenie = 0;
   XSecShape_CCMEC_UBGenie = 0;
   
	//----------------------------------------//

	// detailed xsec uncertainty contributions 

   run_period = 0;
   All_UBGenie = 0;
   AxFFCCQEshape_UBGenie = 0;
   DecayAngMEC_UBGenie = 0;
   NormCCCOH_UBGenie = 0;
   NormNCCOH_UBGenie = 0;
   RPA_CCQE_UBGenie = 0;
   ThetaDelta2NRad_UBGenie = 0;
   Theta_Delta2Npi_UBGenie = 0;
   VecFFCCQEshape_UBGenie = 0;
   XSecShape_CCMEC_UBGenie = 0;
   fluxes = 0;
   reinteractions = 0;
   pi0_MCParticle_Mode = 0;
   pi0_MCParticle_Mom = 0;
   pi0_MCParticle_Phi = 0;
   pi0_MCParticle_CosTheta = 0;
   pi0_MCParticle_StartX = 0;
   pi0_MCParticle_StartY = 0;
   pi0_MCParticle_StartZ = 0;
   pi0_MCParticle_StartContainment = 0;
   pi0_MCParticle_EndX = 0;
   pi0_MCParticle_EndY = 0;
   pi0_MCParticle_EndZ = 0;
   pi0_MCParticle_EndContainment = 0;
   pi0_MCParticle_Pdg = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("T2KWeight", &T2KWeight, &b_T2KWeight);
   fChain->SetBranchAddress("ROOTinoWeight", &ROOTinoWeight, &b_ROOTinoWeight);
   fChain->SetBranchAddress("POTWeight", &POTWeight, &b_POTWeight);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("run_period", &run_period, &b_run_period);
   fChain->SetBranchAddress("All_UBGenie", &All_UBGenie, &b_All_UBGenie);
   fChain->SetBranchAddress("AxFFCCQEshape_UBGenie", &AxFFCCQEshape_UBGenie, &b_AxFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("DecayAngMEC_UBGenie", &DecayAngMEC_UBGenie, &b_DecayAngMEC_UBGenie);
   fChain->SetBranchAddress("NormCCCOH_UBGenie", &NormCCCOH_UBGenie, &b_NormCCCOH_UBGenie);
   fChain->SetBranchAddress("NormNCCOH_UBGenie", &NormNCCOH_UBGenie, &b_NormNCCOH_UBGenie);
   fChain->SetBranchAddress("RPA_CCQE_UBGenie", &RPA_CCQE_UBGenie, &b_RPA_CCQE_UBGenie);
   fChain->SetBranchAddress("ThetaDelta2NRad_UBGenie", &ThetaDelta2NRad_UBGenie, &b_ThetaDelta2NRad_UBGenie);
   fChain->SetBranchAddress("Theta_Delta2Npi_UBGenie", &Theta_Delta2Npi_UBGenie, &b_Theta_Delta2Npi_UBGenie);
   fChain->SetBranchAddress("VecFFCCQEshape_UBGenie", &VecFFCCQEshape_UBGenie, &b_VecFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("XSecShape_CCMEC_UBGenie", &XSecShape_CCMEC_UBGenie, &b_XSecShape_CCMEC_UBGenie);
   fChain->SetBranchAddress("fluxes", &fluxes, &b_fluxes);
   fChain->SetBranchAddress("reinteractions", &reinteractions, &b_reinteractions);
   fChain->SetBranchAddress("True_Ev", &True_Ev, &b_True_Ev);
   fChain->SetBranchAddress("True_Vx", &True_Vx, &b_True_Vx);
   fChain->SetBranchAddress("True_Vy", &True_Vy, &b_True_Vy);
   fChain->SetBranchAddress("True_Vz", &True_Vz, &b_True_Vz);
   fChain->SetBranchAddress("signal", &signal, &b_signal);
   fChain->SetBranchAddress("NCcoh", &NCcoh, &b_NCcoh);
   fChain->SetBranchAddress("NCres", &NCres, &b_NCres);
   fChain->SetBranchAddress("pi0_MCParticle_Mode", &pi0_MCParticle_Mode, &b_pi0_MCParticle_Mode);
   fChain->SetBranchAddress("pi0_MCParticle_Mom", &pi0_MCParticle_Mom, &b_pi0_MCParticle_Mom);
   fChain->SetBranchAddress("pi0_MCParticle_Phi", &pi0_MCParticle_Phi, &b_pi0_MCParticle_Phi);
   fChain->SetBranchAddress("pi0_MCParticle_CosTheta", &pi0_MCParticle_CosTheta, &b_pi0_MCParticle_CosTheta);
   fChain->SetBranchAddress("pi0_MCParticle_StartX", &pi0_MCParticle_StartX, &b_pi0_MCParticle_StartX);
   fChain->SetBranchAddress("pi0_MCParticle_StartY", &pi0_MCParticle_StartY, &b_pi0_MCParticle_StartY);
   fChain->SetBranchAddress("pi0_MCParticle_StartZ", &pi0_MCParticle_StartZ, &b_pi0_MCParticle_StartZ);
   fChain->SetBranchAddress("pi0_MCParticle_StartContainment", &pi0_MCParticle_StartContainment, &b_pi0_MCParticle_StartContainment);
   fChain->SetBranchAddress("pi0_MCParticle_EndX", &pi0_MCParticle_EndX, &b_pi0_MCParticle_EndX);
   fChain->SetBranchAddress("pi0_MCParticle_EndY", &pi0_MCParticle_EndY, &b_pi0_MCParticle_EndY);
   fChain->SetBranchAddress("pi0_MCParticle_EndZ", &pi0_MCParticle_EndZ, &b_pi0_MCParticle_EndZ);
   fChain->SetBranchAddress("pi0_MCParticle_EndContainment", &pi0_MCParticle_EndContainment, &b_pi0_MCParticle_EndContainment);
   fChain->SetBranchAddress("pi0_MCParticle_Pdg", &pi0_MCParticle_Pdg, &b_pi0_MCParticle_Pdg);

   Notify();

}

Bool_t mcc9_10_true_selection::Notify()
{
   return kTRUE;
}

void mcc9_10_true_selection::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mcc9_10_true_selection::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef mcc9_10_true_selection_cxx
