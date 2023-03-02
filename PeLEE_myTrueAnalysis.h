#ifndef PeLEE_myTrueAnalysis_h
#define PeLEE_myTrueAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <vector>
#include <vector>

#include "../../myClasses/Constants.h"

using namespace Constants;

class PeLEE_myTrueAnalysis {

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
   vector<double>  *All_UBGenie;
   vector<double>  *AxFFCCQEshape_UBGenie;
   vector<double>  *DecayAngMEC_UBGenie;
   vector<double>  *NormCCCOH_UBGenie;
   vector<double>  *NormNCCOH_UBGenie;
   vector<double>  *RPA_CCQE_UBGenie;
   vector<double>  *ThetaDelta2NRad_UBGenie;
   vector<double>  *Theta_Delta2Npi_UBGenie;
   vector<double>  *VecFFCCQEshape_UBGenie;
   vector<double>  *XSecShape_CCMEC_UBGenie;

	//----------------------------------------//

	// detailed xsec uncertainty contributions 

	std::vector<double> *AGKYpT1pi_UBGenie;
	std::vector<double> *AGKYxF1pi_UBGenie;
	std::vector<double> *AhtBY_UBGenie;
	std::vector<double> *BhtBY_UBGenie;
	std::vector<double> *CV1uBY_UBGenie;
	std::vector<double> *CV2uBY_UBGenie;
	std::vector<double> *EtaNCEL_UBGenie;
	std::vector<double> *FrAbs_N_UBGenie;
	std::vector<double> *FrAbs_pi_UBGenie;
	std::vector<double> *FrCEx_N_UBGenie;
	std::vector<double> *FrCEx_pi_UBGenie;
	std::vector<double> *FrInel_N_UBGenie;
	std::vector<double> *FrInel_pi_UBGenie;
	std::vector<double> *FrPiProd_N_UBGenie;
	std::vector<double> *FrPiProd_pi_UBGenie;
	std::vector<double> *FracDelta_CCMEC_UBGenie;	
	std::vector<double> *FracPN_CCMEC_UBGenie;
	std::vector<double> *MFP_N_UBGenie;
	std::vector<double> *MFP_pi_UBGenie;
	std::vector<double> *MaCCQE_UBGenie;
	std::vector<double> *MaCCRES_UBGenie;
	std::vector<double> *MaNCEL_UBGenie;
	std::vector<double> *MaNCRES_UBGenie;
	std::vector<double> *MvCCRES_UBGenie;
	std::vector<double> *MvNCRES_UBGenie;
	std::vector<double> *NonRESBGvbarnCC1pi_UBGenie;
	std::vector<double> *NonRESBGvbarnCC2pi_UBGenie;
	std::vector<double> *NonRESBGvbarnNC1pi_UBGenie;
	std::vector<double> *NonRESBGvbarnNC2pi_UBGenie;
	std::vector<double> *NonRESBGvbarpCC1pi_UBGenie;
	std::vector<double> *NonRESBGvbarpCC2pi_UBGenie;
	std::vector<double> *NonRESBGvbarpNC1pi_UBGenie;
	std::vector<double> *NonRESBGvbarpNC2pi_UBGenie;
	std::vector<double> *NonRESBGvnCC1pi_UBGenie;
	std::vector<double> *NonRESBGvnCC2pi_UBGenie;
	std::vector<double> *NonRESBGvnNC1pi_UBGenie;
	std::vector<double> *NonRESBGvnNC2pi_UBGenie;
	std::vector<double> *NonRESBGvpCC1pi_UBGenie;
	std::vector<double> *NonRESBGvpCC2pi_UBGenie;
	std::vector<double> *NonRESBGvpNC1pi_UBGenie;
	std::vector<double> *NonRESBGvpNC2pi_UBGenie;
	std::vector<double> *NormCCMEC_UBGenie;
	std::vector<double> *NormNCMEC_UBGenie;
	std::vector<double> *RDecBR1eta_UBGenie;
	std::vector<double> *RDecBR1gamma_UBGenie;	   
   
   vector<double>  *fluxes;
   vector<double>  *reinteractions;
   Double_t        True_Ev;
   Double_t        True_Vx;
   Double_t        True_Vy;
   Double_t        True_Vz;
   Int_t           NumberMCParticles;
   Int_t           CC1p;
   Int_t           CC1p1pi;
   Int_t           CC2p;
   Int_t           CC2p1pi;
   Int_t           CC3p;
   Int_t           CC3p1pi;
   Int_t           CC3p2pi;
   Int_t           CC3p3pi;
   Int_t           CC4p0pi;
   Int_t           CC4p1pi;
   Int_t           NumberPi0;
   Int_t           NumberNeutrons;
   Int_t           NumberProtons;
   Int_t           NumberMuons;
   Int_t           NumberChargedPions;
   vector<int>     *Muon_MCParticle_Mode;
   vector<double>  *Muon_MCParticle_Mom;
   vector<double>  *Muon_MCParticle_Phi;
   vector<double>  *Muon_MCParticle_CosTheta;
   vector<double>  *Muon_MCParticle_StartX;
   vector<double>  *Muon_MCParticle_StartY;
   vector<double>  *Muon_MCParticle_StartZ;
   vector<int>     *Muon_MCParticle_StartContainment;
   vector<double>  *Muon_MCParticle_EndX;
   vector<double>  *Muon_MCParticle_EndY;
   vector<double>  *Muon_MCParticle_EndZ;
   vector<int>     *Muon_MCParticle_EndContainment;
   vector<int>     *Muon_MCParticle_Pdg;
   vector<int>     *Proton_MCParticle_Mode;
   vector<double>  *Proton_MCParticle_Mom;
   vector<double>  *Proton_MCParticle_Phi;
   vector<double>  *Proton_MCParticle_CosTheta;
   vector<double>  *Proton_MCParticle_Length;
   vector<double>  *Proton_MCParticle_StartX;
   vector<double>  *Proton_MCParticle_StartY;
   vector<double>  *Proton_MCParticle_StartZ;
   vector<int>     *Proton_MCParticle_StartContainment;
   vector<double>  *Proton_MCParticle_EndX;
   vector<double>  *Proton_MCParticle_EndY;
   vector<double>  *Proton_MCParticle_EndZ;
   vector<int>     *Proton_MCParticle_EndContainment;
   vector<int>     *Proton_MCParticle_Pdg;
   vector<double>  *True_A;
   vector<double>  *True_kMiss;
   vector<double>  *True_EMiss;
   vector<double>  *True_PMissMinus;
   vector<double>  *True_PMiss;
   vector<double>  *True_Pt;
   vector<double>  *True_Ptx;
   vector<double>  *True_Pty;
   vector<double>  *True_PL;
   vector<double>  *True_Pn;
   vector<double>  *True_PnPerp;
   vector<double>  *True_PnPerpx;
   vector<double>  *True_PnPerpy;      
   vector<double>  *True_PnPar;      
   vector<double>  *True_DeltaAlphaT;
   vector<double>  *True_DeltaAlpha3Dq;
   vector<double>  *True_DeltaAlpha3DMu;      
   vector<double>  *True_DeltaPhiT;
   vector<double>  *True_DeltaPhi3D;   
   vector<double>  *True_ECal;
   vector<double>  *True_EQE;
   vector<double>  *True_Q2;
   vector<double>  *True_DeltaPhi;
   vector<double>  *True_DeltaTheta;

   // List of branches
   TBranch        *b_Weight;   //!
   TBranch        *b_T2KWeight;   //!
   TBranch        *b_ROOTinoWeight;   //!
   TBranch        *b_POTWeight;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Event;   //!
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
   
	//----------------------------------------//

	// detailed xsec uncertainty contributions 

   TBranch        *b_AGKYpT1pi_UBGenie;
   TBranch        *b_AGKYxF1pi_UBGenie;
   TBranch        *b_AhtBY_UBGenie;
   TBranch        *b_BhtBY_UBGenie;
   TBranch        *b_CV1uBY_UBGenie;
   TBranch        *b_CV2uBY_UBGenie;
   TBranch        *b_EtaNCEL_UBGenie;
   TBranch        *b_FrAbs_N_UBGenie;
   TBranch        *b_FrAbs_pi_UBGenie;
   TBranch        *b_FrCEx_N_UBGenie;
   TBranch        *b_FrCEx_pi_UBGenie;
   TBranch        *b_FrInel_N_UBGenie;
   TBranch        *b_FrInel_pi_UBGenie;
   TBranch        *b_FrPiProd_N_UBGenie;
   TBranch        *b_FrPiProd_pi_UBGenie;
   TBranch        *b_FracDelta_CCMEC_UBGenie;	
   TBranch        *b_FracPN_CCMEC_UBGenie;
   TBranch        *b_MFP_N_UBGenie;
   TBranch        *b_MFP_pi_UBGenie;
   TBranch        *b_MaCCQE_UBGenie;
   TBranch        *b_MaCCRES_UBGenie;
   TBranch        *b_MaNCEL_UBGenie;
   TBranch        *b_MaNCRES_UBGenie;
   TBranch        *b_MvCCRES_UBGenie;
   TBranch        *b_MvNCRES_UBGenie;
   TBranch        *b_NonRESBGvbarnCC1pi_UBGenie;
   TBranch        *b_NonRESBGvbarnCC2pi_UBGenie;
   TBranch        *b_NonRESBGvbarnNC1pi_UBGenie;
   TBranch        *b_NonRESBGvbarnNC2pi_UBGenie;
   TBranch        *b_NonRESBGvbarpCC1pi_UBGenie;
   TBranch        *b_NonRESBGvbarpCC2pi_UBGenie;
   TBranch        *b_NonRESBGvbarpNC1pi_UBGenie;
   TBranch        *b_NonRESBGvbarpNC2pi_UBGenie;
   TBranch        *b_NonRESBGvnCC1pi_UBGenie;
   TBranch        *b_NonRESBGvnCC2pi_UBGenie;
   TBranch        *b_NonRESBGvnNC1pi_UBGenie;
   TBranch        *b_NonRESBGvnNC2pi_UBGenie;
   TBranch        *b_NonRESBGvpCC1pi_UBGenie;
   TBranch        *b_NonRESBGvpCC2pi_UBGenie;
   TBranch        *b_NonRESBGvpNC1pi_UBGenie;
   TBranch        *b_NonRESBGvpNC2pi_UBGenie;
   TBranch        *b_NormCCMEC_UBGenie;
   TBranch        *b_NormNCMEC_UBGenie;
   TBranch        *b_RDecBR1eta_UBGenie;
   TBranch        *b_RDecBR1gamma_UBGenie;		

	//----------------------------------------//      
   
   TBranch        *b_fluxes;   //!
   TBranch        *b_reinteractions;   //!
   TBranch        *b_True_Ev;   //!
   TBranch        *b_True_Vx;   //!
   TBranch        *b_True_Vy;   //!
   TBranch        *b_True_Vz;   //!
   TBranch        *b_NumberMCParticles;   //!
   TBranch        *b_CC1p;   //!
   TBranch        *b_CC1p1pi;   //!
   TBranch        *b_CC2p;   //!
   TBranch        *b_CC2p1pi;   //!
   TBranch        *b_CC3p;   //!
   TBranch        *b_CC3p1pi;   //!
   TBranch        *b_CC3p2pi;   //!
   TBranch        *b_CC3p3pi;   //!
   TBranch        *b_CC4p0pi;   //!
   TBranch        *b_CC4p1pi;   //!
   TBranch        *b_NumberPi0;   //!
   TBranch        *b_NumberNeutrons;   //!
   TBranch        *b_NumberProtons;   //!
   TBranch        *b_NumberMuons;   //!
   TBranch        *b_NumberChargedPions;   //!
   TBranch        *b_Muon_MCParticle_Mode;   //!
   TBranch        *b_Muon_MCParticle_Mom;   //!
   TBranch        *b_Muon_MCParticle_Phi;   //!
   TBranch        *b_Muon_MCParticle_CosTheta;   //!
   TBranch        *b_Muon_MCParticle_StartX;   //!
   TBranch        *b_Muon_MCParticle_StartY;   //!
   TBranch        *b_Muon_MCParticle_StartZ;   //!
   TBranch        *b_Muon_MCParticle_StartContainment;   //!
   TBranch        *b_Muon_MCParticle_EndX;   //!
   TBranch        *b_Muon_MCParticle_EndY;   //!
   TBranch        *b_Muon_MCParticle_EndZ;   //!
   TBranch        *b_Muon_MCParticle_EndContainment;   //!
   TBranch        *b_Muon_MCParticle_Pdg;   //!
   TBranch        *b_Proton_MCParticle_Mode;   //!
   TBranch        *b_Proton_MCParticle_Mom;   //!
   TBranch        *b_Proton_MCParticle_Phi;   //!
   TBranch        *b_Proton_MCParticle_CosTheta;   //!
   TBranch        *b_Proton_MCParticle_Length;   //!
   TBranch        *b_Proton_MCParticle_StartX;   //!
   TBranch        *b_Proton_MCParticle_StartY;   //!
   TBranch        *b_Proton_MCParticle_StartZ;   //!
   TBranch        *b_Proton_MCParticle_StartContainment;   //!
   TBranch        *b_Proton_MCParticle_EndX;   //!
   TBranch        *b_Proton_MCParticle_EndY;   //!
   TBranch        *b_Proton_MCParticle_EndZ;   //!
   TBranch        *b_Proton_MCParticle_EndContainment;   //!
   TBranch        *b_Proton_MCParticle_Pdg;   //!
   TBranch        *b_True_A;   //!
   TBranch        *b_True_kMiss;   //!
   TBranch        *b_True_EMiss;   //!
   TBranch        *b_True_PMissMinus;   //!
   TBranch        *b_True_PMiss;   //!
   TBranch        *b_True_Pt;   //!
   TBranch        *b_True_Ptx;   //!
   TBranch        *b_True_Pty;   //!
   TBranch        *b_True_PL;   //!
   TBranch        *b_True_Pn;   //!
   TBranch        *b_True_PnPerp;   //!
   TBranch        *b_True_PnPerpx;   //!
   TBranch        *b_True_PnPerpy;   //!      
   TBranch        *b_True_PnPar;   //!      
   TBranch        *b_True_DeltaAlphaT;   //!
   TBranch        *b_True_DeltaAlpha3Dq;   //!
   TBranch        *b_True_DeltaAlpha3DMu;   //!      
   TBranch        *b_True_DeltaPhiT;   //!
   TBranch        *b_True_DeltaPhi3D;   //!   
   TBranch        *b_True_ECal;   //!
   TBranch        *b_True_EQE;   //!
   TBranch        *b_True_Q2;   //!
   TBranch        *b_True_DeltaPhi;   //!
   TBranch        *b_True_DeltaTheta;   //!

   PeLEE_myTrueAnalysis(TString WhichSample="",TString Tune="",TString WhichEventWeightLabel="", int UniverseIndex=-1,TTree *tree=0);
   virtual ~PeLEE_myTrueAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif


#ifdef PeLEE_myTrueAnalysis_cxx
PeLEE_myTrueAnalysis::PeLEE_myTrueAnalysis(TString WhichSample, TString Tune, TString WhichEventWeightLabel, int UniverseIndex, TTree *tree) : fChain(0) 
{

   fTune = Tune;
   fWhichSample = WhichSample;
   fEventWeightLabel = WhichEventWeightLabel;
   fUniverseIndex = UniverseIndex;   

//   fPathToFile = "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PeLEETuples/PreTruthSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
	fPathToFile = "/uboone/data/users/apapadop/PeLEETuples/PreTruthSelection_"+fWhichSample+"_"+UBCodeVersion+".root";   
  
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

PeLEE_myTrueAnalysis::~PeLEE_myTrueAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PeLEE_myTrueAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PeLEE_myTrueAnalysis::LoadTree(Long64_t entry)
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

void PeLEE_myTrueAnalysis::Init(TTree *tree)
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

	AGKYpT1pi_UBGenie = 0;
	AGKYxF1pi_UBGenie = 0;
	AhtBY_UBGenie = 0;
	BhtBY_UBGenie = 0;
	CV1uBY_UBGenie = 0;
	CV2uBY_UBGenie = 0;
	EtaNCEL_UBGenie = 0;
	FrAbs_N_UBGenie = 0;
	FrAbs_pi_UBGenie = 0;
	FrCEx_N_UBGenie = 0;
	FrCEx_pi_UBGenie = 0;
	FrInel_N_UBGenie = 0;
	FrInel_pi_UBGenie = 0;
	FrPiProd_N_UBGenie = 0;
	FrPiProd_pi_UBGenie = 0;
	FracDelta_CCMEC_UBGenie = 0;	
	FracPN_CCMEC_UBGenie = 0;
	MFP_N_UBGenie = 0;
	MFP_pi_UBGenie = 0;
	MaCCQE_UBGenie = 0;
	MaCCRES_UBGenie = 0;
	MaNCEL_UBGenie = 0;
	MaNCRES_UBGenie = 0;
	MvCCRES_UBGenie = 0;
	MvNCRES_UBGenie = 0;
	NonRESBGvbarnCC1pi_UBGenie = 0;
	NonRESBGvbarnCC2pi_UBGenie = 0;
	NonRESBGvbarnNC1pi_UBGenie = 0;
	NonRESBGvbarnNC2pi_UBGenie = 0;
	NonRESBGvbarpCC1pi_UBGenie = 0;
	NonRESBGvbarpCC2pi_UBGenie = 0;
	NonRESBGvbarpNC1pi_UBGenie = 0;
	NonRESBGvbarpNC2pi_UBGenie = 0;
	NonRESBGvnCC1pi_UBGenie = 0;
	NonRESBGvnCC2pi_UBGenie = 0;
	NonRESBGvnNC1pi_UBGenie = 0;
	NonRESBGvnNC2pi_UBGenie = 0;
	NonRESBGvpCC1pi_UBGenie = 0;
	NonRESBGvpCC2pi_UBGenie = 0;
	NonRESBGvpNC1pi_UBGenie = 0;
	NonRESBGvpNC2pi_UBGenie = 0;
	NormCCMEC_UBGenie = 0;
	NormNCMEC_UBGenie = 0;
	RDecBR1eta_UBGenie = 0;
	RDecBR1gamma_UBGenie = 0;		

	//----------------------------------------//    
   
   fluxes = 0;
   reinteractions = 0;
   Muon_MCParticle_Mode = 0;
   Muon_MCParticle_Mom = 0;
   Muon_MCParticle_Phi = 0;
   Muon_MCParticle_CosTheta = 0;
   Muon_MCParticle_StartX = 0;
   Muon_MCParticle_StartY = 0;
   Muon_MCParticle_StartZ = 0;
   Muon_MCParticle_StartContainment = 0;
   Muon_MCParticle_EndX = 0;
   Muon_MCParticle_EndY = 0;
   Muon_MCParticle_EndZ = 0;
   Muon_MCParticle_EndContainment = 0;
   Muon_MCParticle_Pdg = 0;
   Proton_MCParticle_Mode = 0;
   Proton_MCParticle_Mom = 0;
   Proton_MCParticle_Phi = 0;
   Proton_MCParticle_CosTheta = 0;
   Proton_MCParticle_Length = 0;
   Proton_MCParticle_StartX = 0;
   Proton_MCParticle_StartY = 0;
   Proton_MCParticle_StartZ = 0;
   Proton_MCParticle_StartContainment = 0;
   Proton_MCParticle_EndX = 0;
   Proton_MCParticle_EndY = 0;
   Proton_MCParticle_EndZ = 0;
   Proton_MCParticle_EndContainment = 0;
   Proton_MCParticle_Pdg = 0;
   True_A = 0;
   True_kMiss = 0;
   True_EMiss = 0;
   True_PMissMinus = 0;
   True_PMiss = 0;
   True_Pt = 0;
   True_Ptx = 0;
   True_Pty = 0;
   True_PL = 0;
   True_Pn = 0;
   True_PnPerp = 0;
   True_PnPerpx = 0;
   True_PnPerpy = 0;      
   True_PnPar = 0;      
   True_DeltaAlphaT = 0;
   True_DeltaAlpha3Dq = 0;
   True_DeltaAlpha3DMu = 0;      
   True_DeltaPhiT = 0;
   True_DeltaPhi3D = 0;   
   True_ECal = 0;
   True_EQE = 0;
   True_Q2 = 0;
   True_DeltaPhi = 0;
   True_DeltaTheta = 0;
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
   
	//----------------------------------------//

	// detailed xsec uncertainties

	if (fWhichSample == "Overlay9_Run1_DecompXSecUnc") {

		fChain->SetBranchAddress("AGKYpT1pi_UBGenie", &AGKYpT1pi_UBGenie, &b_AGKYpT1pi_UBGenie);
		fChain->SetBranchAddress("AGKYxF1pi_UBGenie", &AGKYxF1pi_UBGenie, &b_AGKYxF1pi_UBGenie);
		fChain->SetBranchAddress("AhtBY_UBGenie", &AhtBY_UBGenie, &b_AhtBY_UBGenie);
		fChain->SetBranchAddress("BhtBY_UBGenie", &BhtBY_UBGenie, &b_BhtBY_UBGenie);
		fChain->SetBranchAddress("CV1uBY_UBGenie", &CV1uBY_UBGenie, &b_CV1uBY_UBGenie);
		fChain->SetBranchAddress("CV2uBY_UBGenie", &CV2uBY_UBGenie, &b_CV2uBY_UBGenie);
		fChain->SetBranchAddress("EtaNCEL_UBGenie", &EtaNCEL_UBGenie, &b_EtaNCEL_UBGenie);
		fChain->SetBranchAddress("FrAbs_N_UBGenie", &FrAbs_N_UBGenie, &b_FrAbs_N_UBGenie);
		fChain->SetBranchAddress("FrAbs_pi_UBGenie", &FrAbs_pi_UBGenie, &b_FrAbs_pi_UBGenie);
		fChain->SetBranchAddress("FrCEx_N_UBGenie", &FrCEx_N_UBGenie, &b_FrCEx_N_UBGenie);
		fChain->SetBranchAddress("FrCEx_pi_UBGenie", &FrCEx_pi_UBGenie, &b_FrCEx_pi_UBGenie);
		fChain->SetBranchAddress("FrInel_N_UBGenie", &FrInel_N_UBGenie, &b_FrInel_N_UBGenie);												
		fChain->SetBranchAddress("FrInel_pi_UBGenie", &FrInel_pi_UBGenie, &b_FrInel_pi_UBGenie);
		fChain->SetBranchAddress("FrPiProd_N_UBGenie", &FrPiProd_N_UBGenie, &b_FrPiProd_N_UBGenie);
		fChain->SetBranchAddress("FrPiProd_pi_UBGenie", &FrPiProd_pi_UBGenie, &b_FrPiProd_pi_UBGenie);
		fChain->SetBranchAddress("FracDelta_CCMEC_UBGenie", &FracDelta_CCMEC_UBGenie, &b_FracDelta_CCMEC_UBGenie);	
		fChain->SetBranchAddress("FracPN_CCMEC_UBGenie", &FracPN_CCMEC_UBGenie, &b_FracPN_CCMEC_UBGenie);
		fChain->SetBranchAddress("MFP_N_UBGenie", &MFP_N_UBGenie, &b_MFP_N_UBGenie);
		fChain->SetBranchAddress("MFP_pi_UBGenie", &MFP_pi_UBGenie, &b_MFP_pi_UBGenie);
		fChain->SetBranchAddress("MaCCQE_UBGenie", &MaCCQE_UBGenie, &b_MaCCQE_UBGenie);
		fChain->SetBranchAddress("MaCCRES_UBGenie", &MaCCRES_UBGenie, &b_MaCCRES_UBGenie);
		fChain->SetBranchAddress("MaNCEL_UBGenie", &MaNCEL_UBGenie, &b_MaNCEL_UBGenie);
		fChain->SetBranchAddress("MaNCRES_UBGenie", &MaNCRES_UBGenie, &b_MaNCRES_UBGenie);
		fChain->SetBranchAddress("MvCCRES_UBGenie", &MvCCRES_UBGenie, &b_MvCCRES_UBGenie);
		fChain->SetBranchAddress("MvNCRES_UBGenie", &MvNCRES_UBGenie, &b_MvNCRES_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarnCC1pi_UBGenie", &NonRESBGvbarnCC1pi_UBGenie, &b_NonRESBGvbarnCC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarnCC2pi_UBGenie", &NonRESBGvbarnCC2pi_UBGenie, &b_NonRESBGvbarnCC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarnNC1pi_UBGenie", &NonRESBGvbarnNC1pi_UBGenie, &b_NonRESBGvbarnNC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarnNC2pi_UBGenie", &NonRESBGvbarnNC2pi_UBGenie, &b_NonRESBGvbarnNC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarpCC1pi_UBGenie", &NonRESBGvbarpCC1pi_UBGenie, &b_NonRESBGvbarpCC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarpCC2pi_UBGenie", &NonRESBGvbarpCC2pi_UBGenie, &b_NonRESBGvbarpCC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarpNC1pi_UBGenie", &NonRESBGvbarpNC1pi_UBGenie, &b_NonRESBGvbarpNC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvbarpNC2pi_UBGenie", &NonRESBGvbarpNC2pi_UBGenie, &b_NonRESBGvbarpNC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvnCC1pi_UBGenie", &NonRESBGvnCC1pi_UBGenie, &b_NonRESBGvnCC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvnCC2pi_UBGenie", &NonRESBGvnCC2pi_UBGenie, &b_NonRESBGvnCC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvnNC1pi_UBGenie", &NonRESBGvnNC1pi_UBGenie, &b_NonRESBGvnNC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvnNC2pi_UBGenie", &NonRESBGvnNC2pi_UBGenie, &b_NonRESBGvnNC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvpCC1pi_UBGenie", &NonRESBGvpCC1pi_UBGenie, &b_NonRESBGvpCC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvpCC2pi_UBGenie", &NonRESBGvpCC2pi_UBGenie, &b_NonRESBGvpCC2pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvpNC1pi_UBGenie", &NonRESBGvpNC1pi_UBGenie, &b_NonRESBGvpNC1pi_UBGenie);
		fChain->SetBranchAddress("NonRESBGvpNC2pi_UBGenie", &NonRESBGvpNC2pi_UBGenie, &b_NonRESBGvpNC2pi_UBGenie);
		fChain->SetBranchAddress("NormCCMEC_UBGenie", &NormCCMEC_UBGenie, &b_NormCCMEC_UBGenie);
		fChain->SetBranchAddress("NormNCMEC_UBGenie", &NormNCMEC_UBGenie, &b_NormNCMEC_UBGenie);
		fChain->SetBranchAddress("RDecBR1eta_UBGenie", &RDecBR1eta_UBGenie, &b_RDecBR1eta_UBGenie);
		fChain->SetBranchAddress("RDecBR1gamma_UBGenie", &RDecBR1gamma_UBGenie, &b_RDecBR1gamma_UBGenie);

	}

	//----------------------------------------//      
   
   fChain->SetBranchAddress("fluxes", &fluxes, &b_fluxes);
   fChain->SetBranchAddress("reinteractions", &reinteractions, &b_reinteractions);
   fChain->SetBranchAddress("True_Ev", &True_Ev, &b_True_Ev);
   fChain->SetBranchAddress("True_Vx", &True_Vx, &b_True_Vx);
   fChain->SetBranchAddress("True_Vy", &True_Vy, &b_True_Vy);
   fChain->SetBranchAddress("True_Vz", &True_Vz, &b_True_Vz);
   fChain->SetBranchAddress("NumberMCParticles", &NumberMCParticles, &b_NumberMCParticles);
   fChain->SetBranchAddress("CC1p", &CC1p, &b_CC1p);
   fChain->SetBranchAddress("CC1p1pi", &CC1p1pi, &b_CC1p1pi);
   fChain->SetBranchAddress("CC2p", &CC2p, &b_CC2p);
   fChain->SetBranchAddress("CC2p1pi", &CC2p1pi, &b_CC2p1pi);
   fChain->SetBranchAddress("CC3p", &CC3p, &b_CC3p);
   fChain->SetBranchAddress("CC3p1pi", &CC3p1pi, &b_CC3p1pi);
   fChain->SetBranchAddress("CC3p2pi", &CC3p2pi, &b_CC3p2pi);
   fChain->SetBranchAddress("CC3p3pi", &CC3p3pi, &b_CC3p3pi);
   fChain->SetBranchAddress("CC4p0pi", &CC4p0pi, &b_CC4p0pi);
   fChain->SetBranchAddress("CC4p1pi", &CC4p1pi, &b_CC4p1pi);
   fChain->SetBranchAddress("NumberPi0", &NumberPi0, &b_NumberPi0);
   fChain->SetBranchAddress("NumberNeutrons", &NumberNeutrons, &b_NumberNeutrons);
   fChain->SetBranchAddress("NumberProtons", &NumberProtons, &b_NumberProtons);
   fChain->SetBranchAddress("NumberMuons", &NumberMuons, &b_NumberMuons);
   fChain->SetBranchAddress("NumberChargedPions", &NumberChargedPions, &b_NumberChargedPions);
   fChain->SetBranchAddress("Muon_MCParticle_Mode", &Muon_MCParticle_Mode, &b_Muon_MCParticle_Mode);
   fChain->SetBranchAddress("Muon_MCParticle_Mom", &Muon_MCParticle_Mom, &b_Muon_MCParticle_Mom);
   fChain->SetBranchAddress("Muon_MCParticle_Phi", &Muon_MCParticle_Phi, &b_Muon_MCParticle_Phi);
   fChain->SetBranchAddress("Muon_MCParticle_CosTheta", &Muon_MCParticle_CosTheta, &b_Muon_MCParticle_CosTheta);
   fChain->SetBranchAddress("Muon_MCParticle_StartX", &Muon_MCParticle_StartX, &b_Muon_MCParticle_StartX);
   fChain->SetBranchAddress("Muon_MCParticle_StartY", &Muon_MCParticle_StartY, &b_Muon_MCParticle_StartY);
   fChain->SetBranchAddress("Muon_MCParticle_StartZ", &Muon_MCParticle_StartZ, &b_Muon_MCParticle_StartZ);
   fChain->SetBranchAddress("Muon_MCParticle_StartContainment", &Muon_MCParticle_StartContainment, &b_Muon_MCParticle_StartContainment);
   fChain->SetBranchAddress("Muon_MCParticle_EndX", &Muon_MCParticle_EndX, &b_Muon_MCParticle_EndX);
   fChain->SetBranchAddress("Muon_MCParticle_EndY", &Muon_MCParticle_EndY, &b_Muon_MCParticle_EndY);
   fChain->SetBranchAddress("Muon_MCParticle_EndZ", &Muon_MCParticle_EndZ, &b_Muon_MCParticle_EndZ);
   fChain->SetBranchAddress("Muon_MCParticle_EndContainment", &Muon_MCParticle_EndContainment, &b_Muon_MCParticle_EndContainment);
   fChain->SetBranchAddress("Muon_MCParticle_Pdg", &Muon_MCParticle_Pdg, &b_Muon_MCParticle_Pdg);
   fChain->SetBranchAddress("Proton_MCParticle_Mode", &Proton_MCParticle_Mode, &b_Proton_MCParticle_Mode);
   fChain->SetBranchAddress("Proton_MCParticle_Mom", &Proton_MCParticle_Mom, &b_Proton_MCParticle_Mom);
   fChain->SetBranchAddress("Proton_MCParticle_Phi", &Proton_MCParticle_Phi, &b_Proton_MCParticle_Phi);
   fChain->SetBranchAddress("Proton_MCParticle_CosTheta", &Proton_MCParticle_CosTheta, &b_Proton_MCParticle_CosTheta);
   fChain->SetBranchAddress("Proton_MCParticle_Length", &Proton_MCParticle_Length, &b_Proton_MCParticle_Length);
   fChain->SetBranchAddress("Proton_MCParticle_StartX", &Proton_MCParticle_StartX, &b_Proton_MCParticle_StartX);
   fChain->SetBranchAddress("Proton_MCParticle_StartY", &Proton_MCParticle_StartY, &b_Proton_MCParticle_StartY);
   fChain->SetBranchAddress("Proton_MCParticle_StartZ", &Proton_MCParticle_StartZ, &b_Proton_MCParticle_StartZ);
   fChain->SetBranchAddress("Proton_MCParticle_StartContainment", &Proton_MCParticle_StartContainment, &b_Proton_MCParticle_StartContainment);
   fChain->SetBranchAddress("Proton_MCParticle_EndX", &Proton_MCParticle_EndX, &b_Proton_MCParticle_EndX);
   fChain->SetBranchAddress("Proton_MCParticle_EndY", &Proton_MCParticle_EndY, &b_Proton_MCParticle_EndY);
   fChain->SetBranchAddress("Proton_MCParticle_EndZ", &Proton_MCParticle_EndZ, &b_Proton_MCParticle_EndZ);
   fChain->SetBranchAddress("Proton_MCParticle_EndContainment", &Proton_MCParticle_EndContainment, &b_Proton_MCParticle_EndContainment);
   fChain->SetBranchAddress("Proton_MCParticle_Pdg", &Proton_MCParticle_Pdg, &b_Proton_MCParticle_Pdg);
   fChain->SetBranchAddress("True_A", &True_A, &b_True_A);
   fChain->SetBranchAddress("True_kMiss", &True_kMiss, &b_True_kMiss);
   fChain->SetBranchAddress("True_EMiss", &True_EMiss, &b_True_EMiss);
   fChain->SetBranchAddress("True_PMissMinus", &True_PMissMinus, &b_True_PMissMinus);
   fChain->SetBranchAddress("True_PMiss", &True_PMiss, &b_True_PMiss);
   fChain->SetBranchAddress("True_Pt", &True_Pt, &b_True_Pt);
   fChain->SetBranchAddress("True_Ptx", &True_Ptx, &b_True_Ptx);
   fChain->SetBranchAddress("True_Pty", &True_Pty, &b_True_Pty);
   fChain->SetBranchAddress("True_PL", &True_PL, &b_True_PL);
   fChain->SetBranchAddress("True_Pn", &True_Pn, &b_True_Pn);
   fChain->SetBranchAddress("True_PnPerp", &True_PnPerp, &b_True_PnPerp);
   fChain->SetBranchAddress("True_PnPerpx", &True_PnPerpx, &b_True_PnPerpx);
   fChain->SetBranchAddress("True_PnPerpy", &True_PnPerpy, &b_True_PnPerpy);      
   fChain->SetBranchAddress("True_PnPar", &True_PnPar, &b_True_PnPar);      
   fChain->SetBranchAddress("True_DeltaAlphaT", &True_DeltaAlphaT, &b_True_DeltaAlphaT);
   fChain->SetBranchAddress("True_DeltaAlpha3Dq", &True_DeltaAlpha3Dq, &b_True_DeltaAlpha3Dq);
   fChain->SetBranchAddress("True_DeltaAlpha3DMu", &True_DeltaAlpha3DMu, &b_True_DeltaAlpha3DMu);      
   fChain->SetBranchAddress("True_DeltaPhiT", &True_DeltaPhiT, &b_True_DeltaPhiT);
   fChain->SetBranchAddress("True_DeltaPhi3D", &True_DeltaPhi3D, &b_True_DeltaPhi3D);   
   fChain->SetBranchAddress("True_ECal", &True_ECal, &b_True_ECal);
   fChain->SetBranchAddress("True_EQE", &True_EQE, &b_True_EQE);
   fChain->SetBranchAddress("True_Q2", &True_Q2, &b_True_Q2);
   fChain->SetBranchAddress("True_DeltaPhi", &True_DeltaPhi, &b_True_DeltaPhi);
   fChain->SetBranchAddress("True_DeltaTheta", &True_DeltaTheta, &b_True_DeltaTheta);
   Notify();

}

Bool_t PeLEE_myTrueAnalysis::Notify()
{
   return kTRUE;
}

void PeLEE_myTrueAnalysis::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PeLEE_myTrueAnalysis::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef PeLEE_myTrueAnalysis_cxx
