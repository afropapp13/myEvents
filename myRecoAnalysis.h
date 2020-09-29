#ifndef myRecoAnalysis_h
#define myRecoAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_43/srcs/ubana/ubana/myClasses/Constants.h"

#include <vector>
#include <vector>

using namespace Constants;

class myRecoAnalysis {

private:
	TString fPathToFile;
	TString fWhichSample;
	TString fEventWeightLabel;
	int     fUniverseIndex;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
//   Int_t           PassedSwTrigger;
   double           Weight;
   double           T2KWeight;
   double           ROOTinoWeight;

   vector<double>  *All_UBGenie;
   vector<double>  *AxFFCCQEshape_UBGenie;
   vector<double>  *DecayAngMEC_UBGenie;
   vector<double>  *NormCCCOH_UBGenie;
   vector<double>  *NormNCCOH_UBGenie;
//   vector<double>  *RPA_CCQE_Reduced_UBGenie;
   vector<double>  *RPA_CCQE_UBGenie;
   vector<double>  *ThetaDelta2NRad_UBGenie;
   vector<double>  *Theta_Delta2Npi_UBGenie;
   vector<double>  *VecFFCCQEshape_UBGenie;
   vector<double>  *XSecShape_CCMEC_UBGenie;
   vector<double>  *expskin_FluxUnisim;
   vector<double>  *horncurrent_FluxUnisim;
   vector<double>  *kminus_PrimaryHadronNormalization;
   vector<double>  *kplus_PrimaryHadronFeynmanScaling;
   vector<double>  *kzero_PrimaryHadronSanfordWang;
   vector<double>  *nucleoninexsec_FluxUnisim;
   vector<double>  *nucleonqexsec_FluxUnisim;
   vector<double>  *nucleontotxsec_FluxUnisim;
   vector<double>  *piminus_PrimaryHadronSWCentralSplineVariation;
   vector<double>  *pioninexsec_FluxUnisim;
   vector<double>  *pionqexsec_FluxUnisim;
   vector<double>  *piontotxsec_FluxUnisim;
   vector<double>  *piplus_PrimaryHadronSWCentralSplineVariation;
   vector<double>  *reinteractions_piminus_Geant4;
   vector<double>  *reinteractions_piplus_Geant4;
   vector<double>  *reinteractions_proton_Geant4;
//   vector<double>  *xsr_scc_Fa3_SCC;
//   vector<double>  *xsr_scc_Fv3_SCC;

   
   int             CC1p;
   int             CC1p1pi;
   int             CC2p;
   int             CC2p1pi;
   int             CC3p;
   int             CC3p1pi;
   int             CC3p2pi;   
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

   vector<double>  *CandidateMu_StartX;
   vector<double>  *CandidateMu_StartY;
   vector<double>  *CandidateMu_StartZ;
   vector<double>  *CandidateMu_EndX;
   vector<double>  *CandidateMu_EndY;
   vector<double>  *CandidateMu_EndZ;

   vector<double>  *True_CandidateMu_StartX;
   vector<double>  *True_CandidateMu_StartY;
   vector<double>  *True_CandidateMu_StartZ;
   vector<double>  *True_CandidateMu_EndX;
   vector<double>  *True_CandidateMu_EndY;
   vector<double>  *True_CandidateMu_EndZ;            

   vector<double>  *CandidateMu_P_Range;
   vector<double>  *CandidateMu_P_MCS;   
   vector<double>  *CandidateMu_Phi;
   vector<double>  *CandidateMu_CosTheta;
   vector<double>  *CandidateMu_Length;
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

   vector<double>  *CandidateP_StartX;
   vector<double>  *CandidateP_StartY;
   vector<double>  *CandidateP_StartZ;
   vector<double>  *CandidateP_EndX;
   vector<double>  *CandidateP_EndY;
   vector<double>  *CandidateP_EndZ;

   vector<double>  *True_CandidateP_StartX;
   vector<double>  *True_CandidateP_StartY;
   vector<double>  *True_CandidateP_StartZ;
   vector<double>  *True_CandidateP_EndX;
   vector<double>  *True_CandidateP_EndY;
   vector<double>  *True_CandidateP_EndZ;

   vector<double>  *CandidateP_P_Range;
   vector<double>  *CandidateP_P_MCS;   
   vector<double>  *CandidateP_Phi;
   vector<double>  *CandidateP_CosTheta;
   vector<double>  *CandidateP_Length;
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

   vector<double>  *Reco_DeltaPhi;
   vector<double>  *Reco_DeltaTheta;
   vector<double>  *True_DeltaPhi;
   vector<double>  *True_DeltaTheta;      

   vector<double>  *Reco_kMiss;
   vector<double>  *Reco_PMissMinus;
   vector<double>  *Reco_PMiss;
   vector<double>  *Reco_Pt;
   vector<double>  *Reco_DeltaAlphaT;
   vector<double>  *Reco_DeltaPhiT;
   vector<double>  *Reco_ECal;
   vector<double>  *Reco_EQE;
   vector<double>  *Reco_Q2;

   vector<double>  *True_PMissMinus;
   vector<double>  *True_kMiss;
   vector<double>  *True_PMiss;   
   vector<double>  *True_Pt;
   vector<double>  *True_DeltaAlphaT;
   vector<double>  *True_DeltaPhiT;
   vector<double>  *True_ECal;
   vector<double>  *True_EQE;
   vector<double>  *True_Q2;                     

//   Int_t           NuMuPFParticles;
//   vector<int>     *PFParticle_NuMuDaughters;
//   vector<vector<int> > *PFParticle_NuMuDaughtersPdgCode;

   // List of branches
//   TBranch        *b_PassedSwTrigger;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_T2KWeight;   //!
   TBranch        *b_ROOTinoWeight;   //! 

   TBranch        *b_All_UBGenie;   //!                           
   TBranch        *b_AxFFCCQEshape_UBGenie;   //!                 
   TBranch        *b_DecayAngMEC_UBGenie;   //!                   
   TBranch        *b_NormCCCOH_UBGenie;   //!                     
   TBranch        *b_NormNCCOH_UBGenie;   //!                     
//   TBranch        *b_RPA_CCQE_Reduced_UBGenie;   //!              
   TBranch        *b_RPA_CCQE_UBGenie;   //!                      
   TBranch        *b_ThetaDelta2NRad_UBGenie;   //!               
   TBranch        *b_Theta_Delta2Npi_UBGenie;   //!               
   TBranch        *b_VecFFCCQEshape_UBGenie;   //!                
   TBranch        *b_XSecShape_CCMEC_UBGenie;   //!               
   TBranch        *b_expskin_FluxUnisim;   //!                    
   TBranch        *b_horncurrent_FluxUnisim;   //!                
   TBranch        *b_kminus_PrimaryHadronNormalization;   //!     
   TBranch        *b_kplus_PrimaryHadronFeynmanScaling;   //!     
   TBranch        *b_kzero_PrimaryHadronSanfordWang;   //!        
   TBranch        *b_nucleoninexsec_FluxUnisim;   //!             
   TBranch        *b_nucleonqexsec_FluxUnisim;   //!              
   TBranch        *b_nucleontotxsec_FluxUnisim;   //!             
   TBranch        *b_piminus_PrimaryHadronSWCentralSplineVariation;   //!                                                             
   TBranch        *b_pioninexsec_FluxUnisim;   //!                
   TBranch        *b_pionqexsec_FluxUnisim;   //!                 
   TBranch        *b_piontotxsec_FluxUnisim;   //!                
   TBranch        *b_piplus_PrimaryHadronSWCentralSplineVariation;   //!                                                              
   TBranch        *b_reinteractions_piminus_Geant4;   //!         
   TBranch        *b_reinteractions_piplus_Geant4;   //! 
   TBranch        *b_reinteractions_proton_Geant4;   //!                                                       
//   TBranch        *b_xsr_scc_Fa3_SCC;   //!                                                                                                        
//   TBranch        *b_xsr_scc_Fv3_SCC;   //!
  
   TBranch        *b_CC1p;   //!
   TBranch        *b_CC1p1pi;   //!
   TBranch        *b_CC2p;   //!  
   TBranch        *b_CC2p1pi;   //!
   TBranch        *b_CC3p;   //!  
   TBranch        *b_CC3p1pi;   //!
   TBranch        *b_CC3p2pi;   //!   
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

   TBranch        *b_CandidateMu_StartX;   //!
   TBranch        *b_CandidateMu_StartY;   //!
   TBranch        *b_CandidateMu_StartZ;   //!
   TBranch        *b_CandidateMu_EndX;   //!
   TBranch        *b_CandidateMu_EndY;   //!
   TBranch        *b_CandidateMu_EndZ;   //!

   TBranch        *b_True_CandidateMu_StartX;   //!
   TBranch        *b_True_CandidateMu_StartY;   //!
   TBranch        *b_True_CandidateMu_StartZ;   //!
   TBranch        *b_True_CandidateMu_EndX;   //!
   TBranch        *b_True_CandidateMu_EndY;   //!
   TBranch        *b_True_CandidateMu_EndZ;   //! 
         
   TBranch        *b_CandidateMu_P_Range;   //!
   TBranch        *b_CandidateMu_P_MCS;   //!   
   TBranch        *b_CandidateMu_Phi;   //!
   TBranch        *b_CandidateMu_CosTheta;   //!
   TBranch        *b_CandidateMu_Length;   //!
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

   TBranch        *b_CandidateP_StartX;   //!
   TBranch        *b_CandidateP_StartY;   //!
   TBranch        *b_CandidateP_StartZ;   //!
   TBranch        *b_CandidateP_EndX;   //!
   TBranch        *b_CandidateP_EndY;   //!
   TBranch        *b_CandidateP_EndZ;   //!

   TBranch        *b_True_CandidateP_StartX;   //!
   TBranch        *b_True_CandidateP_StartY;   //!
   TBranch        *b_True_CandidateP_StartZ;   //!
   TBranch        *b_True_CandidateP_EndX;   //!
   TBranch        *b_True_CandidateP_EndY;   //!
   TBranch        *b_True_CandidateP_EndZ;   //!

   TBranch        *b_CandidateP_P_Range;   //!
   TBranch        *b_CandidateP_P_MCS;   //!   
   TBranch        *b_CandidateP_Phi;   //!
   TBranch        *b_CandidateP_CosTheta;   //!
   TBranch        *b_CandidateP_Length;   //!
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

   TBranch        *b_Reco_DeltaPhi;   //!
   TBranch        *b_Reco_DeltaTheta;   //!
   TBranch        *b_True_DeltaPhi;   //!
   TBranch        *b_True_DeltaTheta;   //!      

   TBranch        *b_Reco_kMiss;   //!
   TBranch        *b_Reco_PMissMinus;   //!
   TBranch        *b_Reco_PMiss;   //!   
   TBranch        *b_Reco_Pt;   //!
   TBranch        *b_Reco_DeltaAlphaT;   //!
   TBranch        *b_Reco_DeltaPhiT;   //!
   TBranch        *b_Reco_ECal;   //!
   TBranch        *b_Reco_EQE;   //!
   TBranch        *b_Reco_Q2;   //!

   TBranch        *b_True_kMiss;   //!
   TBranch        *b_True_PMissMinus;   //!
   TBranch        *b_True_PMiss;   //!   
   TBranch        *b_True_Pt;   //!
   TBranch        *b_True_DeltaAlphaT;   //!
   TBranch        *b_True_DeltaPhiT;   //!
   TBranch        *b_True_ECal;   //!
   TBranch        *b_True_EQE;   //!
   TBranch        *b_True_Q2;   //!   
                     

//   TBranch        *b_NuMuPFParticles;   //!
//   TBranch        *b_PFParticle_NuMuDaughters;   //!
//   TBranch        *b_PFParticle_NuMuDaughtersPdgCode;   //!

   myRecoAnalysis(TString WhichSample="",TString WhichEventWeightLabel="", int UniverseIndex=-1, TTree *tree=0);
   virtual ~myRecoAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myRecoAnalysis_cxx
myRecoAnalysis::myRecoAnalysis(TString WhichSample, TString WhichEventWeightLabel, int UniverseIndex, TTree *tree) : fChain(0) 
{

   fWhichSample = WhichSample;

   fEventWeightLabel = WhichEventWeightLabel;
   fUniverseIndex = UniverseIndex;

//   On the gpvm's
   fPathToFile = "/uboone/data/users/apapadop/myEvents/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
//   Locally
//   fPathToFile = "mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fPathToFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(fPathToFile);
      }
      f->GetObject("myPreSelection",tree);

   }
   Init(tree);
}

myRecoAnalysis::~myRecoAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myRecoAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myRecoAnalysis::LoadTree(Long64_t entry)
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

void myRecoAnalysis::Init(TTree *tree)
{

   // Set object pointer

   All_UBGenie = 0;
   AxFFCCQEshape_UBGenie = 0;
   DecayAngMEC_UBGenie = 0;
   NormCCCOH_UBGenie = 0;
   NormNCCOH_UBGenie = 0;
//   RPA_CCQE_Reduced_UBGenie = 0;
   RPA_CCQE_UBGenie = 0;
   ThetaDelta2NRad_UBGenie = 0;
   Theta_Delta2Npi_UBGenie = 0;
   VecFFCCQEshape_UBGenie = 0;
   XSecShape_CCMEC_UBGenie = 0;
   expskin_FluxUnisim = 0;
   horncurrent_FluxUnisim = 0;
   kminus_PrimaryHadronNormalization = 0;
   kplus_PrimaryHadronFeynmanScaling = 0;
   kzero_PrimaryHadronSanfordWang = 0;
   nucleoninexsec_FluxUnisim = 0;
   nucleonqexsec_FluxUnisim = 0;
   nucleontotxsec_FluxUnisim = 0;
   piminus_PrimaryHadronSWCentralSplineVariation = 0;
   pioninexsec_FluxUnisim = 0;
   pionqexsec_FluxUnisim = 0;
   piontotxsec_FluxUnisim = 0;
   piplus_PrimaryHadronSWCentralSplineVariation = 0;
   reinteractions_piminus_Geant4 = 0;
   reinteractions_piplus_Geant4 = 0;
   reinteractions_proton_Geant4 = 0;
//   xsr_scc_Fa3_SCC = 0;
//   xsr_scc_Fv3_SCC = 0;

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

   CandidateMu_StartX = 0;
   CandidateMu_StartY = 0;
   CandidateMu_StartZ = 0;
   CandidateMu_EndX = 0;
   CandidateMu_EndY = 0;
   CandidateMu_EndZ = 0;  

   True_CandidateMu_StartX = 0;
   True_CandidateMu_StartY = 0;
   True_CandidateMu_StartZ = 0;
   True_CandidateMu_EndX = 0;
   True_CandidateMu_EndY = 0;
   True_CandidateMu_EndZ = 0;
          
   CandidateMu_P_Range = 0;
   CandidateMu_P_MCS = 0;   
   CandidateMu_Phi = 0;
   CandidateMu_CosTheta = 0;
   CandidateMu_Length = 0;
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

   CandidateP_StartX = 0;
   CandidateP_StartY = 0;
   CandidateP_StartZ = 0;
   CandidateP_EndX = 0;
   CandidateP_EndY = 0;
   CandidateP_EndZ = 0;

   True_CandidateP_StartX = 0;
   True_CandidateP_StartY = 0;
   True_CandidateP_StartZ = 0;
   True_CandidateP_EndX = 0;
   True_CandidateP_EndY = 0;
   True_CandidateP_EndZ = 0;

   CandidateP_P_Range = 0;
   CandidateP_P_MCS = 0;   
   CandidateP_Phi = 0;
   CandidateP_CosTheta = 0;
   CandidateP_Length = 0;
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

   Reco_DeltaPhi = 0;
   Reco_DeltaTheta = 0;
   True_DeltaPhi = 0;
   True_DeltaTheta = 0;      

   Reco_PMissMinus = 0;
   Reco_kMiss = 0;
   Reco_PMiss = 0;   
   Reco_Pt = 0;
   Reco_DeltaAlphaT = 0;
   Reco_DeltaPhiT = 0;
   Reco_ECal = 0;
   Reco_EQE = 0;
   Reco_Q2 = 0;

   True_kMiss = 0;
   True_PMissMinus = 0;
   True_PMiss = 0;   
   True_Pt = 0;
   True_DeltaAlphaT = 0;
   True_DeltaPhiT = 0;
   True_ECal = 0;
   True_EQE = 0;
   True_Q2 = 0;                     

//   PFParticle_NuMuDaughters = 0;
//   PFParticle_NuMuDaughtersPdgCode = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

//   fChain->SetBranchAddress("PassedSwTrigger", &PassedSwTrigger, &b_PassedSwTrigger);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("T2KWeight", &T2KWeight, &b_T2KWeight);
   fChain->SetBranchAddress("ROOTinoWeight", &ROOTinoWeight, &b_ROOTinoWeight);
   
   fChain->SetBranchAddress("All_UBGenie", &All_UBGenie, &b_All_UBGenie);
   fChain->SetBranchAddress("AxFFCCQEshape_UBGenie", &AxFFCCQEshape_UBGenie, &b_AxFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("DecayAngMEC_UBGenie", &DecayAngMEC_UBGenie, &b_DecayAngMEC_UBGenie);
   fChain->SetBranchAddress("NormCCCOH_UBGenie", &NormCCCOH_UBGenie, &b_NormCCCOH_UBGenie);
   fChain->SetBranchAddress("NormNCCOH_UBGenie", &NormNCCOH_UBGenie, &b_NormNCCOH_UBGenie);
//   fChain->SetBranchAddress("RPA_CCQE_Reduced_UBGenie", &RPA_CCQE_Reduced_UBGenie, &b_RPA_CCQE_Reduced_UBGenie);
   fChain->SetBranchAddress("RPA_CCQE_UBGenie", &RPA_CCQE_UBGenie, &b_RPA_CCQE_UBGenie);
   fChain->SetBranchAddress("ThetaDelta2NRad_UBGenie", &ThetaDelta2NRad_UBGenie, &b_ThetaDelta2NRad_UBGenie);
   fChain->SetBranchAddress("Theta_Delta2Npi_UBGenie", &Theta_Delta2Npi_UBGenie, &b_Theta_Delta2Npi_UBGenie);
   fChain->SetBranchAddress("VecFFCCQEshape_UBGenie", &VecFFCCQEshape_UBGenie, &b_VecFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("XSecShape_CCMEC_UBGenie", &XSecShape_CCMEC_UBGenie, &b_XSecShape_CCMEC_UBGenie);
   fChain->SetBranchAddress("expskin_FluxUnisim", &expskin_FluxUnisim, &b_expskin_FluxUnisim);
   fChain->SetBranchAddress("horncurrent_FluxUnisim", &horncurrent_FluxUnisim, &b_horncurrent_FluxUnisim);
   fChain->SetBranchAddress("kminus_PrimaryHadronNormalization", &kminus_PrimaryHadronNormalization, &b_kminus_PrimaryHadronNormalization);
   fChain->SetBranchAddress("kplus_PrimaryHadronFeynmanScaling", &kplus_PrimaryHadronFeynmanScaling, &b_kplus_PrimaryHadronFeynmanScaling);
   fChain->SetBranchAddress("kzero_PrimaryHadronSanfordWang", &kzero_PrimaryHadronSanfordWang, &b_kzero_PrimaryHadronSanfordWang);
   fChain->SetBranchAddress("nucleoninexsec_FluxUnisim", &nucleoninexsec_FluxUnisim, &b_nucleoninexsec_FluxUnisim);
   fChain->SetBranchAddress("nucleonqexsec_FluxUnisim", &nucleonqexsec_FluxUnisim, &b_nucleonqexsec_FluxUnisim);
   fChain->SetBranchAddress("nucleontotxsec_FluxUnisim", &nucleontotxsec_FluxUnisim, &b_nucleontotxsec_FluxUnisim);
   fChain->SetBranchAddress("piminus_PrimaryHadronSWCentralSplineVariation", &piminus_PrimaryHadronSWCentralSplineVariation, &b_piminus_PrimaryHadronSWCentralSplineVariation);
   fChain->SetBranchAddress("pioninexsec_FluxUnisim", &pioninexsec_FluxUnisim, &b_pioninexsec_FluxUnisim);
   fChain->SetBranchAddress("pionqexsec_FluxUnisim", &pionqexsec_FluxUnisim, &b_pionqexsec_FluxUnisim);
   fChain->SetBranchAddress("piontotxsec_FluxUnisim", &piontotxsec_FluxUnisim, &b_piontotxsec_FluxUnisim);
   fChain->SetBranchAddress("piplus_PrimaryHadronSWCentralSplineVariation", &piplus_PrimaryHadronSWCentralSplineVariation, &b_piplus_PrimaryHadronSWCentralSplineVariation);
   fChain->SetBranchAddress("reinteractions_piminus_Geant4", &reinteractions_piminus_Geant4, &b_reinteractions_piminus_Geant4);
   fChain->SetBranchAddress("reinteractions_piplus_Geant4", &reinteractions_piplus_Geant4, &b_reinteractions_piplus_Geant4);
   fChain->SetBranchAddress("reinteractions_proton_Geant4", &reinteractions_proton_Geant4, &b_reinteractions_proton_Geant4);
//   fChain->SetBranchAddress("xsr_scc_Fa3_SCC", &xsr_scc_Fa3_SCC, &b_xsr_scc_Fa3_SCC);
//   fChain->SetBranchAddress("xsr_scc_Fv3_SCC", &xsr_scc_Fv3_SCC, &b_xsr_scc_Fv3_SCC);
      
   fChain->SetBranchAddress("CC1p", &CC1p, &b_CC1p);
   fChain->SetBranchAddress("CC1p1pi", &CC1p1pi, &b_CC1p1pi);
   fChain->SetBranchAddress("CC2p", &CC2p, &b_CC2p);
   fChain->SetBranchAddress("CC2p1pi", &CC2p1pi, &b_CC2p1pi);
   fChain->SetBranchAddress("CC3p", &CC3p, &b_CC3p);
   fChain->SetBranchAddress("CC3p1pi", &CC3p1pi, &b_CC3p1pi);
   fChain->SetBranchAddress("CC3p2pi", &CC3p2pi, &b_CC3p2pi);  
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

   fChain->SetBranchAddress("CandidateMu_StartX", &CandidateMu_StartX, &b_CandidateMu_StartX);
   fChain->SetBranchAddress("CandidateMu_StartY", &CandidateMu_StartY, &b_CandidateMu_StartY);
   fChain->SetBranchAddress("CandidateMu_StartZ", &CandidateMu_StartZ, &b_CandidateMu_StartZ);
   fChain->SetBranchAddress("CandidateMu_EndX", &CandidateMu_EndX, &b_CandidateMu_EndX);
   fChain->SetBranchAddress("CandidateMu_EndY", &CandidateMu_EndY, &b_CandidateMu_EndY);
   fChain->SetBranchAddress("CandidateMu_EndZ", &CandidateMu_EndZ, &b_CandidateMu_EndZ);

   fChain->SetBranchAddress("True_CandidateMu_StartX", &True_CandidateMu_StartX, &b_True_CandidateMu_StartX);
   fChain->SetBranchAddress("True_CandidateMu_StartY", &True_CandidateMu_StartY, &b_True_CandidateMu_StartY);
   fChain->SetBranchAddress("True_CandidateMu_StartZ", &True_CandidateMu_StartZ, &b_True_CandidateMu_StartZ);
   fChain->SetBranchAddress("True_CandidateMu_EndX", &True_CandidateMu_EndX, &b_True_CandidateMu_EndX);
   fChain->SetBranchAddress("True_CandidateMu_EndY", &True_CandidateMu_EndY, &b_True_CandidateMu_EndY);
   fChain->SetBranchAddress("True_CandidateMu_EndZ", &True_CandidateMu_EndZ, &b_True_CandidateMu_EndZ);    
     
   fChain->SetBranchAddress("CandidateMu_P_Range", &CandidateMu_P_Range, &b_CandidateMu_P_Range);
   fChain->SetBranchAddress("CandidateMu_P_MCS", &CandidateMu_P_MCS, &b_CandidateMu_P_MCS);   
   fChain->SetBranchAddress("CandidateMu_Phi", &CandidateMu_Phi, &b_CandidateMu_Phi);
   fChain->SetBranchAddress("CandidateMu_CosTheta", &CandidateMu_CosTheta, &b_CandidateMu_CosTheta);
   fChain->SetBranchAddress("CandidateMu_Length", &CandidateMu_Length, &b_CandidateMu_Length);
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

   fChain->SetBranchAddress("CandidateP_StartX", &CandidateP_StartX, &b_CandidateP_StartX);
   fChain->SetBranchAddress("CandidateP_StartY", &CandidateP_StartY, &b_CandidateP_StartY);
   fChain->SetBranchAddress("CandidateP_StartZ", &CandidateP_StartZ, &b_CandidateP_StartZ);
   fChain->SetBranchAddress("CandidateP_EndX", &CandidateP_EndX, &b_CandidateP_EndX);
   fChain->SetBranchAddress("CandidateP_EndY", &CandidateP_EndY, &b_CandidateP_EndY);
   fChain->SetBranchAddress("CandidateP_EndZ", &CandidateP_EndZ, &b_CandidateP_EndZ);

   fChain->SetBranchAddress("True_CandidateP_StartX", &True_CandidateP_StartX, &b_True_CandidateP_StartX);
   fChain->SetBranchAddress("True_CandidateP_StartY", &True_CandidateP_StartY, &b_True_CandidateP_StartY);
   fChain->SetBranchAddress("True_CandidateP_StartZ", &True_CandidateP_StartZ, &b_True_CandidateP_StartZ);
   fChain->SetBranchAddress("True_CandidateP_EndX", &True_CandidateP_EndX, &b_True_CandidateP_EndX);
   fChain->SetBranchAddress("True_CandidateP_EndY", &True_CandidateP_EndY, &b_True_CandidateP_EndY);
   fChain->SetBranchAddress("True_CandidateP_EndZ", &True_CandidateP_EndZ, &b_True_CandidateP_EndZ);

   fChain->SetBranchAddress("CandidateP_P_Range", &CandidateP_P_Range, &b_CandidateP_P_Range);
   fChain->SetBranchAddress("CandidateP_P_MCS", &CandidateP_P_MCS, &b_CandidateP_P_MCS);   
   fChain->SetBranchAddress("CandidateP_Phi", &CandidateP_Phi, &b_CandidateP_Phi);
   fChain->SetBranchAddress("CandidateP_CosTheta", &CandidateP_CosTheta, &b_CandidateP_CosTheta);
   fChain->SetBranchAddress("CandidateP_Length", &CandidateP_Length, &b_CandidateP_Length);
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

   fChain->SetBranchAddress("Reco_DeltaPhi", &Reco_DeltaPhi, &b_Reco_DeltaPhi);
   fChain->SetBranchAddress("Reco_DeltaTheta", &Reco_DeltaTheta, &b_Reco_DeltaTheta);
   fChain->SetBranchAddress("True_DeltaPhi", &True_DeltaPhi, &b_True_DeltaPhi);
   fChain->SetBranchAddress("True_DeltaTheta", &True_DeltaTheta, &b_True_DeltaTheta);      

   fChain->SetBranchAddress("Reco_kMiss", &Reco_kMiss, &b_Reco_kMiss);
   fChain->SetBranchAddress("Reco_PMissMinus", &Reco_PMissMinus, &b_Reco_PMissMinus);
   fChain->SetBranchAddress("Reco_PMiss", &Reco_PMiss, &b_Reco_PMiss);   
   fChain->SetBranchAddress("Reco_Pt", &Reco_Pt, &b_Reco_Pt);
   fChain->SetBranchAddress("Reco_DeltaAlphaT", &Reco_DeltaAlphaT, &b_Reco_DeltaAlphaT);
   fChain->SetBranchAddress("Reco_DeltaPhiT", &Reco_DeltaPhiT, &b_Reco_DeltaPhiT);
   fChain->SetBranchAddress("Reco_ECal", &Reco_ECal, &b_Reco_ECal);
   fChain->SetBranchAddress("Reco_EQE", &Reco_EQE, &b_Reco_EQE);
   fChain->SetBranchAddress("Reco_Q2", &Reco_Q2, &b_Reco_Q2);

   fChain->SetBranchAddress("True_kMiss", &True_kMiss, &b_True_kMiss);
   fChain->SetBranchAddress("True_PMissMinus", &True_PMissMinus, &b_True_PMissMinus);
   fChain->SetBranchAddress("True_PMiss", &True_PMiss, &b_True_PMiss);   
   fChain->SetBranchAddress("True_Pt", &True_Pt, &b_True_Pt);
   fChain->SetBranchAddress("True_DeltaAlphaT", &True_DeltaAlphaT, &b_True_DeltaAlphaT);
   fChain->SetBranchAddress("True_DeltaPhiT", &True_DeltaPhiT, &b_True_DeltaPhiT);
   fChain->SetBranchAddress("True_ECal", &True_ECal, &b_True_ECal);
   fChain->SetBranchAddress("True_EQE", &True_EQE, &b_True_EQE);
   fChain->SetBranchAddress("True_Q2", &True_Q2, &b_True_Q2);                     

//   fChain->SetBranchAddress("NuMuPFParticles", &NuMuPFParticles, &b_NuMuPFParticles);
//   fChain->SetBranchAddress("PFParticle_NuMuDaughters", &PFParticle_NuMuDaughters, &b_PFParticle_NuMuDaughters);
//   fChain->SetBranchAddress("PFParticle_NuMuDaughtersPdgCode", &PFParticle_NuMuDaughtersPdgCode, &b_PFParticle_NuMuDaughtersPdgCode);

   Notify();
}

Bool_t myRecoAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myRecoAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myRecoAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myRecoAnalysis_cxx
