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
   TString fEventWeightLabel;
   int     fUniverseIndex;   

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   double          Weight;
   double          T2KWeight;
   double          ROOTinoWeight;   

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
//   vector<double>  *Muon_MCParticle_Length;
   vector<int>     *Muon_MCParticle_StartContainment;
   vector<int>     *Muon_MCParticle_EndContainment;
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
//   vector<double>  *Proton_MCParticle_Length;
   vector<int>     *Proton_MCParticle_StartContainment;
   vector<int>     *Proton_MCParticle_EndContainment;
   vector<int>     *Proton_MCParticle_Pdg;   

   vector<double>  *True_DeltaPhi;
   vector<double>  *True_DeltaTheta;   

   vector<double>  *True_kMiss;
   vector<double>  *True_PMissMinus;
   vector<double>  *True_PMiss;
   
   vector<double>  *True_Pt;
   vector<double>  *True_DeltaAlphaT;
   vector<double>  *True_DeltaPhiT;
   vector<double>  *True_ECal;
   vector<double>  *True_EQE;
   vector<double>  *True_Q2;         

   // List of branches
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
//   TBranch        *b_Muon_MCParticle_Length;   //!
   TBranch        *b_Muon_MCParticle_StartContainment;   //!
   TBranch        *b_Muon_MCParticle_EndContainment;   //!
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
//   TBranch        *b_Proton_MCParticle_Length;   //!
   TBranch        *b_Proton_MCParticle_StartContainment;   //!
   TBranch        *b_Proton_MCParticle_EndContainment;   //!
   TBranch        *b_Proton_MCParticle_Pdg;   //!   
   
   TBranch        *b_True_DeltaPhi;   //!
   TBranch        *b_True_DeltaTheta;   //!      

   TBranch        *b_True_kMiss;   //!
   TBranch        *b_True_PMissMinus;   //!
   TBranch        *b_True_PMiss;   //!
   
   TBranch        *b_True_Pt;   //!
   TBranch        *b_True_DeltaAlphaT;   //!
   TBranch        *b_True_DeltaPhiT;   //!
   TBranch        *b_True_ECal;   //!
   TBranch        *b_True_EQE;   //!
   TBranch        *b_True_Q2;   //!                  

   myTrueAnalysis(TString WhichSample="",TString WhichEventWeightLabel="", int UniverseIndex=-1,TTree *tree=0);
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
myTrueAnalysis::myTrueAnalysis(TString WhichSample, TString WhichEventWeightLabel, int UniverseIndex, TTree *tree) : fChain(0) 
{

   fWhichSample = WhichSample;
   fEventWeightLabel = WhichEventWeightLabel;
   fUniverseIndex = UniverseIndex;   

// LOcally
//   TString PathToFile = "mySamples/"+UBCodeVersion+"/PreTruthSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

// On the gpvm's
   TString PathToFile = "/uboone/data/users/apapadop/myEvents/mySamples/"+UBCodeVersion+"/PreTruthSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

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

   Muon_MCParticle_StartX = 0;
   Muon_MCParticle_StartY = 0;
   Muon_MCParticle_StartZ = 0;
   Muon_MCParticle_EndX = 0;
   Muon_MCParticle_EndY = 0;
   Muon_MCParticle_EndZ = 0;            
   Muon_MCParticle_Mom = 0;
   Muon_MCParticle_Phi = 0;
   Muon_MCParticle_CosTheta = 0;
//   Muon_MCParticle_Length = 0;
   Muon_MCParticle_StartContainment = 0;
   Muon_MCParticle_EndContainment = 0;
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
//   Proton_MCParticle_Length = 0;
   Proton_MCParticle_StartContainment = 0;
   Proton_MCParticle_EndContainment = 0;
   Proton_MCParticle_Pdg = 0;

   True_DeltaPhi = 0;
   True_DeltaTheta = 0;   

   True_kMiss = 0;
   True_PMissMinus = 0;
   True_PMiss = 0;
      
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
//   fChain->SetBranchAddress("Muon_MCParticle_Length", &Muon_MCParticle_Length, &b_Muon_MCParticle_Length);
   fChain->SetBranchAddress("Muon_MCParticle_StartContainment", &Muon_MCParticle_StartContainment, &b_Muon_MCParticle_StartContainment);
   fChain->SetBranchAddress("Muon_MCParticle_EndContainment", &Muon_MCParticle_EndContainment, &b_Muon_MCParticle_EndContainment);
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
//   fChain->SetBranchAddress("Proton_MCParticle_Length", &Proton_MCParticle_Length, &b_Proton_MCParticle_Length);
   fChain->SetBranchAddress("Proton_MCParticle_StartContainment", &Proton_MCParticle_StartContainment, &b_Proton_MCParticle_StartContainment);
   fChain->SetBranchAddress("Proton_MCParticle_EndContainment", &Proton_MCParticle_EndContainment, &b_Proton_MCParticle_EndContainment);
   fChain->SetBranchAddress("Proton_MCParticle_Pdg", &Proton_MCParticle_Pdg, &b_Proton_MCParticle_Pdg);   

   fChain->SetBranchAddress("True_DeltaPhi", &True_DeltaPhi, &b_True_DeltaPhi);
   fChain->SetBranchAddress("True_DeltaTheta", &True_DeltaTheta, &b_True_DeltaTheta);   

   fChain->SetBranchAddress("True_kMiss", &True_kMiss, &b_True_kMiss);
   fChain->SetBranchAddress("True_PMissMinus", &True_PMissMinus, &b_True_PMissMinus);
   fChain->SetBranchAddress("True_PMiss", &True_PMiss, &b_True_PMiss);
   
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
