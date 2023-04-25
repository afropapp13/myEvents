#ifndef PeLEE_myRecoAnalysis_h
#define PeLEE_myRecoAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include "../../myClasses/Constants.h"

#include <vector>
#include <vector>

using namespace Constants;

class PeLEE_myRecoAnalysis {

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

	vector<unsigned short> *AGKYpT1pi_UBGenie;
	vector<unsigned short> *AGKYxF1pi_UBGenie;
	vector<unsigned short> *AhtBY_UBGenie;
	vector<unsigned short> *BhtBY_UBGenie;
	vector<unsigned short> *CV1uBY_UBGenie;
	vector<unsigned short> *CV2uBY_UBGenie;
	vector<unsigned short> *EtaNCEL_UBGenie;
	vector<unsigned short> *FrAbs_N_UBGenie;
	vector<unsigned short> *FrAbs_pi_UBGenie;
	vector<unsigned short> *FrCEx_N_UBGenie;
	vector<unsigned short> *FrCEx_pi_UBGenie;
	vector<unsigned short> *FrInel_N_UBGenie;
	vector<unsigned short> *FrInel_pi_UBGenie;
	vector<unsigned short> *FrPiProd_N_UBGenie;
	vector<unsigned short> *FrPiProd_pi_UBGenie;
	vector<unsigned short> *FracDelta_CCMEC_UBGenie;	
	vector<unsigned short> *FracPN_CCMEC_UBGenie;
	vector<unsigned short> *MFP_N_UBGenie;
	vector<unsigned short> *MFP_pi_UBGenie;
	vector<unsigned short> *MaCCQE_UBGenie;
	vector<unsigned short> *MaCCRES_UBGenie;
	vector<unsigned short> *MaNCEL_UBGenie;
	vector<unsigned short> *MaNCRES_UBGenie;
	vector<unsigned short> *MvCCRES_UBGenie;
	vector<unsigned short> *MvNCRES_UBGenie;
	vector<unsigned short> *NonRESBGvbarnCC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarnCC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarnNC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarnNC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarpCC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarpCC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarpNC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvbarpNC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvnCC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvnCC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvnNC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvnNC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvpCC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvpCC2pi_UBGenie;
	vector<unsigned short> *NonRESBGvpNC1pi_UBGenie;
	vector<unsigned short> *NonRESBGvpNC2pi_UBGenie;
	vector<unsigned short> *NormCCMEC_UBGenie;
	vector<unsigned short> *NormNCMEC_UBGenie;
	vector<unsigned short> *RDecBR1eta_UBGenie;
	vector<unsigned short> *RDecBR1gamma_UBGenie;	
	
	// Unisims
	
	vector<unsigned short> *UnShortAxFFCCQEshape_UBGenie;
	vector<unsigned short> *UnShortDecayAngMEC_UBGenie;
	vector<unsigned short> *UnShortRPA_CCQE_UBGenie;
	vector<unsigned short> *UnShortTheta_Delta2Npi_UBGenie;
	vector<unsigned short> *UnShortVecFFCCQEshape_UBGenie;
	vector<unsigned short> *UnShortXSecShape_CCMEC_UBGenie;		

	//----------------------------------------//   
   
   
   vector<double>  *fluxes;
   vector<double>  *reinteractions;
   Int_t           nue;
   Int_t           NC;
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
   Int_t           MCParticle_Mode;
   Int_t           NumberPi0;
   Int_t           NumberNeutrons;
   Int_t           NumberProtons;
   Int_t           NumberMuons;
   Int_t           NumberChargedPions;
   Double_t        True_Ev;
   Double_t        True_Vx;
   Double_t        True_Vy;
   Double_t        True_Vz;
   Float_t         NuScore;
   Float_t         FlashScore;
   vector<double>  *BeamFlashes_TotalPE;
   vector<double>  *BeamFlashes_Time;
   vector<double>  *CandidateMuP_Distance;
   vector<float>   *Vertex_X;
   vector<float>   *Vertex_Y;
   vector<float>   *Vertex_Z;
   vector<double>  *CandidateMuStartVertexDistance;
   vector<double>  *CandidatePStartVertexDistance;
   vector<double>  *CandidateMuEndVertexDistance;
   vector<double>  *CandidatePEndVertexDistance;
   vector<double>  *CandidateMu_P_Range;
   vector<double>  *CandidateMu_P_MCS;
   vector<double>  *CandidateMu_Phi;
   vector<double>  *CandidateMu_Theta;
   vector<double>  *CandidateMu_CosTheta;
   vector<double>  *CandidateMu_LLR_PID;
   vector<int>     *CandidateMu_StartContainment;
   vector<int>     *CandidateMu_EndContainment;
   vector<double>  *CandidateMu_Length;
   vector<int>     *CandidateMu_MCParticle_Pdg;
   vector<double>  *CandidateMu_MCParticle_Purity;
   vector<double>  *CandidateMu_StartX;
   vector<double>  *CandidateMu_StartY;
   vector<double>  *CandidateMu_StartZ;
   vector<double>  *CandidateMu_EndX;
   vector<double>  *CandidateMu_EndY;
   vector<double>  *CandidateMu_EndZ;
   vector<double>  *True_CandidateMu_P;
   vector<double>  *True_CandidateMu_Phi;
   vector<double>  *True_CandidateMu_Theta;
   vector<double>  *True_CandidateMu_CosTheta;
   vector<double>  *True_CandidateMu_StartX;
   vector<double>  *True_CandidateMu_StartY;
   vector<double>  *True_CandidateMu_StartZ;
   vector<double>  *True_CandidateMu_EndX;
   vector<double>  *True_CandidateMu_EndY;
   vector<double>  *True_CandidateMu_EndZ;
   vector<int>     *True_CandidateMu_StartContainment;
   vector<int>     *True_CandidateMu_EndContainment;
   vector<double>  *CandidateP_P_Range;
   vector<double>  *CandidateP_Phi;
   vector<double>  *CandidateP_Theta;
   vector<double>  *CandidateP_CosTheta;
   vector<double>  *CandidateP_LLR_PID;
   vector<int>     *CandidateP_StartContainment;
   vector<int>     *CandidateP_EndContainment;
   vector<double>  *CandidateP_Length;
   vector<int>     *CandidateP_MCParticle_Pdg;
   vector<double>  *CandidateP_MCParticle_Purity;
   vector<double>  *CandidateP_StartX;
   vector<double>  *CandidateP_StartY;
   vector<double>  *CandidateP_StartZ;
   vector<double>  *CandidateP_EndX;
   vector<double>  *CandidateP_EndY;
   vector<double>  *CandidateP_EndZ;
   vector<double>  *True_CandidateP_P;
   vector<double>  *True_CandidateP_Phi;
   vector<double>  *True_CandidateP_Theta;
   vector<double>  *True_CandidateP_CosTheta;
   vector<double>  *True_CandidateP_StartX;
   vector<double>  *True_CandidateP_StartY;
   vector<double>  *True_CandidateP_StartZ;
   vector<double>  *True_CandidateP_EndX;
   vector<double>  *True_CandidateP_EndY;
   vector<double>  *True_CandidateP_EndZ;
   vector<int>     *True_CandidateP_StartContainment;
   vector<int>     *True_CandidateP_EndContainment;
   vector<double>  *Reco_A;
   vector<double>  *Reco_kMiss;
   vector<double>  *Reco_PMissMinus;
   vector<double>  *Reco_EMiss;
   vector<double>  *Reco_PMiss;
   vector<double>  *Reco_Pt;
   vector<double>  *Reco_Ptx;
   vector<double>  *Reco_Pty;
   vector<double>  *Reco_PL;
   vector<double>  *Reco_Pn;
   vector<double>  *Reco_PnPerp;
   vector<double>  *Reco_PnPerpx;
   vector<double>  *Reco_PnPerpy;      
   vector<double>  *Reco_PnPar;      
   vector<double>  *Reco_DeltaAlphaT;
   vector<double>  *Reco_DeltaAlpha3Dq;
   vector<double>  *Reco_DeltaAlpha3DMu;      
   vector<double>  *Reco_DeltaPhiT;
   vector<double>  *Reco_DeltaPhi3D;   
   vector<double>  *Reco_ECal;
   vector<double>  *Reco_EQE;
   vector<double>  *Reco_Q2;
   vector<double>  *Reco_DeltaPhi;
   vector<double>  *Reco_DeltaTheta;
   vector<double>  *True_A;
   vector<double>  *True_kMiss;
   vector<double>  *True_PMissMinus;
   vector<double>  *True_EMiss;
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
   vector<double>  *StartToStartDistance;
   vector<double>  *EndToEndDistance;

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
   
	// Unisims
	
   TBranch        *b_UnShortAxFFCCQEshape_UBGenie;
   TBranch        *b_UnShortDecayAngMEC_UBGenie;
   TBranch        *b_UnShortRPA_CCQE_UBGenie;
   TBranch        *b_UnShortTheta_Delta2Npi_UBGenie;
   TBranch        *b_UnShortVecFFCCQEshape_UBGenie;
   TBranch        *b_UnShortXSecShape_CCMEC_UBGenie;               			

	//----------------------------------------//   
   
   TBranch        *b_fluxes;   //!
   TBranch        *b_reinteractions;   //!
   TBranch        *b_nue;   //!
   TBranch        *b_NC;   //!
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
   TBranch        *b_MCParticle_Mode;   //!
   TBranch        *b_NumberPi0;   //!
   TBranch        *b_NumberNeutrons;   //!
   TBranch        *b_NumberProtons;   //!
   TBranch        *b_NumberMuons;   //!
   TBranch        *b_NumberChargedPions;   //!
   TBranch        *b_True_Ev;   //!
   TBranch        *b_True_Vx;   //!
   TBranch        *b_True_Vy;   //!
   TBranch        *b_True_Vz;   //!
   TBranch        *b_NuScore;   //!
   TBranch        *b_FlashScore;   //!
   TBranch        *b_BeamFlashes_TotalPE;   //!
   TBranch        *b_BeamFlashes_Time;   //!
   TBranch        *b_CandidateMuP_Distance;   //!
   TBranch        *b_Vertex_X;   //!
   TBranch        *b_Vertex_Y;   //!
   TBranch        *b_Vertex_Z;   //!
   TBranch        *b_CandidateMuStartVertexDistance;   //!
   TBranch        *b_CandidatePStartVertexDistance;   //!
   TBranch        *b_CandidateMuEndVertexDistance;   //!
   TBranch        *b_CandidatePEndVertexDistance;   //!
   TBranch        *b_CandidateMu_P_Range;   //!
   TBranch        *b_CandidateMu_P_MCS;   //!
   TBranch        *b_CandidateMu_Phi;   //!
   TBranch        *b_CandidateMu_Theta;   //!
   TBranch        *b_CandidateMu_CosTheta;   //!
   TBranch        *b_CandidateMu_LLR_PID;   //!
   TBranch        *b_CandidateMu_StartContainment;   //!
   TBranch        *b_CandidateMu_EndContainment;   //!
   TBranch        *b_CandidateMu_Length;   //!
   TBranch        *b_CandidateMu_MCParticle_Pdg;   //!
   TBranch        *b_CandidateMu_MCParticle_Purity;   //!
   TBranch        *b_CandidateMu_StartX;   //!
   TBranch        *b_CandidateMu_StartY;   //!
   TBranch        *b_CandidateMu_StartZ;   //!
   TBranch        *b_CandidateMu_EndX;   //!
   TBranch        *b_CandidateMu_EndY;   //!
   TBranch        *b_CandidateMu_EndZ;   //!
   TBranch        *b_True_CandidateMu_P;   //!
   TBranch        *b_True_CandidateMu_Phi;   //!
   TBranch        *b_True_CandidateMu_Theta;   //!
   TBranch        *b_True_CandidateMu_CosTheta;   //!
   TBranch        *b_True_CandidateMu_StartX;   //!
   TBranch        *b_True_CandidateMu_StartY;   //!
   TBranch        *b_True_CandidateMu_StartZ;   //!
   TBranch        *b_True_CandidateMu_EndX;   //!
   TBranch        *b_True_CandidateMu_EndY;   //!
   TBranch        *b_True_CandidateMu_EndZ;   //!
   TBranch        *b_True_CandidateMu_StartContainment;   //!
   TBranch        *b_True_CandidateMu_EndContainment;   //!
   TBranch        *b_CandidateP_P_Range;   //!
   TBranch        *b_CandidateP_Phi;   //!
   TBranch        *b_CandidateP_Theta;   //!
   TBranch        *b_CandidateP_CosTheta;   //!
   TBranch        *b_CandidateP_LLR_PID;   //!
   TBranch        *b_CandidateP_StartContainment;   //!
   TBranch        *b_CandidateP_EndContainment;   //!
   TBranch        *b_CandidateP_Length;   //!
   TBranch        *b_CandidateP_MCParticle_Pdg;   //!
   TBranch        *b_CandidateP_MCParticle_Purity;   //!
   TBranch        *b_CandidateP_StartX;   //!
   TBranch        *b_CandidateP_StartY;   //!
   TBranch        *b_CandidateP_StartZ;   //!
   TBranch        *b_CandidateP_EndX;   //!
   TBranch        *b_CandidateP_EndY;   //!
   TBranch        *b_CandidateP_EndZ;   //!
   TBranch        *b_True_CandidateP_P;   //!
   TBranch        *b_True_CandidateP_Phi;   //!
   TBranch        *b_True_CandidateP_Theta;   //!
   TBranch        *b_True_CandidateP_CosTheta;   //!
   TBranch        *b_True_CandidateP_StartX;   //!
   TBranch        *b_True_CandidateP_StartY;   //!
   TBranch        *b_True_CandidateP_StartZ;   //!
   TBranch        *b_True_CandidateP_EndX;   //!
   TBranch        *b_True_CandidateP_EndY;   //!
   TBranch        *b_True_CandidateP_EndZ;   //!
   TBranch        *b_True_CandidateP_StartContainment;   //!
   TBranch        *b_True_CandidateP_EndContainment;   //!
   TBranch        *b_Reco_A;   //!
   TBranch        *b_Reco_kMiss;   //!
   TBranch        *b_Reco_PMissMinus;   //!
   TBranch        *b_Reco_EMiss;   //!
   TBranch        *b_Reco_PMiss;   //!
   TBranch        *b_Reco_Pt;   //!
   TBranch        *b_Reco_Ptx;   //!
   TBranch        *b_Reco_Pty;   //!
   TBranch        *b_Reco_PL;   //!
   TBranch        *b_Reco_Pn;   //!
   TBranch        *b_Reco_PnPerp;   //!
   TBranch        *b_Reco_PnPerpx;   //!   
   TBranch        *b_Reco_PnPerpy;   //!   
   TBranch        *b_Reco_PnPar;   //!      
   TBranch        *b_Reco_DeltaAlphaT;   //!
   TBranch        *b_Reco_DeltaAlpha3Dq;   //!
   TBranch        *b_Reco_DeltaAlpha3DMu;   //!      
   TBranch        *b_Reco_DeltaPhiT;   //!
   TBranch        *b_Reco_DeltaPhi3D;   //!   
   TBranch        *b_Reco_ECal;   //!
   TBranch        *b_Reco_EQE;   //!
   TBranch        *b_Reco_Q2;   //!
   TBranch        *b_Reco_DeltaPhi;   //!
   TBranch        *b_Reco_DeltaTheta;   //!
   TBranch        *b_True_A;   //!
   TBranch        *b_True_kMiss;   //!
   TBranch        *b_True_PMissMinus;   //!
   TBranch        *b_True_EMiss;   //!
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
   TBranch        *b_StartToStartDistance;   //!
   TBranch        *b_EndToEndDistance;   //!

   PeLEE_myRecoAnalysis(TString WhichSample="",TString Tune="",TString WhichEventWeightLabel="", int UniverseIndex=-1, TTree *tree=0);
   virtual ~PeLEE_myRecoAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PeLEE_myRecoAnalysis_cxx
PeLEE_myRecoAnalysis::PeLEE_myRecoAnalysis(TString WhichSample, TString Tune, TString WhichEventWeightLabel, int UniverseIndex, TTree *tree) : fChain(0) 
{

   fTune = Tune;
   fWhichSample = WhichSample;

   fEventWeightLabel = WhichEventWeightLabel;
   fUniverseIndex = UniverseIndex;

//   fPathToFile = "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PeLEETuples/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
	fPathToFile = "/uboone/data/users/apapadop/PeLEETuples/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fPathToFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(fPathToFile);
      }
      f->GetObject("myPreSelection",tree);
      fFile = f;

   }
   Init(tree);
}

PeLEE_myRecoAnalysis::~PeLEE_myRecoAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PeLEE_myRecoAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PeLEE_myRecoAnalysis::LoadTree(Long64_t entry)
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

void PeLEE_myRecoAnalysis::Init(TTree *tree)
{

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
	
	// Unisims
	
	UnShortAxFFCCQEshape_UBGenie = 0;
	UnShortDecayAngMEC_UBGenie = 0;
	UnShortRPA_CCQE_UBGenie = 0;
	UnShortTheta_Delta2Npi_UBGenie = 0;
	UnShortVecFFCCQEshape_UBGenie = 0;
	UnShortXSecShape_CCMEC_UBGenie = 0;	 	 	 	 	 	   		

	//----------------------------------------//     
   
   fluxes = 0;
   reinteractions = 0;
   BeamFlashes_TotalPE = 0;
   BeamFlashes_Time = 0;
   CandidateMuP_Distance = 0;
   Vertex_X = 0;
   Vertex_Y = 0;
   Vertex_Z = 0;
   CandidateMuStartVertexDistance = 0;
   CandidatePStartVertexDistance = 0;
   CandidateMuEndVertexDistance = 0;
   CandidatePEndVertexDistance = 0;
   CandidateMu_P_Range = 0;
   CandidateMu_P_MCS = 0;
   CandidateMu_Phi = 0;
   CandidateMu_Theta = 0;
   CandidateMu_CosTheta = 0;
   CandidateMu_LLR_PID = 0;
   CandidateMu_StartContainment = 0;
   CandidateMu_EndContainment = 0;
   CandidateMu_Length = 0;
   CandidateMu_MCParticle_Pdg = 0;
   CandidateMu_MCParticle_Purity = 0;
   CandidateMu_StartX = 0;
   CandidateMu_StartY = 0;
   CandidateMu_StartZ = 0;
   CandidateMu_EndX = 0;
   CandidateMu_EndY = 0;
   CandidateMu_EndZ = 0;
   True_CandidateMu_P = 0;
   True_CandidateMu_Phi = 0;
   True_CandidateMu_Theta = 0;
   True_CandidateMu_CosTheta = 0;
   True_CandidateMu_StartX = 0;
   True_CandidateMu_StartY = 0;
   True_CandidateMu_StartZ = 0;
   True_CandidateMu_EndX = 0;
   True_CandidateMu_EndY = 0;
   True_CandidateMu_EndZ = 0;
   True_CandidateMu_StartContainment = 0;
   True_CandidateMu_EndContainment = 0;
   CandidateP_P_Range = 0;
   CandidateP_Phi = 0;
   CandidateP_Theta = 0;
   CandidateP_CosTheta = 0;
   CandidateP_LLR_PID = 0;
   CandidateP_StartContainment = 0;
   CandidateP_EndContainment = 0;
   CandidateP_Length = 0;
   CandidateP_MCParticle_Pdg = 0;
   CandidateP_MCParticle_Purity = 0;
   CandidateP_StartX = 0;
   CandidateP_StartY = 0;
   CandidateP_StartZ = 0;
   CandidateP_EndX = 0;
   CandidateP_EndY = 0;
   CandidateP_EndZ = 0;
   True_CandidateP_P = 0;
   True_CandidateP_Phi = 0;
   True_CandidateP_Theta = 0;
   True_CandidateP_CosTheta = 0;
   True_CandidateP_StartX = 0;
   True_CandidateP_StartY = 0;
   True_CandidateP_StartZ = 0;
   True_CandidateP_EndX = 0;
   True_CandidateP_EndY = 0;
   True_CandidateP_EndZ = 0;
   True_CandidateP_StartContainment = 0;
   True_CandidateP_EndContainment = 0;
   Reco_A = 0;
   Reco_kMiss = 0;
   Reco_PMissMinus = 0;
   Reco_EMiss = 0;
   Reco_PMiss = 0;
   Reco_Pt = 0;
   Reco_Ptx = 0;
   Reco_Pty = 0;
   Reco_PL = 0;
   Reco_Pn = 0;
   Reco_PnPerp = 0;
   Reco_PnPerpx = 0;
   Reco_PnPerpy = 0;      
   Reco_PnPar = 0;      
   Reco_DeltaAlphaT = 0;
   Reco_DeltaAlpha3Dq = 0;
   Reco_DeltaAlpha3DMu = 0;      
   Reco_DeltaPhiT = 0;
   Reco_DeltaPhi3D = 0;   
   Reco_ECal = 0;
   Reco_EQE = 0;
   Reco_Q2 = 0;
   Reco_DeltaPhi = 0;
   Reco_DeltaTheta = 0;
   True_A = 0;
   True_kMiss = 0;
   True_PMissMinus = 0;
   True_EMiss = 0;
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
   StartToStartDistance = 0;
   EndToEndDistance = 0;
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
		
		// Unisims
		
		fChain->SetBranchAddress("UnShortAxFFCCQEshape_UBGenie", &UnShortAxFFCCQEshape_UBGenie, &b_UnShortAxFFCCQEshape_UBGenie);
		fChain->SetBranchAddress("UnShortDecayAngMEC_UBGenie", &UnShortDecayAngMEC_UBGenie, &b_UnShortDecayAngMEC_UBGenie);
		fChain->SetBranchAddress("UnShortRPA_CCQE_UBGenie", &UnShortRPA_CCQE_UBGenie, &b_UnShortRPA_CCQE_UBGenie);
		fChain->SetBranchAddress("UnShortTheta_Delta2Npi_UBGenie", &UnShortTheta_Delta2Npi_UBGenie, &b_UnShortTheta_Delta2Npi_UBGenie);
		fChain->SetBranchAddress("UnShortVecFFCCQEshape_UBGenie", &UnShortVecFFCCQEshape_UBGenie, &b_UnShortVecFFCCQEshape_UBGenie);
		fChain->SetBranchAddress("UnShortXSecShape_CCMEC_UBGenie", &UnShortXSecShape_CCMEC_UBGenie, &b_UnShortXSecShape_CCMEC_UBGenie);		

	}

	//----------------------------------------//   
   
   fChain->SetBranchAddress("fluxes", &fluxes, &b_fluxes);
   fChain->SetBranchAddress("reinteractions", &reinteractions, &b_reinteractions);
   fChain->SetBranchAddress("nue", &nue, &b_nue);
   fChain->SetBranchAddress("NC", &NC, &b_NC);
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
   fChain->SetBranchAddress("MCParticle_Mode", &MCParticle_Mode, &b_MCParticle_Mode);
   fChain->SetBranchAddress("NumberPi0", &NumberPi0, &b_NumberPi0);
   fChain->SetBranchAddress("NumberNeutrons", &NumberNeutrons, &b_NumberNeutrons);
   fChain->SetBranchAddress("NumberProtons", &NumberProtons, &b_NumberProtons);
   fChain->SetBranchAddress("NumberMuons", &NumberMuons, &b_NumberMuons);
   fChain->SetBranchAddress("NumberChargedPions", &NumberChargedPions, &b_NumberChargedPions);
   fChain->SetBranchAddress("True_Ev", &True_Ev, &b_True_Ev);
   fChain->SetBranchAddress("True_Vx", &True_Vx, &b_True_Vx);
   fChain->SetBranchAddress("True_Vy", &True_Vy, &b_True_Vy);
   fChain->SetBranchAddress("True_Vz", &True_Vz, &b_True_Vz);
   fChain->SetBranchAddress("NuScore", &NuScore, &b_NuScore);
   fChain->SetBranchAddress("FlashScore", &FlashScore, &b_FlashScore);
   fChain->SetBranchAddress("BeamFlashes_TotalPE", &BeamFlashes_TotalPE, &b_BeamFlashes_TotalPE);
   fChain->SetBranchAddress("BeamFlashes_Time", &BeamFlashes_Time, &b_BeamFlashes_Time);
   fChain->SetBranchAddress("CandidateMuP_Distance", &CandidateMuP_Distance, &b_CandidateMuP_Distance);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("CandidateMuStartVertexDistance", &CandidateMuStartVertexDistance, &b_CandidateMuStartVertexDistance);
   fChain->SetBranchAddress("CandidatePStartVertexDistance", &CandidatePStartVertexDistance, &b_CandidatePStartVertexDistance);
   fChain->SetBranchAddress("CandidateMuEndVertexDistance", &CandidateMuEndVertexDistance, &b_CandidateMuEndVertexDistance);
   fChain->SetBranchAddress("CandidatePEndVertexDistance", &CandidatePEndVertexDistance, &b_CandidatePEndVertexDistance);
   fChain->SetBranchAddress("CandidateMu_P_Range", &CandidateMu_P_Range, &b_CandidateMu_P_Range);
   fChain->SetBranchAddress("CandidateMu_P_MCS", &CandidateMu_P_MCS, &b_CandidateMu_P_MCS);
   fChain->SetBranchAddress("CandidateMu_Phi", &CandidateMu_Phi, &b_CandidateMu_Phi);
   fChain->SetBranchAddress("CandidateMu_Theta", &CandidateMu_Theta, &b_CandidateMu_Theta);
   fChain->SetBranchAddress("CandidateMu_CosTheta", &CandidateMu_CosTheta, &b_CandidateMu_CosTheta);
   fChain->SetBranchAddress("CandidateMu_LLR_PID", &CandidateMu_LLR_PID, &b_CandidateMu_LLR_PID);
   fChain->SetBranchAddress("CandidateMu_StartContainment", &CandidateMu_StartContainment, &b_CandidateMu_StartContainment);
   fChain->SetBranchAddress("CandidateMu_EndContainment", &CandidateMu_EndContainment, &b_CandidateMu_EndContainment);
   fChain->SetBranchAddress("CandidateMu_Length", &CandidateMu_Length, &b_CandidateMu_Length);
   fChain->SetBranchAddress("CandidateMu_MCParticle_Pdg", &CandidateMu_MCParticle_Pdg, &b_CandidateMu_MCParticle_Pdg);
   fChain->SetBranchAddress("CandidateMu_MCParticle_Purity", &CandidateMu_MCParticle_Purity, &b_CandidateMu_MCParticle_Purity);
   fChain->SetBranchAddress("CandidateMu_StartX", &CandidateMu_StartX, &b_CandidateMu_StartX);
   fChain->SetBranchAddress("CandidateMu_StartY", &CandidateMu_StartY, &b_CandidateMu_StartY);
   fChain->SetBranchAddress("CandidateMu_StartZ", &CandidateMu_StartZ, &b_CandidateMu_StartZ);
   fChain->SetBranchAddress("CandidateMu_EndX", &CandidateMu_EndX, &b_CandidateMu_EndX);
   fChain->SetBranchAddress("CandidateMu_EndY", &CandidateMu_EndY, &b_CandidateMu_EndY);
   fChain->SetBranchAddress("CandidateMu_EndZ", &CandidateMu_EndZ, &b_CandidateMu_EndZ);
   fChain->SetBranchAddress("True_CandidateMu_P", &True_CandidateMu_P, &b_True_CandidateMu_P);
   fChain->SetBranchAddress("True_CandidateMu_Phi", &True_CandidateMu_Phi, &b_True_CandidateMu_Phi);
   fChain->SetBranchAddress("True_CandidateMu_Theta", &True_CandidateMu_Theta, &b_True_CandidateMu_Theta);
   fChain->SetBranchAddress("True_CandidateMu_CosTheta", &True_CandidateMu_CosTheta, &b_True_CandidateMu_CosTheta);
   fChain->SetBranchAddress("True_CandidateMu_StartX", &True_CandidateMu_StartX, &b_True_CandidateMu_StartX);
   fChain->SetBranchAddress("True_CandidateMu_StartY", &True_CandidateMu_StartY, &b_True_CandidateMu_StartY);
   fChain->SetBranchAddress("True_CandidateMu_StartZ", &True_CandidateMu_StartZ, &b_True_CandidateMu_StartZ);
   fChain->SetBranchAddress("True_CandidateMu_EndX", &True_CandidateMu_EndX, &b_True_CandidateMu_EndX);
   fChain->SetBranchAddress("True_CandidateMu_EndY", &True_CandidateMu_EndY, &b_True_CandidateMu_EndY);
   fChain->SetBranchAddress("True_CandidateMu_EndZ", &True_CandidateMu_EndZ, &b_True_CandidateMu_EndZ);
   fChain->SetBranchAddress("True_CandidateMu_StartContainment", &True_CandidateMu_StartContainment, &b_True_CandidateMu_StartContainment);
   fChain->SetBranchAddress("True_CandidateMu_EndContainment", &True_CandidateMu_EndContainment, &b_True_CandidateMu_EndContainment);
   fChain->SetBranchAddress("CandidateP_P_Range", &CandidateP_P_Range, &b_CandidateP_P_Range);
   fChain->SetBranchAddress("CandidateP_Phi", &CandidateP_Phi, &b_CandidateP_Phi);
   fChain->SetBranchAddress("CandidateP_Theta", &CandidateP_Theta, &b_CandidateP_Theta);
   fChain->SetBranchAddress("CandidateP_CosTheta", &CandidateP_CosTheta, &b_CandidateP_CosTheta);
   fChain->SetBranchAddress("CandidateP_LLR_PID", &CandidateP_LLR_PID, &b_CandidateP_LLR_PID);
   fChain->SetBranchAddress("CandidateP_StartContainment", &CandidateP_StartContainment, &b_CandidateP_StartContainment);
   fChain->SetBranchAddress("CandidateP_EndContainment", &CandidateP_EndContainment, &b_CandidateP_EndContainment);
   fChain->SetBranchAddress("CandidateP_Length", &CandidateP_Length, &b_CandidateP_Length);
   fChain->SetBranchAddress("CandidateP_MCParticle_Pdg", &CandidateP_MCParticle_Pdg, &b_CandidateP_MCParticle_Pdg);
   fChain->SetBranchAddress("CandidateP_MCParticle_Purity", &CandidateP_MCParticle_Purity, &b_CandidateP_MCParticle_Purity);
   fChain->SetBranchAddress("CandidateP_StartX", &CandidateP_StartX, &b_CandidateP_StartX);
   fChain->SetBranchAddress("CandidateP_StartY", &CandidateP_StartY, &b_CandidateP_StartY);
   fChain->SetBranchAddress("CandidateP_StartZ", &CandidateP_StartZ, &b_CandidateP_StartZ);
   fChain->SetBranchAddress("CandidateP_EndX", &CandidateP_EndX, &b_CandidateP_EndX);
   fChain->SetBranchAddress("CandidateP_EndY", &CandidateP_EndY, &b_CandidateP_EndY);
   fChain->SetBranchAddress("CandidateP_EndZ", &CandidateP_EndZ, &b_CandidateP_EndZ);
   fChain->SetBranchAddress("True_CandidateP_P", &True_CandidateP_P, &b_True_CandidateP_P);
   fChain->SetBranchAddress("True_CandidateP_Phi", &True_CandidateP_Phi, &b_True_CandidateP_Phi);
   fChain->SetBranchAddress("True_CandidateP_Theta", &True_CandidateP_Theta, &b_True_CandidateP_Theta);
   fChain->SetBranchAddress("True_CandidateP_CosTheta", &True_CandidateP_CosTheta, &b_True_CandidateP_CosTheta);
   fChain->SetBranchAddress("True_CandidateP_StartX", &True_CandidateP_StartX, &b_True_CandidateP_StartX);
   fChain->SetBranchAddress("True_CandidateP_StartY", &True_CandidateP_StartY, &b_True_CandidateP_StartY);
   fChain->SetBranchAddress("True_CandidateP_StartZ", &True_CandidateP_StartZ, &b_True_CandidateP_StartZ);
   fChain->SetBranchAddress("True_CandidateP_EndX", &True_CandidateP_EndX, &b_True_CandidateP_EndX);
   fChain->SetBranchAddress("True_CandidateP_EndY", &True_CandidateP_EndY, &b_True_CandidateP_EndY);
   fChain->SetBranchAddress("True_CandidateP_EndZ", &True_CandidateP_EndZ, &b_True_CandidateP_EndZ);
   fChain->SetBranchAddress("True_CandidateP_StartContainment", &True_CandidateP_StartContainment, &b_True_CandidateP_StartContainment);
   fChain->SetBranchAddress("True_CandidateP_EndContainment", &True_CandidateP_EndContainment, &b_True_CandidateP_EndContainment);
   fChain->SetBranchAddress("Reco_A", &Reco_A, &b_Reco_A);
   fChain->SetBranchAddress("Reco_kMiss", &Reco_kMiss, &b_Reco_kMiss);
   fChain->SetBranchAddress("Reco_PMissMinus", &Reco_PMissMinus, &b_Reco_PMissMinus);
   fChain->SetBranchAddress("Reco_EMiss", &Reco_EMiss, &b_Reco_EMiss);
   fChain->SetBranchAddress("Reco_PMiss", &Reco_PMiss, &b_Reco_PMiss);
   fChain->SetBranchAddress("Reco_Pt", &Reco_Pt, &b_Reco_Pt);
   fChain->SetBranchAddress("Reco_Ptx", &Reco_Ptx, &b_Reco_Ptx);
   fChain->SetBranchAddress("Reco_Pty", &Reco_Pty, &b_Reco_Pty);
   fChain->SetBranchAddress("Reco_PL", &Reco_PL, &b_Reco_PL);
   fChain->SetBranchAddress("Reco_Pn", &Reco_Pn, &b_Reco_Pn);
   fChain->SetBranchAddress("Reco_PnPerp", &Reco_PnPerp, &b_Reco_PnPerp);
   fChain->SetBranchAddress("Reco_PnPerpx", &Reco_PnPerpx, &b_Reco_PnPerpx);
   fChain->SetBranchAddress("Reco_PnPerpy", &Reco_PnPerpy, &b_Reco_PnPerpy);      
   fChain->SetBranchAddress("Reco_PnPar", &Reco_PnPar, &b_Reco_PnPar);      
   fChain->SetBranchAddress("Reco_DeltaAlphaT", &Reco_DeltaAlphaT, &b_Reco_DeltaAlphaT);
   fChain->SetBranchAddress("Reco_DeltaAlpha3Dq", &Reco_DeltaAlpha3Dq, &b_Reco_DeltaAlpha3Dq);
   fChain->SetBranchAddress("Reco_DeltaAlpha3DMu", &Reco_DeltaAlpha3DMu, &b_Reco_DeltaAlpha3DMu);      
   fChain->SetBranchAddress("Reco_DeltaPhiT", &Reco_DeltaPhiT, &b_Reco_DeltaPhiT);
   fChain->SetBranchAddress("Reco_DeltaPhi3D", &Reco_DeltaPhi3D, &b_Reco_DeltaPhi3D);   
   fChain->SetBranchAddress("Reco_ECal", &Reco_ECal, &b_Reco_ECal);
   fChain->SetBranchAddress("Reco_EQE", &Reco_EQE, &b_Reco_EQE);
   fChain->SetBranchAddress("Reco_Q2", &Reco_Q2, &b_Reco_Q2);
   fChain->SetBranchAddress("Reco_DeltaPhi", &Reco_DeltaPhi, &b_Reco_DeltaPhi);
   fChain->SetBranchAddress("Reco_DeltaTheta", &Reco_DeltaTheta, &b_Reco_DeltaTheta);
   fChain->SetBranchAddress("True_A", &True_A, &b_True_A);
   fChain->SetBranchAddress("True_kMiss", &True_kMiss, &b_True_kMiss);
   fChain->SetBranchAddress("True_PMissMinus", &True_PMissMinus, &b_True_PMissMinus);
   fChain->SetBranchAddress("True_EMiss", &True_EMiss, &b_True_EMiss);
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
   fChain->SetBranchAddress("StartToStartDistance", &StartToStartDistance, &b_StartToStartDistance);
   fChain->SetBranchAddress("EndToEndDistance", &EndToEndDistance, &b_EndToEndDistance);
   Notify();
}

Bool_t PeLEE_myRecoAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PeLEE_myRecoAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PeLEE_myRecoAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PeLEE_myRecoAnalysis_cxx
