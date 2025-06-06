#ifndef mcc9_10_reco_selection_h
#define mcc9_10_reco_selection_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include "../../../generators/constants.h"

#include <vector>

using namespace constants;

class mcc9_10_reco_selection {

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
   Float_t         wc_single_photon_numu_score;
   Float_t         wc_single_photon_other_score;
   Float_t         wc_single_photon_ncpi0_score;
   Float_t         wc_single_photon_nue_score;
   Float_t         wc_numu_score;
   Float_t         wc_nc_pio_score;
   Float_t         wc_kine_pio_vtx_dis;
   Float_t         wc_kine_pio_energy_1;
   Float_t         wc_kine_pio_theta_1;
   Float_t         wc_kine_pio_phi_1;
   Float_t         wc_kine_pio_energy_2;
   Float_t         wc_kine_pio_theta_2;
   Float_t         wc_kine_pio_phi_2;
   Bool_t          wc_match_isFC;
   Int_t           wc_kine_pio_flag;
   vector<int>     *wc_kine_particle_type;
   vector<float>   *wc_kine_energy_particle;
   Int_t           Run;
   Int_t           SubRun;
   Int_t           Event;
   TString         *run_period;
   Int_t           signal;
   Int_t           nc;
   Int_t           numu;
   Int_t           qe;
   Int_t           mec;
   Int_t           res;
   Int_t           dis;
   Int_t           coh;
   Int_t           other;
	int bkg_0pi0_X;
	int bkg_Mpi0_X;
	int bkg_bwds_1pi0_X;		
	int bkg_1n_0p_1pi0_X;	
	int bkg_Nn_0p_1pi0_X;
	int bkg_1p_0n_1pi0_X;	
	int bkg_Np_0n_1pi0_X;	
	int bkg_1pi0_Npipm_X;	
	int bkg_1pi0_Np_Nn_0pipm_X;
	Int_t bkg_1pi0_Np_Nn_Npipm_X;		
	int bkg_1pi0_Nmh_X; // m = mesons (mostly etas), h = heavy particles (Sigmas, Lambdas)	
	int bkg_1pi0_Nl_X; // l = lepton		
	int bkg_other;	   
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
   Int_t           MCParticle_Mode;
   Double_t        True_Ev;
   Double_t        True_Vx;
   Double_t        True_Vy;
   Double_t        True_Vz;
   Float_t         ns_time;
   Float_t         NuScore;
   Float_t         orig_nuscore;
   Int_t           slice_id;
   Float_t         FlashScore;
   Float_t         CosmicIPAll3D;
   Float_t         CosmicDirAll3D;
   Int_t           crtveto;
   Float_t         crthitpe;
   vector<float>   *Vertex_X;
   vector<float>   *Vertex_Y;
   vector<float>   *Vertex_Z;
   vector<int>     *wc_reco_mother;
   vector<vector<float> > *wc_reco_p;
   vector<vector<float> > *wc_reco_start;
   vector<vector<float> > *wc_reco_end;
   vector<int>     *wc_reco_pdg;
   vector<int>     *wc_reco_id;
   vector<double>  *reco_alpha;
   vector<double>  *reco_shower_opening_angle;
   vector<double>  *reco_pi0_p_gammas;
   vector<double>  *reco_pi0_p;
   vector<double>  *reco_pi0_phi;
   vector<double>  *reco_pi0_costheta;
   vector<double>  *reco_pi0_invmass;   
   vector<double>  *reco_g1_p;
   vector<double>  *reco_g1_phi;
   vector<double>  *reco_g1_costheta;
   vector<double>  *reco_g2_p;
   vector<double>  *reco_g2_phi;
   vector<double>  *reco_g2_costheta;
   Int_t           nBlips_saved;
   vector<float>   *Blip_x;
   vector<float>   *Blip_y;
   vector<float>   *Blip_z;
   vector<float>   *Blip_energy;
   vector<float>   *Blip_charge;
   vector<int>     *Blip_nplanes;
   vector<float>   *Blip_proxtrkdist;
   vector<int>     *Blip_proxtrkid;
   vector<bool>    *Blip_touchtrk;
   vector<int>     *Blip_touchtrkid;
   vector<int>     *Blip_pl0_nwires;
   vector<int>     *Blip_pl1_nwires;
   vector<int>     *Blip_pl2_nwires;
   vector<bool>    *Blip_pl0_bydeadwire;
   vector<bool>    *Blip_pl1_bydeadwire;
   vector<bool>    *Blip_pl2_bydeadwire;
   vector<int>     *Blip_true_g4id;
   vector<float>   *Blip_true_energy;
   vector<unsigned int> *pd_generation_v;
   vector<float>   *pd_trk_score_v;
   vector<float>   *pd_trk_llr_pid_score_v;
   Int_t           pd_reco_track_count;
   Int_t           pd_reco_shower_count;

   // List of branches
    TBranch        *b_Weight;   //!
   TBranch        *b_T2KWeight;   //!
   TBranch        *b_ROOTinoWeight;   //!
   TBranch        *b_POTWeight;   //!
   TBranch        *b_wc_single_photon_numu_score;   //!
   TBranch        *b_wc_single_photon_other_score;   //!
   TBranch        *b_wc_single_photon_ncpi0_score;   //!
   TBranch        *b_wc_single_photon_nue_score;   //!
   TBranch        *b_wc_numu_score;   //!
   TBranch        *b_wc_nc_pio_score;   //!
   TBranch        *b_wc_kine_pio_vtx_dis;   //!
   TBranch        *b_wc_kine_pio_energy_1;   //!
   TBranch        *b_wc_kine_pio_theta_1;   //!
   TBranch        *b_wc_kine_pio_phi_1;   //!
   TBranch        *b_wc_kine_pio_energy_2;   //!
   TBranch        *b_wc_kine_pio_theta_2;   //!
   TBranch        *b_wc_kine_pio_phi_2;   //!
   TBranch        *b_wc_match_isFC;   //!
   TBranch        *b_wc_kine_pio_flag;   //!
   TBranch        *b_wc_kine_particle_type;   //!
   TBranch        *b_wc_kine_energy_particle;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_run_period;   //!
   TBranch        *b_signal;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_numu;   //!
   TBranch        *b_qe;   //!
   TBranch        *b_mec;   //!
   TBranch        *b_res;   //!
   TBranch        *b_dis;   //!
   TBranch        *b_coh;   //!
   TBranch        *b_other;   //!
   TBranch        *b_bkg_0pi0_X;
   TBranch        *b_bkg_Mpi0_X;
   TBranch        *b_bkg_bwds_1pi0_X;		
   TBranch        *b_bkg_1n_0p_1pi0_X;	
   TBranch        *b_bkg_Nn_0p_1pi0_X;
   TBranch        *b_bkg_1p_0n_1pi0_X;	
   TBranch        *b_bkg_Np_0n_1pi0_X;	
   TBranch        *b_bkg_1pi0_Npipm_X;	
   TBranch        *b_bkg_1pi0_Np_Nn_0pipm_X;
   TBranch        *b_bkg_1pi0_Np_Nn_Npipm_X;		
   TBranch        *b_bkg_1pi0_Nmh_X; // m = mesons (mostly etas), h = heavy particles (Sigmas, Lambdas)	
   TBranch        *b_bkg_1pi0_Nl_X; // l = lepton		
   TBranch        *b_bkg_other;	   
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
   TBranch        *b_MCParticle_Mode;   //!
   TBranch        *b_True_Ev;   //!
   TBranch        *b_True_Vx;   //!
   TBranch        *b_True_Vy;   //!
   TBranch        *b_True_Vz;   //!
   TBranch        *b_ns_time;   //!
   TBranch        *b_NuScore;   //!
   TBranch        *b_orig_nuscore;   //!
   TBranch        *b_slice_id;   //!
   TBranch        *b_FlashScore;   //!
   TBranch        *b_CosmicIPAll3D;   //!
   TBranch        *b_CosmicDirAll3D;   //!
   TBranch        *b_crtveto;   //!
   TBranch        *b_crthitpe;   //!
   TBranch        *b_Vertex_X;   //!
   TBranch        *b_Vertex_Y;   //!
   TBranch        *b_Vertex_Z;   //!
   TBranch        *b_wc_reco_mother;   //!
   TBranch        *b_wc_reco_p;   //!
   TBranch        *b_wc_reco_start;   //!
   TBranch        *b_wc_reco_end;   //!
   TBranch        *b_wc_reco_pdg;   //!
   TBranch        *b_wc_reco_id;   //!
   TBranch        *b_reco_alpha;   //!
   TBranch        *b_reco_shower_opening_angle;   //!
   TBranch        *b_reco_pi0_p_gammas;   //!
   TBranch        *b_reco_pi0_p;   //!
   TBranch        *b_reco_pi0_phi;   //!
   TBranch        *b_reco_pi0_costheta;   //!
   TBranch        *b_reco_pi0_invmass;   //!   
   TBranch        *b_reco_g1_p;   //!
   TBranch        *b_reco_g1_phi;   //!
   TBranch        *b_reco_g1_costheta;   //!
   TBranch        *b_reco_g2_p;   //!
   TBranch        *b_reco_g2_phi;   //!
   TBranch        *b_reco_g2_costheta;   //!
   TBranch        *b_nBlips_saved;   //!
   TBranch        *b_Blip_x;   //!
   TBranch        *b_Blip_y;   //!
   TBranch        *b_Blip_z;   //!
   TBranch        *b_Blip_energy;   //!
   TBranch        *b_Blip_charge;   //!
   TBranch        *b_Blip_nplanes;   //!
   TBranch        *b_Blip_proxtrkdist;   //!
   TBranch        *b_Blip_proxtrkid;   //!
   TBranch        *b_Blip_touchtrk;   //!
   TBranch        *b_Blip_touchtrkid;   //!
   TBranch        *b_Blip_pl0_nwires;   //!
   TBranch        *b_Blip_pl1_nwires;   //!
   TBranch        *b_Blip_pl2_nwires;   //!
   TBranch        *b_Blip_pl0_bydeadwire;   //!
   TBranch        *b_Blip_pl1_bydeadwire;   //!
   TBranch        *b_Blip_pl2_bydeadwire;   //!
   TBranch        *b_Blip_true_g4id;   //!
   TBranch        *b_Blip_true_energy;   //!
   TBranch        *b_pd_generation_v;   //!
   TBranch        *b_pd_trk_score_v;   //!
   TBranch        *b_pd_trk_llr_pid_score_v;   //!
   TBranch        *b_pd_reco_track_count;   //!
   TBranch        *b_pd_reco_shower_count;   //!

   mcc9_10_reco_selection(TString WhichSample="",TString Tune="",TString WhichEventWeightLabel="", int UniverseIndex=-1, TTree *tree=0);
   virtual ~mcc9_10_reco_selection();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mcc9_10_reco_selection_cxx
mcc9_10_reco_selection::mcc9_10_reco_selection(TString WhichSample, TString Tune, TString WhichEventWeightLabel, int UniverseIndex, TTree *tree) : fChain(0) 
{

   fTune = Tune;
   fWhichSample = WhichSample;

   fEventWeightLabel = WhichEventWeightLabel;
   fUniverseIndex = UniverseIndex;

   fPathToFile = "/exp/uboone/data/users/"+UserID+"/ncpi0/PreSelection_"+fWhichSample+".root";

   // ------------------------ //

   if (fWhichSample == "Overlay9BNBToHonda_Combined" || fWhichSample == "Overlay9BNBToHondaECal_Combined") {

	fPathToFile = "/exp/uboone/data/users/"+UserID+"/ncpi0/PreSelection_Overlay9_Combined.root";

   }

   if (fWhichSample == "OverlayDirt9BNBToHonda_Combined" || fWhichSample == "OverlayDirt9BNBToHondaECal_Combined") {

	fPathToFile = "/exp/uboone/data/users/"+UserID+"/ncpi0/PreSelection_OverlayDirt9_Combined.root";

   }

   if (fWhichSample == "BeamOn9BNBToHonda_Combined" || fWhichSample == "BeamOn9BNBToHondaECal_Combined") {

	fPathToFile = "/exp/uboone/data/users/"+UserID+"/ncpi0/PreSelection_BeamOn9_Combined.root";

   }


   if (fWhichSample == "ExtBNB9BNBToHonda_Combined" || fWhichSample == "ExtBNB9BNBToHondaECal_Combined") {

	fPathToFile = "/exp/uboone/data/users/"+UserID+"/ncpi0/PreSelection_ExtBNB9_Combined.root";

   }


   // ------------------------ //

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

mcc9_10_reco_selection::~mcc9_10_reco_selection()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mcc9_10_reco_selection::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mcc9_10_reco_selection::LoadTree(Long64_t entry)
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

void mcc9_10_reco_selection::Init(TTree *tree)
{

   // Set object pointer
   wc_kine_particle_type = 0;
   wc_kine_energy_particle = 0;
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
   Vertex_X = 0;
   Vertex_Y = 0;
   Vertex_Z = 0;
   wc_reco_mother = 0;
   wc_reco_p = 0;
   wc_reco_start = 0;
   wc_reco_end = 0;
   wc_reco_pdg = 0;
   wc_reco_id = 0;
   reco_alpha = 0;
   reco_shower_opening_angle = 0;
   reco_pi0_p_gammas = 0;
   reco_pi0_p = 0;
   reco_pi0_phi = 0;
   reco_pi0_costheta = 0;
   reco_pi0_invmass = 0;   
   reco_g1_p = 0;
   reco_g1_phi = 0;
   reco_g1_costheta = 0;
   reco_g2_p = 0;
   reco_g2_phi = 0;
   reco_g2_costheta = 0;
   Blip_x = 0;
   Blip_y = 0;
   Blip_z = 0;
   Blip_energy = 0;
   Blip_charge = 0;
   Blip_nplanes = 0;
   Blip_proxtrkdist = 0;
   Blip_proxtrkid = 0;
   Blip_touchtrk = 0;
   Blip_touchtrkid = 0;
   Blip_pl0_nwires = 0;
   Blip_pl1_nwires = 0;
   Blip_pl2_nwires = 0;
   Blip_pl0_bydeadwire = 0;
   Blip_pl1_bydeadwire = 0;
   Blip_pl2_bydeadwire = 0;
   Blip_true_g4id = 0;
   Blip_true_energy = 0;
   pd_generation_v = 0;
   pd_trk_score_v = 0;
   pd_trk_llr_pid_score_v = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("T2KWeight", &T2KWeight, &b_T2KWeight);
   fChain->SetBranchAddress("ROOTinoWeight", &ROOTinoWeight, &b_ROOTinoWeight);
   fChain->SetBranchAddress("POTWeight", &POTWeight, &b_POTWeight);
   fChain->SetBranchAddress("wc_single_photon_numu_score", &wc_single_photon_numu_score, &b_wc_single_photon_numu_score);
   fChain->SetBranchAddress("wc_single_photon_other_score", &wc_single_photon_other_score, &b_wc_single_photon_other_score);
   fChain->SetBranchAddress("wc_single_photon_ncpi0_score", &wc_single_photon_ncpi0_score, &b_wc_single_photon_ncpi0_score);
   fChain->SetBranchAddress("wc_single_photon_nue_score", &wc_single_photon_nue_score, &b_wc_single_photon_nue_score);
   fChain->SetBranchAddress("wc_numu_score", &wc_numu_score, &b_wc_numu_score);
   fChain->SetBranchAddress("wc_nc_pio_score", &wc_nc_pio_score, &b_wc_nc_pio_score);
   fChain->SetBranchAddress("wc_kine_pio_vtx_dis", &wc_kine_pio_vtx_dis, &b_wc_kine_pio_vtx_dis);
   fChain->SetBranchAddress("wc_kine_pio_energy_1", &wc_kine_pio_energy_1, &b_wc_kine_pio_energy_1);
   fChain->SetBranchAddress("wc_kine_pio_theta_1", &wc_kine_pio_theta_1, &b_wc_kine_pio_theta_1);
   fChain->SetBranchAddress("wc_kine_pio_phi_1", &wc_kine_pio_phi_1, &b_wc_kine_pio_phi_1);
   fChain->SetBranchAddress("wc_kine_pio_energy_2", &wc_kine_pio_energy_2, &b_wc_kine_pio_energy_2);
   fChain->SetBranchAddress("wc_kine_pio_theta_2", &wc_kine_pio_theta_2, &b_wc_kine_pio_theta_2);
   fChain->SetBranchAddress("wc_kine_pio_phi_2", &wc_kine_pio_phi_2, &b_wc_kine_pio_phi_2);
   fChain->SetBranchAddress("wc_match_isFC", &wc_match_isFC, &b_wc_match_isFC);
   fChain->SetBranchAddress("wc_kine_pio_flag", &wc_kine_pio_flag, &b_wc_kine_pio_flag);
   fChain->SetBranchAddress("wc_kine_particle_type", &wc_kine_particle_type, &b_wc_kine_particle_type);
   fChain->SetBranchAddress("wc_kine_energy_particle", &wc_kine_energy_particle, &b_wc_kine_energy_particle);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("run_period", &run_period, &b_run_period);
   fChain->SetBranchAddress("signal", &signal, &b_signal);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("numu", &numu, &b_numu);
   fChain->SetBranchAddress("qe", &qe, &b_qe);
   fChain->SetBranchAddress("mec", &mec, &b_mec);
   fChain->SetBranchAddress("res", &res, &b_res);
   fChain->SetBranchAddress("dis", &dis, &b_dis);
   fChain->SetBranchAddress("coh", &coh, &b_coh);
   fChain->SetBranchAddress("other", &other, &b_other);  
	fChain->SetBranchAddress("bkg_0pi0_X",&bkg_0pi0_X,&b_bkg_0pi0_X);
	fChain->SetBranchAddress("bkg_Mpi0_X",&bkg_Mpi0_X,&b_bkg_Mpi0_X);
	fChain->SetBranchAddress("bkg_bwds_1pi0_X",&bkg_bwds_1pi0_X,&b_bkg_bwds_1pi0_X);
	fChain->SetBranchAddress("bkg_1n_0p_1pi0_X",&bkg_1n_0p_1pi0_X,&b_bkg_1n_0p_1pi0_X);			
	fChain->SetBranchAddress("bkg_Nn_0p_1pi0_X",&bkg_Nn_0p_1pi0_X,&b_bkg_Nn_0p_1pi0_X);
	fChain->SetBranchAddress("bkg_1p_0n_1pi0_X",&bkg_1p_0n_1pi0_X,&b_bkg_1p_0n_1pi0_X);
	fChain->SetBranchAddress("bkg_Np_0n_1pi0_X",&bkg_Np_0n_1pi0_X,&b_bkg_Np_0n_1pi0_X);
	fChain->SetBranchAddress("bkg_1pi0_Npipm_X",&bkg_1pi0_Npipm_X,&b_bkg_1pi0_Npipm_X);
	fChain->SetBranchAddress("bkg_1pi0_Np_Nn_0pipm_X",&bkg_1pi0_Np_Nn_0pipm_X,&b_bkg_1pi0_Np_Nn_0pipm_X);
	fChain->SetBranchAddress("bkg_1pi0_Np_Nn_Npipm_X",&bkg_1pi0_Np_Nn_0pipm_X,&b_bkg_1pi0_Np_Nn_0pipm_X);	
	fChain->SetBranchAddress("bkg_1pi0_Nmh_X",&bkg_1pi0_Np_Nn_0pipm_X,&b_bkg_1pi0_Np_Nn_0pipm_X);
	fChain->SetBranchAddress("bkg_1pi0_Nl_X",&bkg_1pi0_Nl_X,&b_bkg_1pi0_Nl_X);
	fChain->SetBranchAddress("bkg_other",&bkg_other,&b_bkg_other);	   
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
   fChain->SetBranchAddress("MCParticle_Mode", &MCParticle_Mode, &b_MCParticle_Mode);
   fChain->SetBranchAddress("True_Ev", &True_Ev, &b_True_Ev);
   fChain->SetBranchAddress("True_Vx", &True_Vx, &b_True_Vx);
   fChain->SetBranchAddress("True_Vy", &True_Vy, &b_True_Vy);
   fChain->SetBranchAddress("True_Vz", &True_Vz, &b_True_Vz);
   fChain->SetBranchAddress("ns_time", &ns_time, &b_ns_time);
   fChain->SetBranchAddress("NuScore", &NuScore, &b_NuScore);
   fChain->SetBranchAddress("orig_nuscore", &orig_nuscore, &b_orig_nuscore);
   fChain->SetBranchAddress("slice_id", &slice_id, &b_slice_id);
   fChain->SetBranchAddress("FlashScore", &FlashScore, &b_FlashScore);
   fChain->SetBranchAddress("CosmicIPAll3D", &CosmicIPAll3D, &b_CosmicIPAll3D);
   fChain->SetBranchAddress("CosmicDirAll3D", &CosmicDirAll3D, &b_CosmicDirAll3D);
   fChain->SetBranchAddress("crtveto", &crtveto, &b_crtveto);
   fChain->SetBranchAddress("crthitpe", &crthitpe, &b_crthitpe);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("wc_reco_mother", &wc_reco_mother, &b_wc_reco_mother);
   fChain->SetBranchAddress("wc_reco_p", &wc_reco_p, &b_wc_reco_p);
   fChain->SetBranchAddress("wc_reco_start", &wc_reco_start, &b_wc_reco_start);
   fChain->SetBranchAddress("wc_reco_end", &wc_reco_end, &b_wc_reco_end);
   fChain->SetBranchAddress("wc_reco_pdg", &wc_reco_pdg, &b_wc_reco_pdg);
   fChain->SetBranchAddress("wc_reco_id", &wc_reco_id, &b_wc_reco_id);
   fChain->SetBranchAddress("reco_alpha", &reco_alpha, &b_reco_alpha);
   fChain->SetBranchAddress("reco_shower_opening_angle", &reco_shower_opening_angle, &b_reco_shower_opening_angle);
   fChain->SetBranchAddress("reco_pi0_p_gammas", &reco_pi0_p_gammas, &b_reco_pi0_p_gammas);
   fChain->SetBranchAddress("reco_pi0_p", &reco_pi0_p, &b_reco_pi0_p);
   fChain->SetBranchAddress("reco_pi0_phi", &reco_pi0_phi, &b_reco_pi0_phi);
   fChain->SetBranchAddress("reco_pi0_costheta", &reco_pi0_costheta, &b_reco_pi0_costheta);
   fChain->SetBranchAddress("reco_pi0_invmass", &reco_pi0_invmass, &b_reco_pi0_invmass);   
   fChain->SetBranchAddress("reco_g1_p", &reco_g1_p, &b_reco_g1_p);
   fChain->SetBranchAddress("reco_g1_phi", &reco_g1_phi, &b_reco_g1_phi);
   fChain->SetBranchAddress("reco_g1_costheta", &reco_g1_costheta, &b_reco_g1_costheta);
   fChain->SetBranchAddress("reco_g2_p", &reco_g2_p, &b_reco_g2_p);
   fChain->SetBranchAddress("reco_g2_phi", &reco_g2_phi, &b_reco_g2_phi);
   fChain->SetBranchAddress("reco_g2_costheta", &reco_g2_costheta, &b_reco_g2_costheta);
   fChain->SetBranchAddress("nBlips_saved", &nBlips_saved, &b_nBlips_saved);
   fChain->SetBranchAddress("Blip_x", &Blip_x, &b_Blip_x);
   fChain->SetBranchAddress("Blip_y", &Blip_y, &b_Blip_y);
   fChain->SetBranchAddress("Blip_z", &Blip_z, &b_Blip_z);
   fChain->SetBranchAddress("Blip_energy", &Blip_energy, &b_Blip_energy);
   fChain->SetBranchAddress("Blip_charge", &Blip_charge, &b_Blip_charge);
   fChain->SetBranchAddress("Blip_nplanes", &Blip_nplanes, &b_Blip_nplanes);
   fChain->SetBranchAddress("Blip_proxtrkdist", &Blip_proxtrkdist, &b_Blip_proxtrkdist);
   fChain->SetBranchAddress("Blip_proxtrkid", &Blip_proxtrkid, &b_Blip_proxtrkid);
   fChain->SetBranchAddress("Blip_touchtrk", &Blip_touchtrk, &b_Blip_touchtrk);
   fChain->SetBranchAddress("Blip_touchtrkid", &Blip_touchtrkid, &b_Blip_touchtrkid);
   fChain->SetBranchAddress("Blip_pl0_nwires", &Blip_pl0_nwires, &b_Blip_pl0_nwires);
   fChain->SetBranchAddress("Blip_pl1_nwires", &Blip_pl1_nwires, &b_Blip_pl1_nwires);
   fChain->SetBranchAddress("Blip_pl2_nwires", &Blip_pl2_nwires, &b_Blip_pl2_nwires);
   fChain->SetBranchAddress("Blip_pl0_bydeadwire", &Blip_pl0_bydeadwire, &b_Blip_pl0_bydeadwire);
   fChain->SetBranchAddress("Blip_pl1_bydeadwire", &Blip_pl1_bydeadwire, &b_Blip_pl1_bydeadwire);
   fChain->SetBranchAddress("Blip_pl2_bydeadwire", &Blip_pl2_bydeadwire, &b_Blip_pl2_bydeadwire);
   fChain->SetBranchAddress("Blip_true_g4id", &Blip_true_g4id, &b_Blip_true_g4id);
   fChain->SetBranchAddress("Blip_true_energy", &Blip_true_energy, &b_Blip_true_energy);
   fChain->SetBranchAddress("pd_generation_v", &pd_generation_v, &b_pd_generation_v);
   fChain->SetBranchAddress("pd_trk_score_v", &pd_trk_score_v, &b_pd_trk_score_v);
   fChain->SetBranchAddress("pd_trk_llr_pid_score_v", &pd_trk_llr_pid_score_v, &b_pd_trk_llr_pid_score_v);
   fChain->SetBranchAddress("pd_reco_track_count", &pd_reco_track_count, &b_pd_reco_track_count);
   fChain->SetBranchAddress("pd_reco_shower_count", &pd_reco_shower_count, &b_pd_reco_shower_count);

   Notify();
}

Bool_t mcc9_10_reco_selection::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mcc9_10_reco_selection::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mcc9_10_reco_selection::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mcc9_10_reco_selection_cxx