////////////////////////////////////////////////////////////////////////
//Class:RecoBenchmarker
//PluginType:analyzer(artv2_05_00)
//File:RecoBenchmarker_module.cc
//
//GeneratedatTueSep1917:11:422017byAdamListerusingcetskelgen
//fromcetlibversionv1_21_00.
////////////////////////////////////////////////////////////////////////

//baseincludes
#include"art/Framework/Core/EDAnalyzer.h"
#include"art/Framework/Core/ModuleMacros.h"
#include"art/Framework/Principal/Event.h"
#include"art/Framework/Principal/Handle.h"
#include"art/Framework/Principal/Run.h"
#include"art/Framework/Principal/SubRun.h"
#include"art/Framework/Services/Optional/TFileService.h"
#include"canvas/Utilities/InputTag.h"
#include"fhiclcpp/ParameterSet.h"
#include"messagefacility/MessageLogger/MessageLogger.h"
#include"canvas/Persistency/Common/FindManyP.h"

//LArSoftincludes
#include"lardataobj/RecoBase/Track.h"
#include"lardataobj/RecoBase/Shower.h"
#include"lardataobj/RecoBase/Cluster.h"
#include"lardataobj/RecoBase/Hit.h"
#include"larreco/RecoAlg/TrackMomentumCalculator.h"
#include"nusimdata/SimulationBase/MCParticle.h"
#include"nusimdata/SimulationBase/MCTruth.h"
#include"lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include"lardataobj/AnalysisBase/Calorimetry.h"
#include"larevt/SpaceChargeServices/SpaceChargeService.h"
#include"lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include"larcoreobj/SimpleTypesAndConstants/geo_types.h"

//ROOTincludes
#include"TFile.h"
#include"TTree.h"
#include"TH2D.h"
#include"TH1D.h"

//localincludes
#include"uboone/RecoBenchmarker/Algos/recoBenchmarkerUtility.h"

#defineisDebug0

#defineisSingleParticle0

namespacerecohelper{
classRecoBenchmarker;
}


classrecohelper::RecoBenchmarker:publicart::EDAnalyzer{

	public:

		explicitRecoBenchmarker(fhicl::ParameterSetconst&p);
		//Thecompiler-generateddestructorisfinefornon-base
		//classeswithoutbarepointersorotherresourceuse.

		//Pluginsshouldnotbecopiedorassigned.
		RecoBenchmarker(RecoBenchmarkerconst&)=delete;
		RecoBenchmarker(RecoBenchmarker&&)=delete;
		RecoBenchmarker&operator=(RecoBenchmarkerconst&)=delete;
		RecoBenchmarker&operator=(RecoBenchmarker&&)=delete;

		//Requiredfunctions.
		voidanalyze(art::Eventconst&e)override;

		//Selectedoptionalfunctions.
		voidbeginJob()override;
		voidendJob()override;

	private:

		//fclinputparameters
		std::stringfTrackLabel;
		std::stringfCalorimetryLabel;
		std::stringfTrackTruthLabel;
		std::stringfShowerLabel;
		std::stringfShowerTruthLabel;
		std::stringfClusterLabel;
		std::stringfHitLabel;
		std::stringfMCTruthLabel;
		std::stringfG4TruthLabel;

		rbutil::recoBenchmarkerUtility_rbutilInstance;

		art::ServiceHandle<art::TFileService>tfs;
		TTree*recoTree;

		//afewparametersforhistoplotting
		longlow_edge=350;
		longhigh_edge=1500;
		longlength_cut=20;//8cm

		//auxiliary

		intfEvent;
		intfRun;
		intfSubRun;
		intfccnc;
		intfinteraction;

		intnimID;
		intrecoNimID;

		//infoontheMCparticle
		std::map<Int_t,Int_t>fparticle_count;
		std::vector<float>flength;
		std::vector<float>fstartT;
		std::vector<float>fstart_x;
		std::vector<float>fstart_y;
		std::vector<float>fstart_z;
		std::vector<float>fend_x;
		std::vector<float>fend_y;
		std::vector<float>fend_z;
		std::vector<int>fn_steps;
		std::vector<int>fpdg;
		std::vector<int>fmother_id;
		std::vector<int>fg4_id;
		std::vector<float>fp0;//initialmomentum
		std::vector<float>fp0x;
		std::vector<float>fp0y;
		std::vector<float>fp0z;
		std::vector<float>fkinE;
		std::vector<float>fcostheta_muon;
		std::vector<bool>fis_leading;

		std::vector<double>fmuon_dqdx;
		std::vector<double>fmuon_residual;
		std::vector<double>fmuon_dedx;
		doublefmuon_range;

		//infocomingfromthetrackingalgorithm-whenthereismctruth
		std::vector<bool>fis_tracked;
		std::vector<bool>fmatch_multiplicity;
		//std::vector<bool>fis_mismatched;//itsaysiftheMCtruthassignmentindifferentplanesisdifferent(possiblehintforwrongorproblematicbacktracking)
		std::vector<float>fcostheta_muon_reco;//MCtruthinfo
		std::vector<int>fpdg_reco;
		std::vector<float>flength_reco;
		std::vector<float>freco_momentum_mcs;//(GeV)MCS
		std::vector<float>freco_momentum_mcs_llhd;//(GeV)MCSLLHD
		std::vector<float>freco_momentum_range;//MeV
		std::vector<float>fpurity;
		std::vector<float>fcompleteness;
		std::vector<double>fnhits;
		std::vector<float>freco_kinE;
		std::vector<float>freco_startx;
		std::vector<float>freco_starty;
		std::vector<float>freco_startz;
		std::vector<float>freco_endx;
		std::vector<float>freco_endy;
		std::vector<float>freco_endz;

		//infocomingfromthetrackingalgorithm-whenthereisNOmctruth
		std::vector<bool>ffake_is_tracked;
		std::vector<bool>ffake_is_mismatched;
		std::vector<float>ffake_costheta_muon_reco;//MCtruthinfo
		std::vector<int>ffake_pdg_reco;
		std::vector<float>ffake_length_reco;
		std::vector<float>ffake_reco_momentum_mcs;//(GeV)MCS
		std::vector<float>ffake_reco_momentum_mcs_llhd;//(GeV)MCSLLHD
		std::vector<float>ffake_reco_momentum_range;//MeV
		std::vector<float>ffake_purity;
		std::vector<float>ffake_completeness;
		std::vector<double>ffake_nhits;
		std::vector<float>ffake_kinE;

		voidAllocateAnalysisHistograms();
		voidFillAnalysisHistograms();
		voidFillCumulativeHistograms();

		//hitliststuff
		std::vector<std::pair<int,float>>allHitPositions;
		TH1D*hitMatchScore;

		std::vector<int>matchIDChecker;

		//trueVertexXZPosition
		std::vector<float>trueVertexXZPosition;

		//matchinginformation
		std::vector<int>isRecoTrackTruthMatched;
		std::vector<int>isRecoShowerTruthMatched;

		//lengthinformation
		std::vector<double>thisRecoLength;
		std::vector<double>thisMcpLength;
		floatthisMatchedLength;

		//angularinformation
		std::vector<double>thisMcpMomentum;
		std::vector<double>nextMcpMomentum;
		std::vector<double>nimMcpMomentum;
		std::vector<double>thisMatchedMomentum;
		std::vector<double>nextRecoMomentum;
		std::vector<double>nimMatchedMomentum;

		//Nim==NeutrinoInducedMuon
		std::vector<float>thisNimMcpAngles;
		floatthisNimMcpAngle;
		std::vector<float>thisNimMcpAnglesXZ;
		floatthisNimMcpAngleXZ;
		std::vector<float>thisNimMatchedMcpAngles;
		floatthisNimMatchedMcpAngle;
		std::vector<float>thisNimMatchedMcpAnglesXZ;
		floatthisNimMatchedMcpAngleXZ;
		floatthisZmatchedAngleYZ;
		floatthisZDirMcpAngleYZ;

		//efficiencyplotsandprerequisites
		TH1D*muMatchedMcpMomentum;
		TH1D*muMcpMomentum;
		TH1D*muMomentumEfficiency;

		TH2D*allMatchedMcpLengthAngle;
		TH2D*allMcpLengthAngleYZ;
		TH2D*allLengthAngleEfficiency;

		TH1D*allMatchedMcpProjectedLength;
		TH1D*allMcpProjectedLength;
		TH1D*allProjectedLengthEfficiency;

		TH1D*allMatchedMcpProjectedAngle;
		TH1D*allMcpProjectedAngle;
		TH1D*allProjectedAngleEfficiency;

		TH2D*mupMatchedMcpAnglePMom;
		TH2D*mupMcpAnglePMom;
		TH2D*mupAngleMomentumEfficiency;

		TH1D*pMatchedMcpProjectedMomentum;
		TH1D*pMcpProjectedMomentum;
		TH1D*pProjectedMomentumEfficiency;

		TH1D*pMatchedMcpProjectedAngle;
		TH1D*pMcpProjectedAngle;
		TH1D*pProjectedAngleEfficiency;

		//cleanlinessplots
		TH1D*showerCleanlinessPrimaryProton;
		TH1D*showerCleanlinessPrimaryMuonOrPion;
		TH1D*showerCleanlinessPrimaryElectron;
		TH1D*showerCleanlinessConv;
		TH1D*showerCleanlinessInelastic;
		TH1D*showerCleanlinessMuIoni;
		TH1D*showerCleanlinessOther;

		TH1D*trackCleanlinessPrimaryProton;
		TH1D*trackCleanlinessPrimaryMuonOrPion;
		TH1D*trackCleanlinessPrimaryElectron;
		TH1D*trackCleanlinessConv;
		TH1D*trackCleanlinessInelastic;
		TH1D*trackCleanlinessMuIoni;
		TH1D*trackCleanlinessOther;

		TH1D*hmuon_pos_res;
		TH1D*hmuon_pos_res_goodprotons;
		TH1D*hmuon_pos_res_badprotons;
		TH1D*hmuon_pos_res_lowprotons;
		TH1D*hproton_pos_res;
		TH1D*hproton_pos_res_goodprotons;
		TH1D*hproton_pos_res_badprotons;
		TH1D*hproton_pos_res_lowprotons;
		TH1D*hmuon_proton_tracked;
		TH1D*hmuon_spectrum;
		TH1D*hmuon_spectrum_all;
		TH1D*hmuon_length;
		TH1D*hmuon_length_all;
		TH1D*hproton_kinE;
		TH1D*hproton_kinE_all;
		TH1D*hproton_p;
		TH1D*hproton_p_all;
		TH1D*hproton_l;
		TH1D*hproton_l_all;
		TH1D*h_pmu_end_not_tracked;
		TH1D*h_pmu_end_tracked;
		TH1D*h_theta_mu_tracked;
		TH1D*h_theta_mu_not_tracked;
		TH1D*h_theta_mu;
		TH2D*h_theta_mu_length;
		TH2D*h_theta_mu_length_all;
		TH2D*h_dqdx_merged;
		TH2D*h_dqdx_not_merged;
		TH2D*h_dqdx_low_protons;
		TH1D*h_dqdx_1d_merged;
		TH1D*h_dqdx_1d_not_merged;
		TH2D*h_dqdx_tailtotot_length_merged;
		TH2D*h_dqdx_tailtotot_length_not_merged;
		TH1D*htail_to_tot_low_protons;
		TH1D*htail_to_tot_merged;
		TH1D*htail_to_tot_not_merged;
		TH2D*h_dqdx_merged_service;
		TH2D*h_dqdx_not_merged_service;
		TH2D*h_dqdx_low_protons_service;

		voidAllocateRecoVectors();
		voidclear_vectors();
};
//____________________________________________________________________________________________________________________________________________________________________________________________________

recohelper::RecoBenchmarker::RecoBenchmarker(fhicl::ParameterSetconst&p)
:
EDAnalyzer(p)//,
//Moreinitializershere.
{
	fTrackLabel=p.get<std::string>("TrackLabel");
	fCalorimetryLabel=p.get<std::string>("CalorimetryLabel");
	fTrackTruthLabel=p.get<std::string>("TrackTruthLabel");
	fShowerLabel=p.get<std::string>("ShowerLabel");
	fMCTruthLabel=p.get<std::string>("MCTruthLabel");
	fG4TruthLabel=p.get<std::string>("G4TruthLabel");
	fShowerTruthLabel=p.get<std::string>("ShowerTruthLabel");
	fClusterLabel=p.get<std::string>("ClusterLabel");
	fHitLabel=p.get<std::string>("HitLabel");
}

voidrecohelper::RecoBenchmarker::beginJob()
{
	recoTree=tfs->make<TTree>("recotree","recotree");

	//definebranches
	recoTree->Branch("Event",&fEvent,"Event/I");
	recoTree->Branch("SubRun",&fSubRun,"SubRun/I");
	recoTree->Branch("Run",&fRun,"Run/I");

	//mctruthinfo
	recoTree->Branch("ccnc",&fccnc);
	recoTree->Branch("interaction",&finteraction);

	//infoontheMCparticle
	recoTree->Branch("length",&flength);
	recoTree->Branch("startT",&fstartT);
	recoTree->Branch("start_x",&fstart_x);
	recoTree->Branch("start_y",&fstart_y);
	recoTree->Branch("start_z",&fstart_z);
	recoTree->Branch("end_x",&fend_x);
	recoTree->Branch("end_y",&fend_y);
	recoTree->Branch("end_z",&fend_z);
	recoTree->Branch("steps",&fn_steps);
	recoTree->Branch("pdg",&fpdg);
	recoTree->Branch("mother_id",&fmother_id);
	recoTree->Branch("g4_id",&fg4_id);//firstindexist
	recoTree->Branch("p0",&fp0);//initialmomentum
	recoTree->Branch("p0x",&fp0x);
	recoTree->Branch("p0y",&fp0y);
	recoTree->Branch("p0z",&fp0z);
	recoTree->Branch("kinE",&fkinE);
	recoTree->Branch("costheta_muon",&fcostheta_muon);
	recoTree->Branch("is_leading",&fis_leading);
	recoTree->Branch("particle_count",&fparticle_count);
	recoTree->Branch("muon_dedx",&fmuon_dedx);
	recoTree->Branch("muon_dqdx",&fmuon_dqdx);
	recoTree->Branch("muon_residual",&fmuon_residual);
	recoTree->Branch("muon_range",&fmuon_range);

	recoTree->Branch("isRecoTrackTruthMatched",&isRecoTrackTruthMatched);
	recoTree->Branch("isRecoShowerTruthMatched",&isRecoShowerTruthMatched);

	//mcpinformation
	recoTree->Branch("thisNimMcpAngles",&thisNimMcpAngles);
	recoTree->Branch("thisNimMcpAnglesXZ",&thisNimMcpAnglesXZ);
	recoTree->Branch("thisMcpLength",&thisMcpLength);

	//matchedmcpinformation
	recoTree->Branch("thisNimMatchedMcpAngles",&thisNimMatchedMcpAngles);
	recoTree->Branch("thisNimMatchedMcpAnglesXZ",&thisNimMatchedMcpAnglesXZ);

	//infocomingfromthetrackingalgorithm-whenthereismctruth
	recoTree->Branch("is_tracked",&fis_tracked);
	recoTree->Branch("match_multiplicity",&fmatch_multiplicity);
	recoTree->Branch("length_reco",&flength_reco);
	recoTree->Branch("reco_momentum_mcs",&freco_momentum_mcs);//(GeV)MCS
	recoTree->Branch("reco_momentum_mcs_llhd",&freco_momentum_mcs_llhd);//(GeV)MCSLLHD
	recoTree->Branch("reco_momentum_range",&freco_momentum_range);//MeV
	recoTree->Branch("purity",&fpurity);
	recoTree->Branch("completeness",&fcompleteness);
	recoTree->Branch("nhits",&fnhits);
	recoTree->Branch("reco_kinE",&freco_kinE);
	recoTree->Branch("reco_startx",&freco_startx);
	recoTree->Branch("reco_starty",&freco_starty);
	recoTree->Branch("reco_startz",&freco_startz);
	recoTree->Branch("reco_endx",&freco_endx);
	recoTree->Branch("reco_endy",&freco_endy);
	recoTree->Branch("reco_endz",&freco_endz);

	AllocateAnalysisHistograms();

	//reconstructedinformation
	recoTree->Branch("thisRecoLength",&thisRecoLength);

	hitMatchScore=tfs->make<TH1D>("hitMatchScore",";hitMatchScoreScore;",2,0,2);

	muMatchedMcpMomentum=tfs->make<TH1D>("muMatchedMcpMomentum",";#muMomentum;",20,0,3.5);
	muMcpMomentum=tfs->make<TH1D>("muMcpMomentum",";#muMomentum;",20,0,3.5);

	allMatchedMcpLengthAngle=tfs->make<TH2D>("allMatchedMcpLengthAngle",";#theta_{YZ}(degrees);Length(cm)",20,0,180,10,0,10);
	allMcpLengthAngleYZ=tfs->make<TH2D>("allMcpLengthAngleYZ",";#theta_{YZ}(degrees);Length(cm)",20,0,180,10,0,10);

	mupMatchedMcpAnglePMom=tfs->make<TH2D>("mupMatchedMcpAnglePMom",";#theta^{#mup}_{XZ}(degrees);P_{p}(Gev)",20,0,180,20,0,1);
	mupMcpAnglePMom=tfs->make<TH2D>("mupMcpAnglePMom",";#theta^{#mup}_{XZ};P_{p}",20,0,180,20,0,1);

	//cleanlinessplots
	showerCleanlinessPrimaryProton=tfs->make<TH1D>("showerCleanlinessPrimaryProton",";cleanliness;",20,0,1);
	showerCleanlinessPrimaryMuonOrPion=tfs->make<TH1D>("showerCleanlinessPrimaryMuonOrPion",";cleanliness;",20,0,1);
	showerCleanlinessPrimaryElectron=tfs->make<TH1D>("showerCleanlinessPrimaryElectron",";cleanliness;",20,0,1);
	showerCleanlinessConv=tfs->make<TH1D>("showerCleanlinessConv",";cleanliness;",20,0,1);
	showerCleanlinessInelastic=tfs->make<TH1D>("showerCleanlinessInelastic",";cleanliness;",20,0,1);
	showerCleanlinessMuIoni=tfs->make<TH1D>("showerCleanlinessMuIoni",";cleanliness;",20,0,1);
	showerCleanlinessOther=tfs->make<TH1D>("showerCleanlinessOther",";cleanliness;",20,0,1);

	trackCleanlinessPrimaryProton=tfs->make<TH1D>("trackCleanlinessPrimaryProton",";cleanliness;",20,0,1);
	trackCleanlinessPrimaryMuonOrPion=tfs->make<TH1D>("trackCleanlinessPrimaryMuonOrPion",";cleanliness;",20,0,1);
	trackCleanlinessPrimaryElectron=tfs->make<TH1D>("trackCleanlinessPrimaryElectron",";cleanliness;",20,0,1);
	trackCleanlinessConv=tfs->make<TH1D>("trackCleanlinessConv",";cleanliness;",20,0,1);
	trackCleanlinessInelastic=tfs->make<TH1D>("trackCleanlinessInelastic",";cleanliness;",20,0,1);
	trackCleanlinessMuIoni=tfs->make<TH1D>("trackCleanlinessMuIoni",";cleanliness;",20,0,1);
	trackCleanlinessOther=tfs->make<TH1D>("trackCleanlinessOther",";cleanliness;",20,0,1);
}
//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::analyze(art::Eventconst&e)
{

	if(e.isRealData()==true)return;

	//autoconst*theDetector=lar::providerFrom<detinfo::DetectorPropertiesService>();
	autoconst*SCE=lar::providerFrom<spacecharge::SpaceChargeService>();
	boolspace_charge=true;

	fEvent=e.id().event();
	fRun=e.run();
	fSubRun=e.subRun();

	//gethandlestoobjectsofinterest

	art::ValidHandle<std::vector<recob::Track>>trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
	std::vector<art::Ptr<recob::Track>>trackList;
	//if(e.getByLabel(fTrackLabel,trackHandle))
	art::fill_ptr_vector(trackList,trackHandle);

	art::FindManyP<simb::MCParticle,anab::BackTrackerMatchingData>mcpsFromTracks(trackHandle,e,fTrackTruthLabel);

	art::FindManyP<anab::Calorimetry>caloFromTracks(trackHandle,e,fCalorimetryLabel);

	art::ValidHandle<std::vector<recob::Shower>>showerHandle = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);

	art::FindManyP<simb::MCParticle,anab::BackTrackerMatchingData>mcpsFromShowers(showerHandle,e,fShowerTruthLabel);

	art::Handle<std::vector<simb::MCParticle>>mcParticleHandle;
	std::vector<art::Ptr<simb::MCParticle>>mcList;
	if(e.getByLabel(fG4TruthLabel,mcParticleHandle))
	art::fill_ptr_vector(mcList,mcParticleHandle);

	art::Handle<std::vector<simb::MCTruth>>mcTruthHandle;
	std::vector<art::Ptr<simb::MCTruth>>mcTruth;
	if(e.getByLabel(fMCTruthLabel,mcTruthHandle))
	art::fill_ptr_vector(mcTruth,mcTruthHandle);

	art::ValidHandle<std::vector<recob::Cluster>>clusterHandle=
	e.getValidHandle<std::vector<recob::Cluster>>(fClusterLabel);

	art::FindManyP<recob::Hit>hitsFromClusters(clusterHandle,e,fClusterLabel);

	//art::ValidHandle<std::vector<recob::Hit>>hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);

	//---------------------------------
	//MCParticles
	//---------------------------------

	//getMCparticlesassociatedtoMCtruth
	art::FindManyP<simb::MCParticle>MCpFromMCtruth(mcTruth,e,fG4TruthLabel);

	//looponallMCtruthframes(mostly1perevent)
	for(unsignedn_truth = 0 ; n_truth < mcTruth.size() ; n_truth++) {

		clear_vectors();

		#ifisSingleParticle==0
		//checkthattheneutrinoisintheactivevolume
		constsimb::MCNeutrinothisNeutrino=mcTruth[n_truth]->GetNeutrino();
		constsimb::MCParticlethisLepton=thisNeutrino.Lepton();
		if(!_rbutilInstance.isInTPC(thisLepton))continue;

		//storetheinteractioninfo
		fccnc=thisNeutrino.CCNC();
		finteraction=thisNeutrino.InteractionType();

		#endif
		//firstloopmuonstofindtrueneutrinoinducedmuon(NIM)
		nimID=-1;
		doublemuon_p=0;
		doublemuon_px=0;
		doublemuon_py=0;
		doublemuon_pz=0;
		doublemuon_endx=0;
		doublemuon_endy=0;
		doublemuon_endz=0;
		intmuon_id=-1;

		for(unsignedi=0;i<MCpFromMCtruth.at(n_truth).size();i++) {

			//for(autoconst&thisMcp:mcList){
			constart::Ptr<simb::MCParticle>&thisMcp=MCpFromMCtruth.at(n_truth).at(i);

			if(std::abs(thisMcp->PdgCode()) == 13 && thisMcp->StatusCode() == 1 && thisMcp->Mother() == 0&& thisMcp->Process() == "primary" && _rbutilInstance.isInTPC(thisMcp) == true) {

				nimID = thisMcp->TrackId();
				nimMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);
				muMcpMomentum->Fill(thisMcp->P());
				trueVertexXZPosition={(float)thisMcp->Vx(),(float)thisMcp->Vz()};
				muon_p=thisMcp->Momentum().Rho();
				muon_px=thisMcp->Momentum().X();
				muon_py=thisMcp->Momentum().Y();
				muon_pz=thisMcp->Momentum().Z();
				muon_endx=thisMcp->EndPosition().X();
				muon_endy=thisMcp->EndPosition().Y();
				muon_endz=thisMcp->EndPosition().Z();
				muon_id=thisMcp->TrackId();
			}
		}

		if(muon_p==0)mf::LogError(__FUNCTION__)<<"Error!Theremustbeatleast1muoninthisanalysis!"<<std::endl;

		//nowotherMCparticles
		for(unsignedi=0;i<MCpFromMCtruth.at(n_truth).size();i++){
		constart::Ptr<simb::MCParticle>&thisMcp=MCpFromMCtruth.at(n_truth).at(i);
	
		//incaseofNC,theneutrinoiskeptandthusweshoulddropit
		//ifit'sanelectron,checkifit'saMichel.Ifyes,keepit.
		if(!(std::abs(thisMcp->PdgCode())==11&&thisMcp->Mother()==muon_id&&abs(thisMcp->Position().X()-muon_endx)<DBL_EPSILON&&
			abs(thisMcp->Position().Y()-muon_endy)<DBL_EPSILON&&abs(thisMcp->Position().Z()-muon_endz)<DBL_EPSILON))
		if(std::abs(thisMcp->PdgCode())==14||thisMcp->Process()!="primary"||thisMcp->StatusCode()!=1||thisMcp->Mother()>0)continue;//weonlywantprimaries

		if((std::abs(thisMcp->PdgCode())==11&&thisMcp->Mother()==muon_id&&abs(thisMcp->Position().X()-muon_endx)<DBL_EPSILON&&
			abs(thisMcp->Position().Y()-muon_endy)<DBL_EPSILON&&abs(thisMcp->Position().Z()-muon_endz)<DBL_EPSILON))
			std::cout<<">>>>>>>>>>>>>>>>>FOUNDAMICHEL!!!!!"<<std::endl;

		#ifisDebug==1
		std::cout<<"----MCParticleInformation----"
		<<"\nTrackID:"<<thisMcp->TrackId()
		<<"\nPdgCode:"<<thisMcp->PdgCode()
		<<"\nProcess:"<<thisMcp->Process()
		<<"\nStatusCode:"<<thisMcp->StatusCode()
		<<"\nMotherPdg:"<<thisMcp->Mother()
		<<"\nPx,Py,Pz:"<<thisMcp->Px()<<","<<thisMcp->Py()<<","<<thisMcp->Pz()
		<<std::endl;
		#endif

		thisMcpLength.push_back(thisMcp->Trajectory().TotalLength());
		thisMcpMomentum=_rbutilInstance.getMomentumVector(thisMcp);

		//alltracks
		if(thisMcp->Process()=="primary"){

			std::vector<double>zDir={0,0,1};

			thisZDirMcpAngleYZ=_rbutilInstance.getAngle(zDir,thisMcpMomentum,_rbutilInstance,"yz");
			allMcpLengthAngleYZ->Fill(thisZDirMcpAngleYZ,thisMcp->Trajectory().TotalLength());
		}

		//specificallyprotons
		if((thisMcp->PdgCode() == 2212) && (thisMcp->Process() == "primary") && (nimID!=-1)) {

			thisNimMcpAngle = _rbutilInstance.getAngle(thisMcpMomentum,nimMcpMomentum,_rbutilInstance,"no");
			thisNimMcpAngles.push_back(thisNimMcpAngle);

			thisNimMcpAngleXZ = _rbutilInstance.getAngle(thisMcpMomentum,nimMcpMomentum,_rbutilInstance,"xz");
			thisNimMcpAnglesXZ.push_back(thisNimMcpAngleXZ);

			mupMcpAnglePMom->Fill(thisNimMcpAngleXZ,thisMcp->P());
		}
	
		//countparticletypes
		if(fparticle_count.find(thisMcp->PdgCode()) == fparticle_count.end())//notfound
			fparticle_count[thisMcp->PdgCode()] = 1;
		else
			fparticle_count[thisMcp->PdgCode()] = fparticle_count[thisMcp->PdgCode()]+1;
	
		flength.push_back(thisMcp->Trajectory().TotalLength());
		fn_steps.push_back(thisMcp->Trajectory().size());
		fstartT.push_back(thisMcp->T());

		//ifspacechargeisON,weshouldcorrecttheMCtruthpositiontoperformtrue-recocomparisons(Giuseppe)
		/*autoscecorr=SCE->GetPosOffsets(geo::Point_t(thisMcp->Position().X(),thisMcp->Position().Y(),thisMcp->Position().Z()));
		doubleg4Ticks=detClocks->TPCG4Time2Tick(thisMcp->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
		doublexOffset=theDetector->ConvertTicksToX(g4Ticks,0,0,0)-scecorr.X();
		doubleyOffset=scecorr.Y();
		doublezOffset=scecorr.Z();*/
		//anatreerecipe:

		doublexOffset=0.7-SCE->GetPosOffsets(thisMcp->Position().X(),thisMcp->Position().Y(),thisMcp->Position().Z())[0];
		doubleyOffset=SCE->GetPosOffsets(thisMcp->Position().X(),thisMcp->Position().Y(),thisMcp->Position().Z())[1];
		doublezOffset=SCE->GetPosOffsets(thisMcp->Position().X(),thisMcp->Position().Y(),thisMcp->Position().Z())[2];

		if(!space_charge) {

			xOffset=0;
			yOffset=0;
			zOffset=0;
		}

		fstart_x.push_back(thisMcp->Position().X()+xOffset);
		fstart_y.push_back(thisMcp->Position().Y()+yOffset);
		fstart_z.push_back(thisMcp->Position().Z()+zOffset);
	
		if(space_charge) {

			xOffset=0.7-SCE->GetPosOffsets(thisMcp->EndPosition().X(),thisMcp->EndPosition().Y(),thisMcp->EndPosition().Z())[0];
			yOffset=SCE->GetPosOffsets(thisMcp->EndPosition().X(),thisMcp->EndPosition().Y(),thisMcp->EndPosition().Z())[1];
			zOffset=SCE->GetPosOffsets(thisMcp->EndPosition().X(),thisMcp->EndPosition().Y(),thisMcp->EndPosition().Z())[2];
		}

		fend_x.push_back(thisMcp->EndPosition().X()+xOffset);
		fend_y.push_back(thisMcp->EndPosition().Y()+yOffset);
		fend_z.push_back(thisMcp->EndPosition().Z()+zOffset);
		fmother_id.push_back(thisMcp->Mother());
		fpdg.push_back(thisMcp->PdgCode());	
		fg4_id.push_back(thisMcp->TrackId());
		fp0.push_back(thisMcp->Momentum().Rho());
		fp0x.push_back(thisMcp->Momentum().X());
		fp0y.push_back(thisMcp->Momentum().Y());
		fp0z.push_back(thisMcp->Momentum().Z());
	
		//looponothertracked_particlesanddecideifthisistheleadingbasedoninitialkineticenergy.ThismustbedonebeforefillingthekinEvectorforthecurrentparticle!
		boolis_leading=true;
		floatcurrent_kinE=thisMcp->E()-thisMcp->Mass();

		for(unsignedjj=0;jj<fkinE.size();jj++) {//looponpreviousparticles

			if(fpdg[jj]==thisMcp->PdgCode()) {
	
				if(fkinE[jj]>current_kinE)
				is_leading=false;
			}
		}

		fis_leading.push_back(is_leading);

		//nowwritecurrentkinE
		fkinE.push_back(thisMcp->E()-thisMcp->Mass());

		if(thisMcp->Momentum().Rho()!=0) {

			fcostheta_muon.push_back(muon_px*thisMcp->Momentum().X()/thisMcp->Momentum().Rho()/muon_p
					+muon_py*thisMcp->Momentum().Y()/thisMcp->Momentum().Rho()/muon_p
					+muon_pz*thisMcp->Momentum().Z()/thisMcp->Momentum().Rho()/muon_p);
		} else {

			fcostheta_muon.push_back(-2);
		}
	
		//adddummyentriesfortherecovariables.Theywillbepossiblyupdatedlater,ifamatchingrecotrackisfound
		AllocateRecoVectors();

	}// End of MCParticles


	//----------------------------
	//Tracks
	//----------------------------

	//looptracksanddotruthmatching
	for(autoconst&thisTrack:trackList) {

	std::vector<art::Ptr<simb::MCParticle>>mcps=mcpsFromTracks.at(thisTrack->ID());
	std::vector<art::Ptr<anab::Calorimetry>>calos=caloFromTracks.at(thisTrack->ID());

	if(mcps.size()>1)mf::LogWarning(__FUNCTION__)<<"Warning!!!Morethan1MCparticleassociatedtothesametrack!"<<std::endl;
	if(calos.size()!=3)mf::LogWarning(__FUNCTION__)<<"Warning!!!Morethan1Calorimetryinfoassociatedtothesametrack!"<<std::endl;

	for(autoconst&thisMcp:mcps){

		//thisifisnecessary,sincesometimesaparticleismatchedtoasecondary(electron,etc)->tobechecked
		if(!(std::abs(thisMcp->PdgCode())==11&&thisMcp->Mother()==muon_id&&abs(thisMcp->Position().X()-muon_endx)<DBL_EPSILON&&
				abs(thisMcp->Position().Y()-muon_endy)<DBL_EPSILON&&abs(thisMcp->Position().Z()-muon_endz)<DBL_EPSILON))
		if(std::abs(thisMcp->PdgCode())==14||thisMcp->Process()!="primary"||thisMcp->StatusCode()!=1||thisMcp->Mother()>0)continue;//weonlywantprimaries

			#ifisDebug==1
			std::cout<<"MATCHEDPARTICLE:"<<std::endl;
			std::cout<<"----MCParticleInformation----"
			<<"\nTrackID:"<<thisMcp->TrackId()
			<<"\nPdgCode:"<<thisMcp->PdgCode()
			<<"\nProcess:"<<thisMcp->Process()
			<<"\nStatusCode:"<<thisMcp->StatusCode()
			<<"\nMotherPdg:"<<thisMcp->Mother()
			<<"\nPx,Py,Pz:"<<thisMcp->Px()<<","<<thisMcp->Py()<<","<<thisMcp->Pz()
			<<std::endl;
			#endif
			autoit_found=std::find(fg4_id.begin(),fg4_id.end(),thisMcp->TrackId());
			if(it_found==fg4_id.end())mf::LogError(__FUNCTION__)<<"ERROR!!!Matchedparticlenotfound!"<<std::endl;
			size_tpos=it_found-fg4_id.begin();

			//saveinformationonrecotrack
			if(fis_tracked[pos]) {

				mf::LogDebug()<<"Probablybrokentrack!"<<std::endl;
				fmatch_multiplicity[pos]=fmatch_multiplicity[pos]+1;
				//decideifkeepingtheoldinfoorthenewone-basedontheminimumtracklengthdifference
				if(abs(thisTrack->Length()-flength[pos])>abs(flength_reco[pos]-flength[pos]))
					continue;
				} else {

					fis_tracked[pos]=true;
					fmatch_multiplicity[pos]=fmatch_multiplicity[pos]+1;
				}

				flength_reco[pos]=thisTrack->Length();
				trkf::TrackMomentumCalculatortrkm;//trackmomentumcalculator
				trkm.SetMinLength(0);//minimumtracklengthformomentumcalculation
				freco_momentum_mcs[pos]=trkm.GetMomentumMultiScatterChi2(thisTrack);
				freco_momentum_mcs_llhd[pos]=trkm.GetMomentumMultiScatterLLHD(thisTrack);
				freco_momentum_range[pos]=trkm.GetTrackMomentum(thisTrack->Length(),fpdg[pos]);//useinfoonpdg
				freco_startx[pos]=thisTrack->Vertex().X();
				freco_starty[pos]=thisTrack->Vertex().Y();
				freco_startz[pos]=thisTrack->Vertex().Z();
				freco_endx[pos]=thisTrack->End().X();
				freco_endy[pos]=thisTrack->End().Y();
				freco_endz[pos]=thisTrack->End().Z();
	
				if(thisMcp->PdgCode() == 13) {//forthemuon,onlywhenthereistheinfo

					if(fmuon_dqdx.size()!=0)mf::LogError(__FUNCTION__)<<"Calorimetryshouldbefilledonlyonce!!!!"<<std::endl;
					
					for(size_tposition=0;position<calos.at(2)->dQdx().size();position++) {

						fmuon_dqdx.push_back(calos.at(2)->dQdx()[position]);//lookonlyatcollection
						fmuon_dedx.push_back(calos.at(2)->dEdx()[position]);//lookonlyatcollection
						fmuon_residual.push_back(calos.at(2)->ResidualRange()[position]);//lookonlyatcollection
					}
					fmuon_range=calos.at(2)->Range();
				}

				autothisMcpCleanliness=mcpsFromTracks.data(thisTrack->ID()).at(0)->cleanliness;
				autothisMcpCompleteness=mcpsFromTracks.data(thisTrack->ID()).at(0)->completeness;
				fpurity[pos]=thisMcpCleanliness;
				fcompleteness[pos]=thisMcpCompleteness;
			}//MCParticle
	}//Tracks

/*
//looptracksanddotruthmatchingtofindneutrinoinducedmuon
recoNimID=-1;
boolisMatched=false;
for(autoconst&thisTrack:(*trackHandle)){

std::vector<art::Ptr<simb::MCParticle>>mcps=mcpsFromTracks.at(thisTrack.ID());

for(autoconst&thisMcp:mcps){

if(thisMcp->PdgCode()==13&&thisMcp->Process()=="primary"&&_rbutilInstance.isInTPC(thisMcp)==true&&isMatched==false){
recoNimID=thisMcp->TrackId();
nimMatchedMomentum=_rbutilInstance.getMomentumVector(thisMcp);
muMatchedMcpMomentum->Fill(thisMcp->P());
isMatched=true;
}

}

}

intit=0;
for(autoconst&thisTrack:(*trackHandle)){

thisRecoLength.push_back(thisTrack.Length());

std::vector<art::Ptr<simb::MCParticle>>mcps=mcpsFromTracks.at(it);
isRecoTrackTruthMatched.push_back(mcps.size());

//checktomakesuretworeconstructedtracksaren'tmatchingtothesamemcp
//i.e.shouldremovebrokentracks
boolisMatched=false;
for(autoconst&thisMcp:mcps){

for(size_ti=0;i<matchIDChecker.size();i++){

if(thisMcp->TrackId()==matchIDChecker.at(i))
isMatched=true;

}
matchIDChecker.push_back(thisMcp->TrackId());

if(isMatched==true)continue;

autothisMcpCleanliness=mcpsFromTracks.data(0).at(0)->cleanliness;
//autothisMcpCompleteness=mcpsFromTracks.data(0).at(0)->completeness;//givesnonsense

//fillingcleanlinessplots
if(thisMcp->Process()=="primary"&&std::abs(thisMcp->PdgCode())==2212)
trackCleanlinessPrimaryProton->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="primary"
&&((std::abs(thisMcp->PdgCode())==13)
||(std::abs(thisMcp->PdgCode())==211)))
trackCleanlinessPrimaryMuonOrPion->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="primary"&&std::abs(thisMcp->PdgCode())==11)
trackCleanlinessPrimaryElectron->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="conv")
trackCleanlinessConv->Fill(thisMcpCleanliness);
elseif((thisMcp->Process()=="neutronInelastic")
||(thisMcp->Process()=="protonInelastic")
||(thisMcp->Process()=="pi+Inelastic"))
trackCleanlinessInelastic->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="muIoni")
trackCleanlinessMuIoni->Fill(thisMcpCleanliness);
elsetrackCleanlinessOther->Fill(thisMcpCleanliness);

if((thisMcp->Process()!="primary")
||(std::abs(thisMcp->PdgCode())==2112)
||(std::abs(thisMcp->PdgCode())==14)
||(std::abs(thisMcp->PdgCode())==12)
||(std::abs(thisMcp->PdgCode())==22)
||(std::abs(thisMcp->PdgCode())==111)
||(std::abs(thisMcp->PdgCode())==11)//removeanythingwhichisshowerlikefornow
||(std::abs(thisMcp->PdgCode())==2212&&thisMcp->P()<0.2)
||(_rbutilInstance.isInTPC(thisMcp))==false)continue;

#ifisDebug==1
std::cout<<"MATCHCLEANLINESS:"<<thisMcpCleanliness<<std::endl;
std::cout<<"MATCHCOMPLETENESS:"<<thisMcpCompleteness<<std::endl;

std::cout<<"----MCParticleInformation----"
<<"\nTrackID:"<<thisMcp->TrackId()
<<"\nPdgCode:"<<thisMcp->PdgCode()
<<"\nProcess:"<<thisMcp->Process()
<<"\nStatusCode:"<<thisMcp->StatusCode()
<<"\nMotherPdg:"<<thisMcp->Mother()
<<"\nPx,Py,Pz:"<<thisMcp->Px()<<","<<thisMcp->Py()<<","<<thisMcp->Pz()
<<std::endl;
#endif

//angularinformation

thisMatchedMomentum=_rbutilInstance.getMomentumVector(thisMcp);
thisMatchedLength=thisMcp->Trajectory().TotalLength();

//alltracks
if(thisMcp->Process()=="primary"){

std::vector<double>zDir={0,0,1};

thisZmatchedAngleYZ=_rbutilInstance.getAngle(zDir,thisMatchedMomentum,_rbutilInstance,"yz");

allMatchedMcpLengthAngle->Fill(thisZmatchedAngleYZ,thisMatchedLength);


}

//specificallyprotons
if(thisMcp->PdgCode()==2212&&thisMcp->Process()=="primary"&&recoNimID!=-1){

thisNimMatchedMcpAngle=_rbutilInstance.getAngle(nimMatchedMomentum,thisMatchedMomentum,_rbutilInstance,"no");
thisNimMatchedMcpAngles.push_back(thisNimMatchedMcpAngle);

thisNimMatchedMcpAngleXZ=_rbutilInstance.getAngle(nimMatchedMomentum,thisMatchedMomentum,_rbutilInstance,"xz");
thisNimMatchedMcpAnglesXZ.push_back(thisNimMatchedMcpAngleXZ);

mupMatchedMcpAnglePMom->Fill(thisNimMatchedMcpAngleXZ,thisMcp->P());
}


}
it++;
}//trackHandle

//------------------------------
//Showers
//------------------------------

it=0;
for(autoconst&thisShower:(*showerHandle)){

//thisistostopLArSoftcomplaining...
std::cout<<"FoundshowerwithID"<<thisShower.ID()<<std::endl;
std::vector<art::Ptr<simb::MCParticle>>mcps=mcpsFromShowers.at(it);

isRecoShowerTruthMatched.push_back(mcps.size());

for(autoconst&thisMcp:mcps){

autothisMcpCleanliness=mcpsFromShowers.data(0).at(0)->cleanliness;
//autothisMcpCompleteness=mcpsFromShowers.data(0).at(0)->completeness;
#ifisDebug==1
std::cout<<"MATCHCLEANLINESS:"<<thisMcpCleanliness<<std::endl;
std::cout<<"MATCHCOMPLETENESS:"<<thisMcpCompleteness<<std::endl;
std::cout<<"----MCParticleInformation----"
<<"\nTrackID:"<<thisMcp->TrackId()
<<"\nPdgCode:"<<thisMcp->PdgCode()
<<"\nProcess:"<<thisMcp->Process()
<<"\nStatusCode:"<<thisMcp->StatusCode()
<<"\nMotherPdg:"<<thisMcp->Mother()
<<"\nPx,Py,Pz:"<<thisMcp->Px()<<","<<thisMcp->Py()<<","<<thisMcp->Pz()
<<std::endl;
#endif

if(thisMcp->Process()=="primary"&&std::abs(thisMcp->PdgCode())==2212)
showerCleanlinessPrimaryProton->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="primary"
&&((std::abs(thisMcp->PdgCode())==13)
||(std::abs(thisMcp->PdgCode())==211)))
showerCleanlinessPrimaryMuonOrPion->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="primary"&&std::abs(thisMcp->PdgCode())==11)
showerCleanlinessPrimaryElectron->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="conv")
showerCleanlinessConv->Fill(thisMcpCleanliness);
elseif((thisMcp->Process()=="neutronInelastic")
||(thisMcp->Process()=="protonInelastic")
||(thisMcp->Process()=="pi+Inelastic"))
showerCleanlinessInelastic->Fill(thisMcpCleanliness);
elseif(thisMcp->Process()=="muIoni")
showerCleanlinessMuIoni->Fill(thisMcpCleanliness);
elseshowerCleanlinessOther->Fill(thisMcpCleanliness);


}
it++;
}//showerHandle

//---------------------------
//Hits
//---------------------------

for(autoconst&thisHit:(*hitHandle)){

if(trueVertexXZPosition.size()==0)break;

//puthitintohitlist
std::pair<int,float>hitPair;
hitPair.first=(int)thisHit.Channel();
hitPair.second=(float)thisHit.PeakTime();

std::vector<float>hitXZpos=_rbutilInstance.getHitXZPosition(thisHit,_rbutilInstance);

boolisHitInRange=_rbutilInstance.isHitNearVertex(trueVertexXZPosition,hitXZpos);

if(isHitInRange==true)
allHitPositions.push_back(hitPair);

}

//---------------------------
//Clusters
//---------------------------

it=0;
for(autoconst&thisCluster:(*clusterHandle)){
if(trueVertexXZPosition.size()==0)break;

std::cout<<thisCluster.ID()<<std::endl;

//getassociatedhits
std::vector<art::Ptr<recob::Hit>>hits=hitsFromClusters.at(it);

for(autoconst&thisHit:hits){

intisMatched=0;
inthitChannel=thisHit->Channel();
inthitPeakTime=thisHit->PeakTime();

std::vector<float>hitXZpos=_rbutilInstance.getHitXZPosition(*thisHit,_rbutilInstance);

boolisHitInRange=_rbutilInstance.isHitNearVertex(trueVertexXZPosition,hitXZpos);

if(isHitInRange==false)continue;

for(size_ti=0;i<allHitPositions.size();i++){

if(hitChannel==allHitPositions.at(i).first&&hitPeakTime==allHitPositions.at(i).second){

isMatched=1;
hitMatchScore->Fill(isMatched);

}

}


}

hitMatchScore->SetBinContent(0,(int)allHitPositions.size()-hitMatchScore->GetBinContent(1));

it++;
}
*/
		recoTree->Fill();
		FillAnalysisHistograms();

	}//MCtruth
}

//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::AllocateAnalysisHistograms() {

	hmuon_pos_res=tfs->make<TH1D>("muon_pos_res","Allrecomuons;Reco-TrueMuonstartposition(cm);",1000,0,50);//displacementbetweenrealmuonstartpositionandreco
	hmuon_pos_res_goodprotons=tfs->make<TH1D>("muon_pos_res_goodprotons","Allrecomuonswithgoodprotonreco(>20MeV);Reco-TrueMuonstartposition(cm);",1000,0,50);//displacementbetweenrealmuonstartpositionandrecoforeventswithallrecoprotons
	hmuon_pos_res_badprotons=tfs->make<TH1D>("muon_pos_res_badprotons","Allrecomuonswithatleastabadprotonreco(>20MeV);Reco-TrueMuonstartposition(cm);",1000,0,50);//displacementbetweenrealmuonstartpositionandrecoforeventswithsomenon-recoprotons
	hmuon_pos_res_lowprotons=tfs->make<TH1D>("muon_pos_res_lowprotons","Allrecomuonswhenthereareprotons<20MeV;Reco-TrueMuonstartposition(cm);",1000,0,50);//displacementbetweenrealmuonstartpositionandrecoforeventswithsomenon-recolowenergyprotons
	hproton_pos_res=tfs->make<TH1D>("proton_pos_res","Allrecoprotons;Reco-TrueProtonstartposition(cm);",1000,0,50);//displacementbetweenrealprotonstartpositionandreco
	hproton_pos_res_goodprotons=tfs->make<TH1D>("proton_pos_res_goodprotons","AllrecoprotonswithE>20MeV;Reco-TrueProtonstartposition(cm);",1000,0,50);//displacementbetweenrealprotonstartpositionandreco
	hproton_pos_res_badprotons=tfs->make<TH1D>("proton_pos_res_badprotons","Allrecoprotonswithatleastabadprotonreco(>20MeV);Reco-TrueProtonstartposition(cm);",1000,0,50);//displacementbetweenrealprotonstartpositionandreco
	hproton_pos_res_lowprotons=tfs->make<TH1D>("proton_pos_res_lowprotons","Allrecoprotonswhenthereareprotons<20MeV;Reco-TrueProtonstartposition(cm);",1000,0,50);//displacementbetweenrealprotonstartpositionandreco
	hmuon_proton_tracked=tfs->make<TH1D>("muon_proton_tracked","RecoProtonStart-RecoMuonStart;Proton-Muon(cm);",1000,0,50);//displacementbetweenrecomuonstartpositionandprotonrecostartposition
	hmuon_spectrum=tfs->make<TH1D>("muon_spectrum","Muonkineticenergy;KineticEnergy(GeV);",1000,0,5);//recomuons
	hmuon_spectrum_all=tfs->make<TH1D>("muon_spectrum_all","Muonkineticenergy;KineticEnergy(GeV);",1000,0,5);//allofthem,notjustreco
	hmuon_length=tfs->make<TH1D>("muon_length","Muonlength;TrueLength(cm);",1000,0,1000);//recomuons
	hmuon_length_all=tfs->make<TH1D>("muon_length_all","Muonlength;TrueLength(cm);",1000,0,1000);//allofthem,notjustreco
	hproton_kinE=tfs->make<TH1D>("proton_kinE","Protonrecoefficiency;KineticEnergy(GeV)",1000,0,2);//recoefficiencyprotonsvskinE
	hproton_kinE_all=tfs->make<TH1D>("proton_kinE_all","Protonrecoefficiency;KineticEnergy(GeV)",1000,0,2);//recoefficiencyprotonsvskinE
	hproton_p=tfs->make<TH1D>("proton_p","Protonrecoefficiency;Momentum(MeV/c)",1000,0,10);//recoefficiencyprotonsvsp
	hproton_p_all=tfs->make<TH1D>("proton_p_all","Protonrecoefficiency;Momentum(MeV/c)",1000,0,10);//recoefficiencyprotonsvsp
	hproton_l=tfs->make<TH1D>("proton_l","Protonrecoefficiency;Truelength(cm)",1000,0,200);//recoefficiencyprotonsvsp
	hproton_l_all=tfs->make<TH1D>("proton_l_all","Protonrecoefficiency;Truelength(cm)",1000,0,200);//recoefficiencyprotonsvsp
	h_pmu_end_not_tracked=tfs->make<TH1D>("pmu_end_not_tracked","NotTrackedprotons;Distance(cm);",1000,0,100);//lateraldistancebetweenprotonendandmuon
	h_pmu_end_tracked=tfs->make<TH1D>("pmu_end_tracked","Trackedprotons;Distance(cm);",1000,0,100);//lateraldistancebetweenprotonendandmuon
	h_theta_mu_tracked=tfs->make<TH1D>("theta_mu_tracked","Cos#thetabetweenmuonandtrackedprotons;cos#theta",1000,-1,1);//costhetabetweenmuonandtrackedprotons
	h_theta_mu_not_tracked=tfs->make<TH1D>("theta_mu_not_tracked","Cos#thetabetweenmuonandnontrackedprotons;cos#theta",1000,-1,1);
	h_theta_mu=tfs->make<TH1D>("theta_mu","Cos#thetabetweenmuonandprotons;cos#theta",1000,-1,1);//truecosthetabetweenmuonandprotons
	h_theta_mu_length=tfs->make<TH2D>("theta_mu_length","Trackingefficiencyvs(length,cos#theta_{p#mu});Cos#theta;l(cm)",1000,-1,1,1000,0,100);
	h_theta_mu_length_all=tfs->make<TH2D>("theta_mu_length_all","Trackingefficiencyvs(length,cos#theta_{p#mu});Cos#theta;l(cm)",1000,-1,1,1000,0,100);
	h_dqdx_merged=tfs->make<TH2D>("dqdx_merged","dq/dxforeventswithatleastamergedproton;Distancefromvertex(cm);dq/dx(ADC)",2500,0,250,1500,0,1500);
	h_dqdx_not_merged=tfs->make<TH2D>("dqdx_not_merged","dq/dxforeventswithnomergedproton;Distancefromvertex(cm);dq/dx(ADC)",2500,0,250,1500,0,1500);
	h_dqdx_low_protons=tfs->make<TH2D>("dqdx_low_protons","dq/dxforeventswithlowEproton;Distancefromvertex(cm);dq/dx(ADC)",2500,0,250,1500,0,1500);
	h_dqdx_1d_merged=tfs->make<TH1D>("dqdx_1d_merged","dq/dxintegratedin(0,8cm)whenprotonsaremerged;dq/dx(ADC);",1500,0,1500);
	h_dqdx_1d_not_merged=tfs->make<TH1D>("dqdx_1d_not_merged","dq/dxintegratedin(0,8cm)whenprotonsarenotmerged;dq/dx(ADC);",1500,0,1500);
	h_dqdx_tailtotot_length_merged=tfs->make<TH2D>("dqdx_tailtotot_length_merged","Tailtototvslengthformergedtracks;IntegrationLength(mm);Tailtotot",1000,0,1000,1000,0,1);
	h_dqdx_tailtotot_length_not_merged=tfs->make<TH2D>("dqdx_tailtotot_length_not_merged","Tailtototvslengthformergedtracks;IntegrationLength(mm);Tailtotot",1000,0,1000,1000,0,1);
	htail_to_tot_low_protons=tfs->make<TH1D>("tailtotot_low_protons","Tailtototw/lowenergyprotons;Tailtotot;",1000,0,1);//tailtototmerged
	htail_to_tot_merged=tfs->make<TH1D>("tailtotot_merged","Tailtototw/mergedprotons;Tailtotot;",100,0,100);//tailtototmerged
	htail_to_tot_not_merged=tfs->make<TH1D>("tailtotot_not_merged","Tailtototw/goodprotons;Tailtotot;",100,0,100);//tailtototnotmerged
	h_dqdx_merged_service=tfs->make<TH2D>("dqdx_merged_service","dqdx_merged_service",2500,0,250,1500,0,1500);
	h_dqdx_not_merged_service=tfs->make<TH2D>("dqdx_not_merged_service","dqdx_not_merged_service",2500,0,250,1500,0,1500);
	h_dqdx_low_protons_service=tfs->make<TH2D>("dqdx_low_protons_service","dqdx_low_protons_service",2500,0,250,1500,0,1500);
}
//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::FillAnalysisHistograms(){

	if(fccnc!=0)return;//commentifyou'rerunningoverastdsamplewithNCinteractionsandyou'reinterestedinthem
	unsignedmuon_pos=-1;
	//checkthatthemuonisreco
	boolreco_muon=false;
	boolis_pion=false;
	boolis_lowmomentum_p=false;
	for(unsignedi=0;i<fpdg.size();i++) {

		if(fpdg[i]==13){//ismuon

			hmuon_length_all->Fill(flength[i]);
			hmuon_spectrum_all->Fill(fkinE[i]);
			if(fis_tracked[i]&&fmuon_dqdx.size()!=0) {//wantthatthemuoniswellmatchedandw/dqdxinfo

				muon_pos=i;
				hmuon_length->Fill(flength[i]);
				hmuon_spectrum->Fill(fkinE[i]);
				reco_muon=true;
				hmuon_pos_res->Fill(sqrt(pow(freco_startx[i]-fstart_x[i],2)+pow(freco_starty[i]-fstart_y[i],2)+pow(freco_startz[i]-fstart_z[i],2)));

				if(sqrt(pow(freco_startx[i]-fstart_x[i],2)+pow(freco_starty[i]-fstart_y[i],2)+pow(freco_startz[i]-fstart_z[i],2))>50) {

					reco_muon=false;//skiptheseevents,mighthavedirectionflipped
				}
			}
		}//recordinfoonthemuon
		if(fpdg[i]==211||fpdg[i]==-211||fpdg[i]==111)is_pion=true;//recordifthereareanypions
		if(fpdg[i]==2212&&fp0[i]<=0.2)is_lowmomentum_p=true;
	}

	longcount_not_tracked=0;
	longcount_tracked=0;

	for(unsignedj=0;j<fpdg.size();j++) {

		if(!reco_muon)break;//selecteventswitharecomuon
		if(fpdg[j]!=2212)continue;//watchonlyprotons
	
		hproton_p_all->Fill(fp0[j]);
		hproton_l_all->Fill(flength[j]);
		hproton_kinE_all->Fill(fkinE[j]);
		h_theta_mu_length_all->Fill(fcostheta_muon[j],flength[j]);

		if(fis_tracked[j]) {//efficiencyplotsforallprotons(butdon'tdivideyet)

			hproton_p->Fill(fp0[j]);
			hproton_l->Fill(flength[j]);
			hproton_kinE->Fill(fkinE[j]);
			h_theta_mu_length->Fill(fcostheta_muon[j],flength[j]);
		}

		//nowstudyalltheprotonswithmomentum>200MeV
		if(fp0[j]>0.2) {

			h_theta_mu->Fill(fcostheta_muon[j]);	

			if(!fis_tracked[j]) {

				count_not_tracked++;
				h_pmu_end_not_tracked->Fill(flength[j]*sqrt(1-pow(fcostheta_muon[j],2)));
				h_theta_mu_not_tracked->Fill(fcostheta_muon[j]);
			} else if(fis_tracked[j]) {

				count_tracked++;
				h_pmu_end_tracked->Fill(flength[j]*sqrt(1-pow(fcostheta_muon[j],2)));
				h_theta_mu_tracked->Fill(fcostheta_muon[j]);
				hmuon_proton_tracked->Fill(sqrt(pow(freco_startx[j]-freco_startx[muon_pos],2)+pow(freco_starty[j]-freco_starty[muon_pos],2)+pow(freco_startz[j]-freco_startz[muon_pos],2)));
				hproton_pos_res->Fill(sqrt(pow(freco_startx[j]-fstart_x[j],2)+pow(freco_starty[j]-fstart_y[j],2)+pow(freco_startz[j]-fstart_z[j],2)));
			}
		}
	}

	if(fmuon_residual.size()!=fmuon_dqdx.size())cout<<"ERRORoncalorimetryvectorsizes!!!"<<endl;

	if(count_not_tracked==0&&count_tracked>0) {//allprotonsaretracked

		for(unsignedjj=1;jj<fmuon_dqdx.size()-1;jj++) {

				if(fmuon_range-fmuon_residual[jj]<8) {//lookat8cmonly

					h_dqdx_1d_not_merged->Fill(fmuon_dqdx[jj]);
				}

				h_dqdx_not_merged->Fill(fmuon_range-fmuon_residual[jj],fmuon_dqdx[jj]);
				h_dqdx_not_merged_service->Fill(fmuon_range-fmuon_residual[jj],fmuon_dqdx[jj]);
		}

		TH1D*h1=NULL;
		h1=h_dqdx_not_merged_service->ProjectionY("h",1,length_cut);
		if(h1)htail_to_tot_not_merged->Fill(h1->Integral(low_edge,high_edge));
		h_dqdx_not_merged_service->Reset();

		hmuon_pos_res_goodprotons->Fill(sqrt(pow(freco_startx[muon_pos]-fstart_x[muon_pos],2)+pow(freco_starty[muon_pos]-fstart_y[muon_pos],2)+pow(freco_startz[muon_pos]-fstart_z[muon_pos],2)));
		for(unsignedj=0;j<fpdg.size();j++) {

			if(!reco_muon)break;//selecteventswitharecomuon
			if(fpdg[j]!=2212&&fp0[j]<=0.2&&!fis_tracked[j])continue;//watchonlyrecoprotonswithp>0.2Gev/c
			hproton_pos_res_goodprotons->Fill(sqrt(pow(freco_startx[j]-fstart_x[j],2)+pow(freco_starty[j]-fstart_y[j],2)+pow(freco_startz[j]-fstart_z[j],2)));
			if(is_lowmomentum_p)hproton_pos_res_lowprotons->Fill(sqrt(pow(freco_startx[j]-fstart_x[j],2)+pow(freco_starty[j]-fstart_y[j],2)+pow(freco_startz[j]-fstart_z[j],2)));
		}
	}

	if(count_not_tracked>0) {//atleast1protonisnottracked

		for(unsignedjj=1;jj<fmuon_dqdx.size()-1;jj++) {

				if(fmuon_range-fmuon_residual[jj]<8)//lookat8cmonly
					h_dqdx_1d_merged->Fill(fmuon_dqdx[jj]);

				h_dqdx_merged->Fill(fmuon_range-fmuon_residual[jj],fmuon_dqdx[jj]);
				h_dqdx_merged_service->Fill(fmuon_range-fmuon_residual[jj],fmuon_dqdx[jj]);
		}

		TH1D*h1=NULL;
		h1=h_dqdx_merged_service->ProjectionY("h",1,length_cut);
		if(h1) htail_to_tot_merged->Fill(h1->Integral(low_edge,high_edge));
		h_dqdx_merged_service->Reset();
		
		hmuon_pos_res_badprotons->Fill(sqrt(pow(freco_startx[muon_pos]-fstart_x[muon_pos],2)+pow(freco_starty[muon_pos]-fstart_y[muon_pos],2)+pow(freco_startz[muon_pos]-fstart_z[muon_pos],2)));
		for(unsignedj=0;j<fpdg.size();j++) {

			if(!reco_muon)break;//selecteventswitharecomuon
			if(fpdg[j]!=2212&&fp0[j]<=0.2&&!fis_tracked[j])continue;//watchonlyrecoprotonswithp>0.2Gev/c
			hproton_pos_res_badprotons->Fill(sqrt(pow(freco_startx[j]-fstart_x[j],2)+pow(freco_starty[j]-fstart_y[j],2)+pow(freco_startz[j]-fstart_z[j],2)));
			if(is_lowmomentum_p)hproton_pos_res_lowprotons->Fill(sqrt(pow(freco_startx[j]-fstart_x[j],2)+pow(freco_starty[j]-fstart_y[j],2)+pow(freco_startz[j]-fstart_z[j],2)));
		}
	}

	if(reco_muon&&count_tracked==0&&count_not_tracked==0&&!is_pion&&is_lowmomentum_p) {//noprotonstrackedbutlowenergyones

		for(unsignedjj=1;jj<fmuon_dqdx.size()-1;jj++) {

				h_dqdx_low_protons->Fill(fmuon_range-fmuon_residual[jj],fmuon_dqdx[jj]);
				h_dqdx_low_protons_service->Fill(fmuon_range-fmuon_residual[jj],fmuon_dqdx[jj]);
		}

		TH1D*h1=NULL;
		h1=h_dqdx_low_protons_service->ProjectionY("h",1,length_cut);
		if(h1) htail_to_tot_low_protons->Fill(h1->Integral(low_edge,high_edge)/h1->Integral(1,high_edge));
		h_dqdx_low_protons_service->Reset();
		
		hmuon_pos_res_lowprotons->Fill(sqrt(pow(freco_startx[muon_pos]-fstart_x[muon_pos],2)+pow(freco_starty[muon_pos]-fstart_y[muon_pos],2)+pow(freco_startz[muon_pos]-fstart_z[muon_pos],2)));
	}


}
//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::FillCumulativeHistograms() {

	for(longii=1;ii<1000;ii++) {

		TH1D*h=h_dqdx_not_merged->ProjectionY("h",1,ii);
		h_dqdx_tailtotot_length_not_merged->Fill(ii,h->Integral(low_edge,high_edge)/h->Integral(1,high_edge));
		h=h_dqdx_merged->ProjectionY("h",1,ii);
		h_dqdx_tailtotot_length_merged->Fill(ii,h->Integral(low_edge,high_edge)/h->Integral(1,high_edge));
	}
}
//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::AllocateRecoVectors() {

	fis_tracked.push_back(0);
	fmatch_multiplicity.push_back(0);
	flength_reco.push_back(-1);
	freco_momentum_mcs.push_back(-1);
	freco_momentum_mcs_llhd.push_back(-1);
	freco_momentum_range.push_back(-1);
	fpurity.push_back(-1);
	fcompleteness.push_back(-1);
	fnhits.push_back(-1);
	freco_kinE.push_back(-1);
	freco_startx.push_back(-1);
	freco_starty.push_back(-1);
	freco_startz.push_back(-1);
	freco_endx.push_back(-1);
	freco_endy.push_back(-1);
	freco_endz.push_back(-1);
}
//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::endJob() {

	TFile&file=tfs->file();
	file.cd();

	//--------------------------
	//Produceeff.plots
	//--------------------------

	//muefficiencies
	muMomentumEfficiency=(TH1D*)muMatchedMcpMomentum->Clone("muMomentumEfficiency");
	muMomentumEfficiency->Divide(muMcpMomentum);

	//mupangleMomentumEfficiencies
	mupAngleMomentumEfficiency=(TH2D*)mupMatchedMcpAnglePMom->Clone("mupAngleMomentumEfficiency");
	mupAngleMomentumEfficiency->Divide(mupMcpAnglePMom);

	//pMomentumEfficiencies
	mupAngleMomentumEfficiency=(TH2D*)mupMatchedMcpAnglePMom->Clone("mupAngleMomentumEfficiency");
	mupAngleMomentumEfficiency->Divide(mupMcpAnglePMom);

	pMatchedMcpProjectedAngle=(TH1D*)mupMatchedMcpAnglePMom->ProjectionY();
	pMcpProjectedAngle=(TH1D*)mupMcpAnglePMom->ProjectionY();

	pProjectedAngleEfficiency=(TH1D*)pMatchedMcpProjectedAngle->Clone("pProjectedAngleEfficiency");
	pProjectedAngleEfficiency->Divide(pMcpProjectedAngle);

	pMatchedMcpProjectedMomentum=(TH1D*)mupMatchedMcpAnglePMom->ProjectionX();
	pMcpProjectedMomentum=(TH1D*)mupMcpAnglePMom->ProjectionX();

	pProjectedMomentumEfficiency=(TH1D*)pMatchedMcpProjectedMomentum->Clone("pProjectedMomentumEfficiency");
	pProjectedMomentumEfficiency->Divide(pMcpProjectedMomentum);

	//trackAngleLengthEfficiencies
	allLengthAngleEfficiency=(TH2D*)allMatchedMcpLengthAngle->Clone("allLengthAngleEfficiency");
	allLengthAngleEfficiency->Divide(allMcpLengthAngleYZ);

	allMatchedMcpProjectedLength=(TH1D*)allMatchedMcpLengthAngle->ProjectionY();
	allMcpProjectedLength=(TH1D*)allMcpLengthAngleYZ->ProjectionY();

	allProjectedLengthEfficiency=(TH1D*)allMatchedMcpProjectedLength->Clone("allProjectedLengthEfficiency");
	allProjectedLengthEfficiency->Divide(allMcpProjectedLength);

	allMatchedMcpProjectedAngle=(TH1D*)allMatchedMcpLengthAngle->ProjectionX();
	allMcpProjectedAngle=(TH1D*)allMcpLengthAngleYZ->ProjectionX();

	allProjectedAngleEfficiency=(TH1D*)allMatchedMcpProjectedAngle->Clone("allProjectedAngleEfficiency");
	allProjectedAngleEfficiency->Divide(allMcpProjectedAngle);

	//aaaaandwrite...
	muMomentumEfficiency->Write();
	mupAngleMomentumEfficiency->Write();
	pMatchedMcpProjectedMomentum->Write();
	pMcpProjectedMomentum->Write();
	pProjectedMomentumEfficiency->Write();
	pMatchedMcpProjectedAngle->Write();
	pMcpProjectedAngle->Write();
	pProjectedAngleEfficiency->Write();
	allLengthAngleEfficiency->Write();
	allMatchedMcpProjectedLength->Write();
	allMcpProjectedLength->Write();
	allProjectedLengthEfficiency->Write();
	allMatchedMcpProjectedAngle->Write();
	allMcpProjectedAngle->Write();
	allProjectedAngleEfficiency->Write();

	FillCumulativeHistograms();
}
//____________________________________________________________________________________________________________________________________________________________________________________________________

voidrecohelper::RecoBenchmarker::clear_vectors() {

	//infoontheMCparticle
	fmuon_dqdx.clear();
	fmuon_dedx.clear();
	fmuon_residual.clear();
	flength.clear();
	fstartT.clear();
	fstart_x.clear();
	fstart_y.clear();
	fstart_z.clear();
	fend_x.clear();
	fend_y.clear();
	fend_z.clear();
	fn_steps.clear();
	fpdg.clear();
	fmother_id.clear();
	fg4_id.clear();
	fp0.clear();//initialmomentum
	fp0x.clear();
	fp0y.clear();
	fp0z.clear();
	fkinE.clear();
	fcostheta_muon.clear();
	fis_leading.clear();
	fparticle_count.clear();

	//infocomingfromthetrackingalgorithm-whenthereismctruth
	fis_tracked.clear();
	fmatch_multiplicity.clear();
	fcostheta_muon_reco.clear();//MCtruthinfo
	freco_momentum_mcs.clear();//(GeV)MCS
	freco_momentum_mcs_llhd.clear();//(GeV)MCSLLHD
	freco_momentum_range.clear();//MeV
	freco_startx.clear();
	freco_starty.clear();
	freco_startz.clear();
	freco_endx.clear();
	freco_endy.clear();
	freco_endz.clear();
	flength_reco.clear();
	fpurity.clear();
	fcompleteness.clear();
	fnhits.clear();
	fkinE.clear();

	//setupvariables

	isRecoTrackTruthMatched.clear();
	isRecoShowerTruthMatched.clear();
	matchIDChecker.clear();
	thisNimMcpAngles.clear();
	thisNimMcpAnglesXZ.clear();
	thisMcpMomentum.clear();
	nimMcpMomentum.clear();
	thisMcpLength.clear();
	thisNimMatchedMcpAngles.clear();
	thisNimMatchedMcpAnglesXZ.clear();
	nimMatchedMomentum.clear();
	thisRecoLength.clear();
	trueVertexXZPosition.clear();
	allHitPositions.clear();

	/*//infocomingfromthetrackingalgorithm-whenthereisNOmctruth
	ffake_is_tracked.clear();
	ffake_is_mismatched.clear();
	ffake_costheta_muon_reco.clear();//MCtruthinfo
	ffake_pdg_reco.clear();
	ffake_length_reco.clear();
	ffake_reco_momentum_mcs.clear();//(GeV)MCS
	ffake_reco_momentum_mcs_llhd.clear();//(GeV)MCSLLHD
	ffake_reco_momentum_range.clear();//MeV
	ffake_purity.clear();
	ffake_completeness.clear();
	ffake_nhits.clear();
	ffake_kinE.clear();
	*/
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
