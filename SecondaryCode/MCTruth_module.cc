//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"

#include "TruncMean.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVectorBase.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/MCCheater/BackTracker.h"

namespace mynamespace {
	class Afro;
}


class mynamespace::Afro : public art::EDAnalyzer {

	public:
		
		explicit Afro(fhicl::ParameterSet const & p);
		Afro(Afro const &) = delete;
		Afro(Afro &&) = delete;
		Afro & operator = (Afro const &) = delete;
		Afro & operator = (Afro &&) = delete;

		void analyze(art::Event const & e) override;
		void beginJob() override;
		void endJob() override;

	private:

		art::ServiceHandle<art::TFileService>tfs;

		// TH1D
		TH1D* Purity;

		// TH2D
		TH2D* Muon_dEdx_vs_ResidualRange;
		TH2D* Proton_dEdx_vs_ResidualRange;
		TH2D* Muon_dQdx_vs_ResidualRange;
		TH2D* Proton_dQdx_vs_ResidualRange;
		
		TH2D* Trunc_Muon_dEdx_vs_ResidualRange;
		TH2D* Trunc_Proton_dEdx_vs_ResidualRange;
		TH2D* Trunc_Muon_dQdx_vs_ResidualRange;
		TH2D* Trunc_Proton_dQdx_vs_ResidualRange;

		bool inFV(double x, double y, double z);
		bool inTPC(double x, double y, double z);

};
// ____________________________________________________________________________________________________________________________________________________________________________________________________

mynamespace::Afro::Afro(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) {

}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

void mynamespace::Afro::beginJob()
{
	// TH1D 
	Purity = tfs->make<TH1D>("Purity",";Purity;",100,0,1);

	// TH2D 
	Muon_dEdx_vs_ResidualRange = tfs->make<TH2D>("Muon_dEdx_vs_ResidualRange",";Muon Residual Range [cm]; dE/dx [MeV/cm]",60,0,300,120,0,30);
	Proton_dEdx_vs_ResidualRange = tfs->make<TH2D>("Proton_dEdx_vs_ResidualRange",";Proton Residual Range [cm]; dE/dx [MeV/cm]",60,0,300,120,0,30);
	Muon_dQdx_vs_ResidualRange = tfs->make<TH2D>("Muon_dQdx_vs_ResidualRange",";Muon Residual Range [cm]; dQ/dx [ADC/cm]",60,0,300,120,0,800);
	Proton_dQdx_vs_ResidualRange = tfs->make<TH2D>("Proton_dQdx_vs_ResidualRange",";Proton Residual Range [cm]; dQ/dx [ADC/cm]",60,0,300,120,0,800);

	Trunc_Muon_dEdx_vs_ResidualRange = tfs->make<TH2D>("Trunc_Muon_dEdx_vs_ResidualRange",";Muon Residual Range [cm]; Truncated dE/dx [MeV/cm]",60,0,300,120,0,30);
	Trunc_Proton_dEdx_vs_ResidualRange = tfs->make<TH2D>("Trunc_Proton_dEdx_vs_ResidualRange",";Proton Residual Range [cm]; Truncated dE/dx [MeV/cm]",60,0,300,120,0,30);
	Trunc_Muon_dQdx_vs_ResidualRange = tfs->make<TH2D>("Trunc_Muon_dQdx_vs_ResidualRange",";Muon Residual Range [cm]; Truncated dQ/dx [ADC/cm]",60,0,300,120,0,800);
	Trunc_Proton_dQdx_vs_ResidualRange = tfs->make<TH2D>("Trunc_Proton_dQdx_vs_ResidualRange",";Proton Residual Range [cm]; Truncated dQ/dx [ADC/cm]",60,0,300,120,0,800);
}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

void mynamespace::Afro::analyze(art::Event const & e)
{

	// Valid Handles
		
	// Hits 
	auto const& hit_handle = e.getValidHandle<std::vector<recob::Hit>>("pandoraCosmicHitRemoval");
	/*auto const& hit_vec(*hit_handle);*/

	// Tracks
	auto const& trk_handle = e.getValidHandle<std::vector<recob::Track>>("pandoraNu");
	auto const& trk_vec(*trk_handle);
// ____________________________________________________________________________________________________________________________________________________________________________________________________

	art::FindManyP<recob::Hit> hits_per_track(trk_handle,e,"pandoraNu");
	std::cout << "\tThere are " << trk_vec.size() << " tracks in this event." << std::endl;

	// Track - MCParticle BackTracking
	// Loop over the tracks

	for (int i_t = 0; i_t < int(trk_vec.size()); ++i_t) {

		if (inTPC(trk_vec.at(i_t).Start().X(),trk_vec.at(i_t).Start().Y(),trk_vec.at(i_t).Start().Z()) == false){ continue; }
		if (inTPC(trk_vec.at(i_t).End()[0],trk_vec.at(i_t).End()[1],trk_vec.at(i_t).End()[2]) == false){ continue; }

		if (inFV(trk_vec.at(i_t).Start().X(),trk_vec.at(i_t).Start().Y(),trk_vec.at(i_t).Start().Z()) == false){ continue; }
		if (inFV(trk_vec.at(i_t).End()[0],trk_vec.at(i_t).End()[1],trk_vec.at(i_t).End()[2]) == false){ continue; }

		art::ServiceHandle<cheat::BackTracker> bt_MCParticle_MCTruth;

		std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i_t);
		std::cout << "\tThere are " << trk_hits_ptrs.size() << " associated hits." << std::endl;

		std::unordered_map<int,double> trkide;
		double maxe = -1, tote = 0;
		art::Ptr< simb::MCParticle > maxp_me; // Pointer for the particle match we will calculate

		art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle,e,"crHitRemovalTruthMatch");
		std::vector<art::Ptr<simb::MCParticle>> particle_vec;
		std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

		// Loop only over the hits
		for (int i_h = 0; i_h < int(trk_hits_ptrs.size()); ++i_h) {

			particle_vec.clear(); match_vec.clear();
			particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);

			//loop over particles
			for(int i_p = 0; i_p < int(particle_vec.size()); ++i_p) {

				trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; // Store energy per track id
				tote += match_vec[i_p]->energy; // Calculate total energy deposited
				if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ) { // Keep track of maximum
					
					maxe = trkide[ particle_vec[i_p]->TrackId() ];
					maxp_me = particle_vec[i_p];

				} // End of the if-statement

			} // End of the loop over particles per hit

		} // End of the loop over the hits

		if (maxp_me.isNull()) continue;
		std::cout << "Final Match (from my loop) is " << maxp_me->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
		<< " \n\tpdg=" << maxp_me->PdgCode()
		<< " trkid=" << maxp_me->TrackId()
		<< " ke=" << maxp_me->E()-maxp_me->Mass()
		<< "\n\tstart (x,y,z)=(" << maxp_me->Vx()
		<< "," << maxp_me->Vy()
		<< "," << maxp_me->Vz()
		<< ")\tend (x,y,z)=(" << maxp_me->EndX()
		<< "," << maxp_me->EndY()
		<< "," << maxp_me->EndZ() << ")" << std::endl;

		double purity = maxe/tote;
		Purity->Fill(purity);

		if (inTPC(maxp_me->Vx(),maxp_me->Vy(),maxp_me->Vz()) == false){ continue; }
		if (inTPC(maxp_me->Vx(),maxp_me->Vy(),maxp_me->Vz()) == false){ continue; }

		if (inFV(maxp_me->Vx(),maxp_me->Vy(),maxp_me->Vz()) == false){ continue; }
		if (inFV(maxp_me->Vx(),maxp_me->Vy(),maxp_me->Vz()) == false){ continue; }

		// Calorimetry
		art::FindManyP<anab::Calorimetry> CalorimetryFromTracks(trk_handle,e,"pandoraNucalo");
		std::vector<art::Ptr<anab::Calorimetry>> Calorimetry_vec = CalorimetryFromTracks.at(i_t);

		if (inTPC(trk_vec.at(i_t).Start().X(),trk_vec.at(i_t).Start().Y(),trk_vec.at(i_t).Start().Z()) == false){ continue; }
		if (inTPC(trk_vec.at(i_t).End()[0],trk_vec.at(i_t).End()[1],trk_vec.at(i_t).End()[2]) == false){ continue; }

		if (inFV(trk_vec.at(i_t).Start().X(),trk_vec.at(i_t).Start().Y(),trk_vec.at(i_t).Start().Z()) == false){ continue; }
		if (inFV(trk_vec.at(i_t).End()[0],trk_vec.at(i_t).End()[1],trk_vec.at(i_t).End()[2]) == false){ continue; }

		if ( fabs(maxp_me->Vx() - trk_vec.at(i_t).Start().X()) > 3 || fabs(maxp_me->Vy() - trk_vec.at(i_t).Start().Y()) > 3 
		  || fabs(maxp_me->Vz() - trk_vec.at(i_t).Start().Z()) > 3 ) { continue; }

		if (purity < 0.8) continue;

		for (int WhichPlane = 0; WhichPlane < 3; WhichPlane ++) {

			art::Ptr<anab::Calorimetry> plane = Calorimetry_vec.at(WhichPlane);
			std::vector< double > dedxPlane = plane->dEdx();
			std::vector< double > dqdxPlane = plane->dQdx();
			std::vector< double > ResidualRangePlane = plane->ResidualRange();

			std::vector<float> dedx_array_trunc_muons;
			std::vector<float> dedx_array_original_muons;
			std::vector<float> dqdx_array_trunc_muons;
			std::vector<float> dqdx_array_original_muons;
			std::vector<float> trkresrg_array_muons;

			std::vector<float> dedx_array_trunc_protons;
			std::vector<float> dedx_array_original_protons;
			std::vector<float> dqdx_array_trunc_protons;
			std::vector<float> dqdx_array_original_protons;
			std::vector<float> trkresrg_array_protons;

			// Loop over the hits on the 3 planes
			for (size_t iHit = 0; iHit < dedxPlane.size(); iHit++){

				if (maxp_me->PdgCode() == 13) {

					std::cout << "maxp_me->PdgCode() = " << maxp_me->PdgCode() << "   dedxPlane.at(iHit) = " << dedxPlane.at(iHit) << std::endl;

					dedx_array_original_muons.push_back(dedxPlane.at(iHit));
					dqdx_array_original_muons.push_back(dqdxPlane.at(iHit));
					trkresrg_array_muons.push_back(ResidualRangePlane.at(iHit));

					Muon_dEdx_vs_ResidualRange->Fill(ResidualRangePlane.at(iHit),dedxPlane.at(iHit));
					Muon_dQdx_vs_ResidualRange->Fill(ResidualRangePlane.at(iHit),dqdxPlane.at(iHit));

				} // End of the if-statement for the muons

				if (maxp_me->PdgCode() == 2212) {

					std::cout << "maxp_me->PdgCode() = " << maxp_me->PdgCode() << "   dedxPlane.at(iHit) = " << dedxPlane.at(iHit) << std::endl;

					dedx_array_original_protons.push_back(dedxPlane.at(iHit));
					dqdx_array_original_protons.push_back(dqdxPlane.at(iHit));
					trkresrg_array_protons.push_back(ResidualRangePlane.at(iHit));

					Proton_dEdx_vs_ResidualRange->Fill(ResidualRangePlane.at(iHit),dedxPlane.at(iHit));
					Proton_dQdx_vs_ResidualRange->Fill(ResidualRangePlane.at(iHit),dqdxPlane.at(iHit));

				} // End of the if-statement for the protons

			} // End of the loop over the hits

			TruncMean truncmean;
			truncmean.CalcTruncMean(trkresrg_array_muons,dedx_array_original_muons,dedx_array_trunc_muons);
			truncmean.CalcTruncMean(trkresrg_array_muons,dqdx_array_original_muons,dqdx_array_trunc_muons);
			truncmean.CalcTruncMean(trkresrg_array_protons,dedx_array_original_protons,dedx_array_trunc_protons);
			truncmean.CalcTruncMean(trkresrg_array_protons,dqdx_array_original_protons,dqdx_array_trunc_protons);

			for (int WhichTrunc = 0 ; WhichTrunc < int(dedx_array_trunc_muons.size()) ; WhichTrunc ++) {

				Trunc_Muon_dEdx_vs_ResidualRange->Fill(trkresrg_array_muons[WhichTrunc],dedx_array_trunc_muons[WhichTrunc]);
				Trunc_Muon_dQdx_vs_ResidualRange->Fill(trkresrg_array_muons[WhichTrunc],dqdx_array_trunc_muons[WhichTrunc]);

			} // End of the loop over the truncated array for the muons

			for (int WhichTrunc = 0 ; WhichTrunc < int(dedx_array_trunc_protons.size()) ; WhichTrunc ++) {

				Trunc_Proton_dEdx_vs_ResidualRange->Fill(trkresrg_array_protons[WhichTrunc],dedx_array_trunc_protons[WhichTrunc]);
				Trunc_Proton_dQdx_vs_ResidualRange->Fill(trkresrg_array_protons[WhichTrunc],dqdx_array_trunc_protons[WhichTrunc]);

			} // End of the loop over the truncated array for the protons
	
		} // End of the loop over the 3-planes

	} // End of the loop over the tracks

} // End of the analysis module
// ____________________________________________________________________________________________________________________________________________________________________________________________________

bool mynamespace::Afro::inFV(double x, double y, double z) {

	double FVx = 256.35;
	double FVy = 233;
	double FVz = 1036.8;
	double borderx = 10.;
	double bordery = 20.;
	double borderz = 10.;

	if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
	else return false;
}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

bool mynamespace::Afro::inTPC(double x, double y, double z) {
  
	if(x<-50 || x>300) return false;
	if(y<-100 || y>100) return false;
	if(z<0.01 || z>1050) return false;
	else return true;
}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

void mynamespace::Afro::endJob()
{
	TFile& file = tfs->file();
	file.cd();
}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

DEFINE_ART_MODULE(mynamespace::Afro)
