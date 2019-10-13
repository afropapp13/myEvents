#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStopwatch.h"

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
using namespace std;

bool ParticleAlreadyMatchedInThisHit(std::vector<int> AlreadyMatched_TrackIDs ,int cTrackID );
TString ToString(int num);
// ____________________________________________________________________________________________________________________________________________________________________________________________________

art::Ptr< simb::MCParticle > BackTrackerTruthMatch
(art::Ptr< simb::MCParticle > maxp_me, art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit, std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs) 
{

	std::unordered_map<int,double> trkide;
	double maxe = -1, tote = 0;

	std::vector<art::Ptr<simb::MCParticle>> particle_vec;
	std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

	double max_dQinTruthMatchedHits = -1., dQinAllHits = 0.;
	std::unordered_map<int,double> trackid_dQinTruthMatchedHits;

	// Loop only over the hits
	for (int i_h = 0; i_h < int(trk_hits_ptrs.size()); ++i_h) {

		float dQinHit = trk_hits_ptrs[i_h]->Integral();
		dQinAllHits += dQinHit;

		particle_vec.clear(); match_vec.clear();
		particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);

		// to avoid from matching the same particle more than once, we introduce a vector of matched TrackId-s
		// and require that the matched particle has not been mathced already for this hit
		std::vector <int> ParticlesMatchedInThisHit;

		//loop over particles
		for(int i_p = 0; i_p < int(particle_vec.size()); ++i_p) {

			float Edep_particle = match_vec[i_p]->energy;  // energy deposited by ionization by this track ID [MeV]

			trkide[ particle_vec[i_p]->TrackId() ] += Edep_particle; //store energy [MeV] deposited by track id
			trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ] += dQinHit; // store the integral on the hit by the track id

			tote += Edep_particle; //calculate total energy deposited

			if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum

				maxe = trkide[ particle_vec[i_p]->TrackId() ];
				maxp_me = particle_vec[i_p];

				if (!ParticleAlreadyMatchedInThisHit( ParticlesMatchedInThisHit , (int)particle_vec[i_p]->TrackId()) ) {

					ParticlesMatchedInThisHit.push_back( (int)particle_vec[i_p]->TrackId() );
					trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ] += dQinHit; // store the integral on the hit by the track id
					max_dQinTruthMatchedHits = trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ];

				} else {

					std::cout << "particle of TrackID "+ToString((int)particle_vec[i_p]->TrackId())+" was already matched in this hit" << std::endl;
				}
		
			} // End of the if-statement

		} // End of the loop over particles per hit

	} // End of the loop over the hits

	double kMin_dQ_inTruthMatchedHits = 0.1;
	if ( (max_dQinTruthMatchedHits / dQinAllHits) < kMin_dQ_inTruthMatchedHits ) { 

		art::Ptr< simb::MCParticle > maxp_me_null;
		maxp_me = maxp_me_null;
	}

	return maxp_me;
}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

bool ParticleAlreadyMatchedInThisHit(std::vector<int> AlreadyMatched_TrackIDs ,int cTrackID )
{
	// to avoid from matching the same particle more than once 
	// we introduce a vector of matched TrackId-s for each hit 
	// and require that the matched particle has not been mathced already for this hit  

	for (auto trk_id:AlreadyMatched_TrackIDs) { 

		if (trk_id == cTrackID) return true; 
	}

	return false; 
}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

TString ToString(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;
}
