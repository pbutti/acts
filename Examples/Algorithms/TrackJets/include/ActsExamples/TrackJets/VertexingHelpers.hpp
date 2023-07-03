#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"

#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"


#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <map>

namespace ActsExamples { 

using HitParticleMap = IndexMultimap<ActsFatras::Barcode>;

// Match fitted tracks to truth tracks

void matchTracks(const TrajectoriesContainer& inputTrajectories,
                 const HitParticleMap& hitParticlesMap,
                 const SimParticleContainer& allTruthParticles,
                 double truthMatchProbMin,
                 std::vector<Acts::BoundTrackParameters>& trackParameters,
                 std::vector<SimParticle>& associatedTruthParticles) {
                 
  
  std::vector<ParticleHitCount> particleHitCounts;
  for (const auto& trajectories : inputTrajectories) {
    for (auto tip : trajectories.tips()) {
      if (!trajectories.hasTrackParameters(tip)) {
        continue;
      }
      
      identifyContributingParticles(hitParticlesMap, trajectories, tip,
                                    particleHitCounts);
      
      ActsFatras::Barcode majorityParticleId =
          particleHitCounts.front().particleId;
      
      size_t nMajorityHits = particleHitCounts.front().hitCount;
      
      auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
          trajectories.multiTrajectory(), tip);
      
      if (nMajorityHits * 1. / trajState.nMeasurements <
          truthMatchProbMin) {
        continue;
      }
      
      auto it = std::find_if(allTruthParticles.begin(),
                             allTruthParticles.end(), [&](const auto& tp) {
                               return tp.particleId() == majorityParticleId;
                             }); 
      
      if (it == allTruthParticles.end()) {
        continue;
      }
      
      const auto& majorityParticle = *it;
      trackParameters.push_back(trajectories.trackParameters(tip));
      associatedTruthParticles.push_back(majorityParticle);
    }
  }
}

// Get Vertex ID and matching fraction

std::pair<int,double> getVtxPrimaryID(const Acts::Vertex<Acts::BoundTrackParameters>& vtx,
                                      const std::vector<Acts::BoundTrackParameters>& trackParameters,
                                      const std::vector<SimParticle>& associatedTruthParticles) {
  
  
  const auto tracks =  vtx.tracks();
  
  SimParticleContainer particleAtVtx;
  std::vector<int> contributingTruthVertices;
  
  for (const auto& trk : tracks) {
    Acts::BoundTrackParameters origTrack = *(trk.originalParams);
    
    // Find associated truth particle now
    // We expect that the vectors `trackParameters` and
    // `associatedTruthParticles` align

    for (std::size_t i = 0; i < trackParameters.size(); ++i) {
      const auto& particle = associatedTruthParticles[i];
      const auto& trackParameter = trackParameters[i].parameters();
      
      if (origTrack.parameters() == trackParameter) {
        int priVtxId = particle.particleId().vertexPrimary();
        
        particleAtVtx.insert(particle);
        contributingTruthVertices.push_back(priVtxId);
      }
    }
  }  // end loop tracks
  
  // Now find true vtx with most matching tracks at reco vtx
  // and check if it contributes more than 50% of all tracks
  std::map<int, int> fmap;
  for (int priVtxId : contributingTruthVertices) {
    fmap[priVtxId]++;
  }
  int maxOccurrence = -1;
  int maxOccurrenceId = -1;
  for (auto it : fmap) {
    if (it.second > maxOccurrence) {
      maxOccurrence = it.second;
      maxOccurrenceId = it.first;
    }
  }

  // Count number of reconstructible tracks on truth vertex
  int nTracksOnTruthVertex = 0;
  for (const auto& particle : associatedTruthParticles) {
    int priVtxId = particle.particleId().vertexPrimary();
    if (priVtxId == maxOccurrenceId) {
      ++nTracksOnTruthVertex;
    }
  }
  
  // Return the max occurrence Id and the track-to-Vtx match fraction
  double trackVtxMatchFraction =
      (double)fmap[maxOccurrenceId] / tracks.size();
  
  return std::pair<int,double>(maxOccurrenceId, trackVtxMatchFraction);
}


double vtxSumPt2(const Acts::Vertex<Acts::BoundTrackParameters>& vtx) {

  double spt2 = 0.;
  for (auto& trk : vtx.tracks()) {
    double trkpt = (*(trk.originalParams)).absoluteMomentum();
    spt2 += trkpt*trkpt;
  }

  return spt2;
}

}

