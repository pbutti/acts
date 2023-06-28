#pragma once


#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "Acts/Definitions/Units.hpp"

#include "fastjet/ClusterSequence.hh"

//#include "TLorentzVector.h"

#include <memory>
#include <vector>

using namespace Acts::UnitLiterals;

namespace ActsExamples {

class TrackJetsAlgorithm final : public IAlgorithm {
 public:

  struct Config {
    /// Input track collection
    std::string inputTrackCollection;
    /// Input trajectories
    std::string inputTrajectories;
    
    /// The simulated particles (for truth studies only, empty otherwise)
    std::string simParticles;
    
    /// Jet cone radius
    double radius = 0.4;
    /// Track mass hypothesis [GeV]
    double trackMass = 0.135 * 1_GeV;


    /// Min TrackJet pT
    double tj_minPt = 10;
    /// Output track jet collection
    std::string outputTrackJets;
  };

  TrackJetsAlgorithm(const Config& config, Acts::Logging::Level level);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;
  
  
  /// Get readonly access to the config parameters
  const Config& config() const {return m_cfg;}
  
  void findJetConstituents(const fastjet::PseudoJet& jet,
                           const SimParticleContainer& simparticles,
                           std::vector<SimParticle>& jetconstituents,
                           float dR = 0.4) const ;

  ActsExamples::jetlabel ClassifyJet(const fastjet::PseudoJet& jet,
                                     const SimParticleContainer& simparticles,
                                     float dR) const;
    
  ActsExamples::hadronlabel defTypeOfHadron(int pdg) const;

  void associateTracksToJets(const TrackParametersContainer& tracks,
                             TrackJetContainer& jetContainer,
                             double DRmax) const ;

 private:
  Config m_cfg;
  std::shared_ptr<fastjet::JetDefinition> m_jet_def;

    ReadDataHandle<TrajectoriesContainer> m_inputTrajectories{this,"InputTrajectories"};
    ReadDataHandle<SimParticleContainer> m_simulatedParticles{this,"SimParticles"};
    
    WriteDataHandle<TrackJetContainer> m_outputTrackJets{this, "outputTrackJets"};
    
};
  
} // namespace ActsExamples
