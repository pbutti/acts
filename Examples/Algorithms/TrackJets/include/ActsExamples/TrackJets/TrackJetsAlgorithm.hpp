#pragma once


#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include "fastjet/ClusterSequence.hh"

#include <memory>
#include <vector>

namespace ActsExamples {

class TrackJetsAlgorithm final : public BareAlgorithm {
 public:

  struct Config {
    /// Input track collection
    std::string inputTrackCollection;
    /// Jet cone radius
    double radius = 0.4;
    /// Output track jet collection
    std::string outputTrackJets;
  };

  TrackJetsAlgorithm(Config config, Acts::Logging::Level level);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;


  /// Get readonly access to the config parameters
  const Config& config() const {return m_cfg;}


 private:
  Config m_cfg;
};
  
} // namespace ActsExamples
