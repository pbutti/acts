#include "ActsExamples/TrackJets/TrackJetsAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>


ActsExamples::TrackJetsAlgorithm::TrackJetsAlgorithm(
    Config config, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TrackJetsAlgorithm", level),
      m_cfg(std::move(config)) {

  if (m_cfg.inputTrackCollection.empty()) {
    throw std::invalid_argument("Missing input track collection");
  }
  if (m_cfg.outputTrackJets.empty()) {
    throw std::invalid_argument("Missing output track jets collection");
  }

}

ActsExamples::ProcessCode ActsExamples::TrackJetsAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  //Read input data

  //For the moment use the 
  const auto& tracks =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackCollection);
  
  ACTS_INFO("Running Fast Jet on tracks");


  std::vector<fastjet::PseudoJet> input_particles;
  input_particles.push_back(fastjet::PseudoJet(0.1,0.2,0.3,0.4));
  input_particles.push_back(fastjet::PseudoJet(1,2,3,4));
  input_particles.push_back(fastjet::PseudoJet(0.2,0.3,0.4,0.5)); 
    
  //Put it in initialize
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, m_cfg.radius);

  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);


  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  double ptmin = 1.0;
  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

  ACTS_INFO("Printing track jets details");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp());
  }

  
  //ctx.eventStore.add(m_cfg.outputTrackJets, std::move(trackJets));
  

  return ActsExamples::ProcessCode::SUCCESS;
  
}
