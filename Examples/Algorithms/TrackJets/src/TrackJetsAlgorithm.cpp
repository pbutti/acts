#include "ActsExamples/TrackJets/TrackJetsAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>


using namespace Acts;

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


  m_jet_def = std::make_shared<fastjet::JetDefinition>(fastjet::antikt_algorithm, m_cfg.radius);
  
}

ActsExamples::ProcessCode ActsExamples::TrackJetsAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {


  //Read input data
  const auto& tracks =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackCollection);

  //Empty set if simParticles is not set in configuration
  const auto& simulatedParticles =
      (!m_cfg.simParticles.empty()) ? ctx.eventStore.get<SimParticleContainer>(m_cfg.simParticles) :
      ::boost::container::flat_set<::ActsFatras::Particle,detail::CompareParticleId>();
  
  ACTS_INFO("Running Fast Jet on tracks");
  

  std::vector<fastjet::PseudoJet> input_particles;

  int index = 0;

  ACTS_INFO("Number of tracks " << tracks.size());
  
  for (auto& track : tracks) {

    //std::cout<<track.momentum()[0]+" "+track.momentum()[1]+" "+track.momentum()[2]<<std::endl;
    
    //Get 4-momentum with track hypothesis
    Vector3 trackp = track.momentum();
    double  trackm = m_cfg.trackMass * 1_GeV;
    double  E = sqrt(track.absoluteMomentum()*track.absoluteMomentum() + trackm*trackm);
    
    fastjet::PseudoJet pj(trackp[0],trackp[1],trackp[2],E);
    pj.set_user_index(index);
    input_particles.push_back(pj);
    index++;
  }
  
  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fastjet::ClusterSequence clust_seq(input_particles, *m_jet_def);
  
  
  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(m_cfg.tj_minPt));
  
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
