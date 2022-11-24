#include "ActsExamples/TrackJets/TrackJetsAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "MCUtils/PIDUtils.h"
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

  std::vector<fastjet::PseudoJet> input_particles;
  
  int index = 0;
  
  ACTS_DEBUG("Number of tracks " << tracks.size());
  
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
  
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    // get the constituents of the jet
    std::vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
    
    //printf("%5u %15.8f %15.8f %15.8f %8u\n",
    //	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
    //	   inclusive_jets[i].perp(), (unsigned int) constituents.size());
    

    fastjet::PseudoJet jet = inclusive_jets[i];
    
    ActsExamples::jetlabel jetType = ClassifyJet(jet,
                                                 simulatedParticles,
                                                 0.4);
    
  }
  
  
  //ctx.eventStore.add(m_cfg.outputTrackJets, std::move(trackJets));
  
  return ActsExamples::ProcessCode::SUCCESS;
  
}

//Taken from https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/PhysicsAnalysis/MCTruthClassifier/Root/MCRecoToTruth.cxx

ActsExamples::jetlabel ActsExamples::TrackJetsAlgorithm::ClassifyJet(const fastjet::PseudoJet& jet,
                                                                     const SimParticleContainer& simparticles,
                                                                     float dR) const {
  //ActsExamples::jetlabel jettype            = ActsExamples::unknown;
  ActsExamples::jetlabel tempjettype        = ActsExamples::unknown;
  ActsExamples::hadronlabel temphadronlabel = ActsExamples::Unknown;

  std::vector<SimParticle> simconstituents;
  simconstituents.clear();
  
  findJetConstituents(jet,
                      simparticles,
                      simconstituents,
                      dR);
  
  ACTS_DEBUG("Found "<<simconstituents.size()<<" constituents");

  
  std::vector<SimParticle>::iterator it;
  for (it = simconstituents.begin(); it != simconstituents.end(); ++it) {

    ACTS_DEBUG("Checking... " << (*it).pdg());
    
    //Skip the leptons
    if (MCUtils::PID::isLepton((*it).pdg()))
      continue;

    //Find an hadron
    if (MCUtils::PID::isHadron((*it).pdg())) {

      //Classify the hadron
      temphadronlabel = defTypeOfHadron((*it).pdg());
      
      //First look for a B-Hadron
      if (temphadronlabel == ActsExamples::BBbarMesonPart  ||
          temphadronlabel == ActsExamples::BottomMesonPart ||
          temphadronlabel == ActsExamples::BottomBaryonPart) {
        tempjettype = ActsExamples::bjet;
      }
      
      //Then look for a C-Hadron. Cs can come from Bs, so keep a b-jet if has been found
      else if (temphadronlabel == ActsExamples::CCbarMesonPart  ||
               temphadronlabel == ActsExamples::CharmedMesonPart  ||
               temphadronlabel == ActsExamples::CharmedBaryonPart) {


        if (tempjettype == ActsExamples::bjet) {/* do nothing */}
        else {tempjettype = ActsExamples::cjet; } 
      }

      //Then look for a S-Hadron. Same logic.
      else if (temphadronlabel == ActsExamples::StrangeMesonPart  ||
               temphadronlabel == ActsExamples::StrangeBaryonPart ||
               temphadronlabel == ActsExamples::LightMesonPart    ||
               temphadronlabel == ActsExamples::LightBaryonPart)   {

        if (tempjettype == ActsExamples::bjet || tempjettype == ActsExamples::cjet) {/* do nothing */}
        else {tempjettype = ActsExamples::ljet;}
      }

      ACTS_DEBUG("Hadron is " << temphadronlabel <<" Jet is " << tempjettype);
      
      
    } //is a hadron
  } // loop on constituents


  return tempjettype;
  
} // classify jet

ActsExamples::hadronlabel ActsExamples::TrackJetsAlgorithm::defTypeOfHadron(int pdg) const {
  
  int q1 = (abs(pdg) / 1000) % 10;
  int q2 = (abs(pdg) / 100) % 10;
  int q3 = (abs(pdg) / 10) % 10;

  // di quark
  // if( q3 == 0 && q2 >=q3 )   cout<<"di quark"<<endl;
  // First two do not have obvious helpers in MCUtils
  if (q1 == 0 && MCUtils::PID::BQUARK == q2 && MCUtils::PID::BQUARK == q3)
    return ActsExamples::BBbarMesonPart;
  else if (q1 == 0 && MCUtils::PID::CQUARK == q3 && MCUtils::PID::CQUARK == q2)
    return ActsExamples::CCbarMesonPart;
  // Now just use the central helper functions
  else if (MCUtils::PID::isBottomMeson(pdg))
    return ActsExamples::BottomMesonPart;
  else if (MCUtils::PID::isCharmMeson(pdg))
    return ActsExamples::CharmedMesonPart;
  else if (MCUtils::PID::isBottomBaryon(pdg))
    return ActsExamples::BottomBaryonPart;
  else if (MCUtils::PID::isCharmBaryon(pdg))
    return ActsExamples::CharmedBaryonPart;
  else if (MCUtils::PID::isStrangeBaryon(pdg))
    return ActsExamples::StrangeBaryonPart;
  else if (MCUtils::PID::isLightBaryon(pdg))
    return ActsExamples::LightBaryonPart;
  else if (MCUtils::PID::isStrangeMeson(pdg))
    return ActsExamples::StrangeMesonPart;
  else if (MCUtils::PID::isLightMeson(pdg))
    return ActsExamples::LightMesonPart;
  else
    return ActsExamples::Unknown;
}


void ActsExamples::TrackJetsAlgorithm::findJetConstituents(const fastjet::PseudoJet& jet,
                                                           const SimParticleContainer& simparticles,
                                                           std::vector<SimParticle>& jetconstituents,
                                                           float dR) const {
  
  if (simparticles.size() < 1)
    return;
  
  for (auto& part : simparticles) {
    
    // use fast-jet to get the deltaR 
    
    fastjet::PseudoJet simpj(part.fourMomentum()(0),
                             part.fourMomentum()(1),
                             part.fourMomentum()(2),
                             part.fourMomentum()(3));
    
    
    if (simpj.delta_R(jet) < dR) {
      
      ACTS_DEBUG("SimParticle Component PdgID: " << part.pdg());
      jetconstituents.push_back(part);
      
    }//DeltaR
    
  }//sim particles
  
} //findJets

