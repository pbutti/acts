#include "ActsExamples/TrackJets/TrackJetsAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "MCUtils/PIDUtils.h"
#include <stdexcept>

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::theta;

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
  
  ACTS_DEBUG("Number of tracks " << tracks.size());
    
  //Empty set if simParticles is not set in configuration
  const auto& simulatedParticles =
      (!m_cfg.simParticles.empty()) ? ctx.eventStore.get<SimParticleContainer>(m_cfg.simParticles) :
      ::boost::container::flat_set<::ActsFatras::Particle,detail::CompareParticleId>();


  //Form jets
  std::vector<fastjet::PseudoJet> input_particles;
  int part_idx = -1;
  for (const auto& particle : simulatedParticles) {

    part_idx++;

    if (!(particle.status() == 1))
      continue;
    if (!particle.isFinal())
      continue;
    if (!particle.isVisible())
      continue;
    if (!particle.isStable())
      continue;

    double particle_pt = std::sqrt(particle.fourMomentum()(0)*particle.fourMomentum()(0) +
                                   particle.fourMomentum()(1)*particle.fourMomentum()(1) );
    
    if (particle_pt < 0.5 * 1_GeV)
      continue;
    
    ACTS_DEBUG("SimParticle Passed to JET building: " << particle.pdg());
    
    //Get 4-momentum of the selected particles
    Acts::Vector4 p4mom = particle.fourMomentum();
    
    fastjet::PseudoJet pj(p4mom(0),p4mom(1),p4mom(2),p4mom(3));
    pj.set_user_index(part_idx);
    input_particles.push_back(pj);
  }
  
  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fastjet::ClusterSequence clust_seq(input_particles, *m_jet_def);
  
  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(m_cfg.tj_minPt));
  
  TrackJetContainer outputTrackJets;
  outputTrackJets.reserve(inclusive_jets.size());
  
  //Prepare each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    
    // get the constituents of the jet
    std::vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
    
    ACTS_DEBUG(i<<" "<<inclusive_jets[i].rap()<<
               " "<<inclusive_jets[i].phi()<<
               " "<<inclusive_jets[i].perp()<<
               " "<<(unsigned int) constituents.size());
    
    fastjet::PseudoJet jet = inclusive_jets[i];
    
    ActsExamples::jetlabel jetType = ClassifyJet(jet,
                                                 simulatedParticles,
                                                 0.4);
    
    Acts::Vector4 pj_fourmom(jet.px(), jet.py(), jet.pz(), jet.E());
    TrackJet tj(pj_fourmom, jetType);
    
    std::vector<int> cons_idxs;
    cons_idxs.reserve(constituents.size());
    
    //Add the particles to the jets
    for (unsigned int j=0; j<constituents.size(); j++) {
      cons_idxs.push_back(constituents[j].user_index());
    }
    
    tj.setConstituents(cons_idxs);
    outputTrackJets.push_back(tj);
  }

  //Associate the tracks to the jets
  associateTracksToJets(tracks, outputTrackJets, 0.4);
  
  ctx.eventStore.add(m_cfg.outputTrackJets, std::move(outputTrackJets));
  
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

/* This method loops over the tracks in the event and try to associate those to jets
   - DeltaRmax < 0.4
   - Associate the track to closest jet. Add the index of the track to the closest jet trk_idxs
*/

void ActsExamples::TrackJetsAlgorithm::associateTracksToJets(const ActsExamples::TrackParametersContainer& tracks,
                                                             ActsExamples::TrackJetContainer& jetContainer,
                                                             double DRmax) const {
  if (tracks.size() < 1)
    return;

  for (size_t itrk=0; itrk<tracks.size(); itrk++) {

    double drMin = DRmax;
    int j_idx = -1;
    
    //Form 4-vectors for DR computation.
    //Mass hypothesis from cfg.
    //TODO remove double call to pseudojet..

    Vector3 trackp = tracks[itrk].momentum();
    double  trackm = m_cfg.trackMass * 1_GeV;
    double  E = sqrt(tracks[itrk].absoluteMomentum()*tracks[itrk].absoluteMomentum() + trackm*trackm);
    
    fastjet::PseudoJet trackj(trackp(0),trackp(1),trackp(2),E);

    //Find the closest jet within dRmax
    //Useless to form the pseudojet every time... TODO FIX 
    
    for (size_t ijet=0; ijet<jetContainer.size(); ijet++) {

      const auto jet = jetContainer[ijet];
      fastjet::PseudoJet jj(jet.getFourMomentum()(0),
                            jet.getFourMomentum()(1),
                            jet.getFourMomentum()(2),
                            jet.getFourMomentum()(3));
      

      if (trackj.delta_R(jj) < drMin) {
        
        j_idx = ijet;
        drMin = trackj.delta_R(jj);
        
      }//dR
    }//loop on jets
    
    //If the matching jet has been found, add the track to the jet.

    ACTS_DEBUG("The found j_idx="<<j_idx);
    
    if (j_idx>=0)
      jetContainer[j_idx].addTrack(itrk);
    
  }//loop on tracks
}//associateTracksToJets

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

