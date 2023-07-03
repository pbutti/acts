#pragma once


#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/Index.hpp"

#include "fastjet/ClusterSequence.hh"

//#include "TLorentzVector.h"

#include <memory>
#include <vector>

using namespace Acts::UnitLiterals;

namespace ActsExamples {

class Lepton {
    
 public:
  
  Lepton(const Acts::Vector3& momentum,
         int idx,
         int pdgID) {
    
    m_mom  = momentum;
    m_pdgId = pdgID;
    m_idx   = idx;
  }

  
  double p() {return m_mom.norm();};
  double pt() {return std::sqrt(m_mom(0)*m_mom(0) + m_mom(1)*m_mom(1));};
  Acts::Vector3 momentum() {return m_mom;};
  
  double getIsolation(){return m_iso;};
  double getPdgId(){return m_pdgId;}

  void setIndex(int idx){m_idx = idx;};
  int  getIndex(){return m_idx;};

  double E() {
    double m = std::abs(m_pdgId) == 11 ? 0.0005 * 1_GeV : 0.1 * 1_GeV;
    return sqrt(m_mom.squaredNorm() + m*m);
  }
  
  void setIsolation(const TrackParametersContainer& tracks,
                    double deltaR = 0.2) {
    
    double iso = 0.;
        
    fastjet::PseudoJet lepton(m_mom(0),m_mom(1),m_mom(2),E());

    //std::cout<<"Running on tracks for isolation"<<std::endl;
    
    for (size_t itrk = 0; itrk < tracks.size(); itrk++) {

      Acts::Vector3 trackp = tracks[itrk].momentum();
      double  trackm = 0.1 * 1_GeV;
      double  E = sqrt(tracks[itrk].absoluteMomentum()*tracks[itrk].absoluteMomentum() + trackm*trackm);

      fastjet::PseudoJet trackj(trackp(0),trackp(1),trackp(2),E);

      if (trackj.delta_R(lepton) < deltaR && itrk != m_idx) {
        //std::cout<<"Dr(Track,Lepton) = "<<trackj.delta_R(lepton)<<std::endl;
        iso += trackp.norm();
      }
    }

    m_iso = iso / m_mom.norm();
    //std::cout<<"iso = "<<iso<<" m_iso="<<m_iso<<std::endl; 
    
  }
  
  
 private:
  
  double m_iso{-1};
  double m_pdgId{0};
  double m_idx{-1};
  Acts::Vector3 m_mom;
  
  };








class TrackJetsAlgorithm final : public IAlgorithm {
 public:

  using HitParticleMap = IndexMultimap<ActsFatras::Barcode>;

  struct Config {
    /// Input track collection
    std::string inputTrackCollection;
    /// Input trajectories
    std::string inputTrajectories;
    
    /// The simulated particles (for truth studies)
    std::string simParticles;

    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;

    /// Input reconstructed vertices
    std::string recoVertices;

    /// Minimum fraction of hits associated to particle to consider
    /// as truth matched.
    double truthMatchProbability = 0.5;
    
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

  
  TrackJetContainer OverlapRemoval (TrackJetContainer& jets,
                                    const TrackParametersContainer& tracks,
                                    const std::vector<SimParticle>& associatedTruthParticles) const;
  
  

 private:
  Config m_cfg;
  std::shared_ptr<fastjet::JetDefinition> m_jet_def;
  
  ReadDataHandle<TrajectoriesContainer> m_inputTrajectories{this,"InputTrajectories"};
  ReadDataHandle<SimParticleContainer> m_simulatedParticles{this,"SimParticles"};
  ReadDataHandle<std::vector<Acts::Vertex<Acts::BoundTrackParameters>>> m_recoVertices{this,"fittedvertices"};
  ReadDataHandle<HitParticleMap> m_inputMeasurementParticlesMap{
    this, "InputMeasurementParticlesMap"};
  
  WriteDataHandle<TrackJetContainer> m_outputTrackJets{this, "outputTrackJets"};
    
};
  
} // namespace ActsExamples
