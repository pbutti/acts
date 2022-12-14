// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TrackJet.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

//This I need to compute the ip2d and ip3d - to do: move this into an algorithm
#include "Acts/Vertexing/ImpactPointEstimator.hpp" 
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"

#include <mutex>
#include <memory>


class TFile;
class TTree;


using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
using PropagatorOptions = Acts::PropagatorOptions<>;
using ImpactPointEstimator =
      Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
using VertexContainer =
    std::vector<Acts::Vertex<Acts::BoundTrackParameters>>;

namespace ActsExamples {


//All the using.
using TrackJetWriter = WriterT<TrackJetContainer>;

/// Write out the Event using track parameters and track jets from both simulation and those estimated from
/// reconstructed seeds into a TTree
///
/// Each entry in the TTree corresponds to one seed for optimum writing
/// speed. The event number is part of the written data.
class RootEventWriter final : public TrackJetWriter {
 public:
  struct Config {
    /// Input jets
    std::string inputJets;
    /// Input estimated track parameters collection.
    std::string inputTrackParameters;
    // Input trajectories
    std::string inputTrajectories;
    // Input Reco Vertices
    std::string recoVertices;
    
    /// Input reconstructed proto tracks collection.
    std::string inputProtoTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// output filename.
    std::string filePath = "events.root";
    /// name of the output tree.
    std::string treeName = "events";
    /// file access mode.
    std::string fileMode = "RECREATE";
    /// the magnetic field provider
    std::shared_ptr<Acts::MagneticFieldProvider> field;
   
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootEventWriter(const Config& config,
                  Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootEventWriter() override;
  
  /// End-of-run hook
  ProcessCode endRun() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }
  
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrackJetContainer& tp) override;  //tp is not used

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number of


  void Clear();
  

  //Vertices//
  std::vector<float> m_vtx_x;
  std::vector<float> m_vtx_y;
  std::vector<float> m_vtx_z;
  std::vector<float> m_vtx_t;
  std::vector<float> m_vtx_sumPt2;
  std::vector<int> m_vtx_isHS;
  std::vector<int> m_vtx_isPU;
  
  
  //Jets//

  //skipping jet_m, jet_q
  
  std::vector<float> m_jet_pt, m_jet_eta, m_jet_phi;
  std::vector<int> m_jet_ncomponents;
  std::vector<std::vector<int>> m_jet_components;
  std::vector<std::vector<int>> m_jet_tracks_idx;
  std::vector<int> m_jet_label;
  std::vector<int> m_jet_isPU;
  std::vector<int> m_jet_isHS;

  
  //Tracks in jets
  std::vector<float> m_tracks_prob;
  std::vector<float> m_trk_d0;        
  std::vector<float> m_trk_z0;
  std::vector<float> m_trk_signed_d0;           
  std::vector<float> m_trk_signed_z0sinTheta;   
  std::vector<float> m_trk_signed_d0sig;        
  std::vector<float> m_trk_signed_z0sinThetasig;
  std::vector<float> m_trk_eta;       
  std::vector<float> m_trk_theta;     
  std::vector<float> m_trk_phi;       
  std::vector<float> m_trk_pt;
  std::vector<float> m_trk_qOverP;
  std::vector<float> m_trk_t;         
  std::vector<float> m_trk_t30;       
  std::vector<float> m_trk_t60;       
  std::vector<float> m_trk_t90;       
  std::vector<float> m_trk_t120;      
  std::vector<float> m_trk_t180;      
  std::vector<float> m_trk_z;         
  std::vector<float> m_trk_var_d0;    
  std::vector<float> m_trk_var_z0;    
  std::vector<float> m_trk_var_phi;   
  std::vector<float> m_trk_var_theta; 
  std::vector<float> m_trk_var_qOverP;
  
  std::vector<float> m_trk_cov_d0z0;       
  std::vector<float> m_trk_cov_d0phi;      
  std::vector<float> m_trk_cov_d0theta;    
  std::vector<float> m_trk_cov_d0qOverP;   
  std::vector<float> m_trk_cov_z0phi;      
  std::vector<float> m_trk_cov_z0theta;    
  std::vector<float> m_trk_cov_z0qOverP;   
  std::vector<float> m_trk_cov_phitheta;   
  std::vector<float> m_trk_cov_phiqOverP;  
  std::vector<float> m_trk_cov_thetaqOverP;
  
  std::vector<int> m_trk_numPix1L;    
  std::vector<int> m_trk_numPix2L;   
  std::vector<int> m_trk_numPix;   
  std::vector<int> m_trk_numSCT;  // short strips
  std::vector<int> m_trk_numLSCT; // long strips
  
  //the index of the jet the tracks belong to
  std::vector<int> m_trk_jet_idx;
  

  //Tools

  std::shared_ptr<Propagator> m_propagator;
  std::shared_ptr<ImpactPointEstimator> m_ipEst;
  
  Acts::GeometryContext gctx_;
  Acts::MagneticFieldContext mctx_;
  
};

}  // namespace ActsExamples
