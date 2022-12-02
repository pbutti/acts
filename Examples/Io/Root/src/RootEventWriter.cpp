// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootEventWriter.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;
using Acts::VectorHelpers::perp;

ActsExamples::RootEventWriter::RootEventWriter(
    const ActsExamples::RootEventWriter::Config& config,
    Acts::Logging::Level level) 
    : TrackJetWriter(config.inputJets,
                           "RootEventWriter", level),
      m_cfg(config) {
  /*
  if (m_cfg.inputJets.empty()) {
    throw std::invalid_argument("Missing jets input collection");
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-simulated-hits map input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  */

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = m_cfg.filePath;
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {

    m_outputTree->Branch("event_nr", &m_eventNr);
    
    // The jets
    m_outputTree->Branch("jet_pt",&m_jet_pt);
    m_outputTree->Branch("jet_eta",&m_jet_eta);
    m_outputTree->Branch("jet_phi",&m_jet_phi);
    m_outputTree->Branch("jet_ncomponents",&m_jet_ncomponents);
    m_outputTree->Branch("jet_tracks_idx",&m_jet_components);
    m_outputTree->Branch("jet_isPU",&m_jet_isPU);
    m_outputTree->Branch("jet_isHS",&m_jet_isHS);
    m_outputTree->Branch("jet_label",&m_jet_label);
    
    
    //Tracks in jets
    m_outputTree->Branch("track_prob",       &m_tracks_prob);
    m_outputTree->Branch("track_d0",         &m_trk_d0);
    m_outputTree->Branch("track_z0",         &m_trk_z0);
    m_outputTree->Branch("track_eta",        &m_trk_eta);
    m_outputTree->Branch("track_theta",      &m_trk_theta);
    m_outputTree->Branch("track_phi",        &m_trk_phi);
    m_outputTree->Branch("track_pt",         &m_trk_pt);
    m_outputTree->Branch("track_qOverP",     &m_trk_qOverP);
    m_outputTree->Branch("track_t",          &m_trk_t);
    m_outputTree->Branch("track_t30",        &m_trk_t30);
    m_outputTree->Branch("track_t60",        &m_trk_t60);
    m_outputTree->Branch("track_t90",        &m_trk_t90);
    m_outputTree->Branch("track_t120",       &m_trk_t120);
    m_outputTree->Branch("track_t180",       &m_trk_t180);
    m_outputTree->Branch("track_z",          &m_trk_z);

    m_outputTree->Branch("track_var_d0",     &m_trk_var_d0);
    m_outputTree->Branch("track_var_z0",     &m_trk_var_z0);
    m_outputTree->Branch("track_var_phi",    &m_trk_var_phi);
    m_outputTree->Branch("track_var_theta",  &m_trk_var_theta);
    m_outputTree->Branch("track_var_qOverP", &m_trk_var_qOverP);

    m_outputTree->Branch("track_cov_d0z0"       ,&m_trk_cov_d0z0);
    m_outputTree->Branch("track_cov_d0phi"      ,&m_trk_cov_d0phi);
    m_outputTree->Branch("track_cov_d0theta"    ,&m_trk_cov_d0theta);
    m_outputTree->Branch("track_cov_d0qOverP"   ,&m_trk_cov_d0qOverP);
    m_outputTree->Branch("track_cov_z0phi"      ,&m_trk_cov_z0phi);
    m_outputTree->Branch("track_cov_z0theta"    ,&m_trk_cov_z0theta);
    m_outputTree->Branch("track_cov_z0qOverP"   ,&m_trk_cov_z0qOverP);
    m_outputTree->Branch("track_cov_phitheta"   ,&m_trk_cov_phitheta);
    m_outputTree->Branch("track_cov_phiqOverP"  ,&m_trk_cov_phiqOverP);
    m_outputTree->Branch("track_cov_tehtaqOverP",&m_trk_cov_tehtaqOverP);

    m_outputTree->Branch("track_numPix1L" ,&m_trk_numPix1L);    
    m_outputTree->Branch("track_numPix2L" ,&m_trk_numPix2L);   
    m_outputTree->Branch("track_numPix"   ,&m_trk_numPix);   
    m_outputTree->Branch("track_numSCT"   ,&m_trk_numSCT);  
    m_outputTree->Branch("track_numLSCT"  ,&m_trk_numLSCT);

    
    
    // The estimated track parameters
    /*
    
    m_outputTree->Branch("loc0", &m_loc0);
    m_outputTree->Branch("loc1", &m_loc1);
    m_outputTree->Branch("phi", &m_phi);
    m_outputTree->Branch("theta", &m_theta);
    m_outputTree->Branch("qop", &m_qop);
    m_outputTree->Branch("time", &m_time);
    m_outputTree->Branch("p", &m_p);
    m_outputTree->Branch("pt", &m_pt);
    m_outputTree->Branch("eta", &m_eta);
    */
   
    // The truth track parameters
    /*
    m_outputTree->Branch("eventNr", &m_eventNr);
    m_outputTree->Branch("t_loc0", &m_t_loc0);
    m_outputTree->Branch("t_loc1", &m_t_loc1);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_qop", &m_t_qop);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("truthMatched", &m_truthMatched);
    */
  }
}

ActsExamples::RootEventWriter::~RootEventWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootEventWriter::endRun() {
  if (m_outputFile != nullptr) {
    m_outputFile->cd();
    m_outputTree->Write();
    //ACTS_INFO("Write estimated parameters from seed to tree '"
    //          << m_cfg.treeName << "' in '" << m_cfg.filePath << "'");
    m_outputFile->Close();
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootEventWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const TrackJetContainer& tp) {
  /*
    using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
    using HitSimHitsMap = IndexMultimap<Index>;
  */

  //Make sure everytghing is clear
  Clear();
  
  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }
  
  // Read input collections
  
  const auto& jets =
      ctx.eventStore.get<TrackJetContainer>(m_cfg.inputJets);
  
  const auto& tracks =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  
  //const auto& particles =
  //    ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  
  
  /*
  const auto& protoTracks =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
  

  
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  const auto& hitSimHitsMap =
      ctx.eventStore.get<HitSimHitsMap>(m_cfg.inputMeasurementSimHitsMap);
  */
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  /*
  // Loop over the estimated track parameters
  for (size_t iparams = 0; iparams < trackParams.size(); ++iparams) {
    // The reference surface of the parameters, i.e. also the reference surface
    // of the first space point
    const auto& surface = trackParams[iparams].referenceSurface();
    // The estimated bound parameters vector
    const auto params = trackParams[iparams].parameters();
    m_loc0 = params[Acts::eBoundLoc0];
    m_loc1 = params[Acts::eBoundLoc1];
    m_phi = params[Acts::eBoundPhi];
    m_theta = params[Acts::eBoundTheta];
    m_qop = params[Acts::eBoundQOverP];
    m_time = params[Acts::eBoundTime];
    m_p = std::abs(1.0 / m_qop);
    m_pt = m_p * std::sin(m_theta);
    m_eta = std::atanh(std::cos(m_theta));

    // Get the proto track from which the track parameters are estimated
    const auto& ptrack = protoTracks[iparams];
    std::vector<ParticleHitCount> particleHitCounts;
    identifyContributingParticles(hitParticlesMap, ptrack, particleHitCounts);
    m_truthMatched = false;
    if (particleHitCounts.size() == 1) {
      m_truthMatched = true;
    }
    // Get the index of the first space point
    const auto& hitIdx = ptrack.front();
    // Get the sim hits via the measurement to sim hits map
    auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
    auto [truthLocal, truthPos4, truthUnitDir] =
        averageSimHits(ctx.geoContext, surface, simHits, indices);
    // Get the truth track parameter at the first space point
    m_t_loc0 = truthLocal[Acts::ePos0];
    m_t_loc1 = truthLocal[Acts::ePos1];
    m_t_phi = phi(truthUnitDir);
    m_t_theta = theta(truthUnitDir);
    m_t_time = truthPos4[Acts::eTime];
    // momemtum averaging makes even less sense than averaging position and
    // direction. use the first momentum or set q/p to zero
    if (not indices.empty()) {
      // we assume that the indices are within valid ranges so we do not
      // need to check their validity again.
      const auto simHitIdx0 = indices.begin()->second;
      const auto& simHit0 = *simHits.nth(simHitIdx0);
      const auto p =
          simHit0.momentum4Before().template segment<3>(Acts::eMom0).norm();
      const auto& particleId = simHit0.particleId();
      // The truth charge has to be retrieved from the sim particle
      auto ip = particles.find(particleId);
      if (ip != particles.end()) {
        const auto& particle = *ip;
        m_t_charge = particle.charge();
        m_t_qop = m_t_charge / p;
      } else {
        ACTS_WARNING("Truth particle with barcode = " << particleId
                                                      << " not found!");
      }
    }
    
    
    m_outputTree->Fill();
    } */

    
  for (size_t ijets = 0; ijets < jets.size(); ++ijets) {
    Acts::Vector4 jet_4mom = jets[ijets].getFourMomentum();

    float jet_theta = theta(jet_4mom);
    m_jet_pt.push_back(perp(jet_4mom));
    //m_jet_eta.push_back(jets[ijets].eta());
    m_jet_eta.push_back(std::atanh(std::cos(jet_theta)));
    m_jet_phi.push_back(phi(jet_4mom));

    m_jet_ncomponents.push_back(jets[ijets].getConstituents().size());
    m_jet_components.push_back(jets[ijets].getConstituents());
    
    m_jet_label.push_back(static_cast<int>(jets[ijets].getLabel()));

    //Get the tracks belonging to this jet
    for (auto& ic : jets[ijets].getTracks()) { 
      
      TrackParameters trk_params = tracks.at(ic);
      const auto params = trk_params.parameters();

      float trk_theta = params[Acts::eBoundTheta];
      float trk_eta   = std::atanh(std::cos(trk_theta));
      float trk_qop   = params[Acts::eBoundQOverP];
      float trk_p     = std::abs(1.0 / trk_qop);
      float trk_pt    = trk_p * std::sin(trk_theta); 
      
      m_tracks_prob.push_back(1.);   // todo
      m_trk_d0.push_back(params[Acts::eBoundLoc0]);             
      m_trk_z0.push_back(params[Acts::eBoundLoc1]);             
      m_trk_eta.push_back(trk_eta);            
      m_trk_theta.push_back(trk_theta);          
      m_trk_phi.push_back(params[Acts::eBoundPhi]);            
      m_trk_pt.push_back(trk_pt);             
      m_trk_qOverP.push_back(trk_p);         
      m_trk_t.push_back(params[Acts::eBoundTime]);              
      m_trk_t30.push_back(1.);            //todo
      m_trk_t60.push_back(1.);            //todo
      m_trk_t90.push_back(1.);            //todo
      m_trk_t120.push_back(1.);           //todo
      m_trk_t180.push_back(1.);           //todo
      m_trk_z.push_back(1.);              //todo
      m_trk_var_d0.push_back(1.);         
      m_trk_var_z0.push_back(1.);         
      m_trk_var_phi.push_back(1.);        
      m_trk_var_theta.push_back(1.);      
      m_trk_var_qOverP.push_back(1.);     
      m_trk_cov_d0z0.push_back(1.);       
      m_trk_cov_d0phi.push_back(1.);      
      m_trk_cov_d0theta.push_back(1.);    
      m_trk_cov_d0qOverP.push_back(1.);   
      m_trk_cov_z0phi.push_back(1.);      
      m_trk_cov_z0theta.push_back(1.);    
      m_trk_cov_z0qOverP.push_back(1.);   
      m_trk_cov_phitheta.push_back(1.);   
      m_trk_cov_phiqOverP.push_back(1.);  
      m_trk_cov_tehtaqOverP.push_back(1.);
            
    } //getTracks
  } // jets
  
  m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}



void ActsExamples::RootEventWriter::Clear() {

  //Jets
  m_jet_pt.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_ncomponents.clear();
  m_jet_components.clear();
  m_jet_tracks_idx.clear();
  m_jet_isPU.clear();
  m_jet_isHS.clear();
  m_jet_label.clear();

  //Tracks
  m_tracks_prob.clear();        
  m_trk_d0.clear();             
  m_trk_z0.clear();             
  m_trk_eta.clear();            
  m_trk_theta.clear();          
  m_trk_phi.clear();            
  m_trk_pt.clear();
  m_trk_qOverP.clear();             
  m_trk_t.clear();              
  m_trk_t30.clear();            
  m_trk_t60.clear();            
  m_trk_t90.clear();            
  m_trk_t120.clear();           
  m_trk_t180.clear();           
  m_trk_z.clear();              
  m_trk_var_d0.clear();         
  m_trk_var_z0.clear();         
  m_trk_var_phi.clear();        
  m_trk_var_theta.clear();      
  m_trk_var_qOverP.clear();     
  m_trk_cov_d0z0.clear();       
  m_trk_cov_d0phi.clear();      
  m_trk_cov_d0theta.clear();    
  m_trk_cov_d0qOverP.clear();   
  m_trk_cov_z0phi.clear();      
  m_trk_cov_z0theta.clear();    
  m_trk_cov_z0qOverP.clear();   
  m_trk_cov_phitheta.clear();   
  m_trk_cov_phiqOverP.clear();  
  m_trk_cov_tehtaqOverP.clear();
  
}
