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


    // The Truth HS Vertex
    
    // The reco vertices
    
    m_outputTree->Branch("recovertex_x",&m_vtx_x);
    m_outputTree->Branch("recovertex_y",&m_vtx_y);
    m_outputTree->Branch("recovertex_z",&m_vtx_z);
    m_outputTree->Branch("recovertex_t",&m_vtx_t);
    m_outputTree->Branch("recovertex_sumPt2",&m_vtx_sumPt2);
    m_outputTree->Branch("recovertex_isHS",&m_vtx_isHS);
    m_outputTree->Branch("recovertex_isPU",&m_vtx_isPU);

    /*
    m_outputTree->Branch("recovertex_cov_xx",&m_vtx_cov_xx);
    m_outputTree->Branch("recovertex_cov_xy",&m_vtx_cov_xy);
    m_outputTree->Branch("recovertex_cov_xz",&m_vtx_cov_xz);
    m_outputTree->Branch("recovertex_cov_xt",&m_vtx_cov_xt);
    m_outputTree->Branch("recovertex_cov_yy",&m_vtx_cov_yy);
    m_outputTree->Branch("recovertex_cov_yz",&m_vtx_cov_yz);
    m_outputTree->Branch("recovertex_cov_yt",&m_vtx_cov_yt);
    m_outputTree->Branch("recovertex_cov_zz",&m_vtx_cov_zz);
    m_outputTree->Branch("recovertex_cov_zt",&m_vtx_cov_zt);
    m_outputTree->Branch("recovertex_cov_tt",&m_vtx_cov_tt);
    */
    
    
    //m_outputTree->Branch("recovertex_tracks_idx",m_vtx_tracks_idx);
    //m_outputTree->Branch("recovertex_tracks_weight",m_vtx_tracks_weight);
    
    
    
    // The jets
    m_outputTree->Branch("jet_pt",&m_jet_pt);
    m_outputTree->Branch("jet_eta",&m_jet_eta);
    m_outputTree->Branch("jet_phi",&m_jet_phi);
    m_outputTree->Branch("jet_ncomponents",&m_jet_ncomponents);
    m_outputTree->Branch("jet_components",&m_jet_components);
    m_outputTree->Branch("jet_tracks_idx",&m_jet_tracks_idx);
    m_outputTree->Branch("jet_isPU",&m_jet_isPU);
    m_outputTree->Branch("jet_isHS",&m_jet_isHS);
    m_outputTree->Branch("jet_label",&m_jet_label);
    
    //Tracks in jets
    m_outputTree->Branch("track_prob",       &m_tracks_prob);
    m_outputTree->Branch("track_d0",         &m_trk_d0);
    m_outputTree->Branch("track_z0",         &m_trk_z0);
    m_outputTree->Branch("track_signedd0",             &m_trk_signed_d0);
    m_outputTree->Branch("track_signedd0sig",          &m_trk_signed_d0sig);
    m_outputTree->Branch("track_signedz0sinTheta",     &m_trk_signed_z0sinTheta);
    m_outputTree->Branch("track_signedz0sinThetasig",  &m_trk_signed_z0sinThetasig);
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
    m_outputTree->Branch("track_cov_tehtaqOverP",&m_trk_cov_thetaqOverP);

    // Number of Innermost Pixel Layer Hits
    m_outputTree->Branch("track_numPix1L" ,&m_trk_numPix1L);
    // Number of Next to Innermost Pixel Layer Hits
    m_outputTree->Branch("track_numPix2L" ,&m_trk_numPix2L);   
    m_outputTree->Branch("track_numPix"   ,&m_trk_numPix);   
    m_outputTree->Branch("track_numSCT"   ,&m_trk_numSCT);  
    m_outputTree->Branch("track_numLSCT"  ,&m_trk_numLSCT);

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

  //Setting up tools

  //Propagator with void navigator
  Acts::EigenStepper<> stepper(m_cfg.field);
  m_propagator = std::make_shared<Propagator>(stepper);
  
  //Setting up the ImpactPointEstimator
  ImpactPointEstimator::Config ipEstCfg(m_cfg.field, m_propagator);

  m_ipEst = std::make_shared<ImpactPointEstimator> (ipEstCfg);


  // Handle inputs

  m_inputJets.initialize(m_cfg.inputJets);
  m_inputTrajectories.initialize(m_cfg.inputTrajectories);
  m_recoVertices.initialize(m_cfg.recoVertices);
  
    
}

ActsExamples::RootEventWriter::~RootEventWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootEventWriter::finalize() {
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

  // Exclusive access to all the writing procedure
  std::lock_guard<std::mutex> lock(m_writeMutex);
  
  //Make sure everything is clear
  Clear();
  
  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }

  
  
  // Read input collections
  const auto& jets = m_inputJets(ctx);
  const auto& inputTrajectories = m_inputTrajectories(ctx);
  const auto& reco_vertices = m_recoVertices(ctx);
  
    
  ACTS_INFO("RootWriter::Number of "<<m_cfg.inputTrajectories << " "<< inputTrajectories.size());

  TrackParametersContainer tracks;
  std::vector<HitInfo> hitInfos;
  
  for (const auto& trajectories : inputTrajectories) {

    std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;

    //The multiTrajectory
    const auto&  mj = trajectories.multiTrajectory();
          
    //Loop over the trajectory entry index
    for (auto tip : trajectories.tips()) {
      if (!trajectories.hasTrackParameters(tip))
        continue;

      //Get the track parameters
      tracks.push_back(trajectories.trackParameters(tip));

      //Collect the trajectory summary info
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, tip);
      
      //Check the hits on track and add them to a vector in line with tracks
      //Kinda idiotic. Define a new track class to keep track parames + track info.
      
      hitInfos.push_back(GetHitInformation(trajState));
      
      
    }//tip loop
  }
  
  ACTS_INFO("RootWriter::Number of tracks" << tracks.size());

  if (tracks.size() != inputTrajectories.size())
    ACTS_WARNING("Valid trajectories size doesn't match expectation..");

  
  
  
  // Get the event number
  m_eventNr = ctx.eventNumber;


  // Get the vertices

  double maxSumPt2 = -1;
  unsigned int HS_idx = 0;

  //If there are 0 vertices in the event, skip the event
  if (reco_vertices.size() < 1)
    return ProcessCode::SUCCESS;
  
  for (size_t ivtx = 0; ivtx < reco_vertices.size(); ivtx++) {
    
    Acts::Vertex<Acts::BoundTrackParameters> recovtx = reco_vertices.at(ivtx);
    
    m_vtx_x.push_back(recovtx.position()[0]);
    m_vtx_y.push_back(recovtx.position()[1]);
    m_vtx_z.push_back(recovtx.position()[2]);
    m_vtx_t.push_back(recovtx.time());
    
    const auto vtx_cov = recovtx.fullCovariance();
    

    //Tracks at the vertex.
    //I will use the fittedParameters
    
    const auto vtxtracks = recovtx.tracks();
    
    double sumPt2 = 0.;
    for (auto& vtrk : vtxtracks) {
      sumPt2+=vtrk.fittedParams.transverseMomentum() * vtrk.fittedParams.transverseMomentum();
      //sumPt2+=vtrk.originalParams.transverseMomentum() * vtrk.originalParams.transverseMomentum();
    }
    
    //Store the vtx index with the maximum sumpt2
    if (sumPt2 > maxSumPt2) {
      maxSumPt2 = sumPt2;
      HS_idx = ivtx;
    }

    m_vtx_sumPt2.push_back(sumPt2);
        
  }

  //Unfortunately I need to make the loop on the vertices twice to assign the flag if the vertex is PU or HS
  
  for (size_t ivtx = 0; ivtx < reco_vertices.size(); ivtx++) {
    
    if (ivtx == HS_idx) {
      m_vtx_isHS.push_back(1);
      m_vtx_isPU.push_back(0);
    } else {
      m_vtx_isHS.push_back(0);
      m_vtx_isPU.push_back(1);
    }
    
  }
    
  for (size_t ijets = 0; ijets < jets.size(); ++ijets) {
    Acts::Vector4 jet_4mom = jets[ijets].getFourMomentum();
    Acts::Vector3 jet_3mom{jet_4mom[0],jet_4mom[1],jet_4mom[2]};
    float jet_theta = theta(jet_3mom);
    m_jet_pt.push_back(perp(jet_4mom));
    //m_jet_eta.push_back(jets[ijets].eta());
    m_jet_eta.push_back(std::atanh(std::cos(jet_theta)));
    m_jet_phi.push_back(phi(jet_4mom));

    m_jet_ncomponents.push_back(jets[ijets].getConstituents().size());
    m_jet_components.push_back(jets[ijets].getConstituents());
    m_jet_tracks_idx.push_back(jets[ijets].getTracks());
    
    m_jet_label.push_back(static_cast<int>(jets[ijets].getLabel()));
    
    m_jet_isPU.push_back(0);
    m_jet_isHS.push_back(1);
    
  } // jets
  

  //Loop on the tracks
  for (size_t itrk = 0; itrk<tracks.size(); itrk++) {

    double signed_d0             = -9999;
    double signed_z0SinTheta     = -9999;
    double signed_d0_err         = 1.;
    double signed_z0SinTheta_err = 1.;

    double covd0         = -999;
    double covz0         = -999;
    double covphi        = -999;
    double covtheta      = -999;
    double covqOverP     = -999;
    double covd0z0       = -999;
    double covd0phi      = -999;
    double covd0theta    = -999;
    double covd0qOverP   = -999;
    double covz0phi      = -999;
    double covz0theta    = -999;
    double covz0qOverP   = -999;
    double covphitheta   = -999;
    double covphiqOverP  = -999;
    double covthetaqOverP= -999;
        

    const auto trk_params = tracks[itrk];
    const auto hit_infos  = hitInfos[itrk];
    const auto params     = trk_params.parameters();
    
    

    //Check if this track belongs to a jet and compute the IPs
    for (size_t ijet = 0; ijet<jets.size(); ++ijet) {
      std::vector<int> jtrks = jets[ijet].getTracks();
      
      if (std::find(jtrks.begin(), jtrks.end(),itrk) != jtrks.end()) {

        Acts::Vector3 jetDir{jets[ijet].getFourMomentum()[0],
          jets[ijet].getFourMomentum()[1],
          jets[ijet].getFourMomentum()[2]};
        
        
        Acts::Result<Acts::ImpactParametersAndSigma> ipAndSigma = m_ipEst->estimateImpactParameters(tracks[itrk], reco_vertices[HS_idx],
                                                                                                    gctx_, mctx_);
        
        Acts::Result<std::pair<double,double>> vszs = m_ipEst->getLifetimesSignOfTrack(tracks[itrk], reco_vertices[HS_idx],
                                                                                       jetDir, gctx_, mctx_);
        
        if (!ipAndSigma.ok() || !vszs.ok())
          continue;

        //This is not unbiased!
        signed_d0             = std::fabs((*ipAndSigma).IPd0)  * (*vszs).first;
        signed_z0SinTheta     = std::fabs((*ipAndSigma).IPz0SinTheta) * (*vszs).second;
        signed_d0_err         = (*ipAndSigma).sigmad0;
        signed_z0SinTheta_err = (*ipAndSigma).sigmaz0SinTheta;

      }//track in jet
    }//loop on jets
    
    
    
    if (trk_params.covariance().has_value()) {

      const auto& cov  = *trk_params.covariance();
      
      float trk_theta = params[Acts::eBoundTheta];
      float trk_eta   = std::atanh(std::cos(trk_theta));
      float trk_qop   = params[Acts::eBoundQOverP];
      float trk_p     = std::abs(1.0 / trk_qop);
      float trk_pt    = trk_p * std::sin(trk_theta); 
      
      m_tracks_prob.push_back(1.);   // todo
      m_trk_d0.push_back(params[Acts::eBoundLoc0]);             
      m_trk_z0.push_back(params[Acts::eBoundLoc1]);
      
      m_trk_signed_d0.push_back(signed_d0);
      m_trk_signed_d0sig.push_back(signed_d0 / signed_d0_err);
      m_trk_signed_z0sinTheta.push_back(signed_z0SinTheta);
      m_trk_signed_z0sinThetasig.push_back(signed_z0SinTheta / signed_z0SinTheta_err);

      m_trk_numPix1L.push_back(hit_infos.nPixInnermost);
      m_trk_numPix2L.push_back(hit_infos.nPixNextToInnermost);
      m_trk_numPix  .push_back(hit_infos.nPix);
      m_trk_numSCT  .push_back(hit_infos.nSStrip);
      m_trk_numLSCT .push_back(hit_infos.nLStrip);
      
      
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
      
      //This give un-initialized warnings
      
      covd0          = cov(Acts::eBoundLoc0,Acts::eBoundLoc0);
      covz0          = cov(Acts::eBoundLoc1,Acts::eBoundLoc1);
      covphi         = cov(Acts::eBoundPhi,Acts::eBoundPhi);
      covtheta       = cov(Acts::eBoundTheta,Acts::eBoundTheta);
      covqOverP      = cov(Acts::eBoundQOverP,Acts::eBoundQOverP);
      covd0z0        = cov(Acts::eBoundLoc0,Acts::eBoundLoc1);
      covd0phi       = cov(Acts::eBoundLoc0,Acts::eBoundPhi);
      covd0theta     = cov(Acts::eBoundLoc0,Acts::eBoundTheta);
      covd0qOverP    = cov(Acts::eBoundLoc0,Acts::eBoundQOverP);
      covz0phi       = cov(Acts::eBoundLoc1,Acts::eBoundPhi);
      covz0theta     = cov(Acts::eBoundLoc1,Acts::eBoundTheta);
      covz0qOverP    = cov(Acts::eBoundLoc1,Acts::eBoundQOverP);
      covphitheta    = cov(Acts::eBoundPhi,Acts::eBoundTheta);
      covphiqOverP   = cov(Acts::eBoundPhi,Acts::eBoundQOverP);
      covthetaqOverP = cov(Acts::eBoundTheta,Acts::eBoundQOverP);
            
      
      m_trk_var_d0.push_back(covd0);
      m_trk_var_z0.push_back(covz0);         
      m_trk_var_phi.push_back(covphi);        
      m_trk_var_theta.push_back(covtheta);      
      m_trk_var_qOverP.push_back(covqOverP);     
      m_trk_cov_d0z0.push_back(covd0z0);       
      m_trk_cov_d0phi.push_back(covd0phi);      
      m_trk_cov_d0theta.push_back(covd0theta);    
      m_trk_cov_d0qOverP.push_back(covd0qOverP);   
      m_trk_cov_z0phi.push_back(covz0phi);      
      m_trk_cov_z0theta.push_back(covz0theta);    
      m_trk_cov_z0qOverP.push_back(covz0qOverP);   
      m_trk_cov_phitheta.push_back(covphitheta);   
      m_trk_cov_phiqOverP.push_back(covphiqOverP);  
      m_trk_cov_thetaqOverP.push_back(covthetaqOverP);
      
    }// cov has value
  } //Tracks
  
  
  m_outputTree->Fill();
  
  return ProcessCode::SUCCESS;
}



ActsExamples::HitInfo ActsExamples::RootEventWriter::GetHitInformation(const Acts::MultiTrajectoryHelpers::TrajectoryState& trajState) const  {
  
  HitInfo hit_info;

  for (size_t imeas =  0; imeas < trajState.nMeasurements; imeas++) {

    int volume = trajState.measurementVolume.at(imeas);
    int layer  = trajState.measurementLayer.at(imeas);
    if (volume == 16 || volume == 17 || volume == 18) {

      hit_info.nPix++;
      if (volume == 17) {
        if (layer==2)
          hit_info.nPixInnermost++;
        else if (layer==4)
          hit_info.nPixNextToInnermost++;
      }
    }
    else if (volume == 23 || volume == 24 || volume == 25) {
        hit_info.nSStrip++;
      }
    else if (volume == 28 || volume == 29 || volume == 30) {
        hit_info.nLStrip++;
      }
  }
  return hit_info;
}

void ActsExamples::RootEventWriter::Clear() {

  //Vertices
  m_vtx_x.clear();
  m_vtx_y.clear();
  m_vtx_z.clear();
  m_vtx_t.clear();
  m_vtx_sumPt2.clear();
  m_vtx_isHS.clear();
  m_vtx_isPU.clear();
  
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
  m_trk_signed_d0.clear();
  m_trk_signed_d0sig.clear();
  m_trk_signed_z0sinTheta.clear();
  m_trk_signed_z0sinThetasig.clear();
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
  m_trk_cov_thetaqOverP.clear();

  m_trk_numPix1L.clear();
  m_trk_numPix2L.clear();
  m_trk_numPix  .clear();
  m_trk_numSCT  .clear();
  m_trk_numLSCT .clear();
  
}
