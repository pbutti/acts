// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>
#include <algorithm>

namespace ActsExamples {

TrackFittingAlgorithm::TrackFittingAlgorithm(Config config,
                                             Acts::Logging::Level level)
    : IAlgorithm("TrackFittingAlgorithm", level), m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurement collection");
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto tracks collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }
  if (!m_cfg.calibrator) {
    throw std::invalid_argument("Missing calibrator");
  }
  if (m_cfg.inputClusters.empty() && m_cfg.calibrator->needsClusters()) {
    throw std::invalid_argument("The configured calibrator needs clusters");
  }

  
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_inputParticles.initialize("particles");
  m_outputTracks.initialize(m_cfg.outputTracks);

  m_outfile = new TFile("outfile.root","RECREATE");
  m_tree = new TTree("tree","tree");
  m_tree->Branch("trackp",            &m_trackp);
  m_tree->Branch("truthp",            &m_truthp);
  m_tree->Branch("tracketa",          &m_tracketa);
  m_tree->Branch("n_time_meas",       &m_n_time_meas);
  m_tree->Branch("chi2min",           &m_chi2min);
  m_tree->Branch("chi2min_prof",      &m_chi2min_prof);
  m_tree->Branch("chi2min_simple",   &m_chi2min_simple);
  m_tree->Branch("deltaChi2",         &m_deltaChi2);
  m_tree->Branch("deltaChi2_prof",    &m_deltaChi2_prof);
  m_tree->Branch("deltaChi2_simple", &m_deltaChi2_simple);
  m_tree->Branch("pididx",            &m_pididx);
  m_tree->Branch("pididx_prof",       &m_pididx_prof);
  m_tree->Branch("pididx_simple"   , &m_pididx_simple);
  m_tree->Branch("t0hat",        &m_t0hat);
  m_tree->Branch("t0hat_simple", &m_t0hat_simple);
  
}

ProcessCode TrackFittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  const auto& particles = m_inputParticles(ctx);

  const ClusterContainer* clusters =
      m_inputClusters.isInitialized() ? &m_inputClusters(ctx) : nullptr;

  ACTS_DEBUG("Input measurements: " << measurements.size());
  ACTS_DEBUG("Input proto tracks: " << protoTracks.size());
  ACTS_DEBUG("Input initial parameters: " << initialParameters.size());
  ACTS_DEBUG("Input clusters: " << (clusters ? clusters->size() : 0));

  // Consistency cross checks
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters "
               << protoTracks.size() << " vs " << initialParameters.size());
    return ProcessCode::ABORT;
  }

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  // Measurement calibrator must be instantiated here, because we need the
  // measurements to construct it. The other extensions are hold by the
  // fit-function-object
  MeasurementCalibratorAdapter calibrator(*(m_cfg.calibrator), measurements,
                                          clusters);

  TrackFitterFunction::GeneralFitterOptions options{
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext, pSurface.get(),
      Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext)};

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  // reserve space for the track containers
  // simply assume 30 states per track for now
  trackContainer->reserve(protoTracks.size());
  trackStateContainer->reserve(protoTracks.size() * 30);

  // Perform the fit for each input track
  std::vector<Acts::SourceLink> trackSourceLinks;

  auto beta = [](double p, double m0) {
    return p / std::sqrt(p * p + m0 * m0);
  };


  // Derivative of the predicted time wrt QoP
  auto dQoPdTp = [](double qOp, double m0, double s) {

    auto p = std::abs(qOp);
    auto q = qOp > 0 ? 1. : -1.;
    return s * m0*m0 / (q*std::sqrt(p*p+m0*m0));
  };
    
  auto thetaToEta =  [](double theta) {
    return -std::log(std::tan(theta / 2.0));
  };
  

  auto computeTimeChi2_v2 = [&](const std::vector<double> arcLs,
                                const std::vector<double> ps,
                                const std::vector<double> qopcovs,
                                const std::vector<double> tmeas,
                                const std::vector<double> tmeascovs,
                                double mass) {
          
    const size_t N = arcLs.size();
    std::vector<double> tpred(N, 0.0);
    
    // The arc length and the predicted time to the previous state
    double s0  = 0 ;
    double t0 = 0 ;
    double Sw = 0.0, Swd = 0.0, Swd2 = 0.0;

    double chi2_nom = 0.;
    
    for (int j = N - 1 ; j>=0; j--) {
      
      //double ds    = arcLs[j] - s;
      double ds    = arcLs[j] - s0;
      double b     = beta(ps[j], mass);
      double dt    = ds / b;
      tpred[j]     = dt - t0;  
      
      double d = tmeas[j] - tpred[j];

      // Charge can be 1 as we don't need it.

      auto dTp = dQoPdTp(1./ps[j], mass, ds);

      
      double cov_tp = dTp*dTp*qopcovs[j];

      //std::cout<<"ds="<<ds<<" m0="<<mass<<" qop="<<1./ps[j]<<" dTp="<<dTp<<" cov_tp="<<cov_tp<<" qopcovs[j]="<<qopcovs[j]<<std::endl;
      
      double w = 1. / (tmeascovs[j] + cov_tp);
      
      Sw   += w;
      Swd  += w * d;
      Swd2 += w * d * d;
      
      //std::cout<<"PF:: "<<arcLs[j]<<" ds="<<ds<<" s="<<s0 <<std::endl;
      //std::cout<<"PF:: "<<tpred[j]<<" dt="<<tpred[j]-t0<<" t0="<<t0 <<std::endl;
      //std::cout<<"PF:: tpred="<<tpred[j]<<" tmeas="<<tmeas[j]<<" d="<<d<<" w="<<w<<std::endl;
      //std::cout<<std::endl;

      // Do not do arc length by arc lenght.
      //s = arcLs[j];
      //t0 = tpred[j];

      chi2_nom+=w*d*d;
    }
    
    double t0_hat = (Sw > 0) ? (Swd / Sw) : 0.0;
    double chi2   = Swd2 - (Swd*Swd)/(Sw > 0 ? Sw : 1.0);
    
    //std::cout<<t0_hat<<"  "<<chi2<<std::endl;


    //std::cout<<"POST FIT residuals"<<std::endl;
    
    double chi2_hat = 0.;
    for (int j = N - 1 ; j>=0; j--) {
      
      double ds    = arcLs[j];
      double b     = beta(ps[j], mass);
      double dt    = ds / b;
      tpred[j]     = dt + t0_hat;  
      
      double d = tmeas[j] - tpred[j];

      // Charge can be 1 as we don't need it.

      auto dTp = dQoPdTp(1./ps[j], mass, ds);

      
      double cov_tp = dTp*dTp*qopcovs[j];

      //std::cout<<"ds="<<ds<<" m0="<<mass<<" qop="<<1./ps[j]<<" dTp="<<dTp<<" cov_tp="<<cov_tp<<" qopcovs[j]="<<qopcovs[j]<<std::endl;
      
      // TODO:: This weight is neglecting the momentum error
      double w = 1. / (tmeascovs[j] + cov_tp);
      
      //std::cout<<"PF:: tpred="<<tpred[j]<<" tmeas="<<tmeas[j]<<" d="<<d<<" w="<<w<<std::endl;
      //std::cout<<std::endl;

      chi2_hat+=d*d*w;
      // Do not do arc length by arc lenght.
      //s = arcLs[j];
      //t0 = tpred[j];
    }

    //std::cout<<"Chi2_hat = "<< chi2_hat<<" chi2="<<chi2<<" chi2_nom="<<chi2_nom<<std::endl;
    
    return std::pair<double,double>(chi2,t0_hat);
    
  };
      
  auto computeTimeChi2 = [](ActsExamples::TrackProxy trk) {
    
    double tchi2_simple = 0;
    double tchi2_profiled = 0;
    
    if (trk.hasReferenceSurface()) {
      
      for (const auto &trackState : trk.trackStatesReversed()) {
        Acts::ConstTrackStateType typeFlags = trackState.typeFlags();
        if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
          
          if (!trackState.hasCalibrated()) {
            throw std::invalid_argument("track state has no calibrated parameters");
          }
          
          auto subspace = trackState.projectorSubspaceHelper();
          
          if (!subspace.contains(Acts::eBoundTime))
            continue;
          
          auto i = subspace.indexOf(Acts::eBoundTime);
          
          if (trackState.hasSmoothed()) {
            
            size_t size = trackState.calibratedSize();
            
            auto smoothed_state = trackState.smoothed();
            auto smoothed_cov   = trackState.smoothedCovariance();
            
            //auto proj_smooth_state = subspace.projectVector(smoothed_state);
            //auto proj_smooth_cov   = subspace.projectMatrix(smoothed_cov);

            auto H =
              trackState.projectorSubspaceHelper().fullProjector()
              .topLeftCorner(
                             trackState.calibratedSize(), Acts::eBoundSize);

            auto proj_smooth_state = H * smoothed_state;
            auto proj_smooth_cov = H*smoothed_cov*H.transpose();
            
            //auto V = state.effectiveCalibratedCovariance();
            
            //auto resCov = V + H * covariance * H.transpose();
                        
            
            auto meas = trackState.effectiveCalibrated();
            auto meas_cov = trackState.effectiveCalibratedCovariance();
            auto res  = meas - proj_smooth_state;
            auto res_cov  = meas_cov - proj_smooth_cov;

            // Check time-only and profiled residuals

            auto ra  = res.head(size-1);
            auto Caa = res_cov.topLeftCorner(size-1,size-1);
            auto Cta = res_cov.bottomLeftCorner(1, size-1);
            auto Ctt = res_cov(size-1, size-1);

            /*
            std::cout<<"###########"<<std::endl;
            std::cout<<"PF::"<<std::endl;
            std::cout<<smoothed_state<<std::endl;
            std::cout<<smoothed_cov<<std::endl;
            std::cout<<"-----------"<<std::endl;
            std::cout<<proj_smooth_state<<std::endl;
            std::cout<<proj_smooth_cov<<std::endl;
            std::cout<<"-----------"<<std::endl;
            std::cout<<meas<<std::endl;
            std::cout<<meas_cov<<std::endl;
            std::cout<<res<<std::endl;
            std::cout<<res_cov<<std::endl;
            std::cout<<"-----------"<<std::endl;
            std::cout<<Caa<<std::endl;
            std::cout<<Cta<<std::endl;
            std::cout<<Ctt<<std::endl;
            std::cout<<"-----------"<<std::endl;
            */
            
            double chi2_simple = res(i) * res(i) / Ctt;
            //std::cout<<res(i)<<" ---  "<<Ctt<<std::endl;
                                                         

            // Use inverse for the moment
            double rt_a = res(i) - (Cta*(Caa.inverse())*ra).value();
            double Ct_a = Ctt - (Cta*(Caa.inverse())*(Cta.transpose())).value();
            
            //std::cout<<rt_a<<" ---  "<<Ct_a<<std::endl;
            double chi2_profiled = rt_a * rt_a / Ct_a;
            
            tchi2_simple+=chi2_simple;
            tchi2_profiled+=chi2_profiled;
          }
        }
      }
    }

    return std::pair<double,double>(tchi2_simple,tchi2_profiled);
  };
  
  
  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // Check if you are not in picking mode
    if (m_cfg.pickTrack > -1 &&
        static_cast<std::size_t>(m_cfg.pickTrack) != itrack) {
      continue;
    }

    if (particles.size() !=1 )
      continue;
    
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];
    const auto& initialParams = initialParameters[itrack];


    //std::cout<<initialParams<<std::endl;

    auto seedQoP = initialParams.parameters()[Acts::eBoundQOverP];

    auto seedP = std::abs(1./seedQoP);

    //std::cout<<"PF::: Processing..."<<std::endl;
    //std::cout<<initialParams<<std::endl;
    
    // initial Params with pion hyp
    const auto& ele_initialParams = Acts::GenericBoundTrackParameters(initialParams.referenceSurface().getSharedPtr(),
                                                                     initialParams.parameters(),
                                                                     initialParams.covariance(),
                                                                     Acts::ParticleHypothesis::electron());
    // initial Params with proton hyp
    const auto& proton_initialParams = Acts::GenericBoundTrackParameters(initialParams.referenceSurface().getSharedPtr(),
                                                                         initialParams.parameters(),
                                                                         initialParams.covariance(),
                                                                         Acts::ParticleHypothesis::proton());
    
    // initial Params with kaon hyp
    const auto& kaon_initialParams = Acts::GenericBoundTrackParameters(initialParams.referenceSurface().getSharedPtr(),
                                                                       initialParams.parameters(),
                                                                       initialParams.covariance(),
                                                                       Acts::ParticleHypothesis::kaon());
    
    // We can have empty tracks which must give empty fit results so the number
    // of entries in input and output containers matches.
    if (protoTrack.empty()) {
      ACTS_WARNING("Empty track " << itrack << " found.");
      continue;
    }

    ACTS_VERBOSE("Initial 4 position: "
                 << initialParams.fourPosition(ctx.geoContext).transpose());
    ACTS_VERBOSE(
        "Initial direction: " << initialParams.direction().transpose());
    ACTS_VERBOSE("Initial momentum: " << initialParams.absoluteMomentum());

    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(protoTrack.size());

    // Fill the source links via their indices from the container
    for (auto measIndex : protoTrack) {
      ConstVariableBoundMeasurementProxy measurement =
          measurements.getMeasurement(measIndex);
      IndexSourceLink sourceLink(measurement.geometryId(), measIndex);
      trackSourceLinks.push_back(Acts::SourceLink(sourceLink));
    }

    ACTS_VERBOSE("Invoke fitter for track " << itrack);
    //std::cout<<"pion track"<<std::endl;
    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, tracks);
    //std::cout<<"ele track"<<std::endl;
    auto result_ele = (*m_cfg.fit)(trackSourceLinks, ele_initialParams, options,
                               calibrator, tracks);
    //std::cout<<"proton track"<<std::endl;
    auto result_proton = (*m_cfg.fit)(trackSourceLinks, proton_initialParams, options,
                               calibrator, tracks);
    //std::cout<<"kaon track"<<std::endl;
    auto result_kaon = (*m_cfg.fit)(trackSourceLinks, kaon_initialParams, options,
                               calibrator, tracks);

    if (result.ok() && result_ele.ok() && result_proton.ok() && result_kaon.ok()) {
      
      // Get the fit output object
      const auto& track        = result.value();
      const auto& track_ele    = result_ele.value();
      const auto& track_proton = result_proton.value();
      const auto& track_kaon   = result_kaon.value();
            
      auto [tchi2_ele,tchi2_ele_profiled] = computeTimeChi2(track_ele);
      auto [tchi2_pi,tchi2_pi_profiled] = computeTimeChi2(track);
      auto [tchi2_proton,tchi2_proton_profiled] = computeTimeChi2(track_proton);
      auto [tchi2_kaon,tchi2_kaon_profiled] = computeTimeChi2(track_kaon);
      
      //std::cout<<tchi2_muon<<" "<<tchi2_pi<<" "<<tchi2_proton<<" "<<tchi2_kaon<<std::endl;
            
      
      if (track.hasReferenceSurface()) {
        ACTS_VERBOSE("Fitted parameters for track " << itrack);
        ACTS_VERBOSE("  " << track.parameters().transpose());
        ACTS_VERBOSE("Measurements: (prototrack->track): "
                     << protoTrack.size() << " -> " << track.nMeasurements());
        
        
        //Compute the time Chi2 with different hypotheses
        
        std::vector<double> arcLs;
        std::vector<double> ps;
        std::vector<double> qopcovs;
        std::vector<double> tmeas;
        std::vector<double> tmeascovs;

        //auto ddtdqop = [](double ds, double qop, double m0) {
                    
        //Assume q = 1
        //return (m0*m0*ds*qop) / (std::sqrt(m0*m0*qop*qop + 1));
        //};


        int n_time_meas = 0;
        double s0 = 0;
        
        for (const auto &trackState : track.trackStatesReversed()) {
          //for (const auto &trackState : track.trackStates()) {
          Acts::ConstTrackStateType typeFlags = trackState.typeFlags();
          if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {

            if (!trackState.hasCalibrated()) {
              throw std::invalid_argument("track state has no calibrated parameters");
            }

            auto subspace = trackState.projectorSubspaceHelper();

            if (!subspace.contains(Acts::eBoundTime))
              continue;
                        
            auto i = subspace.indexOf(Acts::eBoundTime);
            
            if (trackState.hasSmoothed()) {
              auto smoothed_pars  = trackState.smoothed();
              auto smoothed_cov   = trackState.smoothedCovariance();
              
              double qop = smoothed_pars[Acts::eBoundQOverP];
              double eta = thetaToEta(smoothed_pars[Acts::eBoundTheta]);
              
              double p = std::abs(1. / qop);
              double qopcov    = smoothed_cov(Acts::eBoundQOverP,Acts::eBoundQOverP);
                            
              double s = trackState.pathLength();
              n_time_meas++;
              
              arcLs.push_back(s);
              ps.push_back(p);
              qopcovs.push_back(qopcov);
              tmeas.push_back(trackState.effectiveCalibrated()[i]);
              tmeascovs.push_back(trackState.effectiveCalibratedCovariance()(i,i));

              //std::cout<<"PF::"<<n_time_meas<<" "<<s<<" "<<p<<" "<<qopcov<<" "<<eta<<" "<<trackState.effectiveCalibrated()[i]<<" "<<trackState.effectiveCalibratedCovariance()(i,i)<<std::endl;
              
            }// has smoothed
          } // measurement
        } // reversed


        std::vector<double> chi2s;
        std::vector<double> chi2s_prof;
        std::vector<double> chi2s_simple;
        std::vector<double> t0hats;
        std::vector<double> t0hats_simple;
        std::vector<double> deltaChi2;
        std::vector<double> deltaChi2_prof;
        std::vector<double> deltaChi2_simple;


        //auto ele    = computeTimeChi2_v2(arcLs,ps,qopcovs,tmeas,tmeascovs,0.0005);
        //auto mu     = computeTimeChi2_v2(arcLs,ps,qopcovs,tmeas,tmeascovs,0.105);
        auto ele     = computeTimeChi2_v2(arcLs,ps,qopcovs,tmeas,tmeascovs,0.0005); // An electron instead
        auto pion    = computeTimeChi2_v2(arcLs,ps,qopcovs,tmeas,tmeascovs,0.139);
        auto proton  = computeTimeChi2_v2(arcLs,ps,qopcovs,tmeas,tmeascovs,0.930);
        auto kaon    = computeTimeChi2_v2(arcLs,ps,qopcovs,tmeas,tmeascovs,0.494);
        
        t0hats_simple.push_back(ele.second);
        t0hats_simple.push_back(pion.second);
        t0hats_simple.push_back(proton.second);
        t0hats_simple.push_back(kaon.second);
        
        chi2s.push_back(tchi2_ele);
        chi2s_prof.push_back(tchi2_ele_profiled);
        chi2s_simple.push_back(ele.first);
        
        chi2s.push_back(tchi2_pi);
        chi2s_prof.push_back(tchi2_pi_profiled);
        chi2s_simple.push_back(pion.first);
        
        chi2s.push_back(tchi2_proton);
        chi2s_prof.push_back(tchi2_proton_profiled);
        chi2s_simple.push_back(proton.first);
        
        chi2s.push_back(tchi2_kaon);
        chi2s_prof.push_back(tchi2_kaon_profiled);
        chi2s_simple.push_back(kaon.first);
        
        auto it        = std::min_element(chi2s.begin(), chi2s.end());
        auto it_prof   = std::min_element(chi2s_prof.begin(), chi2s_prof.end());
        auto it_simple = std::min_element(chi2s_simple.begin(), chi2s_simple.end());
        
        std::size_t minIndex = std::distance(chi2s.begin(), it);
        std::size_t minIndex_prof = std::distance(chi2s_prof.begin(), it_prof);
        std::size_t minIndex_simple = std::distance(chi2s_simple.begin(), it_simple);
        
        double minChi2 = chi2s[minIndex];
        double minChi2_prof = chi2s_prof[minIndex_prof];
        double minChi2_simple = chi2s_simple[minIndex_simple];
        
        // Compute Δχ² for each hypothesis - time only 
        deltaChi2.reserve(chi2s.size());
        for (double val : chi2s) {
          deltaChi2.push_back(val - minChi2);
        }
        
        // Compute Δχ² for each hypothesis - profiled
        deltaChi2_prof.reserve(chi2s_prof.size());
        for (double val : chi2s_prof) {
          deltaChi2_prof.push_back(val - minChi2_prof);
        }

        // Compute Δχ² for each hypothesis - simple
        deltaChi2_prof.reserve(chi2s_simple.size());
        for (double val : chi2s_simple) {
          deltaChi2_prof.push_back(val - minChi2_simple);
        }
        
        const auto& parameters = track.parameters();
        const auto& covariance = track.covariance();
        double trackp =     std::abs(1. / parameters[Acts::eBoundQOverP]);
        double trackqop = parameters[Acts::eBoundQOverP];
        double trackqopcov = covariance(Acts::eBoundQOverP,Acts::eBoundQOverP);
        double tracktheta = parameters[Acts::eBoundTheta];

        
        //std::cout<<"PF::: Found Track ..."<<std::endl;
        //std::cout<<parameters<<std::endl;
        //std::cout<<"#####################"<<std::endl;
        
        //std::cout<<"CHI2S="<<chi2s[0]<<" "<<chi2s[1]<<" "<<chi2s[2]<<" "<<chi2s[3]<<std::endl;
        //std::cout<<"t0S="<<t0hats[0]<<" "<<t0hats[1]<<" "<<t0hats[2]<<" "<<t0hats[3]<<"    ->"<<parameters[Acts::eBoundTime]<<std::endl;
        
        double tracketa = thetaToEta(tracktheta);
        auto particle = particles.begin();
        m_truthp.push_back(particle->absoluteMomentum());
        m_trackp.push_back(trackp);
        m_trackp.push_back(trackpopcov);
        m_tracketa.push_back(tracketa);
        m_n_time_meas.push_back(n_time_meas);
        
        m_chi2min.push_back(chi2s);
        m_chi2min_prof.push_back(chi2s_prof);
        m_chi2min_simple.push_back(chi2s_simple);
        
        m_deltaChi2.push_back(deltaChi2);
        m_deltaChi2_prof.push_back(deltaChi2_prof);
        m_deltaChi2_simple.push_back(deltaChi2_simple);
        
        m_pididx.push_back(minIndex);
        m_pididx_prof.push_back(minIndex_prof);
        m_pididx_simple.push_back(minIndex_simple);
        
        m_tree->Fill();

        
      } else {
        ACTS_VERBOSE("No fitted parameters for track " << itrack);
      }
    } else {
      if (not result.ok()) 
        ACTS_WARNING("Fit failed for pion track "
                     << itrack << " with error: " << result.error() << ", "
                     << result.error().message());
      else if (not result_ele.ok())
        ACTS_WARNING("Fit failed for ele track "
                     << itrack << " with error: " << result_ele.error() << ", "
                     << result_ele.error().message());
      else if (not result_proton.ok())
        ACTS_WARNING("Fit failed for proton track "
                     << itrack << " with error: " << result_proton.error() << ", "
                     << result_proton.error().message());
      else if (not result_kaon.ok())
        ACTS_WARNING("Fit failed for kaon track "
                     << itrack << " with error: " << result_kaon.error() << ", "
                     << result_kaon.error().message());
    }

    m_truthp.clear();
    m_trackp.clear();
    m_tracketa.clear();
    m_n_time_meas.clear();

    m_chi2min.clear();
    m_chi2min_prof.clear();
    m_chi2min_simple.clear();

    m_deltaChi2.clear();
    m_deltaChi2_prof.clear();
    m_deltaChi2_simple.clear();
    
    m_pididx.clear();
    m_pididx_prof.clear();
    m_pididx_simple.clear();

    m_t0hat.clear();
    m_t0hat_simple.clear();
    
  }

  if (logger().doPrint(Acts::Logging::DEBUG)) {
    std::stringstream ss;
    trackStateContainer->statistics().toStream(ss);
    ACTS_DEBUG(ss.str());
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTracks(ctx, std::move(constTracks));


  
  return ProcessCode::SUCCESS;
}
  
  
  ProcessCode TrackFittingAlgorithm::finalize()  {

    m_outfile->cd();

    m_tree->Write();

    m_outfile->Close();

    
    return ProcessCode::SUCCESS;
  }
  
  
}  // namespace ActsExamples
