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
  m_outputTracks.initialize(m_cfg.outputTracks);

  m_outfile = new TFile("outfile.root","RECREATE");
  m_tree = new TTree("tree","tree");
  m_tree->Branch("trackp",      &m_trackp);
  m_tree->Branch("tracketa",    &m_tracketa);
  m_tree->Branch("n_time_meas", &m_n_time_meas);
  m_tree->Branch("chi2min",     &m_chi2min);
  m_tree->Branch("deltaChi2",   &m_deltaChi2);
  m_tree->Branch("pididx",      &m_pididx);
  
}

ProcessCode TrackFittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

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
  
  auto computeTimeChi2 = [](Acts::TrackProxy trk) {

    double tchi2 = 0;
    
    if (trk.hasReferenceSurface()) {
      
      for (const auto &trackState : track.trackStatesReversed()) {
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
            
            auto smoothed_state = trackState.smoothed();
            auto smoothed_cov = trackState.smoothedCovariance();
            
            //The profiled Chi2
            
            
            
        }
      }
    
  };
  
      
      





  
  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // Check if you are not in picking mode
    if (m_cfg.pickTrack > -1 &&
        static_cast<std::size_t>(m_cfg.pickTrack) != itrack) {
      continue;
    }

    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];
    const auto& initialParams = initialParameters[itrack];

    // initial Params with pion hyp
    const auto pi_initialParams = Acts::GenericBoundTrackParameters(initialParams.referenceSurface().getSharedPtr(),
                                                                    initialParams.parameters(),
                                                                    initialParams.covariance(),
                                                                    Acts::SinglyChargedParticleHypothesis::pion());
    
    // initial Params with proton hyp
    const auto proton_initialParams = Acts::GenericBoundTrackParameters(initialParams.referenceSurface().getSharedPtr(),
                                                                        initialParams.parameters(),
                                                                        initialParams.covariance(),
                                                                        Acts::SinglyChargedParticleHypothesis::proton());
    
    // initial Params with kaon hyp
    const auto kaon_initialParams = Acts::GenericBoundTrackParameters(initialParams.referenceSurface().getSharedPtr(),
                                                                      initialParams.parameters(),
                                                                      initialParams.covariance(),
                                                                      Acts::SinglyChargedParticleHypothesis::kaon());

    
   
    

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
    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, tracks);

    auto result_pi = (*m_cfg.fit)(trackSourceLinks, pi_initialParams, options,
                               calibrator, tracks);

    auto result_proton = (*m_cfg.fit)(trackSourceLinks, proton_initialParams, options,
                               calibrator, tracks);

    auto result_kaon = (*m_cfg.fit)(trackSourceLinks, k_initialParams, options,
                               calibrator, tracks);

    if (result.ok() && result_pi.ok() && result_proton.ok() && result_kaon.ok()) {
      // Get the fit output object
      const auto& track        = result.value();
      const auto& track_pi     = result_pi.value();
      const auto& track_proton = result_proton.value();
      const auto& track_kaon   = result_kaon.value();



                  
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
              double p = std::abs(1. / qop);
              double qopcov    = smoothed_cov(Acts::eBoundQOverP,Acts::eBoundQOverP);

              double s = trackState.pathLength();
              n_time_meas++;
              
              arcLs.push_back(s);
              ps.push_back(p);
              qopcovs.push_back(qopcov);
              tmeas.push_back(trackState.effectiveCalibrated()[i]);
              tmeascovs.push_back(trackState.effectiveCalibratedCovariance()(i,i));
              
            }
          }
        }
                        
        auto computeTimeChi2 = [&](double mass) {

          const size_t N = arcLs.size();
          std::vector<double> tpred(N, 0.0);

          // The arc length and the predicted time to the previous state
          double s  = 0 ;
          double t0 = 0 ;
          double Sw = 0.0, Swd = 0.0, Swd2 = 0.0;
          
          for (int j = N - 1 ; j>=0; j--) {
            
            double ds    = arcLs[j] - s;
            double b     = beta(ps[j], mass);
            double dt    = ds / b;
            tpred[j]     = dt + t0;  
            
            double d = tmeas[j] - tpred[j];
            
            double w = 1. / (tmeascovs[j]);
            
            Sw   += w;
            Swd  += w * d;
            Swd2 += w * d * d;
            
            std::cout<<"PF:: "<<arcLs[j]<<" ds="<<ds<<" s="<<s <<std::endl;
            std::cout<<"PF:: "<<tpred[j]<<" ds="<<tpred[j]-t0<<" t0="<<t0 <<std::endl;
            std::cout<<"PF:: tpred="<<tpred[j]<<" tmeas="<<tmeas[j]<<" d="<<d<<" w="<<w<<std::endl;
            std::cout<<std::endl;
            
            s = arcLs[j];
            t0 = tpred[j];
                        
          }
          
          double t0_hat = (Sw > 0) ? (Swd / Sw) : 0.0;
          double chi2   = Swd2 - (Swd*Swd)/(Sw > 0 ? Sw : 1.0);

          //std::cout<<t0_hat<<"  "<<chi2<<std::endl;
          
          return std::pair<double,double>(chi2,t0_hat);
          
          
        };

        std::vector<double> chi2s;
        std::vector<double> t0hats;
        std::vector<double> deltaChi2;

        auto mu     = computeTimeChi2(0.105);
        auto pion   = computeTimeChi2(0.139);
        auto kaon   = computeTimeChi2(0.494);
        auto proton = computeTimeChi2(0.930);

        chi2s.push_back(mu.first);
        chi2s.push_back(pion.first);
        chi2s.push_back(kaon.first);
        chi2s.push_back(proton.first);
        
        t0hats.push_back(mu.second);
        t0hats.push_back(pion.second);
        t0hats.push_back(kaon.second);
        t0hats.push_back(proton.second);
        
        auto it = std::min_element(chi2s.begin(), chi2s.end());
        std::size_t minIndex = std::distance(chi2s.begin(), it);
        
        double minChi2 = chi2s[minIndex];

        std::cout<<minIndex<<" "<<minChi2<<std::endl;

        // Compute Δχ² for each hypothesis
        
        deltaChi2.reserve(chi2s.size());
        for (double val : chi2s) {
          deltaChi2.push_back(val - minChi2);
        }
        const auto& parameters = track.parameters();
        double trackp =     std::abs(1. / parameters[Acts::eBoundQOverP]);
        double tracktheta = parameters[Acts::eBoundTheta];

        
        
        
        std::cout<<"CHI2S="<<chi2s[0]<<" "<<chi2s[1]<<" "<<chi2s[2]<<" "<<chi2s[3]<<std::endl;
        std::cout<<"t0S="<<t0hats[0]<<" "<<t0hats[1]<<" "<<t0hats[2]<<" "<<t0hats[3]<<"    ->"<<parameters[Acts::eBoundTime]<<std::endl;
        

        
        auto thetaToEta =  [](double theta) {
          return -std::log(std::tan(theta / 2.0));
        };
                
        double tracketa = thetaToEta(tracktheta);
        m_trackp.push_back(trackp);
        m_tracketa.push_back(tracketa);
        m_n_time_meas.push_back(n_time_meas);
        m_chi2min.push_back(chi2s);
        m_deltaChi2.push_back(deltaChi2);
        m_pididx.push_back(minIndex);
        
        m_tree->Fill();

        
      } else {
        ACTS_VERBOSE("No fitted parameters for track " << itrack);
      }
    } else {
      ACTS_WARNING("Fit failed for track "
                   << itrack << " with error: " << result.error() << ", "
                   << result.error().message());
    }

    m_trackp.clear();
    m_tracketa.clear();
    m_n_time_meas.clear();
    m_chi2min.clear();
    m_deltaChi2.clear();
    m_pididx.clear();
    
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
