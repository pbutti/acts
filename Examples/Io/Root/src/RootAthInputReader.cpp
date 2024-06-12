// This file is part of the Acts project.
//
// Copyright (C) 2022-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthInputReader.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <boost/container/static_vector.hpp>
#include <TChain.h>

ActsExamples::RootAthInputReader::RootAthInputReader(
  const ActsExamples::RootAthInputReader::Config& config,
  Acts::Logging::Level level)
  : ActsExamples::IReader(),
  m_cfg(config),
  m_logger(Acts::getDefaultLogger(name(),level)) {

    if(m_cfg.inputfile.empty()) {
      throw std::invalid_argument("Missing input filename"); 
    }
    if (m_cfg.treename.empty()) {
      throw std::invalid_argument("Missing tree name");
    }
  
  m_inputchain = std::make_shared<TChain>(m_cfg.treename.c_str());

  m_outputPixelSpacePoints.initialize(m_cfg.outputPixelSpacePoints);
  m_outputStripSpacePoints.initialize(m_cfg.outputStripSpacePoints);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);



  // Set the branches

  // Set object pointer
  CLhardware = 0;
  CLparticleLink_eventIndex = 0;
  CLparticleLink_barcode = 0;
  CLbarcodesLinked = 0;
  CLparticle_charge = 0;
  CLphis = 0;
  CLetas = 0;
  CLtots = 0;
  CLlocal_cov = 0;
  Part_vParentID = 0;
  Part_vParentBarcode = 0;
  SPtopStripDirection = nullptr;
  SPbottomStripDirection = nullptr;
  SPstripCenterDistance = nullptr;
  SPtopStripCenterPosition = nullptr;
  TRKproperties = 0;
  TRKpattern = 0;
  TRKmeasurementsOnTrack_pixcl_sctcl_index = 0;
  TRKoutliersOnTrack_pixcl_sctcl_index = 0;
  TRKperigee_position = 0;
  TRKperigee_momentum = 0;
  DTTtrajectory_eventindex = 0;
  DTTtrajectory_barcode = 0;
  DTTstTruth_subDetType = 0;
  DTTstTrack_subDetType = 0;
  DTTstCommon_subDetType = 0;


  m_inputchain->SetBranchAddress("run_number", &run_number, &b_run_number);
  m_inputchain->SetBranchAddress("event_number", &event_number, &b_event_number);
  m_inputchain->SetBranchAddress("nSE", &nSE, &b_nSE);
  m_inputchain->SetBranchAddress("SEID", SEID, &b_SEID);
  m_inputchain->SetBranchAddress("nCL", &nCL, &b_nCL);
  m_inputchain->SetBranchAddress("CLindex", CLindex, &b_CLindex);
  m_inputchain->SetBranchAddress("CLhardware", &CLhardware, &b_CLhardware);
  m_inputchain->SetBranchAddress("CLx", CLx, &b_CLx);
  m_inputchain->SetBranchAddress("CLy", CLy, &b_CLy);
  m_inputchain->SetBranchAddress("CLz", CLz, &b_CLz);
  m_inputchain->SetBranchAddress("CLbarrel_endcap", CLbarrel_endcap, &b_CLbarrel_endcap);
  m_inputchain->SetBranchAddress("CLlayer_disk", CLlayer_disk, &b_CLlayer_disk);
  m_inputchain->SetBranchAddress("CLeta_module", CLeta_module, &b_CLeta_module);
  m_inputchain->SetBranchAddress("CLphi_module", CLphi_module, &b_CLphi_module);
  m_inputchain->SetBranchAddress("CLside", CLside, &b_CLside);
  m_inputchain->SetBranchAddress("CLmoduleID", CLmoduleID, &b_CLmoduleID);
  m_inputchain->SetBranchAddress("CLparticleLink_eventIndex", &CLparticleLink_eventIndex, &b_CLparticleLink_eventIndex);
  m_inputchain->SetBranchAddress("CLparticleLink_barcode", &CLparticleLink_barcode, &b_CLparticleLink_barcode);
  m_inputchain->SetBranchAddress("CLbarcodesLinked", &CLbarcodesLinked, &b_CLbarcodesLinked);
  m_inputchain->SetBranchAddress("CLparticle_charge", &CLparticle_charge, &b_CLparticle_charge);
  m_inputchain->SetBranchAddress("CLphis", &CLphis, &b_CLphis);
  m_inputchain->SetBranchAddress("CLetas", &CLetas, &b_CLetas);
  m_inputchain->SetBranchAddress("CLtots", &CLtots, &b_CLtots);
  m_inputchain->SetBranchAddress("CLloc_direction1", CLloc_direction1, &b_CLloc_direction1);
  m_inputchain->SetBranchAddress("CLloc_direction2", CLloc_direction2, &b_CLloc_direction2);
  m_inputchain->SetBranchAddress("CLloc_direction3", CLloc_direction3, &b_CLloc_direction3);
  m_inputchain->SetBranchAddress("CLJan_loc_direction1", CLJan_loc_direction1, &b_CLJan_loc_direction1);
  m_inputchain->SetBranchAddress("CLJan_loc_direction2", CLJan_loc_direction2, &b_CLJan_loc_direction2);
  m_inputchain->SetBranchAddress("CLJan_loc_direction3", CLJan_loc_direction3, &b_CLJan_loc_direction3);
  m_inputchain->SetBranchAddress("CLpixel_count", CLpixel_count, &b_CLpixel_count);
  m_inputchain->SetBranchAddress("CLcharge_count", CLcharge_count, &b_CLcharge_count);
  m_inputchain->SetBranchAddress("CLloc_eta", CLloc_eta, &b_CLloc_eta);
  m_inputchain->SetBranchAddress("CLloc_phi", CLloc_phi, &b_CLloc_phi);
  m_inputchain->SetBranchAddress("CLglob_eta", CLglob_eta, &b_CLglob_eta);
  m_inputchain->SetBranchAddress("CLglob_phi", CLglob_phi, &b_CLglob_phi);
  m_inputchain->SetBranchAddress("CLeta_angle", CLeta_angle, &b_CLeta_angle);
  m_inputchain->SetBranchAddress("CLphi_angle", CLphi_angle, &b_CLphi_angle);
  m_inputchain->SetBranchAddress("CLnorm_x", CLnorm_x, &b_CLnorm_x);
  m_inputchain->SetBranchAddress("CLnorm_y", CLnorm_y, &b_CLnorm_y);
  m_inputchain->SetBranchAddress("CLnorm_z", CLnorm_z, &b_CLnorm_z);
  m_inputchain->SetBranchAddress("CLlocal_cov", &CLlocal_cov, &b_CLlocal_cov);
  m_inputchain->SetBranchAddress("nPartEVT", &nPartEVT, &b_nPartEVT);
  m_inputchain->SetBranchAddress("Part_event_number", Part_event_number, &b_Part_event_number);
  m_inputchain->SetBranchAddress("Part_barcode", Part_barcode, &b_Part_barcode);
  m_inputchain->SetBranchAddress("Part_px", Part_px, &b_Part_px);
  m_inputchain->SetBranchAddress("Part_py", Part_py, &b_Part_py);
  m_inputchain->SetBranchAddress("Part_pz", Part_pz, &b_Part_pz);
  m_inputchain->SetBranchAddress("Part_pt", Part_pt, &b_Part_pt);
  m_inputchain->SetBranchAddress("Part_eta", Part_eta, &b_Part_eta);
  m_inputchain->SetBranchAddress("Part_vx", Part_vx, &b_Part_vx);
  m_inputchain->SetBranchAddress("Part_vy", Part_vy, &b_Part_vy);
  m_inputchain->SetBranchAddress("Part_vz", Part_vz, &b_Part_vz);
  m_inputchain->SetBranchAddress("Part_radius", Part_radius, &b_Part_radius);
  m_inputchain->SetBranchAddress("Part_status", Part_status, &b_Part_status);
  m_inputchain->SetBranchAddress("Part_charge", Part_charge, &b_Part_charge);
  m_inputchain->SetBranchAddress("Part_pdg_id", Part_pdg_id, &b_Part_pdg_id);
  m_inputchain->SetBranchAddress("Part_passed", Part_passed, &b_Part_passed);
  m_inputchain->SetBranchAddress("Part_vProdNin", Part_vProdNin, &b_Part_vProdNin);
  m_inputchain->SetBranchAddress("Part_vProdNout", Part_vProdNout, &b_Part_vProdNout);
  m_inputchain->SetBranchAddress("Part_vProdStatus", Part_vProdStatus, &b_Part_vProdStatus);
  m_inputchain->SetBranchAddress("Part_vProdBarcode", Part_vProdBarcode, &b_Part_vProdBarcode);
  m_inputchain->SetBranchAddress("Part_vParentID", &Part_vParentID, &b_Part_vParentID);
  m_inputchain->SetBranchAddress("Part_vParentBarcode", &Part_vParentBarcode, &b_Part_vParentBarcode);
  m_inputchain->SetBranchAddress("nSP", &nSP, &b_nSP);
  m_inputchain->SetBranchAddress("SPindex", SPindex, &b_SPindex);
  m_inputchain->SetBranchAddress("SPx", SPx, &b_SPx);
  m_inputchain->SetBranchAddress("SPy", SPy, &b_SPy);
  m_inputchain->SetBranchAddress("SPz", SPz, &b_SPz);
  m_inputchain->SetBranchAddress("SPCL1_index", SPCL1_index, &b_SPCL1_index);
  m_inputchain->SetBranchAddress("SPCL2_index", SPCL2_index, &b_SPCL2_index);
  m_inputchain->SetBranchAddress("SPisOverlap", SPisOverlap, &b_SPisOverlap);
  m_inputchain->SetBranchAddress("SPradius", SPradius, &b_SPradius);
  m_inputchain->SetBranchAddress("SPcovr", SPcovr, &b_SPcovr);
  m_inputchain->SetBranchAddress("SPcovz", SPcovz, &b_SPcovz);
  m_inputchain->SetBranchAddress("SPhl_topstrip", SPhl_topstrip, &b_SPhl_topstrip);
  m_inputchain->SetBranchAddress("SPhl_botstrip", SPhl_botstrip, &b_SPhl_botstrip);
  m_inputchain->SetBranchAddress("SPtopStripDirection",
				 &SPtopStripDirection, &b_SPtopStripDirection);
  m_inputchain->SetBranchAddress("SPbottomStripDirection",
				 &SPbottomStripDirection, &b_SPbottomStripDirection);
  m_inputchain->SetBranchAddress("SPstripCenterDistance",
				 &SPstripCenterDistance, &b_SPstripCenterDistance);
  m_inputchain->SetBranchAddress("SPtopStripCenterPosition",
				 &SPtopStripCenterPosition, &b_SPtopStripCenterPosition);
  

  
  m_inputchain->SetBranchAddress("nTRK", &nTRK, &b_nTRK);
  m_inputchain->SetBranchAddress("TRKindex", TRKindex, &b_TRKindex);
  m_inputchain->SetBranchAddress("TRKtrack_fitter", TRKtrack_fitter, &b_TRKtrack_fitter);
  m_inputchain->SetBranchAddress("TRKparticle_hypothesis", TRKparticle_hypothesis, &b_TRKparticle_hypothesis);
  m_inputchain->SetBranchAddress("TRKproperties", &TRKproperties, &b_TRKproperties);
  m_inputchain->SetBranchAddress("TRKpattern", &TRKpattern, &b_TRKpattern);
  m_inputchain->SetBranchAddress("TRKndof", TRKndof, &b_TRKndof);
  m_inputchain->SetBranchAddress("TRKmot", TRKmot, &b_TRKmot);
  m_inputchain->SetBranchAddress("TRKoot", TRKoot, &b_TRKoot);
  m_inputchain->SetBranchAddress("TRKchiSq", TRKchiSq, &b_TRKchiSq);
  m_inputchain->SetBranchAddress("TRKmeasurementsOnTrack_pixcl_sctcl_index", &TRKmeasurementsOnTrack_pixcl_sctcl_index, &b_TRKmeasurementsOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKoutliersOnTrack_pixcl_sctcl_index", &TRKoutliersOnTrack_pixcl_sctcl_index, &b_TRKoutliersOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKcharge", TRKcharge, &b_TRKcharge);
  m_inputchain->SetBranchAddress("TRKperigee_position", &TRKperigee_position, &b_TRKperigee_position);
  m_inputchain->SetBranchAddress("TRKperigee_momentum", &TRKperigee_momentum, &b_TRKperigee_momentum);
  m_inputchain->SetBranchAddress("TTCindex", TTCindex, &b_TTCindex);
  m_inputchain->SetBranchAddress("TTCevent_index", TTCevent_index, &b_TTCevent_index);
  m_inputchain->SetBranchAddress("TTCparticle_link", TTCparticle_link, &b_TTCparticle_link);
  m_inputchain->SetBranchAddress("TTCprobability", TTCprobability, &b_TTCprobability);
  m_inputchain->SetBranchAddress("nDTT", &nDTT, &b_nDTT);
  m_inputchain->SetBranchAddress("DTTindex", DTTindex, &b_DTTindex);
  m_inputchain->SetBranchAddress("DTTsize", DTTsize, &b_DTTsize);
  m_inputchain->SetBranchAddress("DTTtrajectory_eventindex", &DTTtrajectory_eventindex, &b_DTTtrajectory_eventindex);
  m_inputchain->SetBranchAddress("DTTtrajectory_barcode", &DTTtrajectory_barcode, &b_DTTtrajectory_barcode);
  m_inputchain->SetBranchAddress("DTTstTruth_subDetType", &DTTstTruth_subDetType, &b_DTTstTruth_subDetType);
  m_inputchain->SetBranchAddress("DTTstTrack_subDetType", &DTTstTrack_subDetType, &b_DTTstTrack_subDetType);
  m_inputchain->SetBranchAddress("DTTstCommon_subDetType", &DTTstCommon_subDetType, &b_DTTstCommon_subDetType);



  m_inputchain->Add(m_cfg.inputfile.c_str());
  ACTS_DEBUG("Adding file " <<m_cfg.inputfile <<" to tree"<< m_cfg.treename);
  

  m_events = m_inputchain->GetEntries();

  ACTS_DEBUG("End of constructor. In total available events=" << m_events);

  } // constructor


ActsExamples::ProcessCode ActsExamples::RootAthInputReader::read(const ActsExamples::AlgorithmContext& ctx){
  
  // Prepare containers for the hit data using the framework event data types
  //GeometryIdMultimap<Measurement> orderedMeasurements;
  //ClusterContainer clusters;
  //IndexMultimap<Index> measurementSimHitsMap;
  //IndexSourceLinkContainer sourceLinks;
  
  
  ACTS_DEBUG("Starting loop on events");
  auto entry = ctx.eventNumber;
  if (entry >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }
  
  std::lock_guard<std::mutex> lock(m_read_mutex);
  
  m_inputchain->GetEntry(entry);

  //Loop on clusters (measurements) 
  
  for (int im=0; im<5000; im++) {
    
    int bec      = CLbarrel_endcap[im];
    int lydisk   = CLlayer_disk   [im];
    int etamod   = CLeta_module   [im];
    int phimod   = CLphi_module   [im];
    int side     = CLside         [im];
    // ULong64_t moduleID = CLmoduleID     [im];
      
    ACTS_DEBUG(bec<<" "<<lydisk<<" "<<etamod<<" "<<phimod<<" "<<side<<" ");
  
  }
  
  
  // Prepare pixel space points
  SimSpacePointContainer pixelSpacePoints;

  // Prepare Strip space Points
  SimSpacePointContainer stripSpacePoints;

  // Prepare space-point container
  // They contain both pixel and SCT space points
  SimSpacePointContainer spacePoints;


  ACTS_DEBUG("Found "<< nSP <<" space points");
  
  // Bookkeeping of strip space points
  int stripsp = 0;

  //Loop on space points
  for (int isp=0; isp<nSP;isp++) {
    
    Acts::Vector3 globalPos{SPx[isp],SPy[isp],SPz[isp]};
    double sp_covr = SPcovr[isp];
    double sp_covz = SPcovz[isp];

    //PIX=1  STRIP = 2
    int type = SPCL2_index[isp] == -1 ?  1 : 2 ;
    
    ACTS_DEBUG("SP:: "<<type<< " " <<globalPos << " " << sp_covr<< " "<<sp_covz);

    
    boost::container::static_vector<Acts::SourceLink,2> sLinks;
    IndexSourceLink emptysl(0,isp);
    sLinks.emplace_back(emptysl);
  
    
    if (type == 1) { // pixel

      pixelSpacePoints.emplace_back(globalPos,std::nullopt,
        sp_covr, sp_covz, std::nullopt, sLinks);

    } else  { // Strip 

      
      float hl_topstrip = SPhl_topstrip[isp];
      float hl_botstrip = SPhl_botstrip[isp];

      ACTS_DEBUG("Checking SP:: "<<isp<< " strip sp:"<<stripsp);  
      
    
      Acts::Vector3 topStripDirection{
        ((*SPtopStripDirection)[stripsp])[0],
        ((*SPtopStripDirection)[stripsp])[1],
        ((*SPtopStripDirection)[stripsp])[2]};

      ACTS_DEBUG("Top Strip Direction "<<topStripDirection);

      
      Acts::Vector3 botStripDirection{
        ((*SPbottomStripDirection)[stripsp])[0],
        ((*SPbottomStripDirection)[stripsp])[1],
        ((*SPbottomStripDirection)[stripsp])[2]
      };

      ACTS_DEBUG("Bot Strip Direction "<<botStripDirection);
      
      Acts::Vector3 stripCenterDistance{
        ((*SPstripCenterDistance)[stripsp])[0],
        ((*SPstripCenterDistance)[stripsp])[1],
        ((*SPstripCenterDistance)[stripsp])[2]
      };


      
      Acts::Vector3 topStripCenterPosition{
        ((*SPtopStripCenterPosition)[stripsp])[0],
        ((*SPtopStripCenterPosition)[stripsp])[1],
        ((*SPtopStripCenterPosition)[stripsp])[2]
      };
      
      stripSpacePoints.emplace_back(globalPos, std::nullopt,
        sp_covr, sp_covz, std::nullopt,sLinks,
        hl_topstrip, hl_botstrip,
        topStripDirection, botStripDirection,
        stripCenterDistance, topStripCenterPosition
        );
      
      

      stripsp++;
    }
    
    spacePoints.emplace_back(globalPos, std::nullopt,
			    sp_covr, sp_covz, std::nullopt,sLinks);

    
  }
  

  ACTS_DEBUG("Created"<<pixelSpacePoints.size()<<" "<<
      m_cfg.outputPixelSpacePoints<< "space points");

  ACTS_DEBUG("Created"<<stripSpacePoints.size()<<" "<<
      m_cfg.outputStripSpacePoints<< "space points");

  ACTS_DEBUG("Created "<<spacePoints.size() <<" " <<
	     m_cfg.outputSpacePoints << "space points");
  
  m_outputPixelSpacePoints(ctx,std::move(pixelSpacePoints));
  m_outputStripSpacePoints(ctx,std::move(stripSpacePoints));
  m_outputSpacePoints(ctx, std::move(spacePoints));
  
  return ProcessCode::SUCCESS;
}

