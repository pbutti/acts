// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/SpacePointContainer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ostream>
#include "TFile.h"
#include "TTree.h"

namespace Acts {
template <std::size_t>
class GridBinFinder;
}  // namespace Acts

namespace ActsExamples {
struct AlgorithmContext;

/// Construct track seeds from space points.
class SeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    ///
    /// We allow multiple space point collections to allow different parts of
    /// the detector to use different algorithms for space point construction,
    /// e.g. single-hit space points for pixel-like detectors or double-hit
    /// space points for strip-like detectors.
    std::vector<std::string> inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;

    Acts::SeedFilterConfig seedFilterConfig;

    Acts::SeedFinderConfig<typename Acts::SpacePointContainer<
        ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
        Acts::detail::RefHolder>::SpacePointProxyType>
        seedFinderConfig;
    Acts::CylindricalSpacePointGridConfig gridConfig;
    Acts::CylindricalSpacePointGridOptions gridOptions;

    Acts::SeedFinderOptions seedFinderOptions;

    // allow for different values of rMax in gridConfig and seedFinderConfig
    bool allowSeparateRMax = false;

    // vector containing the map of z bins in the top and bottom layers
    std::vector<std::pair<int, int>> zBinNeighborsTop;
    std::vector<std::pair<int, int>> zBinNeighborsBottom;
    // number of phiBin neighbors at each side of the current bin that will be
    // used to search for SPs
    int numPhiNeighbors = 1;

    // Connect custom selections on the space points or to the doublet
    // compatibility
    bool useExtraCuts = false;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  using SpacePointProxy_t = typename Acts::SpacePointContainer<
      ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
      Acts::detail::RefHolder>::SpacePointProxyType;

  Acts::SeedFinder<SpacePointProxy_t,
                   Acts::CylindricalSpacePointGrid<SpacePointProxy_t>>
      m_seedFinder;
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_bottomBinFinder{nullptr};
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_topBinFinder{nullptr};

  Config m_cfg;

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  
  static inline Acts::simple_doubletMap loadModuleMap() {
    
    std::unique_ptr<TFile> mapfile = std::make_unique<TFile>("/mnt/hdd1/pibutti/data/GNN_ModuleMaps/ModuleMap_closest_all_dr_rel24_89809evts_double.doublets.tol1e-10.root");
  
    TTree* ttree = static_cast<TTree*>(mapfile->Get("TreeModuleDoublet"));
    
    double        detamax_12;
    double        detamin_12;
    double        dphimax_12;
    double        dphimin_12;
    double        phiSlopemax_12;
    double        phiSlopemin_12;
    double        z0max_12;
    double        z0min_12;
    uint64_t      Module1;
    uint64_t      Module2;
    double        z0sum_12;
    double        dphisum_12;
    double        phiSlopesum_12;
    double        detasum_12;
    double        z0sumSq_12;
    double        dphisumSq_12;
    double        phiSlopesumSq_12;
    double        detasumSq_12;
    
    
    TBranch        *b_detamax_12;
    TBranch        *b_detamin_12;
    TBranch        *b_dphimax_12;
    TBranch        *b_dphimin_12;
    TBranch        *b_phiSlopemax_12;
    TBranch        *b_phiSlopemin_12;
    TBranch        *b_z0max_12;
    TBranch        *b_z0min_12;
    TBranch        *b_Module1;
    TBranch        *b_Module2;
    TBranch        *b_z0sum_12;
    TBranch        *b_dphisum_12;
    TBranch        *b_phiSlopesum_12;
    TBranch        *b_detasum_12;
    TBranch        *b_z0sumSq_12;
    TBranch        *b_dphisumSq_12;
    TBranch        *b_phiSlopesumSq_12;
    TBranch        *b_detasumSq_12;
    
    ttree->SetBranchAddress("detamax_12", &detamax_12, &b_detamax_12);
    ttree->SetBranchAddress("detamin_12", &detamin_12, &b_detamin_12);
    ttree->SetBranchAddress("dphimax_12", &dphimax_12, &b_dphimax_12);
    ttree->SetBranchAddress("dphimin_12", &dphimin_12, &b_dphimin_12);
    ttree->SetBranchAddress("phiSlopemax_12", &phiSlopemax_12, &b_phiSlopemax_12);
    ttree->SetBranchAddress("phiSlopemin_12", &phiSlopemin_12, &b_phiSlopemin_12);
    ttree->SetBranchAddress("z0max_12", &z0max_12, &b_z0max_12);
    ttree->SetBranchAddress("z0min_12", &z0min_12, &b_z0min_12);
    ttree->SetBranchAddress("Module1", &Module1, &b_Module1);
    ttree->SetBranchAddress("Module2", &Module2, &b_Module2);
    ttree->SetBranchAddress("z0sum_12", &z0sum_12, &b_z0sum_12);
    ttree->SetBranchAddress("dphisum_12", &dphisum_12, &b_dphisum_12);
    ttree->SetBranchAddress("phiSlopesum_12", &phiSlopesum_12, &b_phiSlopesum_12);
    ttree->SetBranchAddress("detasum_12", &detasum_12, &b_detasum_12);
    ttree->SetBranchAddress("z0sumSq_12", &z0sumSq_12, &b_z0sumSq_12);
    ttree->SetBranchAddress("dphisumSq_12", &dphisumSq_12, &b_dphisumSq_12);
    ttree->SetBranchAddress("phiSlopesumSq_12", &phiSlopesumSq_12, &b_phiSlopesumSq_12);
    ttree->SetBranchAddress("detasumSq_12", &detasumSq_12, &b_detasumSq_12);
    
    
    Long64_t nentries = ttree->GetEntriesFast();

    Acts::simple_doubletMap doubletmap;
    Acts::commutativePairHash hasher;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      ttree->GetEntry(jentry);
      
      //std::cout<<Module1<<" "<<Module2<<" hash="<<hasher(std::make_pair(Module1,Module2))<<std::endl;
      uint64_t hash = hasher(std::make_pair(Module1,Module2));
      doubletmap[hash] = true;
      
    }
    
    std::cout<<"PF: SeedingAlgorithm--Size of doublet map "<<doubletmap.size()<<std::endl;
    return doubletmap;
    
  }
  
  static inline bool itkFastTrackingCuts(float bottomRadius, float cotTheta) {
    static float rMin = 45.;
    static float cotThetaMax = 1.5;

    if (bottomRadius < rMin &&
        (cotTheta > cotThetaMax || cotTheta < -cotThetaMax)) {
      return false;
    }
    return true;
  }

  static inline bool itkFastTrackingSPselect(const SpacePointProxy_t& sp) {
    // At small r we remove points beyond |z| > 200.
    float r = sp.radius();
    float zabs = std::abs(sp.z());
    if (zabs > 200. && r < 45.) {
      return false;
    }

    /// Remove space points beyond eta=4 if their z is
    /// larger than the max seed z0 (150.)
    float cotTheta = 27.2899;  // corresponds to eta=4
    if ((zabs - 150.) > cotTheta * r) {
      return false;
    }
    return true;
  }
};

}  // namespace ActsExamples
