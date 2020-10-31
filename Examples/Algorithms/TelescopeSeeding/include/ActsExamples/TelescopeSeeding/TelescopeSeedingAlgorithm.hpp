#pragma once

#include "Acts/Seeding/Seed.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "Acts/LdmxEDM/LdmxSpacePoint.hpp"

#include "TH1F.h"
#include "TH2F.h"

#include <memory>

namespace ActsExamples {
    
    class TelescopeSeedingAlgorithm final : public BareAlgorithm {
        
    public:
        
        struct Config {
            /// Output collection of clusters.
            std::string outputSeeds;
            /// Output prototracks.
            std::string outputProtoTracks;
            /// input space points.
            std::string inputSpacePoints;
        }; ///end configuration
        
        TelescopeSeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);
        
        ProcessCode execute(const AlgorithmContext& ctx) const final override;
        
    private:
        Config m_cfg;
        
        //TH1F* h_p;
        //TH1F* h_tanTheta;
        //TH1F* h_d0;
        
    };
} // namespace ActsExamples
