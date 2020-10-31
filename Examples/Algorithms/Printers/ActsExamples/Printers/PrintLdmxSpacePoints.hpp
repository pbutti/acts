#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

namespace ActsExamples {

    /// Print hits within some geometric region-of-interest.
    class PrintLdmxSpacePoints : public BareAlgorithm {
    public:
        struct Config {
            /// Input collection.
            std::string inputCollection;
        };

        PrintLdmxSpacePoints(const Config& cfg, Acts::Logging::Level level);

        ProcessCode execute(const AlgorithmContext& ctx) const;

    private:
        Config m_cfg;
    };

}  // namespace ActsExamples
