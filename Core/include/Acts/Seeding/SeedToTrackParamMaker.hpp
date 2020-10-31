#pragma once

#include "Acts/Seeding/Seed.hpp"

namespace Acts {
    
    class SeedToTrackParamMaker {
    public:
        
        SeedToTrackParamMaker(){};
        
        template <typename external_spacepoint_t>
        bool KarimakiFit(const std::vector<external_spacepoint_t*>&sp, std::array<double,9>& data);

        /// This resembles the method used in ATLAS for the seed fitting
        /// L811 https://acode-browser.usatlas.bnl.gov/lxr/source/athena/InnerDetector/InDetRecTools/SiTrackMakerTool_xk/src/SiTrackMaker_xk.cxx
        template <typename external_spacepoint_t>
        bool FitSeedAtlas(const Seed<external_spacepoint_t>& seed, std::array<double, 9>& data, const Transform3D& Tp, const double& bFieldZ);
        
        template <typename external_spacepoint_t>
        bool FitSeedAtlas(const std::vector<external_spacepoint_t>& sp, std::array<double,9>& data, const Transform3D& Tp, const double& bFieldZ);
        
        /// This is a simple Line and Circle fit based on Taubin Circle fit
        template <typename external_spacepoint_t>
        bool FitSeedLinCircle(const Seed<external_spacepoint_t>& seed, std::vector<double>& data);
        
        
        /// This is a simple Line and Parabola fit (from HPS reconstruction by Robert Johnson)
        template <typename external_spacepoint_t>
        bool FitSeedLinPar(const Seed<external_spacepoint_t>& seed, std::vector<double>& data);
        
    private:
        
        double m_pTmin;
        
        
        template <typename external_spacepoint_t>
        bool LinearParabolaFit(std::vector<external_spacepoint_t>& sPs);

    };
    
}

#include "Acts/Seeding/SeedToTrackParamMaker.ipp"
