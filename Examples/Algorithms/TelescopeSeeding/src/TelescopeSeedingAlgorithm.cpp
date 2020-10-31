#include "ActsExamples/TelescopeSeeding/TelescopeSeedingAlgorithm.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "Acts/Seeding/SeedToTrackParamMaker.hpp"

ActsExamples::TelescopeSeedingAlgorithm::TelescopeSeedingAlgorithm(
    ActsExamples::TelescopeSeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("TelescopeSeedingAlgorithm",lvl),
      m_cfg(std::move(cfg)) {
    
    ACTS_INFO("Constructed TelescopeSeeding Algorithm");
    ACTS_INFO("Input SpacePoint collection:" << m_cfg.inputSpacePoints);
    

}

ActsExamples::ProcessCode ActsExamples::TelescopeSeedingAlgorithm::execute( 
    const AlgorithmContext& ctx) const {
    
    //Configure the tools --- is it necessary to do the configuration in execute?

    Acts::SeedfinderConfig<LdmxSpacePoint> config;
    
    //Tagger r max
    config.rMax = 1000.;
    config.deltaRMin = 3.; 
    config.deltaRMax = 220.; 
    config.collisionRegionMin = -50; 
    config.collisionRegionMax =  50; 
    config.zMin = -300;
    config.zMax = 300.;
    config.maxSeedsPerSpM = 5;
    
    //More or less the max angle is something of the order of 50 / 600 (assuming that the track hits all the layers)
    //Theta for the seeder is like ATLAS eta, so it's 90-lambda.
    //Max lamba is of the order of ~0.1 so cotThetaMax will be 1./tan(pi/2 - 0.1) ~ 1.4. 
    config.cotThetaMax = 1.5;
    
    //cotThetaMax and deltaRMax matter to choose the binning in z. The bin size is given by cotThetaMax*deltaRMax
    
    config.sigmaScattering = 2.25;
    config.minPt = 500.;
    config.bFieldInZ = 1.5e-3;  // in kT
    config.beamPos = {0, 0}; // units?
    config.impactMax = 20.;
    
    //setup the spacepoint grid configuration 
    
    // setup spacepoint grid config
    Acts::SpacePointGridConfig gridConf;
    gridConf.bFieldInZ = config.bFieldInZ;
    gridConf.minPt = config.minPt;
    gridConf.rMax = config.rMax;
    gridConf.zMax = config.zMax;
    gridConf.zMin = config.zMin;
    gridConf.deltaRMax = config.deltaRMax;
    gridConf.cotThetaMax = config.cotThetaMax;
    
    
    auto bottomBinFinder = std::make_shared<Acts::BinFinder<LdmxSpacePoint> > (
        Acts::BinFinder<LdmxSpacePoint>());

    auto topBinFinder = std::make_shared<Acts::BinFinder<LdmxSpacePoint> > (
        Acts::BinFinder<LdmxSpacePoint>());

    //The seed finder needs a seed filter instance
    
    //In the seed finder there is the correction for the beam axis (?), which you could ignore if you set the penalty for high impact parameters. So removing that in the seeder config.
    
    Acts::SeedFilterConfig seedFilter_cfg;
    seedFilter_cfg.impactWeightFactor = 0.;
    
    //For the moment no experiment dependent cuts are assigned to the filter
    config.seedFilter = std::make_unique<Acts::SeedFilter<LdmxSpacePoint>>(
        Acts::SeedFilter<LdmxSpacePoint>(seedFilter_cfg));
    
    Acts::Seedfinder<LdmxSpacePoint> seedFinder(config);
    
    
    //Retrieve the SpacePoints
    using ldmx_sps = std::vector<LdmxSpacePoint>;
    const auto& allSps = ctx.eventStore.get<ldmx_sps>(m_cfg.inputSpacePoints);
    
    //Select the ones to be used for seeding
    //I'll use layer 3,5,7
    
    std::vector<const LdmxSpacePoint*> spVec;

    for (size_t isp = 0; isp < allSps.size(); isp++) {
        
        if (allSps[isp].layer()==3 || allSps[isp].layer()==7 || allSps[isp].layer()==9) 
            spVec.push_back(&(allSps[isp]));
    }
    
    //crosscheck
    for(size_t isp = 0; isp < spVec.size(); isp++) {
        std::cout<<"spVec ["<<isp<<"] layer "<<spVec[isp]->layer()<<std::endl;
        std::cout<<"(x,y,z) = ("<<spVec[isp]->x()<<","<<spVec[isp]->y()<<","<<spVec[isp]->z()<<")"<<std::endl;
    }
    
    // create grid with bin sizes according to the configured geometry
    std::unique_ptr<Acts::SpacePointGrid<LdmxSpacePoint>> grid =
        Acts::SpacePointGridCreator::createGrid<LdmxSpacePoint>(gridConf);
    
    ACTS_DEBUG("Formed grid");
    
    // covariance tool, sets covariances per spacepoint as required  --- for the moment 0 covariance
    auto ct = [=](const LdmxSpacePoint& sp, float, float,
                  float) -> Acts::Vector2D {
        return {0., 0.};
    };

    ACTS_DEBUG("Defined covariance tool");

    // create the space point group
    auto spGroup = Acts::BinnedSPGroup<LdmxSpacePoint>(
        spVec.begin(), spVec.end(), ct, bottomBinFinder, topBinFinder,
        std::move(grid), config);

    ACTS_DEBUG("BinnedSPGroup defined");
    
    // seed vector
    std::vector<std::vector<Acts::Seed<LdmxSpacePoint>>> seedVector;
    
    // find the seeds
    auto groupIt = spGroup.begin();
    auto endOfGroups = spGroup.end();
    
    for (; !(groupIt == endOfGroups); ++groupIt) {
        
        seedVector.push_back(seedFinder.createSeedsForGroup(
                                 groupIt.bottom(), groupIt.middle(), groupIt.top()));
    }
        
    
    //The seed to track parameters fitter
    Acts::SeedToTrackParamMaker seedToTrackMaker;
    
    int numSeeds = 0;
    for (auto& outVec : seedVector) {
        numSeeds += outVec.size();
        
        for (size_t iSeed = 0; iSeed < outVec.size(); iSeed++) {
            std::array<double,9> outData;
            Acts::Transform3D tP;
            
            tP.setIdentity();
            tP(0,0)=0;
            tP(0,1)=0;
            tP(0,2)=1;

            tP(1,0)=0;
            tP(1,1)=1;
            tP(1,2)=0;

            tP(2,0)=1;
            tP(2,1)=0;
            tP(2,2)=0;
            
            tP.translation().x() = -213.1;
            tP.translation().y() = 0.;
            tP.translation().z() = 0.;
            std::cout<<"PF:: CHECK CHECK \n"<<tP.translation()<<std::endl;
            seedToTrackMaker.FitSeedAtlas(outVec[iSeed], outData, tP,0.0015);
            std::cout<<outData[0]<<" "<<outData[1]<<" "<<outData[2]<<" "<<outData[3]<<" "<<outData[4]<<std::endl;
            
            seedToTrackMaker.KarimakiFit(outVec[iSeed].sp(),outData);
            //h_p->Fill(outData[4]);
            
            
        }
    }
    
    std::cout<<spVec.size() << " hits, " << seedVector.size() << " regions, "
             << numSeeds << " seeds"<<std::endl;
    
    
    return ActsExamples::ProcessCode::SUCCESS;
}




