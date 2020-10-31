#include "PrintLdmxSpacePoints.hpp"
#include "Acts/LdmxEDM/LdmxSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <vector>

ActsExamples::PrintLdmxSpacePoints::PrintLdmxSpacePoints(const ActsExamples::PrintLdmxSpacePoints::Config& cfg,
                                                         Acts::Logging::Level level) 
    : BareAlgorithm("PrintLdmxSpacePoints", level), m_cfg(cfg){}

ActsExamples::ProcessCode ActsExamples::PrintLdmxSpacePoints::execute(const ActsExamples::AlgorithmContext& ctx) const { 
    using ldmxsps = std::vector<LdmxSpacePoint>; 
    const auto& ldmxspc = ctx.eventStore.get<ldmxsps>(m_cfg.inputCollection);
    
    std::cout<<"Printing:: " << m_cfg.inputCollection<< "with size " << ldmxspc.size() << std::endl;
     
    for (size_t ih = 0; ih < ldmxspc.size(); ih++) {
        
        LdmxSpacePoint ldmxsp = ldmxspc[ih];
        //ACTS_INFO("Ldmx Sp " << ldmxsp.id()<< " x=" <<ldmxsp.x()<< " y=" <<ldmxsp.y()<< " z=" <<ldmxsp.z()<< " \n" << " t="<< ldmxsp.t()<<" ly="<<ldmxsp.layer());

        std::cout<<"Ldmx Sp " << ldmxsp.id()<< " x=" <<ldmxsp.x()<< " y=" <<ldmxsp.y()<< " z=" <<ldmxsp.z()<< " \n" << " t="<< ldmxsp.t()<<" ly="<<ldmxsp.layer() <<std::endl;
    }

    return ProcessCode::SUCCESS;

}
    
