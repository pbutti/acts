#pragma once


#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include "TString.h"
#include "TArrayD.h"

#include <mutex>
#include <vector>

using namespace std;

class TChain;
class TBranch;


namespace ActsExamples {

    /// @class RootVertexAndTracksReader
    ///
    /// @brief Reads in ldmx hit information from a root file
    /// and fills it in the event store

    class RootLdmxSimHitReader final : public IReader {

        /// @brief The nested configuration struct
    public:

        struct Config {
            std::string outputCollection = "LdmxSimHitCollection";
            std::string treeName = "LDMX_Events";     ///< name of the output tree
            std::vector<std::string> fileList;  ///< The name of the input file
            unsigned int batchSize = 1;         ///!< Batch
        };


        /// Constructor
        /// @param cfg The Configuration struct
        /// @param lvl Message level declaration
        RootLdmxSimHitReader(Config cfg, Acts::Logging::Level lvl);

        /// Destructor
        ~RootLdmxSimHitReader() final override;

        /// Framework name() method
        std::string name() const final override;

        /// Return the available events range.
        std::pair<size_t, size_t> availableEvents() const final override;

        /// Read out data from the input stream
        ///
        /// @param context The algorithm context
        ProcessCode read(
            const ActsExamples::AlgorithmContext& context) final override;

    private:
        /// The config class;
        Config m_cfg;
        /// mutex used to protect multi-threaded reads
        std::mutex m_read_mutex;
        /// The number of events
        size_t m_events = 0;
        /// The input tree
        TChain* fChain = nullptr;

        int m_eventNr = 0;

        /*
        
        // Fixed size dimensions of array or collections stored in the TTree if any.
        static constexpr int kMaxeventNumber = 1;
        static constexpr int kMaxrun = 1;
        static constexpr int kMaxweight = 1;
        static constexpr int kMaxisRealData = 1;
        static constexpr int kMaxintParameters_ = 1;
        static constexpr int kMaxfloatParameters_ = 1;
        static constexpr int kMaxstringParameters_ = 1;
        static constexpr int kMaxSimParticles_v12 = 92;
        static constexpr int kMaxSimParticles_v12_second_energy = 92;
        static constexpr int kMaxSimParticles_v12_second_pdgID = 92;
        static constexpr int kMaxSimParticles_v12_second_genStatus = 92;
        static constexpr int kMaxSimParticles_v12_second_time = 92;
        static constexpr int kMaxSimParticles_v12_second_x = 92;
        static constexpr int kMaxSimParticles_v12_second_y = 92;
        static constexpr int kMaxSimParticles_v12_second_z = 92;
        static constexpr int kMaxSimParticles_v12_second_endX = 92;
        static constexpr int kMaxSimParticles_v12_second_endY = 92;
        static constexpr int kMaxSimParticles_v12_second_endZ = 92;
        static constexpr int kMaxSimParticles_v12_second_px = 92;
        static constexpr int kMaxSimParticles_v12_second_py = 92;
        static constexpr int kMaxSimParticles_v12_second_pz = 92;
        static constexpr int kMaxSimParticles_v12_second_endpx = 92;
        static constexpr int kMaxSimParticles_v12_second_endpy = 92;
        static constexpr int kMaxSimParticles_v12_second_endpz = 92;
        static constexpr int kMaxSimParticles_v12_second_mass = 92;
        static constexpr int kMaxSimParticles_v12_second_charge = 92;
        static constexpr int kMaxSimParticles_v12_second_daughters = 92;
        static constexpr int kMaxSimParticles_v12_second_parents = 92;
        static constexpr int kMaxSimParticles_v12_second_processType = 92;
        static constexpr int kMaxSimParticles_v12_second_vertexVolume = 92;
        static constexpr int kMaxTaggerSimHits_v12 = 177;
        static constexpr int kMaxTaggerSimHits_v12_id = 177;
        static constexpr int kMaxTaggerSimHits_v12_layerID = 177;
        static constexpr int kMaxTaggerSimHits_v12_moduleID = 177;
        static constexpr int kMaxTaggerSimHits_v12_edep = 177;
        static constexpr int kMaxTaggerSimHits_v12_time = 177;
        static constexpr int kMaxTaggerSimHits_v12_px = 177;
        static constexpr int kMaxTaggerSimHits_v12_py = 177;
        static constexpr int kMaxTaggerSimHits_v12_pz = 177;
        static constexpr int kMaxTaggerSimHits_v12_energy = 177;
        static constexpr int kMaxTaggerSimHits_v12_x = 177;
        static constexpr int kMaxTaggerSimHits_v12_y = 177;
        static constexpr int kMaxTaggerSimHits_v12_z = 177;
        static constexpr int kMaxTaggerSimHits_v12_pathLength = 177;
        static constexpr int kMaxTaggerSimHits_v12_trackID = 177;
        static constexpr int kMaxTaggerSimHits_v12_pdgID = 177;
        static constexpr int kMaxTriggerPadUpSimHits_v12 = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_id = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_edep = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_x = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_y = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_z = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_time = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_trackIDContribs = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_incidentIDContribs = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_pdgCodeContribs = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_edepContribs = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_timeContribs = 69;
        static constexpr int kMaxTriggerPadUpSimHits_v12_nContribs = 69;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12 = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_id = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_edep = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_x = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_y = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_z = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_time = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_trackIDContribs = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_incidentIDContribs = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_pdgCodeContribs = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_edepContribs = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_timeContribs = 73;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_nContribs = 73;
        static constexpr int kMaxTargetSimHits_v12 = 76;
        static constexpr int kMaxTargetSimHits_v12_id = 76;
        static constexpr int kMaxTargetSimHits_v12_edep = 76;
        static constexpr int kMaxTargetSimHits_v12_x = 76;
        static constexpr int kMaxTargetSimHits_v12_y = 76;
        static constexpr int kMaxTargetSimHits_v12_z = 76;
        static constexpr int kMaxTargetSimHits_v12_time = 76;
        static constexpr int kMaxTargetSimHits_v12_trackIDContribs = 76;
        static constexpr int kMaxTargetSimHits_v12_incidentIDContribs = 76;
        static constexpr int kMaxTargetSimHits_v12_pdgCodeContribs = 76;
        static constexpr int kMaxTargetSimHits_v12_edepContribs = 76;
        static constexpr int kMaxTargetSimHits_v12_timeContribs = 76;
        static constexpr int kMaxTargetSimHits_v12_nContribs = 76;
        static constexpr int kMaxTriggerPadDownSimHits_v12 = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_id = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_edep = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_x = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_y = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_z = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_time = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_trackIDContribs = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_incidentIDContribs = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_pdgCodeContribs = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_edepContribs = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_timeContribs = 70;
        static constexpr int kMaxTriggerPadDownSimHits_v12_nContribs = 70;
        static constexpr int kMaxRecoilSimHits_v12 = 184;
        static constexpr int kMaxRecoilSimHits_v12_id = 184;
        static constexpr int kMaxRecoilSimHits_v12_layerID = 184;
        static constexpr int kMaxRecoilSimHits_v12_moduleID = 184;
        static constexpr int kMaxRecoilSimHits_v12_edep = 184;
        static constexpr int kMaxRecoilSimHits_v12_time = 184;
        static constexpr int kMaxRecoilSimHits_v12_px = 184;
        static constexpr int kMaxRecoilSimHits_v12_py = 184;
        static constexpr int kMaxRecoilSimHits_v12_pz = 184;
        static constexpr int kMaxRecoilSimHits_v12_energy = 184;
        static constexpr int kMaxRecoilSimHits_v12_x = 184;
        static constexpr int kMaxRecoilSimHits_v12_y = 184;
        static constexpr int kMaxRecoilSimHits_v12_z = 184;
        static constexpr int kMaxRecoilSimHits_v12_pathLength = 184;
        static constexpr int kMaxRecoilSimHits_v12_trackID = 184;
        static constexpr int kMaxRecoilSimHits_v12_pdgID = 184;
        static constexpr int kMaxEcalSimHits_v12 = 162;
        static constexpr int kMaxEcalSimHits_v12_id = 162;
        static constexpr int kMaxEcalSimHits_v12_edep = 162;
        static constexpr int kMaxEcalSimHits_v12_x = 162;
        static constexpr int kMaxEcalSimHits_v12_y = 162;
        static constexpr int kMaxEcalSimHits_v12_z = 162;
        static constexpr int kMaxEcalSimHits_v12_time = 162;
        static constexpr int kMaxEcalSimHits_v12_trackIDContribs = 162;
        static constexpr int kMaxEcalSimHits_v12_incidentIDContribs = 162;
        static constexpr int kMaxEcalSimHits_v12_pdgCodeContribs = 162;
        static constexpr int kMaxEcalSimHits_v12_edepContribs = 162;
        static constexpr int kMaxEcalSimHits_v12_timeContribs = 162;
        static constexpr int kMaxEcalSimHits_v12_nContribs = 162;
        static constexpr int kMaxHcalSimHits_v12 = 1862;
        static constexpr int kMaxHcalSimHits_v12_id = 1862;
        static constexpr int kMaxHcalSimHits_v12_edep = 1862;
        static constexpr int kMaxHcalSimHits_v12_x = 1862;
        static constexpr int kMaxHcalSimHits_v12_y = 1862;
        static constexpr int kMaxHcalSimHits_v12_z = 1862;
        static constexpr int kMaxHcalSimHits_v12_time = 1862;
        static constexpr int kMaxHcalSimHits_v12_trackIDContribs = 1862;
        static constexpr int kMaxHcalSimHits_v12_incidentIDContribs = 1862;
        static constexpr int kMaxHcalSimHits_v12_pdgCodeContribs = 1862;
        static constexpr int kMaxHcalSimHits_v12_edepContribs = 1862;
        static constexpr int kMaxHcalSimHits_v12_timeContribs = 1862;
        static constexpr int kMaxHcalSimHits_v12_nContribs = 1862;
        static constexpr int kMaxMagnetScoringPlaneHits_v12 = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_id = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_layerID = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_moduleID = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_edep = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_time = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_px = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_py = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_pz = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_energy = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_x = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_y = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_z = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_pathLength = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_trackID = 3208;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_pdgID = 3208;
        static constexpr int kMaxTargetScoringPlaneHits_v12 = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_id = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_layerID = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_moduleID = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_edep = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_time = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_px = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_py = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_pz = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_energy = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_x = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_y = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_z = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_pathLength = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_trackID = 3351;
        static constexpr int kMaxTargetScoringPlaneHits_v12_pdgID = 3351;
        static constexpr int kMaxTrackerScoringPlaneHits_v12 = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_id = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_layerID = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_moduleID = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_edep = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_time = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_px = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_py = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_pz = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_energy = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_x = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_y = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_z = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_pathLength = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_trackID = 2730;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_pdgID = 2730;
        static constexpr int kMaxEcalScoringPlaneHits_v12 = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_id = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_layerID = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_moduleID = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_edep = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_time = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_px = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_py = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_pz = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_energy = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_x = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_y = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_z = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_pathLength = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_trackID = 391;
        static constexpr int kMaxEcalScoringPlaneHits_v12_pdgID = 391;
        static constexpr int kMaxHcalScoringPlaneHits_v12 = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_id = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_layerID = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_moduleID = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_edep = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_time = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_px = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_py = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_pz = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_energy = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_x = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_y = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_z = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_pathLength = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_trackID = 3430;
        static constexpr int kMaxHcalScoringPlaneHits_v12_pdgID = 3430;

        // Declaration of leaf types
        //ldmx::EventHeader *EventHeader;
        int           eventNumber_;
        int           run_;
        int           timestamp__fSec;
        int           timestamp__fNanoSec;
        double        weight_;
        bool          isRealData_;
        int           intParameters__;
        std::string          intParameters__first[kMaxintParameters_];
        int           intParameters__second[kMaxintParameters_];   //[intParameters__]
        int           floatParameters__;
        std::string          floatParameters__first[kMaxfloatParameters_];
        float         floatParameters__second[kMaxfloatParameters_];   //[floatParameters__]
        int           stringParameters__;
        std::string          stringParameters__first[kMaxstringParameters_];
        std::string          stringParameters__second[kMaxstringParameters_];
        int           SimParticles_v12_;
        int           SimParticles_v12_first[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_energy_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        int           SimParticles_v12_second_pdgID_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        int           SimParticles_v12_second_genStatus_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_time_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_x_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_y_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_z_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endX_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endY_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endZ_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_px_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_py_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_pz_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endpx_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endpy_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endpz_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_mass_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_charge_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        vector<int>     SimParticles_v12_second_daughters_[kMaxSimParticles_v12];
        vector<int>     SimParticles_v12_second_parents_[kMaxSimParticles_v12];
        int           SimParticles_v12_second_processType_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        std::string          SimParticles_v12_second_vertexVolume_[kMaxSimParticles_v12];
        int           TaggerSimHits_v12_;
        int           TaggerSimHits_v12_id_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_layerID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_moduleID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_edep_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_time_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_px_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_py_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_pz_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_energy_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_x_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_y_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_z_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_pathLength_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_trackID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_pdgID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TriggerPadUpSimHits_v12_;
        int           TriggerPadUpSimHits_v12_id_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_edep_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_x_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_y_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_z_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_time_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        vector<int>     TriggerPadUpSimHits_v12_trackIDContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<int>     TriggerPadUpSimHits_v12_incidentIDContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<int>     TriggerPadUpSimHits_v12_pdgCodeContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<float>   TriggerPadUpSimHits_v12_edepContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<float>   TriggerPadUpSimHits_v12_timeContribs_[kMaxTriggerPadUpSimHits_v12];
        unsigned int          TriggerPadUpSimHits_v12_nContribs_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        int           TriggerPadTaggerSimHits_v12_;
        int           TriggerPadTaggerSimHits_v12_id_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_edep_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_x_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_y_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_z_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_time_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        vector<int>     TriggerPadTaggerSimHits_v12_trackIDContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<int>     TriggerPadTaggerSimHits_v12_incidentIDContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<int>     TriggerPadTaggerSimHits_v12_pdgCodeContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<float>   TriggerPadTaggerSimHits_v12_edepContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<float>   TriggerPadTaggerSimHits_v12_timeContribs_[kMaxTriggerPadTaggerSimHits_v12];
        unsigned int          TriggerPadTaggerSimHits_v12_nContribs_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        int           TargetSimHits_v12_;
        int           TargetSimHits_v12_id_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_edep_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_x_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_y_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_z_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_time_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        vector<int>     TargetSimHits_v12_trackIDContribs_[kMaxTargetSimHits_v12];
        vector<int>     TargetSimHits_v12_incidentIDContribs_[kMaxTargetSimHits_v12];
        vector<int>     TargetSimHits_v12_pdgCodeContribs_[kMaxTargetSimHits_v12];
        vector<float>   TargetSimHits_v12_edepContribs_[kMaxTargetSimHits_v12];
        vector<float>   TargetSimHits_v12_timeContribs_[kMaxTargetSimHits_v12];
        unsigned int          TargetSimHits_v12_nContribs_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        int           TriggerPadDownSimHits_v12_;
        int           TriggerPadDownSimHits_v12_id_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_edep_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_x_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_y_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_z_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_time_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        vector<int>     TriggerPadDownSimHits_v12_trackIDContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<int>     TriggerPadDownSimHits_v12_incidentIDContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<int>     TriggerPadDownSimHits_v12_pdgCodeContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<float>   TriggerPadDownSimHits_v12_edepContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<float>   TriggerPadDownSimHits_v12_timeContribs_[kMaxTriggerPadDownSimHits_v12];
        unsigned int          TriggerPadDownSimHits_v12_nContribs_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        int           RecoilSimHits_v12_;
        int           RecoilSimHits_v12_id_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_layerID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_moduleID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_edep_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_time_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_px_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_py_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_pz_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_energy_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_x_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_y_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_z_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_pathLength_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_trackID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_pdgID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           EcalSimHits_v12_;
        int           EcalSimHits_v12_id_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_edep_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_x_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_y_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_z_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_time_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        vector<int>     EcalSimHits_v12_trackIDContribs_[kMaxEcalSimHits_v12];
        vector<int>     EcalSimHits_v12_incidentIDContribs_[kMaxEcalSimHits_v12];
        vector<int>     EcalSimHits_v12_pdgCodeContribs_[kMaxEcalSimHits_v12];
        vector<float>   EcalSimHits_v12_edepContribs_[kMaxEcalSimHits_v12];
        vector<float>   EcalSimHits_v12_timeContribs_[kMaxEcalSimHits_v12];
        unsigned int          EcalSimHits_v12_nContribs_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        int           HcalSimHits_v12_;
        int           HcalSimHits_v12_id_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_edep_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_x_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_y_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_z_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_time_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        vector<int>     HcalSimHits_v12_trackIDContribs_[kMaxHcalSimHits_v12];
        vector<int>     HcalSimHits_v12_incidentIDContribs_[kMaxHcalSimHits_v12];
        vector<int>     HcalSimHits_v12_pdgCodeContribs_[kMaxHcalSimHits_v12];
        vector<float>   HcalSimHits_v12_edepContribs_[kMaxHcalSimHits_v12];
        vector<float>   HcalSimHits_v12_timeContribs_[kMaxHcalSimHits_v12];
        unsigned int          HcalSimHits_v12_nContribs_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        int           MagnetScoringPlaneHits_v12_;
        int           MagnetScoringPlaneHits_v12_id_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_layerID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_moduleID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_edep_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_time_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_px_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_py_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_pz_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_energy_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_x_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_y_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_z_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_pathLength_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_trackID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_pdgID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_;
        int           TargetScoringPlaneHits_v12_id_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_layerID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_moduleID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_edep_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_time_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_px_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_py_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_pz_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_energy_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_x_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_y_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_z_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_pathLength_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_trackID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_pdgID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_;
        int           TrackerScoringPlaneHits_v12_id_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_layerID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_moduleID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_edep_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_time_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_px_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_py_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_pz_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_energy_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_x_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_y_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_z_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_pathLength_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_trackID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_pdgID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_;
        int           EcalScoringPlaneHits_v12_id_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_layerID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_moduleID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_edep_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_time_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_px_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_py_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_pz_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_energy_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_x_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_y_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_z_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_pathLength_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_trackID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_pdgID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_;
        int           HcalScoringPlaneHits_v12_id_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_layerID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_moduleID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_edep_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_time_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_px_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_py_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_pz_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_energy_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_x_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_y_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_z_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_pathLength_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_trackID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_pdgID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]

        // List of branches
        TBranch        *b_EventHeader_eventNumber_;   //!
        TBranch        *b_EventHeader_run_;   //!
        TBranch        *b_EventHeader_timestamp__fSec;   //!
        TBranch        *b_EventHeader_timestamp__fNanoSec;   //!
        TBranch        *b_EventHeader_weight_;   //!
        TBranch        *b_EventHeader_isRealData_;   //!
        TBranch        *b_EventHeader_intParameters__;   //!
        TBranch        *b_intParameters__first;   //!
        TBranch        *b_intParameters__second;   //!
        TBranch        *b_EventHeader_floatParameters__;   //!
        TBranch        *b_floatParameters__first;   //!
        TBranch        *b_floatParameters__second;   //!
        TBranch        *b_EventHeader_stringParameters__;   //!
        TBranch        *b_stringParameters__first;   //!
        TBranch        *b_stringParameters__second;   //!
        TBranch        *b_SimParticles_v12_;   //!
        TBranch        *b_SimParticles_v12_first;   //!
        TBranch        *b_SimParticles_v12_second_energy_;   //!
        TBranch        *b_SimParticles_v12_second_pdgID_;   //!
        TBranch        *b_SimParticles_v12_second_genStatus_;   //!
        TBranch        *b_SimParticles_v12_second_time_;   //!
        TBranch        *b_SimParticles_v12_second_x_;   //!
        TBranch        *b_SimParticles_v12_second_y_;   //!
        TBranch        *b_SimParticles_v12_second_z_;   //!
        TBranch        *b_SimParticles_v12_second_endX_;   //!
        TBranch        *b_SimParticles_v12_second_endY_;   //!
        TBranch        *b_SimParticles_v12_second_endZ_;   //!
        TBranch        *b_SimParticles_v12_second_px_;   //!
        TBranch        *b_SimParticles_v12_second_py_;   //!
        TBranch        *b_SimParticles_v12_second_pz_;   //!
        TBranch        *b_SimParticles_v12_second_endpx_;   //!
        TBranch        *b_SimParticles_v12_second_endpy_;   //!
        TBranch        *b_SimParticles_v12_second_endpz_;   //!
        TBranch        *b_SimParticles_v12_second_mass_;   //!
        TBranch        *b_SimParticles_v12_second_charge_;   //!
        TBranch        *b_SimParticles_v12_second_daughters_;   //!
        TBranch        *b_SimParticles_v12_second_parents_;   //!
        TBranch        *b_SimParticles_v12_second_processType_;   //!
        TBranch        *b_SimParticles_v12_second_vertexVolume_;   //!
        TBranch        *b_TaggerSimHits_v12_;   //!
        TBranch        *b_TaggerSimHits_v12_id_;   //!
        TBranch        *b_TaggerSimHits_v12_layerID_;   //!
        TBranch        *b_TaggerSimHits_v12_moduleID_;   //!
        TBranch        *b_TaggerSimHits_v12_edep_;   //!
        TBranch        *b_TaggerSimHits_v12_time_;   //!
        TBranch        *b_TaggerSimHits_v12_px_;   //!
        TBranch        *b_TaggerSimHits_v12_py_;   //!
        TBranch        *b_TaggerSimHits_v12_pz_;   //!
        TBranch        *b_TaggerSimHits_v12_energy_;   //!
        TBranch        *b_TaggerSimHits_v12_x_;   //!
        TBranch        *b_TaggerSimHits_v12_y_;   //!
        TBranch        *b_TaggerSimHits_v12_z_;   //!
        TBranch        *b_TaggerSimHits_v12_pathLength_;   //!
        TBranch        *b_TaggerSimHits_v12_trackID_;   //!
        TBranch        *b_TaggerSimHits_v12_pdgID_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_id_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_edep_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_x_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_y_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_z_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_time_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_edepContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_timeContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_nContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_id_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_edep_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_x_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_y_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_z_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_time_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_edepContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_timeContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_nContribs_;   //!
        TBranch        *b_TargetSimHits_v12_;   //!
        TBranch        *b_TargetSimHits_v12_id_;   //!
        TBranch        *b_TargetSimHits_v12_edep_;   //!
        TBranch        *b_TargetSimHits_v12_x_;   //!
        TBranch        *b_TargetSimHits_v12_y_;   //!
        TBranch        *b_TargetSimHits_v12_z_;   //!
        TBranch        *b_TargetSimHits_v12_time_;   //!
        TBranch        *b_TargetSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TargetSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TargetSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TargetSimHits_v12_edepContribs_;   //!
        TBranch        *b_TargetSimHits_v12_timeContribs_;   //!
        TBranch        *b_TargetSimHits_v12_nContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_id_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_edep_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_x_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_y_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_z_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_time_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_edepContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_timeContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_nContribs_;   //!
        TBranch        *b_RecoilSimHits_v12_;   //!
        TBranch        *b_RecoilSimHits_v12_id_;   //!
        TBranch        *b_RecoilSimHits_v12_layerID_;   //!
        TBranch        *b_RecoilSimHits_v12_moduleID_;   //!
        TBranch        *b_RecoilSimHits_v12_edep_;   //!
        TBranch        *b_RecoilSimHits_v12_time_;   //!
        TBranch        *b_RecoilSimHits_v12_px_;   //!
        TBranch        *b_RecoilSimHits_v12_py_;   //!
        TBranch        *b_RecoilSimHits_v12_pz_;   //!
        TBranch        *b_RecoilSimHits_v12_energy_;   //!
        TBranch        *b_RecoilSimHits_v12_x_;   //!
        TBranch        *b_RecoilSimHits_v12_y_;   //!
        TBranch        *b_RecoilSimHits_v12_z_;   //!
        TBranch        *b_RecoilSimHits_v12_pathLength_;   //!
        TBranch        *b_RecoilSimHits_v12_trackID_;   //!
        TBranch        *b_RecoilSimHits_v12_pdgID_;   //!
        TBranch        *b_EcalSimHits_v12_;   //!
        TBranch        *b_EcalSimHits_v12_id_;   //!
        TBranch        *b_EcalSimHits_v12_edep_;   //!
        TBranch        *b_EcalSimHits_v12_x_;   //!
        TBranch        *b_EcalSimHits_v12_y_;   //!
        TBranch        *b_EcalSimHits_v12_z_;   //!
        TBranch        *b_EcalSimHits_v12_time_;   //!
        TBranch        *b_EcalSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_EcalSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_EcalSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_EcalSimHits_v12_edepContribs_;   //!
        TBranch        *b_EcalSimHits_v12_timeContribs_;   //!
        TBranch        *b_EcalSimHits_v12_nContribs_;   //!
        TBranch        *b_HcalSimHits_v12_;   //!
        TBranch        *b_HcalSimHits_v12_id_;   //!
        TBranch        *b_HcalSimHits_v12_edep_;   //!
        TBranch        *b_HcalSimHits_v12_x_;   //!
        TBranch        *b_HcalSimHits_v12_y_;   //!
        TBranch        *b_HcalSimHits_v12_z_;   //!
        TBranch        *b_HcalSimHits_v12_time_;   //!
        TBranch        *b_HcalSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_HcalSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_HcalSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_HcalSimHits_v12_edepContribs_;   //!
        TBranch        *b_HcalSimHits_v12_timeContribs_;   //!
        TBranch        *b_HcalSimHits_v12_nContribs_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_id_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_time_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_px_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_py_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_x_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_y_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_z_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_id_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_time_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_px_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_py_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_x_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_y_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_z_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_id_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_time_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_px_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_py_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_x_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_y_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_z_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_id_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_time_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_px_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_py_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_x_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_y_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_z_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_id_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_time_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_px_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_py_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_x_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_y_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_z_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_pdgID_;   //!
        
        */
        // Fixed size dimensions of array or collections stored in the TTree if any.
        static constexpr int kMaxSimParticles_v12 = 54;
        static constexpr int kMaxSimParticles_v12_second_energy = 54;
        static constexpr int kMaxSimParticles_v12_second_pdgID = 54;
        static constexpr int kMaxSimParticles_v12_second_genStatus = 54;
        static constexpr int kMaxSimParticles_v12_second_time = 54;
        static constexpr int kMaxSimParticles_v12_second_x = 54;
        static constexpr int kMaxSimParticles_v12_second_y = 54;
        static constexpr int kMaxSimParticles_v12_second_z = 54;
        static constexpr int kMaxSimParticles_v12_second_endX = 54;
        static constexpr int kMaxSimParticles_v12_second_endY = 54;
        static constexpr int kMaxSimParticles_v12_second_endZ = 54;
        static constexpr int kMaxSimParticles_v12_second_px = 54;
        static constexpr int kMaxSimParticles_v12_second_py = 54;
        static constexpr int kMaxSimParticles_v12_second_pz = 54;
        static constexpr int kMaxSimParticles_v12_second_endpx = 54;
        static constexpr int kMaxSimParticles_v12_second_endpy = 54;
        static constexpr int kMaxSimParticles_v12_second_endpz = 54;
        static constexpr int kMaxSimParticles_v12_second_mass = 54;
        static constexpr int kMaxSimParticles_v12_second_charge = 54;
        static constexpr int kMaxSimParticles_v12_second_daughters = 54;
        static constexpr int kMaxSimParticles_v12_second_parents = 54;
        static constexpr int kMaxSimParticles_v12_second_processType = 54;
        static constexpr int kMaxSimParticles_v12_second_vertexVolume = 54;
        static constexpr int kMaxTaggerSimHits_v12 = 243;
        static constexpr int kMaxTaggerSimHits_v12_id = 243;
        static constexpr int kMaxTaggerSimHits_v12_layerID = 243;
        static constexpr int kMaxTaggerSimHits_v12_moduleID = 243;
        static constexpr int kMaxTaggerSimHits_v12_edep = 243;
        static constexpr int kMaxTaggerSimHits_v12_time = 243;
        static constexpr int kMaxTaggerSimHits_v12_px = 243;
        static constexpr int kMaxTaggerSimHits_v12_py = 243;
        static constexpr int kMaxTaggerSimHits_v12_pz = 243;
        static constexpr int kMaxTaggerSimHits_v12_energy = 243;
        static constexpr int kMaxTaggerSimHits_v12_x = 243;
        static constexpr int kMaxTaggerSimHits_v12_y = 243;
        static constexpr int kMaxTaggerSimHits_v12_z = 243;
        static constexpr int kMaxTaggerSimHits_v12_pathLength = 243;
        static constexpr int kMaxTaggerSimHits_v12_trackID = 243;
        static constexpr int kMaxTaggerSimHits_v12_pdgID = 243;
        static constexpr int kMaxTriggerPadUpSimHits_v12 = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_id = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_edep = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_x = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_y = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_z = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_time = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_trackIDContribs = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_incidentIDContribs = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_pdgCodeContribs = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_edepContribs = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_timeContribs = 75;
        static constexpr int kMaxTriggerPadUpSimHits_v12_nContribs = 75;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12 = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_id = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_edep = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_x = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_y = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_z = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_time = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_trackIDContribs = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_incidentIDContribs = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_pdgCodeContribs = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_edepContribs = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_timeContribs = 123;
        static constexpr int kMaxTriggerPadTaggerSimHits_v12_nContribs = 123;
        static constexpr int kMaxTargetSimHits_v12 = 118;
        static constexpr int kMaxTargetSimHits_v12_id = 118;
        static constexpr int kMaxTargetSimHits_v12_edep = 118;
        static constexpr int kMaxTargetSimHits_v12_x = 118;
        static constexpr int kMaxTargetSimHits_v12_y = 118;
        static constexpr int kMaxTargetSimHits_v12_z = 118;
        static constexpr int kMaxTargetSimHits_v12_time = 118;
        static constexpr int kMaxTargetSimHits_v12_trackIDContribs = 118;
        static constexpr int kMaxTargetSimHits_v12_incidentIDContribs = 118;
        static constexpr int kMaxTargetSimHits_v12_pdgCodeContribs = 118;
        static constexpr int kMaxTargetSimHits_v12_edepContribs = 118;
        static constexpr int kMaxTargetSimHits_v12_timeContribs = 118;
        static constexpr int kMaxTargetSimHits_v12_nContribs = 118;
        static constexpr int kMaxTriggerPadDownSimHits_v12 = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_id = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_edep = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_x = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_y = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_z = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_time = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_trackIDContribs = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_incidentIDContribs = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_pdgCodeContribs = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_edepContribs = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_timeContribs = 107;
        static constexpr int kMaxTriggerPadDownSimHits_v12_nContribs = 107;
        static constexpr int kMaxRecoilSimHits_v12 = 195;
        static constexpr int kMaxRecoilSimHits_v12_id = 195;
        static constexpr int kMaxRecoilSimHits_v12_layerID = 195;
        static constexpr int kMaxRecoilSimHits_v12_moduleID = 195;
        static constexpr int kMaxRecoilSimHits_v12_edep = 195;
        static constexpr int kMaxRecoilSimHits_v12_time = 195;
        static constexpr int kMaxRecoilSimHits_v12_px = 195;
        static constexpr int kMaxRecoilSimHits_v12_py = 195;
        static constexpr int kMaxRecoilSimHits_v12_pz = 195;
        static constexpr int kMaxRecoilSimHits_v12_energy = 195;
        static constexpr int kMaxRecoilSimHits_v12_x = 195;
        static constexpr int kMaxRecoilSimHits_v12_y = 195;
        static constexpr int kMaxRecoilSimHits_v12_z = 195;
        static constexpr int kMaxRecoilSimHits_v12_pathLength = 195;
        static constexpr int kMaxRecoilSimHits_v12_trackID = 195;
        static constexpr int kMaxRecoilSimHits_v12_pdgID = 195;
        static constexpr int kMaxEcalSimHits_v12 = 143;
        static constexpr int kMaxEcalSimHits_v12_id = 143;
        static constexpr int kMaxEcalSimHits_v12_edep = 143;
        static constexpr int kMaxEcalSimHits_v12_x = 143;
        static constexpr int kMaxEcalSimHits_v12_y = 143;
        static constexpr int kMaxEcalSimHits_v12_z = 143;
        static constexpr int kMaxEcalSimHits_v12_time = 143;
        static constexpr int kMaxEcalSimHits_v12_trackIDContribs = 143;
        static constexpr int kMaxEcalSimHits_v12_incidentIDContribs = 143;
        static constexpr int kMaxEcalSimHits_v12_pdgCodeContribs = 143;
        static constexpr int kMaxEcalSimHits_v12_edepContribs = 143;
        static constexpr int kMaxEcalSimHits_v12_timeContribs = 143;
        static constexpr int kMaxEcalSimHits_v12_nContribs = 143;
        static constexpr int kMaxHcalSimHits_v12 = 1071;
        static constexpr int kMaxHcalSimHits_v12_id = 1071;
        static constexpr int kMaxHcalSimHits_v12_edep = 1071;
        static constexpr int kMaxHcalSimHits_v12_x = 1071;
        static constexpr int kMaxHcalSimHits_v12_y = 1071;
        static constexpr int kMaxHcalSimHits_v12_z = 1071;
        static constexpr int kMaxHcalSimHits_v12_time = 1071;
        static constexpr int kMaxHcalSimHits_v12_trackIDContribs = 1071;
        static constexpr int kMaxHcalSimHits_v12_incidentIDContribs = 1071;
        static constexpr int kMaxHcalSimHits_v12_pdgCodeContribs = 1071;
        static constexpr int kMaxHcalSimHits_v12_edepContribs = 1071;
        static constexpr int kMaxHcalSimHits_v12_timeContribs = 1071;
        static constexpr int kMaxHcalSimHits_v12_nContribs = 1071;
        static constexpr int kMaxMagnetScoringPlaneHits_v12 = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_id = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_layerID = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_moduleID = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_edep = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_time = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_px = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_py = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_pz = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_energy = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_x = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_y = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_z = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_pathLength = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_trackID = 61375;
        static constexpr int kMaxMagnetScoringPlaneHits_v12_pdgID = 61375;
        static constexpr int kMaxTargetScoringPlaneHits_v12 = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_id = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_layerID = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_moduleID = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_edep = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_time = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_px = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_py = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_pz = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_energy = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_x = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_y = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_z = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_pathLength = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_trackID = 55339;
        static constexpr int kMaxTargetScoringPlaneHits_v12_pdgID = 55339;
        static constexpr int kMaxTrackerScoringPlaneHits_v12 = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_id = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_layerID = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_moduleID = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_edep = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_time = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_px = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_py = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_pz = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_energy = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_x = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_y = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_z = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_pathLength = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_trackID = 17556;
        static constexpr int kMaxTrackerScoringPlaneHits_v12_pdgID = 17556;
        static constexpr int kMaxEcalScoringPlaneHits_v12 = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_id = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_layerID = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_moduleID = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_edep = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_time = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_px = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_py = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_pz = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_energy = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_x = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_y = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_z = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_pathLength = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_trackID = 268;
        static constexpr int kMaxEcalScoringPlaneHits_v12_pdgID = 268;
        static constexpr int kMaxHcalScoringPlaneHits_v12 = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_id = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_layerID = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_moduleID = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_edep = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_time = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_px = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_py = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_pz = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_energy = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_x = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_y = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_z = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_pathLength = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_trackID = 1872;
        static constexpr int kMaxHcalScoringPlaneHits_v12_pdgID = 1872;
        static constexpr int kMaxeventNumber = 1;
        static constexpr int kMaxrun = 1;
        static constexpr int kMaxweight = 1;
        static constexpr int kMaxisRealData = 1;
        static constexpr int kMaxintParameters_ = 1;
        static constexpr int kMaxfloatParameters_ = 1;
        static constexpr int kMaxstringParameters_ = 1;
        static constexpr int kMaxchannelIDs = 1;
        static constexpr int kMaxsamples_ = 1440;
        static constexpr int kMaxsamples__word = 1440;
        static constexpr int kMaxnumSamplesPerDigi = 1;
        static constexpr int kMaxsampleOfInterest = 1;
        static constexpr int kMaxEcalRecHits_reco = 144;
        static constexpr int kMaxEcalRecHits_reco_id = 144;
        static constexpr int kMaxEcalRecHits_reco_amplitude = 144;
        static constexpr int kMaxEcalRecHits_reco_energy = 144;
        static constexpr int kMaxEcalRecHits_reco_time = 144;
        static constexpr int kMaxEcalRecHits_reco_xpos = 144;
        static constexpr int kMaxEcalRecHits_reco_ypos = 144;
        static constexpr int kMaxEcalRecHits_reco_zpos = 144;
        static constexpr int kMaxEcalRecHits_reco_isNoise = 144;
        static constexpr int kMaxnReadoutHits = 1;
        static constexpr int kMaxdeepestLayerHit = 1;
        static constexpr int kMaxsummedDet = 1;
        static constexpr int kMaxsummedTightIso = 1;
        static constexpr int kMaxmaxCellDep = 1;
        static constexpr int kMaxshowerRMS = 1;
        static constexpr int kMaxxStd = 1;
        static constexpr int kMaxyStd = 1;
        static constexpr int kMaxavgLayerHit = 1;
        static constexpr int kMaxstdLayerHit = 1;
        static constexpr int kMaxecalBackEnergy = 1;
        static constexpr int kMaxelectronContainmentEnergy = 1;
        static constexpr int kMaxphotonContainmentEnergy = 1;
        static constexpr int kMaxoutsideContainmentEnergy = 1;
        static constexpr int kMaxoutsideContainmentNHits = 1;
        static constexpr int kMaxoutsideContainmentXStd = 1;
        static constexpr int kMaxoutsideContainmentYStd = 1;
        static constexpr int kMaxdiscValue = 1;
        static constexpr int kMaxrecoilPx = 1;
        static constexpr int kMaxrecoilPy = 1;
        static constexpr int kMaxrecoilPz = 1;
        static constexpr int kMaxrecoilX = 1;
        static constexpr int kMaxrecoilY = 1;
        static constexpr int kMaxecalLayerEdepReadout = 1;
        static constexpr int kMaxHcalRecHits_reco = 358;
        static constexpr int kMaxHcalRecHits_reco_id = 358;
        static constexpr int kMaxHcalRecHits_reco_amplitude = 358;
        static constexpr int kMaxHcalRecHits_reco_energy = 358;
        static constexpr int kMaxHcalRecHits_reco_time = 358;
        static constexpr int kMaxHcalRecHits_reco_xpos = 358;
        static constexpr int kMaxHcalRecHits_reco_ypos = 358;
        static constexpr int kMaxHcalRecHits_reco_zpos = 358;
        static constexpr int kMaxHcalRecHits_reco_isNoise = 358;
        static constexpr int kMaxHcalRecHits_reco_pe = 358;
        static constexpr int kMaxHcalRecHits_reco_minpe = 358;
        static constexpr int kMaxmaxPEHit__id = 1;
        static constexpr int kMaxmaxPEHit__amplitude = 1;
        static constexpr int kMaxmaxPEHit__energy = 1;
        static constexpr int kMaxmaxPEHit__time = 1;
        static constexpr int kMaxmaxPEHit__xpos = 1;
        static constexpr int kMaxmaxPEHit__ypos = 1;
        static constexpr int kMaxmaxPEHit__zpos = 1;
        static constexpr int kMaxmaxPEHit__isNoise = 1;
        static constexpr int kMaxmaxPEHit__pe = 1;
        static constexpr int kMaxmaxPEHit__minpe = 1;
        static constexpr int kMaxpassesVeto = 1;
        static constexpr int kMaxtrigScintDigisTag_reco = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_id = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_amplitude = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_energy = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_time = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_xpos = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_ypos = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_zpos = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_isNoise = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_pe = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_minpe = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_barID = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_moduleID = 27;
        static constexpr int kMaxtrigScintDigisTag_reco_beamEfrac = 27;
        static constexpr int kMaxtrigScintDigisUp_reco = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_id = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_amplitude = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_energy = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_time = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_xpos = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_ypos = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_zpos = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_isNoise = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_pe = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_minpe = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_barID = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_moduleID = 18;
        static constexpr int kMaxtrigScintDigisUp_reco_beamEfrac = 18;
        static constexpr int kMaxtrigScintDigisDn_reco = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_id = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_amplitude = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_energy = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_time = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_xpos = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_ypos = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_zpos = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_isNoise = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_pe = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_minpe = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_barID = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_moduleID = 18;
        static constexpr int kMaxtrigScintDigisDn_reco_beamEfrac = 18;
        static constexpr int kMaxTriggerPadTaggerClusters_reco = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_hitIDs = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_energy = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_nHits = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_PE = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_seed = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_centroid = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_centroidX = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_centroidY = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_centroidZ = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_beamEfrac = 14;
        static constexpr int kMaxTriggerPadTaggerClusters_reco_time = 14;
        static constexpr int kMaxTriggerPadUpClusters_reco = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_hitIDs = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_energy = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_nHits = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_PE = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_seed = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_centroid = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_centroidX = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_centroidY = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_centroidZ = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_beamEfrac = 12;
        static constexpr int kMaxTriggerPadUpClusters_reco_time = 12;
        static constexpr int kMaxTriggerPadDownClusters_reco = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_hitIDs = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_energy = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_nHits = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_PE = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_seed = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_centroid = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_centroidX = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_centroidY = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_centroidZ = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_beamEfrac = 11;
        static constexpr int kMaxTriggerPadDownClusters_reco_time = 11;
        static constexpr int kMaxTriggerPadTracks_reco = 2;
        static constexpr int kMaxTriggerPadTracks_reco_centroid = 2;
        static constexpr int kMaxTriggerPadTracks_reco_centroidX = 2;
        static constexpr int kMaxTriggerPadTracks_reco_centroidY = 2;
        static constexpr int kMaxTriggerPadTracks_reco_centroidZ = 2;
        static constexpr int kMaxTriggerPadTracks_reco_residual = 2;
        static constexpr int kMaxTriggerPadTracks_reco_residualX = 2;
        static constexpr int kMaxTriggerPadTracks_reco_residualY = 2;
        static constexpr int kMaxTriggerPadTracks_reco_residualZ = 2;
        static constexpr int kMaxTriggerPadTracks_reco_nClusters = 2;
        static constexpr int kMaxTriggerPadTracks_reco_constituents = 2;
        static constexpr int kMaxTriggerPadTracks_reco_beamEfrac = 2;
        static constexpr int kMaxTriggerPadTracks_reco_px = 2;
        static constexpr int kMaxTriggerPadTracks_reco_py = 2;
        static constexpr int kMaxTriggerPadTracks_reco_pz = 2;
        static constexpr int kMaxTriggerPadTracks_reco_x = 2;
        static constexpr int kMaxTriggerPadTracks_reco_y = 2;
        static constexpr int kMaxTriggerPadTracks_reco_z = 2;
        static constexpr int kMaxSiStripHits_reco = 193;
        static constexpr int kMaxSiStripHits_reco_adcValues = 193;
        static constexpr int kMaxSiStripHits_reco_time = 193;
        static constexpr int kMaxSiStripHits_reco_simTrackerHits = 193;
        static constexpr int kMaxname = 1;
        static constexpr int kMaxpass = 1;
        static constexpr int kMaxvariables = 1;
        static constexpr int kMaxFindableTracks_reco = 12;
        static constexpr int kMaxFindableTracks_reco_particleTrackID = 12;
        static constexpr int kMaxFindableTracks_reco_is4sFindable = 12;
        static constexpr int kMaxFindableTracks_reco_is3s1aFindable = 12;
        static constexpr int kMaxFindableTracks_reco_is2s2aFindable = 12;
        static constexpr int kMaxFindableTracks_reco_is2aFindable = 12;
        static constexpr int kMaxFindableTracks_reco_is2sFindable = 12;
        static constexpr int kMaxFindableTracks_reco_is3sFindable = 12;
        

        // Declaration of leaf types
        int           SimParticles_v12_;
        int           SimParticles_v12_first[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_energy_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        int           SimParticles_v12_second_pdgID_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        int           SimParticles_v12_second_genStatus_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_time_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_x_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_y_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_z_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endX_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endY_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endZ_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_px_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_py_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_pz_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endpx_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endpy_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_endpz_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_mass_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        double        SimParticles_v12_second_charge_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        vector<int>     SimParticles_v12_second_daughters_[kMaxSimParticles_v12];
        vector<int>     SimParticles_v12_second_parents_[kMaxSimParticles_v12];
        int           SimParticles_v12_second_processType_[kMaxSimParticles_v12];   //[SimParticles_v12_]
        string          SimParticles_v12_second_vertexVolume_[kMaxSimParticles_v12];
        int           TaggerSimHits_v12_;
        int           TaggerSimHits_v12_id_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_layerID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_moduleID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_edep_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_time_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_px_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_py_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_pz_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_energy_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_x_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_y_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_z_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        float         TaggerSimHits_v12_pathLength_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_trackID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TaggerSimHits_v12_pdgID_[kMaxTaggerSimHits_v12];   //[TaggerSimHits_v12_]
        int           TriggerPadUpSimHits_v12_;
        int           TriggerPadUpSimHits_v12_id_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_edep_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_x_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_y_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_z_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        float         TriggerPadUpSimHits_v12_time_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        vector<int>     TriggerPadUpSimHits_v12_trackIDContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<int>     TriggerPadUpSimHits_v12_incidentIDContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<int>     TriggerPadUpSimHits_v12_pdgCodeContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<float>   TriggerPadUpSimHits_v12_edepContribs_[kMaxTriggerPadUpSimHits_v12];
        vector<float>   TriggerPadUpSimHits_v12_timeContribs_[kMaxTriggerPadUpSimHits_v12];
        unsigned int          TriggerPadUpSimHits_v12_nContribs_[kMaxTriggerPadUpSimHits_v12];   //[TriggerPadUpSimHits_v12_]
        int           TriggerPadTaggerSimHits_v12_;
        int           TriggerPadTaggerSimHits_v12_id_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_edep_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_x_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_y_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_z_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        float         TriggerPadTaggerSimHits_v12_time_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        vector<int>     TriggerPadTaggerSimHits_v12_trackIDContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<int>     TriggerPadTaggerSimHits_v12_incidentIDContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<int>     TriggerPadTaggerSimHits_v12_pdgCodeContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<float>   TriggerPadTaggerSimHits_v12_edepContribs_[kMaxTriggerPadTaggerSimHits_v12];
        vector<float>   TriggerPadTaggerSimHits_v12_timeContribs_[kMaxTriggerPadTaggerSimHits_v12];
        unsigned int          TriggerPadTaggerSimHits_v12_nContribs_[kMaxTriggerPadTaggerSimHits_v12];   //[TriggerPadTaggerSimHits_v12_]
        int           TargetSimHits_v12_;
        int           TargetSimHits_v12_id_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_edep_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_x_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_y_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_z_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        float         TargetSimHits_v12_time_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        vector<int>     TargetSimHits_v12_trackIDContribs_[kMaxTargetSimHits_v12];
        vector<int>     TargetSimHits_v12_incidentIDContribs_[kMaxTargetSimHits_v12];
        vector<int>     TargetSimHits_v12_pdgCodeContribs_[kMaxTargetSimHits_v12];
        vector<float>   TargetSimHits_v12_edepContribs_[kMaxTargetSimHits_v12];
        vector<float>   TargetSimHits_v12_timeContribs_[kMaxTargetSimHits_v12];
        unsigned int          TargetSimHits_v12_nContribs_[kMaxTargetSimHits_v12];   //[TargetSimHits_v12_]
        int           TriggerPadDownSimHits_v12_;
        int           TriggerPadDownSimHits_v12_id_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_edep_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_x_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_y_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_z_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        float         TriggerPadDownSimHits_v12_time_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        vector<int>     TriggerPadDownSimHits_v12_trackIDContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<int>     TriggerPadDownSimHits_v12_incidentIDContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<int>     TriggerPadDownSimHits_v12_pdgCodeContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<float>   TriggerPadDownSimHits_v12_edepContribs_[kMaxTriggerPadDownSimHits_v12];
        vector<float>   TriggerPadDownSimHits_v12_timeContribs_[kMaxTriggerPadDownSimHits_v12];
        unsigned int          TriggerPadDownSimHits_v12_nContribs_[kMaxTriggerPadDownSimHits_v12];   //[TriggerPadDownSimHits_v12_]
        int           RecoilSimHits_v12_;
        int           RecoilSimHits_v12_id_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_layerID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_moduleID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_edep_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_time_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_px_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_py_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_pz_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_energy_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_x_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_y_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_z_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        float         RecoilSimHits_v12_pathLength_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_trackID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           RecoilSimHits_v12_pdgID_[kMaxRecoilSimHits_v12];   //[RecoilSimHits_v12_]
        int           EcalSimHits_v12_;
        int           EcalSimHits_v12_id_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_edep_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_x_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_y_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_z_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        float         EcalSimHits_v12_time_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        vector<int>     EcalSimHits_v12_trackIDContribs_[kMaxEcalSimHits_v12];
        vector<int>     EcalSimHits_v12_incidentIDContribs_[kMaxEcalSimHits_v12];
        vector<int>     EcalSimHits_v12_pdgCodeContribs_[kMaxEcalSimHits_v12];
        vector<float>   EcalSimHits_v12_edepContribs_[kMaxEcalSimHits_v12];
        vector<float>   EcalSimHits_v12_timeContribs_[kMaxEcalSimHits_v12];
        unsigned int          EcalSimHits_v12_nContribs_[kMaxEcalSimHits_v12];   //[EcalSimHits_v12_]
        int           HcalSimHits_v12_;
        int           HcalSimHits_v12_id_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_edep_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_x_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_y_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_z_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        float         HcalSimHits_v12_time_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        vector<int>     HcalSimHits_v12_trackIDContribs_[kMaxHcalSimHits_v12];
        vector<int>     HcalSimHits_v12_incidentIDContribs_[kMaxHcalSimHits_v12];
        vector<int>     HcalSimHits_v12_pdgCodeContribs_[kMaxHcalSimHits_v12];
        vector<float>   HcalSimHits_v12_edepContribs_[kMaxHcalSimHits_v12];
        vector<float>   HcalSimHits_v12_timeContribs_[kMaxHcalSimHits_v12];
        unsigned int          HcalSimHits_v12_nContribs_[kMaxHcalSimHits_v12];   //[HcalSimHits_v12_]
        int           MagnetScoringPlaneHits_v12_;
        int           MagnetScoringPlaneHits_v12_id_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_layerID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_moduleID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_edep_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_time_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_px_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_py_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_pz_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_energy_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_x_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_y_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_z_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        float         MagnetScoringPlaneHits_v12_pathLength_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_trackID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           MagnetScoringPlaneHits_v12_pdgID_[kMaxMagnetScoringPlaneHits_v12];   //[MagnetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_;
        int           TargetScoringPlaneHits_v12_id_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_layerID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_moduleID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_edep_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_time_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_px_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_py_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_pz_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_energy_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_x_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_y_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_z_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        float         TargetScoringPlaneHits_v12_pathLength_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_trackID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TargetScoringPlaneHits_v12_pdgID_[kMaxTargetScoringPlaneHits_v12];   //[TargetScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_;
        int           TrackerScoringPlaneHits_v12_id_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_layerID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_moduleID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_edep_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_time_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_px_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_py_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_pz_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_energy_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_x_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_y_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_z_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        float         TrackerScoringPlaneHits_v12_pathLength_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_trackID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           TrackerScoringPlaneHits_v12_pdgID_[kMaxTrackerScoringPlaneHits_v12];   //[TrackerScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_;
        int           EcalScoringPlaneHits_v12_id_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_layerID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_moduleID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_edep_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_time_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_px_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_py_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_pz_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_energy_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_x_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_y_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_z_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        float         EcalScoringPlaneHits_v12_pathLength_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_trackID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           EcalScoringPlaneHits_v12_pdgID_[kMaxEcalScoringPlaneHits_v12];   //[EcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_;
        int           HcalScoringPlaneHits_v12_id_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_layerID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_moduleID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_edep_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_time_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_px_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_py_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_pz_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_energy_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_x_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_y_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_z_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        float         HcalScoringPlaneHits_v12_pathLength_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_trackID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        int           HcalScoringPlaneHits_v12_pdgID_[kMaxHcalScoringPlaneHits_v12];   //[HcalScoringPlaneHits_v12_]
        //ldmx::EventHeader *EventHeader;
        int           eventNumber_;
        int           run_;
        int           timestamp__fSec;
        int           timestamp__fNanoSec;
        double        weight_;
        bool          isRealData_;
        int           intParameters__;
        string          intParameters__first[kMaxintParameters_];
        int           intParameters__second[kMaxintParameters_];   //[intParameters__]
        int           floatParameters__;
        string          floatParameters__first[kMaxfloatParameters_];
        float         floatParameters__second[kMaxfloatParameters_];   //[floatParameters__]
        int           stringParameters__;
        string          stringParameters__first[kMaxstringParameters_];
        string          stringParameters__second[kMaxstringParameters_];
        //ldmx::HgcrocDigiCollection *EcalDigis_reco;
        vector<unsigned int> channelIDs_;
        int           samples__;
        unsigned int          samples__word_[kMaxsamples_];   //[samples__]
        unsigned int          numSamplesPerDigi_;
        unsigned int          sampleOfInterest_;
        int           EcalRecHits_reco_;
        int           EcalRecHits_reco_id_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        float         EcalRecHits_reco_amplitude_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        float         EcalRecHits_reco_energy_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        float         EcalRecHits_reco_time_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        float         EcalRecHits_reco_xpos_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        float         EcalRecHits_reco_ypos_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        float         EcalRecHits_reco_zpos_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        bool          EcalRecHits_reco_isNoise_[kMaxEcalRecHits_reco];   //[EcalRecHits_reco_]
        //ldmx::EcalVetoResult *EcalVeto_reco;
        bool          passesVeto_;
        int           nReadoutHits_;
        int           deepestLayerHit_;
        float         summedDet_;
        float         summedTightIso_;
        float         maxCellDep_;
        float         showerRMS_;
        float         xStd_;
        float         yStd_;
        float         avgLayerHit_;
        float         stdLayerHit_;
        float         ecalBackEnergy_;
        vector<float>   electronContainmentEnergy_;
        vector<float>   photonContainmentEnergy_;
        vector<float>   outsideContainmentEnergy_;
        vector<int>     outsideContainmentNHits_;
        vector<float>   outsideContainmentXStd_;
        vector<float>   outsideContainmentYStd_;
        float         discValue_;
        double        recoilPx_;
        double        recoilPy_;
        double        recoilPz_;
        float         recoilX_;
        float         recoilY_;
        vector<float>   ecalLayerEdepReadout_;
        int           HcalRecHits_reco_;
        int           HcalRecHits_reco_id_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_amplitude_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_energy_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_time_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_xpos_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_ypos_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_zpos_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        bool          HcalRecHits_reco_isNoise_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_pe_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        float         HcalRecHits_reco_minpe_[kMaxHcalRecHits_reco];   //[HcalRecHits_reco_]
        //ldmx::HcalVetoResult *HcalVeto_reco;
        int           maxPEHit__id_;
        float         maxPEHit__amplitude_;
        float         maxPEHit__energy_;
        float         maxPEHit__time_;
        float         maxPEHit__xpos_;
        float         maxPEHit__ypos_;
        float         maxPEHit__zpos_;
        bool          maxPEHit__isNoise_;
        float         maxPEHit__pe_;
        float         maxPEHit__minpe_;
        int           trigScintDigisTag_reco_;
        int           trigScintDigisTag_reco_id_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_amplitude_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_energy_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_time_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_xpos_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_ypos_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_zpos_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        bool          trigScintDigisTag_reco_isNoise_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_pe_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_minpe_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        int           trigScintDigisTag_reco_barID_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        int           trigScintDigisTag_reco_moduleID_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        float         trigScintDigisTag_reco_beamEfrac_[kMaxtrigScintDigisTag_reco];   //[trigScintDigisTag_reco_]
        int           trigScintDigisUp_reco_;
        int           trigScintDigisUp_reco_id_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_amplitude_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_energy_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_time_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_xpos_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_ypos_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_zpos_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        bool          trigScintDigisUp_reco_isNoise_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_pe_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_minpe_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        int           trigScintDigisUp_reco_barID_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        int           trigScintDigisUp_reco_moduleID_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        float         trigScintDigisUp_reco_beamEfrac_[kMaxtrigScintDigisUp_reco];   //[trigScintDigisUp_reco_]
        int           trigScintDigisDn_reco_;
        int           trigScintDigisDn_reco_id_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_amplitude_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_energy_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_time_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_xpos_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_ypos_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_zpos_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        bool          trigScintDigisDn_reco_isNoise_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_pe_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_minpe_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        int           trigScintDigisDn_reco_barID_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        int           trigScintDigisDn_reco_moduleID_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        float         trigScintDigisDn_reco_beamEfrac_[kMaxtrigScintDigisDn_reco];   //[trigScintDigisDn_reco_]
        int           TriggerPadTaggerClusters_reco_;
        vector<unsigned int> TriggerPadTaggerClusters_reco_hitIDs_[kMaxTriggerPadTaggerClusters_reco];
        double        TriggerPadTaggerClusters_reco_energy_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        int           TriggerPadTaggerClusters_reco_nHits_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        int           TriggerPadTaggerClusters_reco_PE_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        int           TriggerPadTaggerClusters_reco_seed_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        double        TriggerPadTaggerClusters_reco_centroid_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        double        TriggerPadTaggerClusters_reco_centroidX_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        double        TriggerPadTaggerClusters_reco_centroidY_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        double        TriggerPadTaggerClusters_reco_centroidZ_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        float         TriggerPadTaggerClusters_reco_beamEfrac_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        float         TriggerPadTaggerClusters_reco_time_[kMaxTriggerPadTaggerClusters_reco];   //[TriggerPadTaggerClusters_reco_]
        int           TriggerPadUpClusters_reco_;
        vector<unsigned int> TriggerPadUpClusters_reco_hitIDs_[kMaxTriggerPadUpClusters_reco];
        double        TriggerPadUpClusters_reco_energy_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        int           TriggerPadUpClusters_reco_nHits_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        int           TriggerPadUpClusters_reco_PE_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        int           TriggerPadUpClusters_reco_seed_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        double        TriggerPadUpClusters_reco_centroid_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        double        TriggerPadUpClusters_reco_centroidX_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        double        TriggerPadUpClusters_reco_centroidY_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        double        TriggerPadUpClusters_reco_centroidZ_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        float         TriggerPadUpClusters_reco_beamEfrac_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        float         TriggerPadUpClusters_reco_time_[kMaxTriggerPadUpClusters_reco];   //[TriggerPadUpClusters_reco_]
        int           TriggerPadDownClusters_reco_;
        vector<unsigned int> TriggerPadDownClusters_reco_hitIDs_[kMaxTriggerPadDownClusters_reco];
        double        TriggerPadDownClusters_reco_energy_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        int           TriggerPadDownClusters_reco_nHits_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        int           TriggerPadDownClusters_reco_PE_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        int           TriggerPadDownClusters_reco_seed_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        double        TriggerPadDownClusters_reco_centroid_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        double        TriggerPadDownClusters_reco_centroidX_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        double        TriggerPadDownClusters_reco_centroidY_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        double        TriggerPadDownClusters_reco_centroidZ_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        float         TriggerPadDownClusters_reco_beamEfrac_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        float         TriggerPadDownClusters_reco_time_[kMaxTriggerPadDownClusters_reco];   //[TriggerPadDownClusters_reco_]
        int           TriggerPadTracks_reco_;
        float         TriggerPadTracks_reco_centroid_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_centroidX_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_centroidY_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_centroidZ_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_residual_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_residualX_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_residualY_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_residualZ_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        int           TriggerPadTracks_reco_nClusters_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        //vector<ldmx::TrigScintCluster> TriggerPadTracks_reco_constituents_[kMaxTriggerPadTracks_reco];
        float         TriggerPadTracks_reco_beamEfrac_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_px_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_py_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_pz_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_x_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_y_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        float         TriggerPadTracks_reco_z_[kMaxTriggerPadTracks_reco];   //[TriggerPadTracks_reco_]
        int           SiStripHits_reco_;
        vector<short>   SiStripHits_reco_adcValues_[kMaxSiStripHits_reco];
        float         SiStripHits_reco_time_[kMaxSiStripHits_reco];   //[SiStripHits_reco_]
        //vector<ldmx::SimTrackerHit> SiStripHits_reco_simTrackerHits_[kMaxSiStripHits_reco];
        //ldmx::TriggerResult *Trigger_reco;
        TString         name_;
        bool          pass_;
        TArrayD         variables_;
        int           FindableTracks_reco_;
        int           FindableTracks_reco_particleTrackID_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        bool          FindableTracks_reco_is4sFindable_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        bool          FindableTracks_reco_is3s1aFindable_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        bool          FindableTracks_reco_is2s2aFindable_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        bool          FindableTracks_reco_is2aFindable_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        bool          FindableTracks_reco_is2sFindable_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        bool          FindableTracks_reco_is3sFindable_[kMaxFindableTracks_reco];   //[FindableTracks_reco_]
        //ldmx::TrackerVetoResult *TrackerVeto_reco;
        
        // List of branches
        TBranch        *b_SimParticles_v12_;   //!
        TBranch        *b_SimParticles_v12_first;   //!
        TBranch        *b_SimParticles_v12_second_energy_;   //!
        TBranch        *b_SimParticles_v12_second_pdgID_;   //!
        TBranch        *b_SimParticles_v12_second_genStatus_;   //!
        TBranch        *b_SimParticles_v12_second_time_;   //!
        TBranch        *b_SimParticles_v12_second_x_;   //!
        TBranch        *b_SimParticles_v12_second_y_;   //!
        TBranch        *b_SimParticles_v12_second_z_;   //!
        TBranch        *b_SimParticles_v12_second_endX_;   //!
        TBranch        *b_SimParticles_v12_second_endY_;   //!
        TBranch        *b_SimParticles_v12_second_endZ_;   //!
        TBranch        *b_SimParticles_v12_second_px_;   //!
        TBranch        *b_SimParticles_v12_second_py_;   //!
        TBranch        *b_SimParticles_v12_second_pz_;   //!
        TBranch        *b_SimParticles_v12_second_endpx_;   //!
        TBranch        *b_SimParticles_v12_second_endpy_;   //!
        TBranch        *b_SimParticles_v12_second_endpz_;   //!
        TBranch        *b_SimParticles_v12_second_mass_;   //!
        TBranch        *b_SimParticles_v12_second_charge_;   //!
        TBranch        *b_SimParticles_v12_second_daughters_;   //!
        TBranch        *b_SimParticles_v12_second_parents_;   //!
        TBranch        *b_SimParticles_v12_second_processType_;   //!
        TBranch        *b_SimParticles_v12_second_vertexVolume_;   //!
        TBranch        *b_TaggerSimHits_v12_;   //!
        TBranch        *b_TaggerSimHits_v12_id_;   //!
        TBranch        *b_TaggerSimHits_v12_layerID_;   //!
        TBranch        *b_TaggerSimHits_v12_moduleID_;   //!
        TBranch        *b_TaggerSimHits_v12_edep_;   //!
        TBranch        *b_TaggerSimHits_v12_time_;   //!
        TBranch        *b_TaggerSimHits_v12_px_;   //!
        TBranch        *b_TaggerSimHits_v12_py_;   //!
        TBranch        *b_TaggerSimHits_v12_pz_;   //!
        TBranch        *b_TaggerSimHits_v12_energy_;   //!
        TBranch        *b_TaggerSimHits_v12_x_;   //!
        TBranch        *b_TaggerSimHits_v12_y_;   //!
        TBranch        *b_TaggerSimHits_v12_z_;   //!
        TBranch        *b_TaggerSimHits_v12_pathLength_;   //!
        TBranch        *b_TaggerSimHits_v12_trackID_;   //!
        TBranch        *b_TaggerSimHits_v12_pdgID_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_id_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_edep_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_x_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_y_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_z_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_time_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_edepContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_timeContribs_;   //!
        TBranch        *b_TriggerPadUpSimHits_v12_nContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_id_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_edep_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_x_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_y_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_z_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_time_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_edepContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_timeContribs_;   //!
        TBranch        *b_TriggerPadTaggerSimHits_v12_nContribs_;   //!
        TBranch        *b_TargetSimHits_v12_;   //!
        TBranch        *b_TargetSimHits_v12_id_;   //!
        TBranch        *b_TargetSimHits_v12_edep_;   //!
        TBranch        *b_TargetSimHits_v12_x_;   //!
        TBranch        *b_TargetSimHits_v12_y_;   //!
        TBranch        *b_TargetSimHits_v12_z_;   //!
        TBranch        *b_TargetSimHits_v12_time_;   //!
        TBranch        *b_TargetSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TargetSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TargetSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TargetSimHits_v12_edepContribs_;   //!
        TBranch        *b_TargetSimHits_v12_timeContribs_;   //!
        TBranch        *b_TargetSimHits_v12_nContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_id_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_edep_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_x_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_y_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_z_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_time_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_edepContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_timeContribs_;   //!
        TBranch        *b_TriggerPadDownSimHits_v12_nContribs_;   //!
        TBranch        *b_RecoilSimHits_v12_;   //!
        TBranch        *b_RecoilSimHits_v12_id_;   //!
        TBranch        *b_RecoilSimHits_v12_layerID_;   //!
        TBranch        *b_RecoilSimHits_v12_moduleID_;   //!
        TBranch        *b_RecoilSimHits_v12_edep_;   //!
        TBranch        *b_RecoilSimHits_v12_time_;   //!
        TBranch        *b_RecoilSimHits_v12_px_;   //!
        TBranch        *b_RecoilSimHits_v12_py_;   //!
        TBranch        *b_RecoilSimHits_v12_pz_;   //!
        TBranch        *b_RecoilSimHits_v12_energy_;   //!
        TBranch        *b_RecoilSimHits_v12_x_;   //!
        TBranch        *b_RecoilSimHits_v12_y_;   //!
        TBranch        *b_RecoilSimHits_v12_z_;   //!
        TBranch        *b_RecoilSimHits_v12_pathLength_;   //!
        TBranch        *b_RecoilSimHits_v12_trackID_;   //!
        TBranch        *b_RecoilSimHits_v12_pdgID_;   //!
        TBranch        *b_EcalSimHits_v12_;   //!
        TBranch        *b_EcalSimHits_v12_id_;   //!
        TBranch        *b_EcalSimHits_v12_edep_;   //!
        TBranch        *b_EcalSimHits_v12_x_;   //!
        TBranch        *b_EcalSimHits_v12_y_;   //!
        TBranch        *b_EcalSimHits_v12_z_;   //!
        TBranch        *b_EcalSimHits_v12_time_;   //!
        TBranch        *b_EcalSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_EcalSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_EcalSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_EcalSimHits_v12_edepContribs_;   //!
        TBranch        *b_EcalSimHits_v12_timeContribs_;   //!
        TBranch        *b_EcalSimHits_v12_nContribs_;   //!
        TBranch        *b_HcalSimHits_v12_;   //!
        TBranch        *b_HcalSimHits_v12_id_;   //!
        TBranch        *b_HcalSimHits_v12_edep_;   //!
        TBranch        *b_HcalSimHits_v12_x_;   //!
        TBranch        *b_HcalSimHits_v12_y_;   //!
        TBranch        *b_HcalSimHits_v12_z_;   //!
        TBranch        *b_HcalSimHits_v12_time_;   //!
        TBranch        *b_HcalSimHits_v12_trackIDContribs_;   //!
        TBranch        *b_HcalSimHits_v12_incidentIDContribs_;   //!
        TBranch        *b_HcalSimHits_v12_pdgCodeContribs_;   //!
        TBranch        *b_HcalSimHits_v12_edepContribs_;   //!
        TBranch        *b_HcalSimHits_v12_timeContribs_;   //!
        TBranch        *b_HcalSimHits_v12_nContribs_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_id_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_time_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_px_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_py_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_x_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_y_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_z_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_MagnetScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_id_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_time_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_px_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_py_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_x_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_y_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_z_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_TargetScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_id_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_time_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_px_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_py_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_x_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_y_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_z_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_TrackerScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_id_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_time_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_px_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_py_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_x_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_y_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_z_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_EcalScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_id_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_layerID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_moduleID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_edep_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_time_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_px_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_py_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_pz_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_energy_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_x_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_y_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_z_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_pathLength_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_trackID_;   //!
        TBranch        *b_HcalScoringPlaneHits_v12_pdgID_;   //!
        TBranch        *b_EventHeader_eventNumber_;   //!
        TBranch        *b_EventHeader_run_;   //!
        TBranch        *b_EventHeader_timestamp__fSec;   //!
        TBranch        *b_EventHeader_timestamp__fNanoSec;   //!
        TBranch        *b_EventHeader_weight_;   //!
        TBranch        *b_EventHeader_isRealData_;   //!
        TBranch        *b_EventHeader_intParameters__;   //!
        TBranch        *b_intParameters__first;   //!
        TBranch        *b_intParameters__second;   //!
        TBranch        *b_EventHeader_floatParameters__;   //!
        TBranch        *b_floatParameters__first;   //!
        TBranch        *b_floatParameters__second;   //!
        TBranch        *b_EventHeader_stringParameters__;   //!
        TBranch        *b_stringParameters__first;   //!
        TBranch        *b_stringParameters__second;   //!
        TBranch        *b_EcalDigis_reco_channelIDs_;   //!
        TBranch        *b_EcalDigis_reco_samples__;   //!
        TBranch        *b_samples__word_;   //!
        TBranch        *b_EcalDigis_reco_numSamplesPerDigi_;   //!
        TBranch        *b_EcalDigis_reco_sampleOfInterest_;   //!
        TBranch        *b_EcalRecHits_reco_;   //!
        TBranch        *b_EcalRecHits_reco_id_;   //!
        TBranch        *b_EcalRecHits_reco_amplitude_;   //!
        TBranch        *b_EcalRecHits_reco_energy_;   //!
        TBranch        *b_EcalRecHits_reco_time_;   //!
        TBranch        *b_EcalRecHits_reco_xpos_;   //!
        TBranch        *b_EcalRecHits_reco_ypos_;   //!
        TBranch        *b_EcalRecHits_reco_zpos_;   //!
        TBranch        *b_EcalRecHits_reco_isNoise_;   //!
        TBranch        *b_EcalVeto_reco_passesVeto_;   //!
        TBranch        *b_EcalVeto_reco_nReadoutHits_;   //!
        TBranch        *b_EcalVeto_reco_deepestLayerHit_;   //!
        TBranch        *b_EcalVeto_reco_summedDet_;   //!
        TBranch        *b_EcalVeto_reco_summedTightIso_;   //!
        TBranch        *b_EcalVeto_reco_maxCellDep_;   //!
        TBranch        *b_EcalVeto_reco_showerRMS_;   //!
        TBranch        *b_EcalVeto_reco_xStd_;   //!
        TBranch        *b_EcalVeto_reco_yStd_;   //!
        TBranch        *b_EcalVeto_reco_avgLayerHit_;   //!
        TBranch        *b_EcalVeto_reco_stdLayerHit_;   //!
        TBranch        *b_EcalVeto_reco_ecalBackEnergy_;   //!
        TBranch        *b_EcalVeto_reco_electronContainmentEnergy_;   //!
        TBranch        *b_EcalVeto_reco_photonContainmentEnergy_;   //!
        TBranch        *b_EcalVeto_reco_outsideContainmentEnergy_;   //!
        TBranch        *b_EcalVeto_reco_outsideContainmentNHits_;   //!
        TBranch        *b_EcalVeto_reco_outsideContainmentXStd_;   //!
        TBranch        *b_EcalVeto_reco_outsideContainmentYStd_;   //!
        TBranch        *b_EcalVeto_reco_discValue_;   //!
        TBranch        *b_EcalVeto_reco_recoilPx_;   //!
        TBranch        *b_EcalVeto_reco_recoilPy_;   //!
        TBranch        *b_EcalVeto_reco_recoilPz_;   //!
        TBranch        *b_EcalVeto_reco_recoilX_;   //!
        TBranch        *b_EcalVeto_reco_recoilY_;   //!
        TBranch        *b_EcalVeto_reco_ecalLayerEdepReadout_;   //!
        TBranch        *b_HcalRecHits_reco_;   //!
        TBranch        *b_HcalRecHits_reco_id_;   //!
        TBranch        *b_HcalRecHits_reco_amplitude_;   //!
        TBranch        *b_HcalRecHits_reco_energy_;   //!
        TBranch        *b_HcalRecHits_reco_time_;   //!
        TBranch        *b_HcalRecHits_reco_xpos_;   //!
        TBranch        *b_HcalRecHits_reco_ypos_;   //!
        TBranch        *b_HcalRecHits_reco_zpos_;   //!
        TBranch        *b_HcalRecHits_reco_isNoise_;   //!
        TBranch        *b_HcalRecHits_reco_pe_;   //!
        TBranch        *b_HcalRecHits_reco_minpe_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__id_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__amplitude_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__energy_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__time_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__xpos_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__ypos_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__zpos_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__isNoise_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__pe_;   //!
        TBranch        *b_HcalVeto_reco_maxPEHit__minpe_;   //!
        TBranch        *b_HcalVeto_reco_passesVeto_;   //!
        TBranch        *b_trigScintDigisTag_reco_;   //!
        TBranch        *b_trigScintDigisTag_reco_id_;   //!
        TBranch        *b_trigScintDigisTag_reco_amplitude_;   //!
        TBranch        *b_trigScintDigisTag_reco_energy_;   //!
        TBranch        *b_trigScintDigisTag_reco_time_;   //!
        TBranch        *b_trigScintDigisTag_reco_xpos_;   //!
        TBranch        *b_trigScintDigisTag_reco_ypos_;   //!
        TBranch        *b_trigScintDigisTag_reco_zpos_;   //!
        TBranch        *b_trigScintDigisTag_reco_isNoise_;   //!
        TBranch        *b_trigScintDigisTag_reco_pe_;   //!
        TBranch        *b_trigScintDigisTag_reco_minpe_;   //!
        TBranch        *b_trigScintDigisTag_reco_barID_;   //!
        TBranch        *b_trigScintDigisTag_reco_moduleID_;   //!
        TBranch        *b_trigScintDigisTag_reco_beamEfrac_;   //!
        TBranch        *b_trigScintDigisUp_reco_;   //!
        TBranch        *b_trigScintDigisUp_reco_id_;   //!
        TBranch        *b_trigScintDigisUp_reco_amplitude_;   //!
        TBranch        *b_trigScintDigisUp_reco_energy_;   //!
        TBranch        *b_trigScintDigisUp_reco_time_;   //!
        TBranch        *b_trigScintDigisUp_reco_xpos_;   //!
        TBranch        *b_trigScintDigisUp_reco_ypos_;   //!
        TBranch        *b_trigScintDigisUp_reco_zpos_;   //!
        TBranch        *b_trigScintDigisUp_reco_isNoise_;   //!
        TBranch        *b_trigScintDigisUp_reco_pe_;   //!
        TBranch        *b_trigScintDigisUp_reco_minpe_;   //!
        TBranch        *b_trigScintDigisUp_reco_barID_;   //!
        TBranch        *b_trigScintDigisUp_reco_moduleID_;   //!
        TBranch        *b_trigScintDigisUp_reco_beamEfrac_;   //!
        TBranch        *b_trigScintDigisDn_reco_;   //!
        TBranch        *b_trigScintDigisDn_reco_id_;   //!
        TBranch        *b_trigScintDigisDn_reco_amplitude_;   //!
        TBranch        *b_trigScintDigisDn_reco_energy_;   //!
        TBranch        *b_trigScintDigisDn_reco_time_;   //!
        TBranch        *b_trigScintDigisDn_reco_xpos_;   //!
        TBranch        *b_trigScintDigisDn_reco_ypos_;   //!
        TBranch        *b_trigScintDigisDn_reco_zpos_;   //!
        TBranch        *b_trigScintDigisDn_reco_isNoise_;   //!
        TBranch        *b_trigScintDigisDn_reco_pe_;   //!
        TBranch        *b_trigScintDigisDn_reco_minpe_;   //!
        TBranch        *b_trigScintDigisDn_reco_barID_;   //!
        TBranch        *b_trigScintDigisDn_reco_moduleID_;   //!
        TBranch        *b_trigScintDigisDn_reco_beamEfrac_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_hitIDs_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_energy_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_nHits_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_PE_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_seed_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_centroid_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_centroidX_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_centroidY_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_centroidZ_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_beamEfrac_;   //!
        TBranch        *b_TriggerPadTaggerClusters_reco_time_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_hitIDs_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_energy_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_nHits_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_PE_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_seed_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_centroid_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_centroidX_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_centroidY_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_centroidZ_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_beamEfrac_;   //!
        TBranch        *b_TriggerPadUpClusters_reco_time_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_hitIDs_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_energy_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_nHits_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_PE_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_seed_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_centroid_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_centroidX_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_centroidY_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_centroidZ_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_beamEfrac_;   //!
        TBranch        *b_TriggerPadDownClusters_reco_time_;   //!
        TBranch        *b_TriggerPadTracks_reco_;   //!
        TBranch        *b_TriggerPadTracks_reco_centroid_;   //!
        TBranch        *b_TriggerPadTracks_reco_centroidX_;   //!
        TBranch        *b_TriggerPadTracks_reco_centroidY_;   //!
        TBranch        *b_TriggerPadTracks_reco_centroidZ_;   //!
        TBranch        *b_TriggerPadTracks_reco_residual_;   //!
        TBranch        *b_TriggerPadTracks_reco_residualX_;   //!
        TBranch        *b_TriggerPadTracks_reco_residualY_;   //!
        TBranch        *b_TriggerPadTracks_reco_residualZ_;   //!
        TBranch        *b_TriggerPadTracks_reco_nClusters_;   //!
        TBranch        *b_TriggerPadTracks_reco_beamEfrac_;   //!
        TBranch        *b_TriggerPadTracks_reco_px_;   //!
        TBranch        *b_TriggerPadTracks_reco_py_;   //!
        TBranch        *b_TriggerPadTracks_reco_pz_;   //!
        TBranch        *b_TriggerPadTracks_reco_x_;   //!
        TBranch        *b_TriggerPadTracks_reco_y_;   //!
        TBranch        *b_TriggerPadTracks_reco_z_;   //!
        TBranch        *b_SiStripHits_reco_;   //!
        TBranch        *b_SiStripHits_reco_adcValues_;   //!
        TBranch        *b_SiStripHits_reco_time_;   //!
        TBranch        *b_Trigger_reco_name_;   //!
        TBranch        *b_Trigger_reco_pass_;   //!
        TBranch        *b_Trigger_reco_variables_;   //!
        TBranch        *b_FindableTracks_reco_;   //!
        TBranch        *b_FindableTracks_reco_particleTrackID_;   //!
        TBranch        *b_FindableTracks_reco_is4sFindable_;   //!
        TBranch        *b_FindableTracks_reco_is3s1aFindable_;   //!
        TBranch        *b_FindableTracks_reco_is2s2aFindable_;   //!
        TBranch        *b_FindableTracks_reco_is2aFindable_;   //!
        TBranch        *b_FindableTracks_reco_is2sFindable_;   //!
        TBranch        *b_FindableTracks_reco_is3sFindable_;   //!
        TBranch        *b_TrackerVeto_reco_passesVeto_;   //!
           
        std::unique_ptr<const Acts::Logger> m_logger;

        const Acts::Logger& logger() const {return *m_logger;}
        
    }; // RootLdmxSimHitReader

} //namespace ActsExamples
