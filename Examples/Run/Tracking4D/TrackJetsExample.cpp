
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleSelectorOptions.hpp"
#include "ActsExamples/Options/ParticleSmearingOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/Pythia8Options.hpp"

//Algorithms
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TrackJets/TrackJetsAlgorithm.hpp"


#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;


int main(int argc, char* argv[]) {


  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  Options::addParticleSelectorOptions(desc);
  Options::addParticleSmearingOptions(desc);
  Options::addMagneticFieldOptions(desc);

  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // Setup the magnetic field - not used
  auto magneticField = Options::readMagneticField(vars);

  Sequencer sequencer(Options::readSequencerConfig(vars));

  // setup event generator
  EventGenerator::Config evgen = Options::readPythia8Options(vars,logLevel);
  evgen.outputParticles = "particles_generated";
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen,logLevel));

  
  
  // pre-select particles
  ParticleSelector::Config selectParticles =
      Options::readParticleSelectorConfig(vars);
  selectParticles.inputParticles = evgen.outputParticles;
  selectParticles.outputParticles = "particles_selected";
  // smearing only works with charge particles for now
  selectParticles.removeNeutral = true;
  selectParticles.absEtaMax = 2.5;
  selectParticles.rhoMax =4 * 1_mm;
  selectParticles.ptMin = 500 * 1_MeV;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectParticles, logLevel));


  // Run the particle smearing
  auto particleSmearingCfg = setupParticleSmearing(
      vars, sequencer, rnd, selectParticles.outputParticles);

  
  
  // Pass the smeared particles to the track jets algorithm
  TrackJetsAlgorithm::Config trackJetsConfig;
  trackJetsConfig.inputTrackCollection = particleSmearingCfg.outputTrackParameters;
  trackJetsConfig.radius = 0.4;
  trackJetsConfig.outputTrackJets = "TrackJetsAntikt04";

  sequencer.addAlgorithm(std::make_shared<TrackJetsAlgorithm>(
      trackJetsConfig,logLevel));

  
      
      
  return sequencer.run();
                          
  
}
