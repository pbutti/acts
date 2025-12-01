#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples

u = acts.UnitConstants


def runTruthTrackingKalman(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    inputHitsPath: Optional[Path] = None,
    decorators=[],
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=1,
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addDigitization,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )
    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        addKalmanTracks,
    )

    s = s or acts.examples.Sequencer(
        events=100000, numThreads=1, logLevel=acts.logging.INFO
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    logger = acts.logging.getLogger("Truth tracking example")

    if inputParticlePath is None:
        addParticleGun(
            s,
            #ParticleConfig(num=1, pdg=acts.PdgParticle.eProton, randomizeCharge=True),
            #ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
            #ParticleConfig(num=1, pdg=acts.PdgParticle.ePionPlus, randomizeCharge=True),
            ParticleConfig(num=1, pdg=acts.PdgParticle.eKaonPlus, randomizeCharge=True),
            #EtaConfig(-3.0, 3.0, uniform=True),
            EtaConfig(-1.5, 1.5, uniform=True),
            MomentumConfig(0.1 * u.GeV, 2.5 * u.GeV, transverse=False),
            PhiConfig(0.0, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=1,
            rnd=rnd,
            printParticles=False
        )
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        s.addReader(
            acts.examples.RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_generated",
            )
        )
        s.addWhiteboardAlias("particles", "particles_generated")

    if inputHitsPath is None:

        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
            enableInteractions=True,
            pMin = 0.1
        )
    else:
        logger.info("Reading hits from %s", inputHitsPath.resolve())
        assert inputHitsPath.exists()
        s.addReader(
            acts.examples.RootSimHitReader(
                level=acts.logging.INFO,
                filePath=str(inputHitsPath.resolve()),
                outputSimHits="simhits",
            )
        )
    
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.3 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.pion,
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold,
        reverseFilteringCovarianceScaling,
        logLevel=acts.logging.INFO,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected-tracks")

    s.addWriter(
        acts.examples.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf.root"),
        )
    )

#    s.addWriter(
#        acts.examples.RootParticleWriter(
#            level=acts.logging.INFO,
#            inputParticles="particles_selected",
#            filePath=str(outputDir / "particles.root"),
#        )
#    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"

    ## GenericDetector
    # detector = acts.examples.GenericDetector()
    # trackingGeometry = detector.trackingGeometry()
    # digiConfigFile = (
    #     srcdir
    #     / "Examples/Configs/generic-digi-smearing-config.json"
    # )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTrackingKalman(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
