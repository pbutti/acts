// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <random>
#include <utility>

#include "FitterTestsCommon.hpp"

using namespace Acts;
using namespace Acts::detail::Test;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using KalmanFitter =
    Acts::KalmanFitter<ConstantFieldPropagator, VectorMultiTrajectory>;

static const auto pion = Acts::ParticleHypothesis::pion();

KalmanUpdater kfUpdater;
KalmanSmoother kfSmoother;

// Construct initial track parameters.
Acts::BoundTrackParameters makeParameters(
    const Acts::ParticleHypothesis& particleHypothesis = pion) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-3_m, 0., 0., 42_ns);
  return Acts::BoundTrackParameters::createCurvilinear(
      mPos4, 0_degree, 90_degree, 1_e / 1_GeV, cov, particleHypothesis);
}

// Instantiate the tester
const FitterTester tester;

// reconstruction propagator and fitter
auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
const auto kfZeroPropagator =
    makeConstantFieldPropagator<ConstantFieldStepper>(tester.geometry, 0_T);
const auto kfZero = KalmanFitter(kfZeroPropagator, std::move(kfLogger));

std::default_random_engine rng(42);

auto makeDefaultKalmanFitterOptions() {
  KalmanFitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibratorStrict<VectorMultiTrajectory>>();
  extensions.updater.connect<&KalmanUpdater::operator()<VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother
      .connect<&KalmanSmoother::operator()<VectorMultiTrajectory>>(&kfSmoother);
  extensions.surfaceAccessor.connect<
      &Acts::detail::Test::TestSourceLink::SurfaceAccessor::operator()>(
      &tester.surfaceAccessor);

  return KalmanFitterOptions(
      tester.geoCtx, tester.magCtx, tester.calCtx, extensions,
      PropagatorPlainOptions(tester.geoCtx, tester.magCtx));
}

/// Checks that the calibrator is never handed a track state whose "best
/// available" parameters have not been computed yet, see
/// https://github.com/acts-project/acts/issues/5777
struct CalibratorContractChecker {
  mutable std::size_t nCalls = 0;

  template <typename trajectory_t>
  void operator()(const GeometryContext& gctx, const CalibrationContext& cctx,
                  const SourceLink& sourceLink,
                  typename trajectory_t::TrackStateProxy trackState) const {
    // filtering and smoothing have not happened yet, so the corresponding
    // components must not be allocated
    BOOST_CHECK(!trackState.hasFiltered());
    BOOST_CHECK(!trackState.hasSmoothed());
    // consequently the best available parameters are the predicted ones. we
    // compare the storage addresses to make sure no copy is involved
    BOOST_REQUIRE(trackState.hasPredicted());
    BOOST_CHECK_EQUAL(trackState.parameters().data(),
                      trackState.predicted().data());
    BOOST_CHECK_EQUAL(trackState.covariance().data(),
                      trackState.predictedCovariance().data());
    BOOST_CHECK(trackState.parameters().allFinite());
    BOOST_CHECK(trackState.covariance().allFinite());

    ++nCalls;

    testSourceLinkCalibrator<trajectory_t>(gctx, cctx, sourceLink, trackState);
  }
};

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

// The calibrator must see valid parameters. Before
// https://github.com/acts-project/acts/issues/5777 the filtered parameters
// were allocated up front, so `parameters()` returned uninitialized memory.
BOOST_AUTO_TEST_CASE(CalibratorSeesValidParameters) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  CalibratorContractChecker checker;
  kfOptions.extensions.calibrator
      .connect<&CalibratorContractChecker::operator()<VectorMultiTrajectory>>(
          &checker);

  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldNoSurfaceForward(kfZero, kfOptions, start, rng,
                                        expected_reversed, expected_smoothed,
                                        true);

  BOOST_CHECK_GT(checker.nCalls, 0u);
}

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldNoSurfaceForward(kfZero, kfOptions, start, rng,
                                        expected_reversed, expected_smoothed,
                                        true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  // regular smoothing
  kfOptions.reverseFiltering = false;
  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldWithSurfaceForward(kfZero, kfOptions, start, rng,
                                          expected_reversed, expected_smoothed,
                                          true);

  // reverse filtering instead of smoothing
  kfOptions.reverseFiltering = true;
  kfOptions.reverseFilteringCovarianceScaling = 100.0;
  expected_reversed = true;
  expected_smoothed = false;
  tester.test_ZeroFieldWithSurfaceForward(kfZero, kfOptions, start, rng,
                                          expected_reversed, expected_smoothed,
                                          true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  // regular smoothing
  kfOptions.reverseFiltering = false;
  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldWithSurfaceBackward(kfZero, kfOptions, start, rng,
                                           expected_reversed, expected_smoothed,
                                           true);

  // reverse filtering instead of smoothing
  kfOptions.reverseFiltering = true;
  kfOptions.reverseFilteringCovarianceScaling = 100.0;
  expected_reversed = true;
  expected_smoothed = false;
  tester.test_ZeroFieldWithSurfaceBackward(kfZero, kfOptions, start, rng,
                                           expected_reversed, expected_smoothed,
                                           true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldWithSurfaceAtExit(kfZero, kfOptions, start, rng,
                                         expected_reversed, expected_smoothed,
                                         true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldShuffled(kfZero, kfOptions, start, rng,
                                expected_reversed, expected_smoothed, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldWithHole(kfZero, kfOptions, start, rng,
                                expected_reversed, expected_smoothed, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  auto start = makeParameters();

  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  auto kfOptions = makeDefaultKalmanFitterOptions();

  TestOutlierFinder tof{5_mm};
  kfOptions.extensions.outlierFinder
      .connect<&TestOutlierFinder::operator()<VectorMultiTrajectory>>(&tof);

  bool expected_reversed = false;
  bool expected_smoothed = true;
  tester.test_ZeroFieldWithOutliers(kfZero, kfOptions, start, rng,
                                    expected_reversed, expected_smoothed, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithReverseFiltering) {
  auto start = makeParameters();

  auto test = [&](double threshold, bool reverse, bool expected_reversed,
                  bool expected_smoothed) {
    auto kfOptions = makeDefaultKalmanFitterOptions();

    TestReverseFilteringLogic trfl{threshold};
    kfOptions.extensions.reverseFilteringLogic
        .connect<&TestReverseFilteringLogic::operator()<VectorMultiTrajectory>>(
            &trfl);

    kfOptions.reverseFiltering = reverse;
    kfOptions.reverseFilteringCovarianceScaling = 100.0;

    tester.test_ZeroFieldWithReverseFiltering(kfZero, kfOptions, start, rng,
                                              expected_reversed,
                                              expected_smoothed, true);
  };

  // Track of 1 GeV with a threshold set at 0.1 GeV, reversed filtering should
  // not be used
  test(0.1_GeV, false, false, true);

  // Track of 1 GeV with a threshold set at 10 GeV, reversed filtering should
  // be used
  test(10._GeV, false, true, false);

  // Track of 1 GeV with a threshold set at 10 GeV, reversed filtering should
  // be used
  test(0.1_GeV, true, true, false);
}

// TODO this is not really Kalman fitter specific. is probably better tested
// with a synthetic trajectory.
BOOST_AUTO_TEST_CASE(GlobalCovariance) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_GlobalCovariance(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(OutputTrackParticleHypothesis) {
  auto start = makeParameters(Acts::ParticleHypothesis::muon());
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ParticleHypothesis(kfZero, kfOptions, start, rng);
}

namespace {

/// Measurements of one straight track through the cubic test detector, in
/// track order, from a fixed seed
std::vector<Acts::SourceLink> makeOrderedSourceLinks(
    const Acts::BoundTrackParameters& start) {
  std::default_random_engine seqRng(4242);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, seqRng);
  return FitterTester::prepareSourceLinks(measurements.sourceLinks);
}

}  // namespace

// Ordered input on a track that crosses every surface once: the sequential
// assignment gives the same fit as the default assignment by surface, with
// and without reverse filtering.
BOOST_AUTO_TEST_CASE(SequentialAssignmentOrderedEqualsBySurface) {
  auto start = makeParameters();
  const auto sourceLinks = makeOrderedSourceLinks(start);
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), FitterTester::nMeasurements);

  for (const bool reverse : {false, true}) {
    BOOST_TEST_INFO("reverse filtering " << reverse);
    auto kfOptions = makeDefaultKalmanFitterOptions();
    kfOptions.referenceSurface = &start.referenceSurface();
    kfOptions.reverseFiltering = reverse;
    BOOST_CHECK(kfOptions.measurementAssignment ==
                KalmanFitterMeasurementAssignment::BySurface);

    Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                Acts::VectorMultiTrajectory{}};
    auto bySurface = kfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                                kfOptions, tracks);
    BOOST_REQUIRE(bySurface.ok());

    kfOptions.measurementAssignment =
        KalmanFitterMeasurementAssignment::Sequential;
    auto sequential = kfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                                 kfOptions, tracks);
    BOOST_REQUIRE(sequential.ok());

    BOOST_CHECK_EQUAL(sequential->nMeasurements(), sourceLinks.size());
    BOOST_CHECK_EQUAL(sequential->nMeasurements(), bySurface->nMeasurements());
    BOOST_CHECK_EQUAL(sequential->nHoles(), bySurface->nHoles());
    BOOST_REQUIRE(sequential->hasReferenceSurface());
    BOOST_REQUIRE(bySurface->hasReferenceSurface());
    BOOST_CHECK(
        sequential->parameters().isApprox(bySurface->parameters(), 1e-12));
    BOOST_CHECK(
        sequential->covariance().isApprox(bySurface->covariance(), 1e-12));
  }
}

// A measurement moved to the front of the input: the default assignment
// ignores the order, the sequential one passes over the misplaced measurement
// and its surface becomes a hole.
BOOST_AUTO_TEST_CASE(SequentialAssignmentSkipsMisplacedMeasurement) {
  auto start = makeParameters();
  auto sourceLinks = makeOrderedSourceLinks(start);
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), FitterTester::nMeasurements);
  // order: m3, m0, m1, m2, m4, m5
  std::rotate(sourceLinks.begin(), std::next(sourceLinks.begin(), 3),
              std::next(sourceLinks.begin(), 4));

  auto kfOptions = makeDefaultKalmanFitterOptions();
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  auto bySurface = kfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                              kfOptions, tracks);
  BOOST_REQUIRE(bySurface.ok());
  BOOST_CHECK_EQUAL(bySurface->nMeasurements(), sourceLinks.size());
  BOOST_CHECK_EQUAL(bySurface->nHoles(), 0u);

  kfOptions.measurementAssignment =
      KalmanFitterMeasurementAssignment::Sequential;
  auto sequential = kfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                               kfOptions, tracks);
  BOOST_REQUIRE(sequential.ok());
  BOOST_CHECK_EQUAL(sequential->nMeasurements(), sourceLinks.size() - 1u);
  BOOST_CHECK_EQUAL(sequential->nHoles(), 1u);

  // No measurement may be passed over: m3 waits at the cursor until its
  // surface is reached, after that m0 blocks the cursor for good
  kfOptions.maxSkippedMeasurements = 0;
  auto strict = kfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                           kfOptions, tracks);
  BOOST_REQUIRE(strict.ok());
  BOOST_CHECK_EQUAL(strict->nMeasurements(), 1u);
}

// Two input measurements on the same surface of a track that crosses it once:
// the default assignment keeps the first one; the sequential one uses one on
// the single visit and passes over the other. Neither applies one twice.
BOOST_AUTO_TEST_CASE(SequentialAssignmentRepeatedSurface) {
  auto start = makeParameters();
  auto sourceLinks = makeOrderedSourceLinks(start);
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), FitterTester::nMeasurements);
  // a second measurement on the surface of m2, right after it
  sourceLinks.insert(std::next(sourceLinks.begin(), 3), sourceLinks[2]);

  for (const auto assignment :
       {KalmanFitterMeasurementAssignment::BySurface,
        KalmanFitterMeasurementAssignment::Sequential}) {
    auto kfOptions = makeDefaultKalmanFitterOptions();
    kfOptions.measurementAssignment = assignment;
    Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                Acts::VectorMultiTrajectory{}};
    auto res = kfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                          kfOptions, tracks);
    BOOST_REQUIRE(res.ok());
    BOOST_CHECK_EQUAL(res->nMeasurements(), FitterTester::nMeasurements);
    BOOST_CHECK_EQUAL(res->nHoles(), 0u);
  }
}

// Sequential assignment with a chi2 cut. The true measurement on one surface
// replaced by a far one: without a cut the sequential assignment uses it; with
// the cut it is rejected at the cursor and kept as an outlier, so the fit
// equals the one without that measurement, also when no measurement may be
// passed over. A cut that every true measurement passes changes nothing.
BOOST_AUTO_TEST_CASE(SequentialAssignmentChi2CutMakesOutlier) {
  auto start = makeParameters();
  std::default_random_engine seqRng(4244);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, seqRng);
  BOOST_REQUIRE_EQUAL(measurements.sourceLinks.size(),
                      FitterTester::nMeasurements);
  constexpr std::size_t kReplaced = 4;

  const auto clean = FitterTester::prepareSourceLinks(measurements.sourceLinks);

  auto displacedLinks = measurements.sourceLinks;
  displacedLinks[kReplaced].parameters[0] += 20_mm;
  const auto displaced = FitterTester::prepareSourceLinks(displacedLinks);

  auto missingLinks = measurements.sourceLinks;
  missingLinks.erase(std::next(missingLinks.begin(), kReplaced));
  const auto missing = FitterTester::prepareSourceLinks(missingLinks);

  for (const bool reverse : {false, true}) {
    BOOST_TEST_INFO("reverse filtering " << reverse);
    auto kfOptions = makeDefaultKalmanFitterOptions();
    kfOptions.referenceSurface = &start.referenceSurface();
    kfOptions.reverseFiltering = reverse;
    Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                Acts::VectorMultiTrajectory{}};

    auto cleanReference =
        kfZero.fit(clean.begin(), clean.end(), start, kfOptions, tracks);
    BOOST_REQUIRE(cleanReference.ok());
    auto reference =
        kfZero.fit(missing.begin(), missing.end(), start, kfOptions, tracks);
    BOOST_REQUIRE(reference.ok());
    BOOST_CHECK_EQUAL(reference->nMeasurements(),
                      FitterTester::nMeasurements - 1u);
    BOOST_CHECK_EQUAL(reference->nHoles(), 1u);

    kfOptions.measurementAssignment =
        KalmanFitterMeasurementAssignment::Sequential;
    auto noCut = kfZero.fit(displaced.begin(), displaced.end(), start,
                            kfOptions, tracks);
    BOOST_REQUIRE(noCut.ok());
    BOOST_CHECK_EQUAL(noCut->nMeasurements(), FitterTester::nMeasurements);
    BOOST_CHECK_GT(noCut->chi2(), 100. * reference->chi2());

    kfOptions.maxPredictedChi2 = 25.;
    auto cleanCut =
        kfZero.fit(clean.begin(), clean.end(), start, kfOptions, tracks);
    BOOST_REQUIRE(cleanCut.ok());
    BOOST_CHECK_EQUAL(cleanCut->nMeasurements(), FitterTester::nMeasurements);
    BOOST_REQUIRE(cleanCut->hasReferenceSurface());
    BOOST_CHECK(
        cleanCut->parameters().isApprox(cleanReference->parameters(), 1e-12));

    for (const std::size_t maxSkip :
         {std::numeric_limits<std::size_t>::max(), std::size_t{0}}) {
      BOOST_TEST_INFO("maxSkippedMeasurements " << maxSkip);
      kfOptions.maxSkippedMeasurements = maxSkip;
      auto cut = kfZero.fit(displaced.begin(), displaced.end(), start,
                            kfOptions, tracks);
      BOOST_REQUIRE(cut.ok());
      BOOST_CHECK_EQUAL(cut->nMeasurements(), FitterTester::nMeasurements - 1u);
      BOOST_CHECK_EQUAL(cut->nOutliers(), 1u);
      BOOST_CHECK_EQUAL(cut->nHoles(), 0u);
      BOOST_CHECK_CLOSE(cut->chi2(), reference->chi2(), 1e-6);
      BOOST_REQUIRE(cut->hasReferenceSurface());
      BOOST_CHECK(cut->parameters().isApprox(reference->parameters(), 1e-10));
      BOOST_CHECK(cut->covariance().isApprox(reference->covariance(), 1e-10));
      // The outlier state carries the rejected measurement
      std::size_t nOutlierStates = 0;
      for (const auto& state : cut->trackStatesReversed()) {
        if (state.typeFlags().isOutlier()) {
          ++nOutlierStates;
          BOOST_CHECK(state.hasCalibrated());
          BOOST_CHECK_EQUAL(state.getUncalibratedSourceLink()
                                .template get<TestSourceLink>()
                                .index(),
                            displacedLinks[kReplaced].index());
        }
      }
      BOOST_CHECK_EQUAL(nOutlierStates, 1u);
    }
  }
}

// Sequential assignment with a chi2 cut: a displaced decoy right before the
// true measurement on the same surface. Without a cut the decoy is used and
// the true one passed over; with the cut the decoy is rejected and the next
// measurement in the window, the true one, is used: the fit equals the clean
// one.
BOOST_AUTO_TEST_CASE(SequentialAssignmentChi2CutTriesNextInWindow) {
  auto start = makeParameters();
  std::default_random_engine seqRng(4245);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, seqRng);
  BOOST_REQUIRE_EQUAL(measurements.sourceLinks.size(),
                      FitterTester::nMeasurements);
  constexpr std::size_t kDecoyed = 2;
  const auto clean = FitterTester::prepareSourceLinks(measurements.sourceLinks);

  auto decoyLinks = measurements.sourceLinks;
  auto decoy = decoyLinks[kDecoyed];
  decoy.parameters[0] += 20_mm;
  decoyLinks.insert(std::next(decoyLinks.begin(), kDecoyed), decoy);
  const auto decoyed = FitterTester::prepareSourceLinks(decoyLinks);

  auto kfOptions = makeDefaultKalmanFitterOptions();
  kfOptions.referenceSurface = &start.referenceSurface();
  kfOptions.measurementAssignment =
      KalmanFitterMeasurementAssignment::Sequential;
  kfOptions.maxSkippedMeasurements = 1;
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  auto reference =
      kfZero.fit(clean.begin(), clean.end(), start, kfOptions, tracks);
  BOOST_REQUIRE(reference.ok());
  BOOST_REQUIRE_EQUAL(reference->nMeasurements(), FitterTester::nMeasurements);

  auto noCut =
      kfZero.fit(decoyed.begin(), decoyed.end(), start, kfOptions, tracks);
  BOOST_REQUIRE(noCut.ok());
  BOOST_CHECK_EQUAL(noCut->nMeasurements(), FitterTester::nMeasurements);
  BOOST_CHECK_GT(noCut->chi2(), 100. * reference->chi2());

  kfOptions.maxPredictedChi2 = 25.;
  auto cut =
      kfZero.fit(decoyed.begin(), decoyed.end(), start, kfOptions, tracks);
  BOOST_REQUIRE(cut.ok());
  BOOST_CHECK_EQUAL(cut->nMeasurements(), FitterTester::nMeasurements);
  BOOST_CHECK_EQUAL(cut->nHoles(), 0u);
  BOOST_CHECK_CLOSE(cut->chi2(), reference->chi2(), 1e-6);
  BOOST_REQUIRE(cut->hasReferenceSurface());
  BOOST_CHECK(cut->parameters().isApprox(reference->parameters(), 1e-10));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
