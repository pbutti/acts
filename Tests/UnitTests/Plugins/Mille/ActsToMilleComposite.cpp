// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"
#include "ActsPlugins/Mille/ActsToMille.hpp"
#include "ActsTests/CommonHelpers/AlignmentHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <map>
#include <string>

#include <Mille/MilleDecoder.h>
#include <Mille/MilleFactory.h>

using namespace Acts;
using namespace ActsTests;
using namespace ActsTests::AlignmentUtils;
using namespace Acts::detail::Test;
using namespace Acts::UnitConstants;

/// Dump a telescope track with all aligned surfaces grouped into one tilted
/// structure, and check that each measurement carries the sum of the surface
/// derivatives chained through the composite Jacobians.
BOOST_AUTO_TEST_CASE(CompositeStructureToMille) {
  aliTestUtils utils;
  TelescopeDetector detector(utils.geoCtx);
  const auto geometry = detector();

  auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
  const auto kfZeroPropagator =
      makeConstantFieldPropagator(geometry, 0_T, std::move(kfLogger));
  auto kfZero = KalmanFitterType(kfZeroPropagator);
  auto alignLogger = getDefaultLogger("Alignment", Logging::INFO);
  const auto alignZero =
      ActsAlignment::Alignment(std::move(kfZero), std::move(alignLogger));

  const auto& trajectories = createTrajectories(geometry, 1, utils);

  auto extensions = getExtensions(utils);
  TestSourceLink::SurfaceAccessor surfaceAccessor{*geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);
  KalmanFitterOptions kfOptions(
      utils.geoCtx, utils.magCtx, utils.calCtx, extensions,
      PropagatorPlainOptions(utils.geoCtx, utils.magCtx));

  // align all surfaces but layer 8, as a single structure
  const Transform3 frame = Translation3(5., -3., 100.) *
                           AngleAxis3(0.2, Vector3(0., 1., 1.).normalized());
  std::unordered_map<const Surface*, std::size_t> idxedAlignSurfaces;
  ActsPlugins::ActsToMille::CompositeMap composites;
  for (auto& det : detector.detectorStore) {
    const auto& surface = det->surface();
    if (surface.geometryId().layer() != 8) {
      idxedAlignSurfaces.emplace(&surface, idxedAlignSurfaces.size());
      composites.emplace(
          &surface,
          ActsPlugins::ActsToMille::makeCompositeLink(
              0, frame, surface.localToGlobalTransform(utils.geoCtx)));
    }
  }

  const auto& inputTraj = trajectories.front();
  kfOptions.referenceSurface = &(*inputTraj.startParameters).referenceSurface();
  auto evaluateRes = alignZero.evaluateTrackAlignmentState(
      kfOptions.geoContext, inputTraj.sourceLinks, *inputTraj.startParameters,
      kfOptions, idxedAlignSurfaces, ActsAlignment::AlignmentMask::All);
  BOOST_REQUIRE(evaluateRes.ok());
  const auto& alignState = evaluateRes.value();
  BOOST_REQUIRE_GT(alignState.alignedSurfaces.size(), 1u);

  const std::string fname = "myCompositeRecord.root";
  auto milleRecord = Mille::spawnMilleRecord(fname, true);
  ActsPlugins::ActsToMille::dumpToMille(alignState, *milleRecord, false,
                                        &composites);
  milleRecord->flushOutputFile();
  milleRecord.reset();

  auto milleReader = Mille::spawnMilleReader(fname);
  BOOST_REQUIRE(milleReader->open(fname));
  Mille::MilleDecoder decoder;
  std::vector<Mille::MilleMeasurement> measurements;
  BOOST_REQUIRE(decoder.decode(*milleReader, measurements) ==
                Mille::MilleDecoder::ReadResult::OK);

  // the surface measurements come first, then the Kalman correlation
  // pseudo-measurements, which carry no global derivatives
  BOOST_REQUIRE_GE(measurements.size(), alignState.measurementDim);
  for (std::size_t iMeas = 0; iMeas < measurements.size(); ++iMeas) {
    const auto& measurement = measurements[iMeas];
    if (iMeas >= alignState.measurementDim) {
      BOOST_CHECK(measurement.globalLabels.empty());
      continue;
    }

    AlignmentRowVector expected = AlignmentRowVector::Zero();
    for (const auto& [surface, indices] : alignState.alignedSurfaces) {
      expected +=
          alignState.alignmentToResidualDerivative.block<1, eAlignmentSize>(
              iMeas, eAlignmentSize * indices.second) *
          composites.at(surface).jacobian;
    }

    // each structure label at most once, zero derivatives are not written
    std::map<int, double> written;
    for (std::size_t i = 0; i < measurement.globalLabels.size(); ++i) {
      BOOST_CHECK(written
                      .emplace(measurement.globalLabels[i],
                               measurement.globalDerivatives[i])
                      .second);
    }
    for (std::size_t iPar = 0; iPar < eAlignmentSize; ++iPar) {
      auto it = written.find(static_cast<int>(iPar + 1));
      const double value = it == written.end() ? 0. : it->second;
      CHECK_CLOSE_ABS(value, expected(iPar), 1e-9);
    }
    BOOST_CHECK_LE(written.size(), eAlignmentSize);
  }
}
