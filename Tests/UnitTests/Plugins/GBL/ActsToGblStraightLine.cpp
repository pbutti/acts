// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/GBL/ActsToGbl.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

using namespace Acts;
using namespace ActsPlugins::ActsToGbl;

namespace {

/// Closure test of the Acts -> GBL conversion on a straight line.
///
/// The reference trajectory runs along @p theta / @p phi through five planes
/// perpendicular to it, whose local axes are the curvilinear U and V. The
/// "true" track is a *different* straight line, offset and tilted with respect
/// to the reference, which makes its measurements a linear function of arc
/// length.
///
/// A straight-line GBL model can absorb that deviation exactly, so a correct
/// set of jacobians must yield a vanishing chi2 and recover the injected
/// offsets and slopes. An error in the curvilinear basis change, in the
/// parameter ordering or in the slope convention cannot absorb it and shows up
/// immediately as a large chi2.
void checkStraightLineClosure(double theta, double phi) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;
  auto logger = getDefaultLogger("ActsToGbl", Logging::WARNING);

  const Vector3 dir(std::sin(theta) * std::cos(phi),
                    std::sin(theta) * std::sin(phi), std::cos(theta));
  const RotationMatrix3 frame = CurvilinearSurface(dir).referenceFrame();

  const std::vector<double> arcLength = {0., 100., 200., 300., 400.};

  // the true track, as offset + slope along the curvilinear U and V
  const double offsetU = 0.35;
  const double slopeU = 2.0e-3;
  const double offsetV = -0.20;
  const double slopeV = -1.1e-3;

  const double sigma = 0.010;
  const double qOverP = 1. / 1000.;

  std::vector<std::shared_ptr<PlaneSurface>> surfaces;
  std::vector<GblTrackState> states;

  for (std::size_t i = 0; i < arcLength.size(); ++i) {
    Transform3 transform = Transform3::Identity();
    transform.linear() = frame;
    transform.translation() = arcLength[i] * dir;
    surfaces.push_back(Surface::makeShared<PlaneSurface>(
        transform, std::make_shared<RectangleBounds>(5000., 5000.)));

    const Vector3 position = arcLength[i] * dir;

    GblTrackState state;
    state.surface = surfaces.back().get();
    state.boundParameters[eBoundPhi] = phi;
    state.boundParameters[eBoundTheta] = theta;
    state.boundParameters[eBoundQOverP] = qOverP;
    state.freeParameters.segment<3>(eFreePos0) = position;
    state.freeParameters.segment<3>(eFreeDir0) = dir;
    state.freeParameters[eFreeQOverP] = qOverP;

    // straight-line bound-to-bound jacobian from the previous plane
    if (i > 0) {
      const double ds = arcLength[i] - arcLength[i - 1];
      FreeMatrix freeTransport = FreeMatrix::Identity();
      freeTransport.block<3, 3>(eFreePos0, eFreeDir0) =
          ds * SquareMatrix3::Identity();
      state.jacobian =
          state.surface->freeToBoundJacobian(gctx, position, dir) *
          freeTransport *
          surfaces[i - 1]->boundToFreeJacobian(gctx, arcLength[i - 1] * dir,
                                               dir);
    }

    state.isMeasurement = true;
    state.measurementDim = 2;
    state.measurementIndices[0] = eBoundLoc0;
    state.measurementIndices[1] = eBoundLoc1;
    state.calibrated[0] = offsetU + slopeU * arcLength[i];
    state.calibrated[1] = offsetV + slopeV * arcLength[i];
    state.calibratedCovariance(0, 0) = sigma * sigma;
    state.calibratedCovariance(1, 1) = sigma * sigma;

    states.push_back(state);
  }

  GblConversionOptions options;
  options.addScatterers = false;
  options.addGlobals = false;

  auto result = buildGblPoints(gctx, mctx, states, {}, options, *logger);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->points.size(), arcLength.size());

  // straight line: no curvature parameter in the fit
  gbl::GblTrajectory trajectory(result->points, false);
  BOOST_REQUIRE(trajectory.isValid());

  double chi2 = 0.;
  double lostWeight = 0.;
  int ndf = 0;
  trajectory.fit(chi2, ndf, lostWeight);

  // 5 planes x 2D measurements - 4 straight-line parameters
  BOOST_CHECK_EQUAL(ndf, 6);
  BOOST_CHECK_SMALL(chi2, 1e-6);

  Eigen::VectorXd parameters;
  Eigen::MatrixXd covariance;
  trajectory.getResults(1, parameters, covariance);

  // GBL local parameters are (q/p, u'_1, u'_2, u_1, u_2)
  BOOST_CHECK_SMALL(parameters(1) - slopeU, 1e-12);
  BOOST_CHECK_SMALL(parameters(2) - slopeV, 1e-12);
  BOOST_CHECK_SMALL(parameters(3) - offsetU, 1e-9);
  BOOST_CHECK_SMALL(parameters(4) - offsetV, 1e-9);
}

/// Closure test at non-normal incidence.
///
/// The planes here are tilted with respect to the track, with alternating tilt
/// sign, and the measurements are the *real* intersections of the true line
/// with them. This is what exercises the path corrections in the curvilinear
/// <-> bound coordinate changes: the curvilinear plane and the sensor pass
/// through the same point but are tilted relative to each other, so a varied
/// track crosses them at points separated along the trajectory.
///
/// Dropping those corrections leaves the jacobian wrong at first order by
/// sin^2(incidence) - 144 units of chi2 at a 30 degree tilt, against the 6e-2
/// that the (purely second-order) nonlinearity of this test construction
/// produces. The deviation is kept small here so the remaining nonlinearity is
/// far below the tolerance.
void checkTiltedClosure(double tiltAngle) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;
  auto logger = getDefaultLogger("ActsToGbl", Logging::WARNING);

  const double theta = 1.1;
  const double phi = 0.4;
  const Vector3 dir(std::sin(theta) * std::cos(phi),
                    std::sin(theta) * std::sin(phi), std::cos(theta));
  const RotationMatrix3 curvFrame = CurvilinearSurface(dir).referenceFrame();
  const Vector3 axisU = curvFrame.col(0);
  const Vector3 axisV = curvFrame.col(1);

  // small deviation: keep the test construction in its linear regime
  const double offsetU = 3.5e-3;
  const double slopeU = 2.0e-5;
  const double offsetV = -2.0e-3;
  const double slopeV = -1.1e-5;

  const Vector3 trueOrigin = offsetU * axisU + offsetV * axisV;
  const Vector3 trueDir = (dir + slopeU * axisU + slopeV * axisV).normalized();

  const std::vector<double> arcLength = {0., 100., 200., 300., 400.};
  const double sigma = 0.010;
  const double qOverP = 1. / 1000.;

  FreeVector pathDerivatives = FreeVector::Zero();
  pathDerivatives.head<3>() = dir;

  std::vector<std::shared_ptr<PlaneSurface>> surfaces;
  std::vector<GblTrackState> states;

  for (std::size_t i = 0; i < arcLength.size(); ++i) {
    // alternate the tilt so the planes are not all parallel
    const double sign = (i % 2 == 0) ? 1. : -1.;
    const Vector3 normal =
        (std::cos(tiltAngle) * dir + sign * std::sin(tiltAngle) * axisU)
            .normalized();
    const Vector3 seed =
        (std::abs(normal.z()) < 0.9 ? Vector3::UnitZ() : Vector3::UnitX());
    RotationMatrix3 rotation;
    rotation.col(0) = seed.cross(normal).normalized();
    rotation.col(1) = normal.cross(rotation.col(0));
    rotation.col(2) = normal;

    const Vector3 position = arcLength[i] * dir;
    Transform3 transform = Transform3::Identity();
    transform.linear() = rotation;
    transform.translation() = position;
    surfaces.push_back(Surface::makeShared<PlaneSurface>(
        transform, std::make_shared<RectangleBounds>(10000., 10000.)));

    GblTrackState state;
    state.surface = surfaces.back().get();
    state.boundParameters[eBoundPhi] = phi;
    state.boundParameters[eBoundTheta] = theta;
    state.boundParameters[eBoundQOverP] = qOverP;
    state.freeParameters.segment<3>(eFreePos0) = position;
    state.freeParameters.segment<3>(eFreeDir0) = dir;
    state.freeParameters[eFreeQOverP] = qOverP;

    if (i > 0) {
      const double ds = arcLength[i] - arcLength[i - 1];
      FreeMatrix freeTransport = FreeMatrix::Identity();
      freeTransport.block<3, 3>(eFreePos0, eFreeDir0) =
          ds * SquareMatrix3::Identity();
      state.jacobian = detail::boundToBoundTransportJacobian(
          gctx, state.freeParameters,
          surfaces[i - 1]->boundToFreeJacobian(gctx, arcLength[i - 1] * dir,
                                               dir),
          freeTransport, pathDerivatives, *state.surface);
    }

    // real intersection of the true line with this tilted plane
    const double s = normal.dot(position - trueOrigin) / normal.dot(trueDir);
    const Vector3 delta = trueOrigin + s * trueDir - position;

    state.isMeasurement = true;
    state.measurementDim = 2;
    state.measurementIndices[0] = eBoundLoc0;
    state.measurementIndices[1] = eBoundLoc1;
    state.calibrated[0] = delta.dot(rotation.col(0));
    state.calibrated[1] = delta.dot(rotation.col(1));
    state.calibratedCovariance(0, 0) = sigma * sigma;
    state.calibratedCovariance(1, 1) = sigma * sigma;

    states.push_back(state);
  }

  GblConversionOptions options;
  options.addScatterers = false;
  options.addGlobals = false;

  auto result = buildGblPoints(gctx, mctx, states, {}, options, *logger);
  BOOST_REQUIRE(result.ok());

  gbl::GblTrajectory trajectory(result->points, false);
  BOOST_REQUIRE(trajectory.isValid());

  double chi2 = 0.;
  double lostWeight = 0.;
  int ndf = 0;
  trajectory.fit(chi2, ndf, lostWeight);

  BOOST_CHECK_EQUAL(ndf, 6);
  BOOST_CHECK_SMALL(chi2, 1e-8);

  Eigen::VectorXd parameters;
  Eigen::MatrixXd covariance;
  trajectory.getResults(1, parameters, covariance);
  BOOST_CHECK_SMALL(parameters(3) - offsetU, 1e-9);
  BOOST_CHECK_SMALL(parameters(4) - offsetV, 1e-9);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ActsToGblTests)

BOOST_AUTO_TEST_CASE(StraightLineClosureCentral) {
  checkStraightLineClosure(std::numbers::pi / 2., 0.);
}

BOOST_AUTO_TEST_CASE(StraightLineClosureForward) {
  checkStraightLineClosure(0.7, 0.4);
  checkStraightLineClosure(0.05, 0.4);
}

/// Acts::CurvilinearSurface falls back to U = X x T for directions within
/// s_curvilinearProjTolerance of the z axis. There dT/dphi = sin(theta) U and
/// dT/dtheta = -V no longer hold, so a conversion that hard-codes them is
/// silently wrong. These two cases sit inside that region.
BOOST_AUTO_TEST_CASE(StraightLineClosureFrameSwitch) {
  checkStraightLineClosure(1.0e-3, 0.4);
  checkStraightLineClosure(std::numbers::pi - 1.0e-3, 0.4);
}

/// Non-normal incidence exercises the path corrections in the curvilinear <->
/// bound coordinate changes. Every sensor Acts crosses at an angle depends on
/// them, so this is the case that matters for a real barrel geometry.
BOOST_AUTO_TEST_CASE(TiltedPlaneClosure) {
  for (double tiltDegrees : {5., 15., 30., 45., 60.}) {
    checkTiltedClosure(tiltDegrees * std::numbers::pi / 180.);
  }
}

BOOST_AUTO_TEST_SUITE_END()
