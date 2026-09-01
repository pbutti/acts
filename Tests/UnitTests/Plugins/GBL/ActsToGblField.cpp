// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsPlugins/GBL/ActsToGbl.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsPlugins::ActsToGbl;

namespace {

/// Outcome of one closure run.
struct ClosureOutcome {
  double chi2 = 0.;
  int ndf = 0;
  /// Fitted GBL local parameters at the first point, (q/p, u'_1, u'_2, u_1,
  /// u_2)
  Eigen::VectorXd parameters;
  /// The same quantities read off the two propagated trajectories directly
  Eigen::Matrix<double, 5, 1> expected = Eigen::Matrix<double, 5, 1>::Zero();
};

/// @brief Closure test of the Acts -> GBL conversion in a magnetic field.
///
/// Unlike the straight-line tests, every ingredient here comes from the real
/// propagator: the reference trajectory is a helix produced by @c EigenStepper,
/// and the per-point jacobians are the stepper's own bound-to-bound transport
/// jacobians -- exactly what a Kalman fit would store on its track states.
///
/// The "true" track is a second helix launched from slightly perturbed start
/// parameters on the first plane. Its intersections with the (tilted) planes
/// are the measurements. Because both tracks obey the same equation of motion,
/// a correct set of point-to-point jacobians absorbs the difference to first
/// order and the fit must return a vanishing chi2.
///
/// This is the configuration that exercises the Lorentz term of the path
/// derivative, @c dT/ds = (q/p) T x B. That term is zero in a field-free setup,
/// so the straight-line closure tests cannot see it at all. Here it enters
/// through the (dir, pos) block of the path correction @c (1 + D sigma_S), and
/// dropping it (@c useField == false) leaves the jacobians inconsistent with
/// the transport the stepper actually performed.
///
/// @param tiltAngle Tilt of the planes wrt the reference direction
/// @param momentum Momentum of the reference track
/// @param useField Whether to hand the field to the conversion
ClosureOutcome runFieldClosure(double tiltAngle, double momentum,
                               bool useField) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;
  auto logger = getDefaultLogger("ActsToGbl", Logging::WARNING);

  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2._T));

  const double theta = 1.1;
  const double phi = 0.3;
  const Vector3 dir0 = makeDirectionFromPhiTheta(phi, theta);
  const double qOverP = 1. / momentum;

  // curvilinear frame at the start; the planes are tilted about its U axis
  const Vector3 axisU0 = CurvilinearSurface(dir0).referenceFrame().col(0);

  const std::vector<double> nominal = {0., 100._mm, 200._mm, 300._mm, 400._mm};

  std::vector<std::shared_ptr<PlaneSurface>> surfaces;
  for (std::size_t i = 0; i < nominal.size(); ++i) {
    // alternate the tilt so the planes are not all parallel
    const double sign = (i % 2 == 0) ? 1. : -1.;
    const Vector3 normal =
        (std::cos(tiltAngle) * dir0 + sign * std::sin(tiltAngle) * axisU0)
            .normalized();
    const Vector3 seed =
        (std::abs(normal.z()) < 0.9 ? Vector3::UnitZ() : Vector3::UnitX());
    RotationMatrix3 rotation;
    rotation.col(0) = seed.cross(normal).normalized();
    rotation.col(1) = normal.cross(rotation.col(0));
    rotation.col(2) = normal;

    Transform3 transform = Transform3::Identity();
    transform.linear() = rotation;
    transform.translation() = nominal[i] * dir0;
    surfaces.push_back(Surface::makeShared<PlaneSurface>(
        transform, std::make_shared<RectangleBounds>(10000., 10000.)));
  }

  using Stepper = EigenStepper<>;
  using PropagatorType = Propagator<Stepper, VoidNavigator>;
  PropagatorType propagator(Stepper(bField), VoidNavigator(),
                            getDefaultLogger("Propagator", Logging::WARNING));

  /// Bound parameters and the incoming per-leg jacobian on every plane.
  struct Leg {
    BoundVector pars = BoundVector::Zero();
    BoundMatrix jacobian = BoundMatrix::Identity();
  };

  // Propagate one set of start parameters through the whole plane stack,
  // collecting the stepper's bound-to-bound jacobian for each leg.
  auto runThrough =
      [&](const BoundVector& startPars) -> std::optional<std::vector<Leg>> {
    std::vector<Leg> legs;
    legs.push_back({startPars, BoundMatrix::Identity()});

    // a covariance has to be present for the propagator to transport (and
    // therefore report) the jacobian at all
    BoundTrackParameters current(surfaces[0], startPars,
                                 BoundMatrix::Identity(),
                                 ParticleHypothesis::pion());

    for (std::size_t i = 1; i < surfaces.size(); ++i) {
      PropagatorType::Options<> options(gctx, mctx);
      // VoidNavigator means nothing bounds the step size, so the stepper can
      // overshoot the target plane badly (the default initial step is 10 m).
      // A real fit always has a navigator doing this.
      options.stepping.maxStepSize = 10._mm;
      auto res = propagator.propagate(current, *surfaces[i], options);
      if (!res.ok() || !res->endParameters.has_value() ||
          !res->transportJacobian.has_value()) {
        return std::nullopt;
      }
      current = *res->endParameters;
      legs.push_back({current.parameters(), *res->transportJacobian});
    }
    return legs;
  };

  BoundVector refStart = BoundVector::Zero();
  refStart[eBoundPhi] = phi;
  refStart[eBoundTheta] = theta;
  refStart[eBoundQOverP] = qOverP;

  // the true track: a second helix, offset and tilted wrt the reference
  BoundVector trueStart = refStart;
  trueStart[eBoundLoc0] += 0.05_mm;
  trueStart[eBoundLoc1] += -0.03_mm;
  trueStart[eBoundPhi] += 2.0e-5;
  trueStart[eBoundTheta] += -1.5e-5;

  auto refLegs = runThrough(refStart);
  auto trueLegs = runThrough(trueStart);
  BOOST_REQUIRE(refLegs.has_value());
  BOOST_REQUIRE(trueLegs.has_value());

  const double sigma = 0.010_mm;

  std::vector<GblTrackState> states;
  for (std::size_t i = 0; i < surfaces.size(); ++i) {
    GblTrackState state;
    state.surface = surfaces[i].get();
    state.jacobian = (*refLegs)[i].jacobian;
    state.boundParameters = (*refLegs)[i].pars;
    state.freeParameters = transformBoundToFreeParameters(*surfaces[i], gctx,
                                                          (*refLegs)[i].pars);

    state.isMeasurement = true;
    state.measurementDim = 2;
    state.measurementIndices[0] = eBoundLoc0;
    state.measurementIndices[1] = eBoundLoc1;
    state.calibrated[0] = (*trueLegs)[i].pars[eBoundLoc0];
    state.calibrated[1] = (*trueLegs)[i].pars[eBoundLoc1];
    state.calibratedCovariance(0, 0) = sigma * sigma;
    state.calibratedCovariance(1, 1) = sigma * sigma;

    states.push_back(state);
  }

  GblConversionOptions options;
  options.addScatterers = false;
  options.addGlobals = false;
  options.particleHypothesis = ParticleHypothesis::pion();
  options.magneticField = useField ? bField : nullptr;

  auto result = buildGblPoints(gctx, mctx, states, {}, options, *logger);
  BOOST_REQUIRE(result.ok());

  // curved trajectory: q/p is part of the local fit
  gbl::GblTrajectory trajectory(result->points, true);
  BOOST_REQUIRE(trajectory.isValid());

  ClosureOutcome outcome;
  double lostWeight = 0.;
  trajectory.fit(outcome.chi2, outcome.ndf, lostWeight);

  Eigen::MatrixXd covariance;
  // compare at an interior point, so the truth propagation below is a
  // well-defined forward propagation
  constexpr std::size_t compareAt = 2;
  trajectory.getResults(static_cast<int>(compareAt) + 1, outcome.parameters,
                        covariance);

  // Truth for the comparison, taken from the propagator rather than from an
  // analytic expansion. GBL's local parameters at a point are the deviation of
  // the true trajectory from the reference, measured *on the curvilinear plane*
  // through the reference point -- not on the sensor. Propagating the true
  // track onto exactly that plane therefore reads them off directly: its local
  // coordinates there are the offsets, and its direction gives the slopes.
  //
  // Doing this analytically instead would need the along-track slide between
  // the sensor crossing and the curvilinear crossing, which is precisely the
  // term the path correction supplies -- so it would assume the answer.
  const FreeVector refFree = transformBoundToFreeParameters(
      *surfaces[compareAt], gctx, (*refLegs)[compareAt].pars);
  const Vector3 refPos = refFree.segment<3>(eFreePos0);
  const Vector3 refDir = refFree.segment<3>(eFreeDir0);

  const CurvilinearSurface curvilinear(refPos, refDir);
  const RotationMatrix3 frame = curvilinear.referenceFrame();
  auto curvilinearSurface = curvilinear.surface();

  BoundTrackParameters trueStartPars(surfaces[0], trueStart,
                                     BoundMatrix::Identity(),
                                     ParticleHypothesis::pion());
  PropagatorType::Options<> truthOptions(gctx, mctx);
  truthOptions.stepping.maxStepSize = 10._mm;
  auto truthRes =
      propagator.propagate(trueStartPars, *curvilinearSurface, truthOptions);
  BOOST_REQUIRE(truthRes.ok());
  BOOST_REQUIRE(truthRes->endParameters.has_value());

  const auto& truthPars = *truthRes->endParameters;
  const Vector3 deltaDir = truthPars.direction() - refDir;

  outcome.expected[0] = truthPars.parameters()[eBoundQOverP] - refFree[eFreeQOverP];
  outcome.expected[1] = deltaDir.dot(frame.col(0));
  outcome.expected[2] = deltaDir.dot(frame.col(1));
  outcome.expected[3] = truthPars.parameters()[eBoundLoc0];
  outcome.expected[4] = truthPars.parameters()[eBoundLoc1];

  return outcome;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ActsToGblFieldTests)

/// Check the fitted local parameters against the propagated truth.
void checkAgainstTruth(const ClosureOutcome& outcome) {
  // 5 planes x 2D measurements - 5 helix parameters
  BOOST_CHECK_EQUAL(outcome.ndf, 5);
  BOOST_CHECK_SMALL(outcome.chi2, 1e-6);
  BOOST_CHECK_SMALL(outcome.parameters(0) - outcome.expected[0], 1e-6);
  BOOST_CHECK_SMALL(outcome.parameters(1) - outcome.expected[1], 1e-6);
  BOOST_CHECK_SMALL(outcome.parameters(2) - outcome.expected[2], 1e-6);
  BOOST_CHECK_SMALL(outcome.parameters(3) - outcome.expected[3], 1e-5);
  BOOST_CHECK_SMALL(outcome.parameters(4) - outcome.expected[4], 1e-5);
}

/// Closure across a range of incidence angles.
BOOST_AUTO_TEST_CASE(FieldClosureTilted) {
  for (double tiltDegrees : {5., 15., 30., 45.}) {
    checkAgainstTruth(
        runFieldClosure(tiltDegrees * std::numbers::pi / 180., 1._GeV, true));
  }
}

/// Closure across a momentum range. The Lorentz term scales with q/p, so a
/// low-momentum track is where it matters most.
BOOST_AUTO_TEST_CASE(FieldClosureMomenta) {
  const double tilt = 30. * std::numbers::pi / 180.;
  for (double momentum : {0.5_GeV, 1._GeV, 10._GeV}) {
    checkAgainstTruth(runFieldClosure(tilt, momentum, true));
  }
}

/// Leaving the field unset falls back to the linear track model, i.e. it drops
/// the Lorentz term dT/ds = (q/p) T x B from the path derivative D.
///
/// That omission is a pure shear of the GBL local parameters: it leaves the
/// offsets untouched and only shifts the slopes. The reason is structural -- a
/// plane's sigma_S has only a position block (Surface::freeToPathDerivative),
/// while the dropped term sits in D's direction block, so the product D sigma_S
/// can only populate the (slope, offset) corner of the path correction.
///
/// Two consequences, both asserted here:
///
///  - chi2 and the offsets are *completely* insensitive to it, bit for bit.
///    GBL measurements constrain offsets only, and the shear leaves the offset
///    rows of the jacobian exactly alone. So no chi2-based closure test can
///    ever detect a wrong path derivative, however tilted the planes or strong
///    the field. This is why the field-free closure tests in
///    ActsToGblStraightLine.cpp cannot cover this code path.
///  - the slope along U is wrong by three orders of magnitude more than the
///    field-aware conversion. That propagates into everything that constrains
///    slope space: scatterer kink precisions, the alignment global derivatives,
///    and the fitted direction itself.
///
/// With B along z the dropped term (q/p) T x B is parallel to the curvilinear
/// U axis (which is z x T), so the whole error lands in u'_1 and u'_2 is
/// untouched -- a useful cross-check that the effect seen here is really the
/// Lorentz term and not a generic tolerance failure.
BOOST_AUTO_TEST_CASE(FieldTermMovesOnlyTheSlopes) {
  const double tilt = 30. * std::numbers::pi / 180.;

  for (double momentum : {0.5_GeV, 1._GeV}) {
    auto withField = runFieldClosure(tilt, momentum, true);
    auto withoutField = runFieldClosure(tilt, momentum, false);

    // chi2, the offsets and q/p are blind to the missing term
    BOOST_CHECK_CLOSE(withField.chi2, withoutField.chi2, 1e-6);
    BOOST_CHECK_SMALL(withField.parameters(0) - withoutField.parameters(0),
                      1e-14);
    BOOST_CHECK_SMALL(withField.parameters(3) - withoutField.parameters(3),
                      1e-12);
    BOOST_CHECK_SMALL(withField.parameters(4) - withoutField.parameters(4),
                      1e-12);
    // B is along z, so is the slope transverse to it
    BOOST_CHECK_SMALL(withField.parameters(2) - withoutField.parameters(2),
                      1e-12);

    // the slope along U is not blind: only the field-aware conversion
    // reproduces the propagated truth
    const double devWith =
        std::abs(withField.parameters(1) - withField.expected[1]);
    const double devWithout =
        std::abs(withoutField.parameters(1) - withoutField.expected[1]);
    BOOST_TEST_MESSAGE("p = " << momentum << ": |dev u'_1| with field = "
                              << devWith << ", without = " << devWithout);
    BOOST_CHECK_SMALL(devWith, 1e-6);
    BOOST_CHECK_GT(devWithout, 1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
