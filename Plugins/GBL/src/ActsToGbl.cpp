// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/GBL/ActsToGbl.hpp"

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include <cmath>

namespace ActsPlugins::ActsToGbl {

namespace {

/// GBL orders its local parameters as (q/p, u'_1, u'_2, u_1, u_2), while Acts
/// curvilinear parameters are (loc0, loc1, phi, theta, qOverP). The two are
/// related by a fixed, purely local transformation.
///
/// The offsets need no conversion: Acts' curvilinear loc0/loc1 are already the
/// offsets along the curvilinear U and V axes (the bound-to-free jacobian sets
/// its position block to [U V]).
///
/// The slopes do. GBL's u'_1, u'_2 are the components of a direction variation
/// along U and V, so they follow from projecting dT onto the frame:
///
///     u'_1 = U . (dT/dphi dphi + dT/dtheta dtheta)
///     u'_2 = V . (dT/dphi dphi + dT/dtheta dtheta)
///
/// For the usual frame U = Z x T, V = T x U this reduces to the familiar
/// u'_1 = sin(theta) dphi and u'_2 = -dtheta, because dT/dphi = sin(theta) U
/// and dT/dtheta = -V. That shortcut is deliberately *not* hard-coded here:
/// Acts::CurvilinearSurface falls back to U = X x T for directions within
/// s_curvilinearProjTolerance of the z axis (|eta| >~ 6.4), where those two
/// identities no longer hold. Taking the projections from the jacobian Acts
/// itself builds keeps the conversion correct in both branches.
using Matrix5 = Eigen::Matrix<double, 5, 5>;
using Matrix2 = Eigen::Matrix2d;

/// The 2x2 block relating (dphi, dtheta) to the GBL slopes (u'_1, u'_2),
/// taken from the curvilinear bound-to-free jacobian so that it follows
/// whichever reference frame Acts picked for this direction.
Matrix2 slopeBlock(const Acts::Vector3& direction) {
  const Acts::CurvilinearSurface surface(direction);
  const Acts::RotationMatrix3 frame = surface.referenceFrame();
  const Acts::BoundToFreeMatrix boundToFree = surface.boundToFreeJacobian();

  const Acts::Vector3 dTdPhi =
      boundToFree.block<3, 1>(Acts::eFreeDir0, Acts::eBoundPhi);
  const Acts::Vector3 dTdTheta =
      boundToFree.block<3, 1>(Acts::eFreeDir0, Acts::eBoundTheta);

  Matrix2 block;
  block << frame.col(0).dot(dTdPhi), frame.col(0).dot(dTdTheta),
      frame.col(1).dot(dTdPhi), frame.col(1).dot(dTdTheta);
  return block;
}

/// Acts curvilinear -> GBL
Matrix5 actsToGblBasis(const Matrix2& slopes) {
  Matrix5 p = Matrix5::Zero();
  p(0, Acts::eBoundQOverP) = 1.;                    // q/p
  p.block<2, 2>(1, Acts::eBoundPhi) = slopes;       // u'_1, u'_2
  p(3, Acts::eBoundLoc0) = 1.;                      // u_1
  p(4, Acts::eBoundLoc1) = 1.;                      // u_2
  return p;
}

/// GBL -> Acts curvilinear
Matrix5 gblToActsBasis(const Matrix2& slopes) {
  Matrix5 p = Matrix5::Zero();
  p(Acts::eBoundLoc0, 3) = 1.;
  p(Acts::eBoundLoc1, 4) = 1.;
  p.block<2, 2>(Acts::eBoundPhi, 1) = slopes.inverse();
  p(Acts::eBoundQOverP, 0) = 1.;
  return p;
}

/// @brief Path derivative of the free parameters, d(free)/d(path length).
///
/// Uses the linear track model, i.e. the change of direction along the path is
/// neglected. This is the same approximation the Acts alignment kernel makes
/// when it builds @c alignmentToBoundDerivative, so the local fit and the
/// alignment derivatives stay consistent with each other.
Acts::FreeVector pathDerivative(const Acts::Vector3& direction) {
  Acts::FreeVector derivative = Acts::FreeVector::Zero();
  derivative.head<3>() = direction;
  return derivative;
}

/// @brief Coordinate change from curvilinear parameters onto a surface.
///
/// Not simply the surface's free-to-bound projection applied to the
/// curvilinear-to-free jacobian: the curvilinear plane and the surface pass
/// through the same point but are tilted with respect to each other, so a
/// varied track crosses them at points separated along the trajectory. The
/// path correction (I + dFree/ds * ds/dFree) accounts for that slide. Without
/// it the round trip curvilinear -> bound -> curvilinear loses
/// sin^2(incidence), i.e. everything except normal incidence is wrong.
///
/// @see Acts::detail::boundToBoundTransportJacobian, which applies the same
/// factor at its destination surface.
Acts::BoundMatrix curvilinearToBound(const Acts::GeometryContext& gctx,
                                     const Acts::Surface& surface,
                                     const Acts::Vector3& position,
                                     const Acts::Vector3& direction) {
  const Acts::FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(gctx, position, direction);
  return surface.freeToBoundJacobian(gctx, position, direction) *
         (Acts::FreeMatrix::Identity() + pathDerivative(direction) * freeToPath) *
         Acts::CurvilinearSurface(direction).boundToFreeJacobian();
}

/// @brief Coordinate change from a surface onto the curvilinear frame.
///
/// The mirror image of @c curvilinearToBound: here the curvilinear plane is
/// the destination, so it carries the path correction. Its normal is the track
/// direction, which is why the correction reduces to the -direction outer
/// product used by @c Acts::detail::boundToCurvilinearTransportJacobian.
Acts::BoundMatrix boundToCurvilinear(const Acts::GeometryContext& gctx,
                                     const Acts::Surface& surface,
                                     const Acts::Vector3& position,
                                     const Acts::Vector3& direction) {
  Acts::FreeToBoundMatrix freeToCurvilinear =
      Acts::CurvilinearSurface(direction).freeToBoundJacobian();
  freeToCurvilinear.topLeftCorner<6, 3>() +=
      (freeToCurvilinear * pathDerivative(direction)) *
      (-1.0 * direction).transpose();
  return freeToCurvilinear *
         surface.boundToFreeJacobian(gctx, position, direction);
}

/// Transform the fitter's bound-to-bound jacobian into the curvilinear,
/// GBL-ordered 5x5 jacobian between two consecutive points.
///
/// The bound jacobian goes from the bound frame of @p prev to that of @p curr.
/// Composing it with the two coordinate changes above turns it into a map from
/// curvilinear parameters at @p prev onto curvilinear parameters at @p curr,
/// which is then rebased into GBL's parameter ordering.
Matrix5 toGblJacobian(const Acts::GeometryContext& gctx,
                      const GblTrackState& prev, const GblTrackState& curr) {
  const Acts::Vector3 prevDir = prev.direction();
  const Acts::Vector3 currDir = curr.direction();

  const Acts::BoundMatrix curvJacobian =
      boundToCurvilinear(gctx, *curr.surface, curr.position(), currDir) *
      curr.jacobian *
      curvilinearToBound(gctx, *prev.surface, prev.position(), prevDir);

  // drop the time row/column and change basis into GBL's ordering. The
  // slope block is direction dependent, hence evaluated at either end.
  const Matrix5 curv5 = curvJacobian.topLeftCorner<5, 5>();
  return actsToGblBasis(slopeBlock(currDir)) * curv5 *
         gblToActsBasis(slopeBlock(prevDir));
}

}  // namespace

Acts::Result<GblConversionResult> buildGblPoints(
    const Acts::GeometryContext& gctx,
    const std::vector<GblTrackState>& states,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces,
    const GblConversionOptions& opts, const Acts::Logger& logger) {
  ACTS_VERBOSE("Converting " << states.size() << " track states to GBL points");

  if (states.size() < 2) {
    return Acts::Result<GblConversionResult>::failure(GblError::TooFewStates);
  }

  GblConversionResult result;
  result.points.reserve(states.size());

  std::size_t nMeasurements = 0;

  for (std::size_t iState = 0; iState < states.size(); ++iState) {
    const GblTrackState& state = states[iState];

    // (1) jacobian from the previous point. GBL expects the identity for the
    // very first point of the trajectory.
    const Matrix5 jacobian =
        (iState == 0) ? Matrix5::Identity()
                      : toGblJacobian(gctx, states[iState - 1], state);

    gbl::GblPoint point(jacobian);

    const Acts::Vector3 position = state.position();
    const Acts::Vector3 direction = state.direction();

    // (2) the measurement, if any
    if (state.isMeasurement) {
      const std::size_t measdim = state.measurementDim;

      // all measured components must be bound local coordinates for the
      // projection below to be meaningful
      for (std::size_t i = 0; i < measdim; ++i) {
        const auto idx = state.measurementIndices[i];
        if (idx != Acts::eBoundLoc0 && idx != Acts::eBoundLoc1) {
          ACTS_WARNING("Measurement component "
                       << i << " on surface " << state.surface->geometryId()
                       << " is not a bound local coordinate - skipping track");
          return Acts::Result<GblConversionResult>::failure(
              GblError::UnsupportedMeasurement);
        }
      }

      // curvilinear (U, V) axes and the measurement axes, both in global
      // coordinates
      const Acts::RotationMatrix3 curvFrame =
          Acts::CurvilinearSurface(direction).referenceFrame();
      const Acts::RotationMatrix3 measFrame =
          state.surface->referenceFrame(gctx, position, direction);

      // projection measurement -> local (curvilinear) and its inverse
      const Eigen::Matrix2d proM2l =
          curvFrame.leftCols<2>().transpose() * measFrame.leftCols<2>();
      const double det = proM2l.determinant();
      if (std::abs(det) < 1e-12) {
        ACTS_WARNING("Singular measurement projection on surface "
                     << state.surface->geometryId() << " (det = " << det
                     << ") - skipping track");
        return Acts::Result<GblConversionResult>::failure(
            GblError::SingularProjection);
      }
      const Eigen::Matrix2d proL2mFull = proM2l.inverse();

      // keep only the rows that are actually measured
      Acts::DynamicMatrix proL2m =
          Acts::DynamicMatrix::Zero(static_cast<int>(measdim), 2);
      Acts::DynamicVector residual =
          Acts::DynamicVector::Zero(static_cast<int>(measdim));
      for (std::size_t i = 0; i < measdim; ++i) {
        const auto idx = state.measurementIndices[i];
        proL2m.row(static_cast<int>(i)) = proL2mFull.row(idx);
        residual[static_cast<int>(i)] =
            state.calibrated[static_cast<int>(i)] - state.boundParameters[idx];
      }

      const Acts::DynamicMatrix precision =
          state.calibratedCovariance
              .topLeftCorner(static_cast<int>(measdim),
                             static_cast<int>(measdim))
              .inverse();

      point.addMeasurement(proL2m, residual, precision, opts.minPrecision);
      ++nMeasurements;

      // (3) the alignment global derivatives, if this surface is alignable
      if (opts.addGlobals) {
        if (auto it = idxedAlignSurfaces.find(state.surface);
            it != idxedAlignSurfaces.end()) {
          // linear track model, as in the Acts alignment kernel
          Acts::FreeVector pathDerivative = Acts::FreeVector::Zero();
          pathDerivative.head<3>() = direction;

          Acts::AlignmentToBoundMatrix alignToBound =
              state.surface->alignmentToBoundDerivative(gctx, position,
                                                        direction,
                                                        pathDerivative);
          ActsAlignment::detail::resetAlignmentDerivative(alignToBound,
                                                          opts.alignMask);

          // residual = measurement - prediction, hence the minus sign. This is
          // the same convention ActsToMille uses, so the labels and signs of
          // the resulting Millepede constants are directly comparable.
          Acts::DynamicMatrix globalDerivatives = Acts::DynamicMatrix::Zero(
              static_cast<int>(measdim), Acts::eAlignmentSize);
          std::vector<int> labels(Acts::eAlignmentSize, 0);

          for (std::size_t i = 0; i < measdim; ++i) {
            globalDerivatives.row(static_cast<int>(i)) =
                -alignToBound.row(state.measurementIndices[i]);
          }
          for (std::size_t iDof = 0; iDof < Acts::eAlignmentSize; ++iDof) {
            // Millepede labels are 1-based, matching ActsToMille
            labels[iDof] = static_cast<int>(
                it->second * Acts::eAlignmentSize + iDof + 1);
          }

          point.addGlobals(labels, globalDerivatives);
          ++result.nGlobalPoints;
        }
      }
    }

    // (4) the scatterer. Acts does not keep the material slab on the track
    // state (the Kalman fitter folds it into the covariance and drops it), so
    // it is re-evaluated from the surface material here.
    if (opts.addScatterers && state.surface->surfaceMaterial() != nullptr) {
      auto slabRes = Acts::detail::evaluateMaterialSlab(
          gctx, *state.surface, Acts::Direction::Forward(), position, direction,
          Acts::MaterialUpdateMode::FullUpdate);
      if (slabRes.ok() && !slabRes->isVacuum()) {
        const Acts::detail::PointwiseMaterialEffects effects =
            Acts::detail::computeMaterialEffects(
                *slabRes, opts.particleHypothesis, direction,
                static_cast<float>(state.qOverP()),
                /*multipleScattering=*/true, /*energyLoss=*/false,
                /*covTransport=*/true);
        // varianceTheta is theta0^2; in the curvilinear frame both projected
        // kink angles have that same variance.
        const double theta0Squared = effects.varianceTheta;
        if (theta0Squared >
            opts.minScatteringTheta * opts.minScatteringTheta) {
          const Eigen::Vector2d scatPrecision(1. / theta0Squared,
                                              1. / theta0Squared);
          point.addScatterer(Eigen::Vector2d::Zero(), scatPrecision);
          ++result.nScatterers;
        }
      }
    }

    result.points.push_back(std::move(point));

    // GBL labels are 1-based indices into the point list
    if (state.isMeasurement) {
      result.measurementSurfaces.emplace(
          static_cast<unsigned int>(result.points.size()), state.surface);
    }
  }

  if (nMeasurements == 0) {
    return Acts::Result<GblConversionResult>::failure(
        GblError::NoMeasurementOnTrack);
  }

  ACTS_VERBOSE("Built " << result.points.size() << " GBL points, "
                        << nMeasurements << " with measurements, "
                        << result.nGlobalPoints << " with globals, "
                        << result.nScatterers << " with scatterers");

  return result;
}

}  // namespace ActsPlugins::ActsToGbl
