// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsPlugins/GBL/GblError.hpp"

#include <cstddef>
#include <map>
#include <unordered_map>
#include <vector>

#include "GblPoint.h"
#include "GblTrajectory.h"

namespace ActsPlugins::ActsToGbl {

/// @brief Per-track-state input to the GBL conversion.
///
/// This is a plain struct, decoupled from the Acts track and
/// trajectory types. It lets the (Eigen-heavy) conversion itself live in a
/// single translation unit instead of being instantiated for every track
/// container type, mirroring what the alignment kernel does with
/// @c finaliseTrackAlignState.
struct GblTrackState {
  /// The surface this state lives on
  const Acts::Surface* surface = nullptr;

  /// Bound-to-bound jacobian from the *previous* state along the track.
  /// Identity for the first state.
  Acts::BoundMatrix jacobian = Acts::BoundMatrix::Identity();

  /// Reference (smoothed) bound parameters on @c surface
  /// Using the smoothed reference is only an approximation at second-order.
  /// Should be fine for small correction wrt a reference trajectory
  /// However the right thing would be propagating a single set of parameters
  Acts::BoundVector boundParameters = Acts::BoundVector::Zero();

  /// The same reference parameters in free coordinates. Kept explicitly
  /// because position/direction/qOverP are needed for the curvilinear frame,
  /// the alignment derivatives and the material lookup.
  Acts::FreeVector freeParameters = Acts::FreeVector::Zero();

  /// Whether this state carries a measurement. States that do not are
  /// converted into pure scattering points.
  bool isMeasurement = false;

  /// Dimension of the measurement (0 if @c isMeasurement is false)
  std::size_t measurementDim = 0;

  /// Which bound parameters are measured. Only the first @c measurementDim
  /// entries are meaningful.
  Acts::BoundSubspaceIndices measurementIndices{};

  /// Calibrated measurement, first @c measurementDim entries used
  /// TODO:: Not sure of this. This should be the measurementDim, not eBoundSize
  Acts::Vector<Acts::eBoundSize> calibrated =
      Acts::Vector<Acts::eBoundSize>::Zero();

  /// Calibrated measurement covariance, leading @c measurementDim block used
  Acts::SquareMatrix<Acts::eBoundSize> calibratedCovariance =
      Acts::SquareMatrix<Acts::eBoundSize>::Zero();

  /// Convenience accessors on the reference parameters
  Acts::Vector3 position() const {
    return freeParameters.segment<3>(Acts::eFreePos0);
  }
  Acts::Vector3 direction() const {
    return freeParameters.segment<3>(Acts::eFreeDir0);
  }
  double qOverP() const { return freeParameters[Acts::eFreeQOverP]; }
};

/// @brief Options steering the Acts -> GBL conversion
struct GblConversionOptions {
  /// Add thin scatterers derived from the surface material
  bool addScatterers = true;

  /// Add Millepede global derivatives for alignable surfaces
  bool addGlobals = true;

  /// Which alignment degrees of freedom to expose as globals
  ActsAlignment::AlignmentMask alignMask = ActsAlignment::AlignmentMask::All;

  /// Particle hypothesis used for the multiple-scattering estimate.
  /// Only relevant when @c addScatterers is enabled.
  /// Default is pion
  Acts::ParticleHypothesis particleHypothesis = Acts::ParticleHypothesis::pion();

  /// Minimal measurement precision accepted by GBL
  double minPrecision = 0.;

  /// Scattering angles below this value (rad) are treated as "no material"
  /// and the scatterer is skipped, to avoid injecting infinite precision.
  double minScatteringTheta = 1e-12;
};

/// @brief Outcome of the Acts -> GBL conversion
struct GblConversionResult {
  /// The GBL points, ordered along the trajectory
  std::vector<gbl::GblPoint> points;

  /// Maps the (1-based) GBL label of every point that carries a measurement
  /// back to its surface
  std::map<unsigned int, const Acts::Surface*> measurementSurfaces;

  /// Number of points that received alignment global derivatives
  std::size_t nGlobalPoints = 0;

  /// Number of points that received a scatterer
  std::size_t nScatterers = 0;
};

/// @brief Build the GBL point list from the extracted track states.
///
/// The points are expressed in the curvilinear system with local parameters
/// (q/p, u'_1, u'_2, u_1, u_2), which is what @c gbl::GblPoint expects.
///
/// @param gctx The geometry context the track was fitted in
/// @param states The track states, ordered along the trajectory
/// @param idxedAlignSurfaces Alignable surfaces and their global index. Used
///        to build the Millepede labels, with the exact same convention as
///        @c ActsPlugins::ActsToMille (label = eAlignmentSize * index + dof + 1)
/// @param opts Conversion options
/// @param logger Logger instance
///
/// @return The GBL points, or an error
Acts::Result<GblConversionResult> buildGblPoints(
    const Acts::GeometryContext& gctx,
    const std::vector<GblTrackState>& states,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces,
    const GblConversionOptions& opts, const Acts::Logger& logger);

/// @brief Extract the GBL inputs from a fitted Acts track.
///
/// @note The states are returned ordered along the trajectory, i.e. the
/// reverse of the native @c trackStatesReversed() storage order. Each state's
/// jacobian is the one the fitter stored, i.e. the transport from its
/// predecessor in *filtering* order. This assumes the fit ran in the
/// along-track direction, which holds for the standard forward KF/CKF.
/// TODO:: CHECK WHAT HAPPENS IN THE TWO WAY TRACKS
///
/// @param gctx The geometry context the track was fitted in
/// @param track The fitted track
///
/// @return The extracted states, ordered along the trajectory
template <typename track_proxy_t>
std::vector<GblTrackState> extractGblTrackStates(
    const Acts::GeometryContext& gctx, const track_proxy_t& track) {
  std::vector<GblTrackState> states;
  states.reserve(track.nTrackStates());

  for (const auto& ts : track.trackStatesReversed()) {
    // without smoothed parameters we have no reference trajectory to
    // linearise around
    if (!ts.hasSmoothed()) {
      continue;
    }
    const auto typeFlags = ts.typeFlags();
    const bool isMeasurement = typeFlags.isMeasurement();
    // pure hole / passive states carry nothing we can use
    // TODO:: If hole add as scatter?!
    
    if (!isMeasurement && !typeFlags.isMaterial()) {
      continue;
    }

    GblTrackState state;
    state.surface = &ts.referenceSurface();
    state.jacobian = ts.jacobian();
    state.boundParameters = ts.smoothed();
    state.freeParameters = Acts::MultiTrajectoryHelpers::freeSmoothed(gctx, ts);
    state.isMeasurement = isMeasurement;

    if (isMeasurement) {
      const std::size_t measdim = ts.calibratedSize();
      state.measurementDim = measdim;
      state.measurementIndices = ts.projectorSubspaceIndices();
      state.calibrated.head(measdim) = ts.effectiveCalibrated();
      state.calibratedCovariance.topLeftCorner(measdim, measdim) =
          ts.effectiveCalibratedCovariance();
    }

    states.push_back(std::move(state));
  }

  // trackStatesReversed() walks from the last to the first state; GBL needs
  // the points ordered along the trajectory.
  std::reverse(states.begin(), states.end());
  return states;
}

/// @brief Convenience wrapper: convert a fitted Acts track into GBL points.
///
/// @see extractGblTrackStates and buildGblPoints
template <typename track_proxy_t>
Acts::Result<GblConversionResult> convertTrackToGbl(
    const Acts::GeometryContext& gctx, const track_proxy_t& track,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces,
    const GblConversionOptions& opts, const Acts::Logger& logger) {
  GblConversionOptions localOpts = opts;
  localOpts.particleHypothesis = track.particleHypothesis();
  return buildGblPoints(gctx, extractGblTrackStates(gctx, track),
                        idxedAlignSurfaces, localOpts, logger);
}

}  // namespace ActsPlugins::ActsToGbl
