// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>
#include <unordered_map>

#include "Mille/IMilleReader.h"
#include "Mille/MilleDecoder.h"
#include "Mille/MilleRecord.h"

namespace ActsPlugins::ActsToMille {

/// The MilleRecord is Millepede's interface for
/// writing out alignment fit inputs. It can be instantiated
/// using Mille::spawnMilleRecord(desired_file_name),
/// provided by Mille/MilleFactory.h
using Mille::MilleRecord;

/// Link from an aligned surface to the composite structure (e.g. a stave or
/// a layer) whose alignment parameters replace the surface ones in the Mille
/// record.
struct CompositeLink {
  /// Index of the structure. Its Mille labels are 6 * index + dof + 1, with
  /// dof the local-frame parameter (dT0, dT1, dT2, dW0, dW1, dW2).
  std::size_t structureIndex = 0;
  /// Jacobian d(ACTS alignment parameters of the surface) / d(local-frame
  /// alignment parameters of the structure)
  Acts::AlignmentMatrix jacobian = Acts::AlignmentMatrix::Identity();
};

/// Composite structure of each aligned surface
using CompositeMap = std::unordered_map<const Acts::Surface*, CompositeLink>;

/// @brief Link a surface to a composite structure.
/// @param structureIndex: Index of the structure, defining its Mille labels.
/// @param structureTransform: Local-to-global transform of the structure frame,
/// in which its alignment parameters are defined.
/// @param surfaceTransform: Local-to-global transform of the surface.
/// A surface can be its own structure, which gives Mille labels for
/// local-frame sensor alignment parameters.
CompositeLink makeCompositeLink(std::size_t structureIndex,
                                const Acts::Transform3& structureTransform,
                                const Acts::Transform3& surfaceTransform);

/// @brief Dump a Kalman track encoded as a TrackAlignmentState into
/// a Mille record.
/// @param state: Alignment state to dump.
/// @param record: Mille record to write to.
/// Note: Not very efficient - we have to "un-fit" the kalman track.
/// Used for R&D, recommending the GBL track model (under development)
/// for production use.
/// @param removeUnconstrainedTrackPar If enabled, will remove
/// poorly constrained parameters from the (local) track fits.
/// @param composites: [optional]: If given, write the derivatives w.r.t. the
/// alignment parameters of the composite structures instead of those of the
/// surfaces, summing the contributions of all surfaces of a structure.
/// Aligned surfaces missing from the map are not written.
void dumpToMille(const ActsAlignment::detail::TrackAlignmentState& state,
                 MilleRecord& record, bool removeUnconstrainedTrackPar = true,
                 const CompositeMap* composites = nullptr);

/// @brief read one record (= track or (constrained) track pair) from
/// a Mille binary into the equivalent matrices of a TrackAlignmentState.
/// Allows to use Mille to collect tracks across multiple events and
/// align them with the ACTS solver, and to validate the outputs of dumpToMille.
/// @param reader: A Mille Reader, connected to a valid input file.
/// @param targetState: The TrackAlignmentState to populate.
/// @param idxedAlignSurfaces: [optional]: Indexed alignment surfaces from the geometry. If passed,
/// the internal `alignedSurfaces` member of the state will be configured to
/// link back to the correct surfaces.
/// @return a ReadResult enum with 3 possible states to indicate the outcome- ok / end-of-file / read-error.
/// The targetState will only be modified if the result is 'ok'.
Mille::MilleDecoder::ReadResult unpackMilleRecord(
    Mille::IMilleReader& reader,
    ActsAlignment::detail::TrackAlignmentState& targetState,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces);

/// Writes an alignment outcome into a text file in the format
/// used by Millepede. Allows the constants to be processed
/// in the same visualisation code use for MP-II results,
/// or to be read back as initial values for a fit with MP-II
void dumpAsMillepedeRes(const ActsAlignment::AlignmentResult& result,
                        std::ostream& out);

}  // namespace ActsPlugins::ActsToMille
