// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <system_error>

namespace ActsPlugins {

/// Errors raised while converting an Acts track into a GBL trajectory
enum class GblError {
  TooFewStates = 1,
  NoMeasurementOnTrack = 2,
  UnsupportedMeasurement = 3,
  SingularProjection = 4,
  MissingSmoothedState = 5,
};

namespace detail {
class GblErrorCategory : public std::error_category {
 public:
  const char* name() const noexcept final { return "GblError"; }
  std::string message(int c) const final {
    switch (static_cast<GblError>(c)) {
      case GblError::TooFewStates:
        return "Too few track states to build a GBL trajectory";
      case GblError::NoMeasurementOnTrack:
        return "No measurement states on the track";
      case GblError::UnsupportedMeasurement:
        return "Measurement is not expressed in bound local coordinates";
      case GblError::SingularProjection:
        return "Local-to-measurement projection is singular";
      case GblError::MissingSmoothedState:
        return "Track state has no smoothed parameters";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

extern inline const detail::GblErrorCategory& GblErrorCategory() {
  static detail::GblErrorCategory c;
  return c;
}

inline std::error_code make_error_code(ActsPlugins::GblError e) {
  return {static_cast<int>(e), ActsPlugins::GblErrorCategory()};
}

}  // namespace ActsPlugins

namespace std {
template <>
struct is_error_code_enum<ActsPlugins::GblError> : std::true_type {};
}  // namespace std
