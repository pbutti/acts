// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/BluePrintBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"

#include "TGeoManager.h"

namespace ActsPlugins {
std::shared_ptr<TGeoDetectorElement> BlueprintBuilder::defaultElementFactory(
    const TGeoDetectorElement::Identifier& identifier,
    const TGeoNode& detElement, const std::string& axes, double lengthScale) {
  return std::make_shared<TGeoDetectorElement>(identifier, detElement, axes,
                                               lengthScale);
}

std::shared_ptr<TGeoDetectorElement> BlueprintBuilder::createDetectorElement(
    const TGeoDetectorElement::Identifier& identifier,
    const TGeoNode& detElement, const std::string& axes) const {
  auto elem =
      m_cfg.elementFactory(identifier, detElement, axes, m_cfg.lengthScale);

  return elem;
}

}  // namespace ActsPlugins
