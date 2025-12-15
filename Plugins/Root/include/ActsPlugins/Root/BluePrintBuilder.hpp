// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/TGeoDetectorElement.hpp"

#include <functional>
#include <memory>
#include <regex>
#include <string>
#include <utility>

#include "TGeoManager.h"

namespace Acts {
namespace Experimental {
class CylinderContainerBlueprintNode;
}  // namespace Experimental
}  // namespace Acts

namespace ActsPlugins {

class BlueprintBuilder {
 public:
  // Here I could use a derived class to return a whatever detector element that
  // inherits from this TGeoDetectorElement
  using ElementFactory = std::function<std::shared_ptr<TGeoDetectorElement>(
      const TGeoDetectorElement::Identifier& identifier,
      const TGeoNode& detElement, const std::string& axes, double lengthScale)>;

  static std::shared_ptr<TGeoDetectorElement> defaultElementFactory(
      const TGeoDetectorElement::Identifier& identifier,
      const TGeoNode& detElement, const std::string& axes, double lengthScale);

  struct Config {
    ElementFactory elementFactory = defaultElementFactory;
    // TODO:: This should be directly the TGeo world volume
    std::string path = "";
    double lengthScale = 1.0;
    TGeoNode* world;
  };

  explicit BlueprintBuilder(const Config& cfg,
                            std::unique_ptr<const Acts::Logger> logger_ =
                                Acts::getDefaultLogger("BlueprintBuilder",
                                                       Acts::Logging::INFO))
      : m_cfg(cfg), m_logger{std::move(logger_)} {
    if (m_cfg.path.empty()) {
      throw std::invalid_argument(
          "BlueprintBuilder:: path to TGeo root file is null");
    }
  }

  std::shared_ptr<TGeoDetectorElement> createDetectorElement(
      const TGeoDetectorElement::Identifier& identifier,
      const TGeoNode& detElement, const std::string& axes) const;

  // std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
  //       const TGeoNode& tGeoNode, const std::string& axes,
  //       std::optional<std::string> layerAxes = std::nullopt) const;

  // std::shared_ptr<Acts::Experimental::StaticBlueprintNode> makeBeampipe()
  // const;

  //[[deprecated("Consider using .layerHelper() to produce the layers")]]

  // LayerHelper layerHelper();

  const TGeoNode* findDetElementByName(
     const TGeoNode* parent, const std::string& name);

  const TGeoNode* findDetElementByName(
     const std::string& name);

  //static std::vector<TGeoNode> findDetElementByNamePattern(
  //   const TGeoNode& parent, const std::regex& pattern);

  //static std::vector<TGeoNode> findDetElementByNamePattern(
  //  const std::regex& pattern);

  const Acts::Logger& logger() const { return *m_logger; }

 private:
  // static std::vector<TGeoNode> resolveSensitives(
  //   const TGeoNode& detElement);

  TGeoNode* world() const;
  
  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger;
};
}  // namespace ActsPlugins
