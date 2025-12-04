// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/GdmlDetectorConstruction.hpp"

#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"

#include <utility>

#include <G4GDMLParser.hh>

namespace ActsExamples {

GdmlDetectorConstruction::GdmlDetectorConstruction(
    std::string path, const Geant4ConstructionOptions& options)
    : G4VUserDetectorConstruction(),
      m_path(std::move(path)),
      m_options(options) {}

G4VPhysicalVolume* GdmlDetectorConstruction::Construct() {
  if (m_world == nullptr) {
    std::cout << "PF:::GdmlDetectorConstruction::Construct()" << std::endl;
    G4GDMLParser parser;
    // TODO how to handle errors
    std::cout << "PF::m_path:" << m_path << std::endl;
    parser.Read(m_path);
    m_world = parser.GetWorldVolume();
    std::cout << "PF::m_world" << m_world << std::endl;

    // Create regions
    for (const auto& regionCreator : m_options.regionCreators) {
      std::cout << "PF::Region Creator->buildRegion()" << std::endl;
      regionCreator->buildRegion();
    }
  }
  return m_world;
}

}  // namespace ActsExamples
