// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "ActsExamples/TrackJets/TrackJetsAlgorithm.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>



namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTrackJets(Context& ctx) {
  auto mex = ctx.get("examples");

  /*
  {
    using Alg = ActsExamples::TrackJetsAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg,IAlgorithm,std::shared_ptr<Alg>>(
            mex,"TrackJetsAlgorithm")
        .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
             py::arg("config"), py::arg("level"))
        .def_property_readonly("config",&Alg::config);


    auto c = py::class_<Config>(alg,"Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrackCollection);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(simParticles);
    ACTS_PYTHON_MEMBER(radius);
    ACTS_PYTHON_MEMBER(trackMass);
    ACTS_PYTHON_MEMBER(tj_minPt);
    ACTS_PYTHON_MEMBER(outputTrackJets);
    ACTS_PYTHON_STRUCT_END();
    
    }
  */
  
  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackJetsAlgorithm, mex,"TrackJetsAlgorithm",
      inputTrackCollection, inputTrajectories, simParticles,inputMeasurementParticlesMap, recoVertices,
      truthMatchProbability, radius, trackMass, tj_minPt, outputTrackJets);
}

} // namespace Acts:Python
