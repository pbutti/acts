// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AlignmentMillePede/ActsSolverFromMille.hpp"
#include "ActsExamples/AlignmentMillePede/MillePedeAlignmentSandbox.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsAlignmentMillePede, m) {
  {
    auto [alg, c] = declareAlgorithm<MillePedeAlignmentSandbox, IAlgorithm>(
        m, "MillePedeAlignmentSandbox");
    ACTS_PYTHON_STRUCT(c, milleOutput, inputMeasurements, inputTracks,
                       trackingGeometry, magneticField, fixModules,
                       discardUnconstrainedTrackPar, outFileInternalSolving,
                       outFileDecomposition, structures, outFileStructures);

    using AlignmentStructure = MillePedeAlignmentSandbox::AlignmentStructure;
    auto s = py::class_<AlignmentStructure>(alg, "AlignmentStructure")
                 .def(py::init<>())
                 .def(py::init([](const GeometryIdentifier& selector,
                                  const std::optional<Transform3>& transform) {
                        return AlignmentStructure{selector, transform};
                      }),
                      py::arg("selector"), py::arg("transform") = py::none());
    ACTS_PYTHON_STRUCT(s, selector, transform);
  }
  ACTS_PYTHON_DECLARE_ALGORITHM(ActsSolverFromMille, m, "ActsSolverFromMille",
                                milleInput, trackingGeometry, magneticField,
                                fixModules, outFile);
}
