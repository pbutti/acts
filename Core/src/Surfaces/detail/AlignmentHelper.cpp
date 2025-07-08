// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/detail/AlignmentHelper.hpp"

#include <utility>

namespace Acts {

detail::RotationToAxes detail::rotationToLocalAxesDerivative(
    const RotationMatrix3& compositeRotation,
    const RotationMatrix3& relRotation) {
  // Suppose the local axes of the composite have small rotation first around
  // its original local x axis by alpha, then around its original local y by
  // beta, then last around its original local z by gamma, the new rotation
  // matrix of the composite is then compositeRotation*deltaRotation, where
  // deltaRotation has the following form:
  //  | cbeta*cgamma  salpha*sbeta*cgamma-calpha*sgamma calpha*sbeta*cgamma +
  //  salpha*sgamma|, | cbeta*sgamma  salpha*sbeta*sgamma+calpha*cgamma
  //  calpha*sbeta*sgamma-salpha*cgamma  |, | -sbeta        salpha*cbeta
  //  calpha*cbeta                       | where prefix 's' means sin and 'c'
  //  means cos, then:
  //  1) the derivatives of new local x axis of the composite
  //  w.r.t. (alpha, beta, gamma) is rotToCompositeLocalXAxis =
  //  compositeRotation* |0  0  0|,
  //                     |0  0  1|
  //                     |0 -1  0|
  //  2) the derivatives of new local y axis of the composite
  //  w.r.t. (alpha, beta, gamma) is rotToCompositeLocalYAxis =
  //  compositeRotation* |0  0 -1|,
  //                     |0  0  0|
  //                     |1  0  0|
  //  3) the derivatives of new local z axis of the composite
  //  w.r.t. (alpha, beta, gamma) is rotToCompositeLocalZAxis =
  //  compositeRotation* | 0  1  0|,
  //                     |-1  0  0|
  //                     | 0  0  0|

  // The object rotation is objectRotation = compositeRotation*relRotation, then
  // 1) the derivate of the new local
  // x axis of the object w.r.t. (alpha, beta, gamma) is
  // rotToCompositeLocalXAxis * relRotation(0,0) +
  // rotToCompositeLocalYAxis*relRotation(1,0) +
  // rotToCompositeLocalZAxis*relRotation(2,0),
  // 2) the derivate of the new local
  // y axis of the object w.r.t. (alpha, beta, gamma) is
  // rotToCompositeLocalXAxis * relRotation(0,1) +
  // rotToCompositeLocalYAxis*relRotation(1,1) +
  // rotToCompositeLocalZAxis*relRotation(2,1),
  // 3) the derivate of the new local
  // z axis of the object w.r.t. (alpha, beta, gamma) is
  // rotToCompositeLocalXAxis * relRotation(0,2) +
  // rotToCompositeLocalYAxis*relRotation(1,2) +
  // rotToCompositeLocalZAxis*relRotation(2,2),

  // Derivative of local x axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToCompositeLocalXAxis = RotationMatrix3::Zero();
  rotToCompositeLocalXAxis.col(0) = compositeRotation * Vector3(0, 0, 0);
  rotToCompositeLocalXAxis.col(1) = compositeRotation * Vector3(0, 0, -1);
  rotToCompositeLocalXAxis.col(2) = compositeRotation * Vector3(0, 1, 0);
  // Derivative of local y axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToCompositeLocalYAxis = RotationMatrix3::Zero();
  rotToCompositeLocalYAxis.col(0) = compositeRotation * Vector3(0, 0, 1);
  rotToCompositeLocalYAxis.col(1) = compositeRotation * Vector3(0, 0, 0);
  rotToCompositeLocalYAxis.col(2) = compositeRotation * Vector3(-1, 0, 0);
  // Derivative of local z axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3 rotToCompositeLocalZAxis = RotationMatrix3::Zero();
  rotToCompositeLocalZAxis.col(0) = compositeRotation * Vector3(0, -1, 0);
  rotToCompositeLocalZAxis.col(1) = compositeRotation * Vector3(1, 0, 0);
  rotToCompositeLocalZAxis.col(2) = compositeRotation * Vector3(0, 0, 0);

  RotationMatrix3 rotToLocalXAxis = RotationMatrix3::Zero();
  RotationMatrix3 rotToLocalYAxis = RotationMatrix3::Zero();
  RotationMatrix3 rotToLocalZAxis = RotationMatrix3::Zero();
  rotToLocalXAxis = rotToCompositeLocalXAxis * relRotation(0, 0) +
                    rotToCompositeLocalYAxis * relRotation(1, 0) +
                    rotToCompositeLocalZAxis * relRotation(2, 0);
  rotToLocalYAxis = rotToCompositeLocalXAxis * relRotation(0, 1) +
                    rotToCompositeLocalYAxis * relRotation(1, 1) +
                    rotToCompositeLocalZAxis * relRotation(2, 1);
  rotToLocalZAxis = rotToCompositeLocalXAxis * relRotation(0, 2) +
                    rotToCompositeLocalYAxis * relRotation(1, 2) +
                    rotToCompositeLocalZAxis * relRotation(2, 2);

  return {std::move(rotToLocalXAxis), std::move(rotToLocalYAxis),
          std::move(rotToLocalZAxis)};
}

  // compositeTransform is the global to composite      G2Composite
  // componentTransform is the composite to component   Composite2Component
  ActsMatrix<6,6> FrameJacobian(const Transform3& compositeTransform,
				const Transform3& componentTransform) {
    
    
    ActsMatrix<6,6> jacobian = ActsMatrix<6,6>::Zero();
    
    
    // Set the top-left 3x3 block of mat to R(g2l)
    jacobian.block<3, 3>(0, 0) = componentTransform.linear();
    jacobian.block<3, 3>(3, 3) = componentTransform.linear();

    return jacobian;
    
  }
  
  ActsMatrix<3,3> dtdA(const ActsVector3& T0,
		       const RotationMatrix3& R) {
    
    Acts::Vector3 xs{1,0,0};
    Acts::Vector3 ys{1,0,0};
    Acts::Vector3 zs{1,0,0};
    
    Acts::Vector3 xshat = R * xs;
    Acts::Vector3 yshat = R * ys;
    Acts::Vector3 zshat = R * zs;
    
    ActsMatrix<3,3> jacPosRot = ActsMatrix<3,3>::Zero();
    
    Acts::Vector3 Row0 = T0.cross(xs);
    Acts::Vector3 Row1 = T0.cross(ys);
    Acts::Vector3 Row2 = T0.cross(zs);
    
    jacPosRot.row(0) = Row0.transpose();
    jacPosRot.row(1) = Row1.transpose();
    jacPosRot.row(2) = Row2.transpose();
  }
  
}  // namespace Acts
