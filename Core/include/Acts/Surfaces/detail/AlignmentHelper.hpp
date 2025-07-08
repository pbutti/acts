// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <tuple>

namespace Acts::detail {

// The container for derivative of local frame axis w.r.t. its
// rotation parameters. The first element is for x axis, second for y axis and
// last for z axis
using RotationToAxes =
    std::tuple<RotationMatrix3, RotationMatrix3, RotationMatrix3>;

/// @brief Evaluate the derivative of local frame axes vector w.r.t.
/// its rotation around local x/y/z axis
/// @Todo: add parameter for rotation axis order
///
/// @param compositeRotation The rotation that help places the composite object being rotated
/// @param relRotation The relative rotation of the surface with respect to the composite object being rotated
///
/// @return Derivative of local frame x/y/z axis vector w.r.t. its
/// rotation angles (extrinsic Euler angles) around local x/y/z axis
RotationToAxes rotationToLocalAxesDerivative(
    const RotationMatrix3& compositeRotation,
    const RotationMatrix3& relRotation = RotationMatrix3::Identity());


class FrameToFrameJacobian { 

  // Skew-symmetric rotation derivatives 
  ActsMatrix<3,3> R_Dalpha = ActsMatrix<3,3>::Zero;
  ActsMatrix<3,3> R_Dbeta  = ActsMatrix<3,3>::Zero;
  ActsMatrix<3,3> R_Dgamma = ActsMatrix<3,3>::Zero;

  double sign = 1.;
  
  R_Dalpha(1,2) = sign;
  R_Dalpha(2,1) = -sign;
  
  R_Dbeta(0,2) = -sign;
  R_Dbeta(2,0) = sign;

  R_Dgamma(0,1) = sign;
  R_Dgamma(1,0) = -sign;
    

/// @brief Evaluate the da/dA jacobians
/// where "a" are the alignment parameters (tx,ty,tz,rx,ry,rz) of the sub-component
/// and "A" are the alignment parameters (Tx,Ty,Tz,Rx,Ry,Rz) of the composite structure
///
/// @param compositeTransform The local-to-global transformation of the composite structure
/// @param componentTransform The local-to-global transformation of the sub-component detector element

/// 
/// The Jacobian is a 6x6 matrix:
///
/// / dtx/dTx dtx/dTy dtx/dTz dtx/dRx dtx/dRy dtx/dRz \        /       |                \
/// | dty/dTx dty/dTy dty/dTz dty/dRx dty/dRy dty/dRz |        |  R^-1 |  (T x a_s) . b |
/// | dtz/dTx dtz/dTy dtz/dTz dtz/dRx dtz/dRy dtz/dRz |   =    |_______|________________|
/// | drx/dTx drx/dTy drx/dTz drx/dRx drx/dRy drx/dRz |        |       |                |
/// | dry/dTx dry/dTy dry/dTz dry/dRx dry/dRy dry/dRz |        |  0    |      R^-1      |
/// \ drz/dTx drz/dTy drz/dTz drz/dRx drz/dRy drz/dRz /        \       |                /
///
/// where |T,R| is local to composite transform
/// a_s are the vectors x_s = R * (1 0 0 ), y_s = R*(0 1 0), z_s = R*(0 0 1), b are the composite object axes

ActsMatrix<6,6> FrameJacobian(Transform3 compositeTransform,
			      Transform3 componentTransform);


// @brief Evaluate the jacobian between the translation of component wrt rotation of composite
// @
ActsMatrix<3,3> dtdA(const ActsVector3& T0,
		     const RotationMatrix3& R);



ActsMatrix<3,3> dtdA_cms(const ActsVector3 & T0,
			 const RotationMatrix3 RotDet,
			 const RotationMAtrix3 RotRot) {

  ActsVector3 dEulerA;
  ActsVector

    }
};

}  // namespace Acts::detail
