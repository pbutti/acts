#pragma once

#include <vector>
#include <algorithm>


#include <iostream>
#include <hwy/highway.h>
//#include "hwy/print.h"

#define _DEBUG_

namespace hn = hwy::HWY_NAMESPACE;

namespace Acts {



  template<typename external_spacepoint_t, bool isBottomCandidate> 
  std::vector<uint8_t> optimizeSpacePointSearch(
			const std::vector<const external_spacepoint_t*>& otherSPs,
			const external_spacepoint_t& mediumSP,
			const SeedFinderConfig<external_spacepoint_t>& config
						) {

    
    float rM = mediumSP.radius();
    float xM = mediumSP.x();
    float yM = mediumSP.y();
    float zM = mediumSP.z();
    float uIP = -1. / rM;
    float cosPhiM = -xM / uIP;
    float sinPhiM = -yM / uIP;
    float varianceRM = mediumSP.varianceR();
    float varianceZM = mediumSP.varianceZ();

    float deltaRMaxSP = config.deltaRMaxTopSP;
    if constexpr (isBottomCandidate) {
      deltaRMaxSP = config.deltaRMaxBottomSP;
    }
    
    float deltaRMinSP = config.deltaRMinTopSP;
    if constexpr (isBottomCandidate) {
      deltaRMinSP = config.deltaRMinBottomSP;
    }
    
    
    float impactMax = config.impactMax;
    if constexpr (isBottomCandidate) {
      impactMax = -impactMax;
    }
    
    float cotThetaMax = config.cotThetaMax;
    float deltaZMax   = config.deltaZMax;
    float collisionRegionMin = config.collisionRegionMin;
    float collisionRegionMax = config.collisionRegionMax;

  #ifdef _DEBUG_
  // Debug: Print input parameters
    std::cout << "Debug: optimizeSpacePointSearch " 
	      << (isBottomCandidate ? "[Bottom Candidate]" : "[Top Candidate]") << std::endl;
    std::cout << "  Reference Radius (rM): " << rM << std::endl;
    std::cout << "  Delta R Min: " << deltaRMinSP << std::endl;
    std::cout << "  Delta R Max: " << deltaRMaxSP << std::endl;
    std::cout << "  Total Space Points: " << otherSPs.size() << std::endl;
#endif
    
    // Highway SIMD vector configuration for float
    const hn::ScalableTag<float> d;
        
    const auto lanes = hn::Lanes(d);
    
    // Output mask. Use uint8_t for better memory layout.
    // Alternatively consider <bool>
    std::vector<uint8_t> valid_mask;
    valid_mask.reserve(otherSPs.size());
    size_t paddingSize = 0;
    
    // Prepare vectorized processing
    for (size_t i = 0; i < otherSPs.size(); i += lanes) {    
      // Prepare buffers for vectorized processing
      alignas(64) float radiusBuffer[lanes];
      alignas(64) float zBuffer[lanes];
      alignas(64) float xBuffer[lanes];
      alignas(64) float yBuffer[lanes];
      
      // Load radii into buffer
      for (size_t j = 0; j < lanes; ++j) {
	if (i + j < otherSPs.size()) {
	  radiusBuffer[j] = otherSPs[i + j]->radius();
	  
	  // I load them all in the first pass
	  // TODO consider masked loads
	  // Vec<D> MaskedLoad(M mask, D d, const T* p): equivalent to MaskedLoadOr(Zero(d), mask, d, p), but potentially slightly more efficient.
	  zBuffer[j] = otherSPs[i + j]->z();
	  xBuffer[j] = otherSPs[i + j]->x();
	  yBuffer[j] = otherSPs[i + j]->y();
	  
	} else {
	  // Pad with sentinel values to avoid processing invalid elements
	  radiusBuffer[j] = isBottomCandidate 
	    ? std::numeric_limits<float>::max() 
	    : std::numeric_limits<float>::lowest();
	  
	  
	  zBuffer[j] = 0.f;
	  xBuffer[j] = 0.f;
	  yBuffer[j] = 0.f;
	  
	  
	  
	  // Last iteration
	  paddingSize = otherSPs.size() - i;
	}
      } // buffer preparation
      
      // Load the buffer into the scalable tag
      auto radiusVector = hn::Load(d, radiusBuffer);
      
      // Set the reference vector
      auto rMs = hn::Set(d, rM);
      
#ifdef _DEBUG_
      std::cout << "  Reference Radius Vector: " << std::endl;
      
      auto dbg_buffer = std::make_unique<int[]>(lanes);
      
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(rMs,idbg) << std::endl;
      }
      
      std::cout << "  Radius Vector: " << std::endl;
      
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(radiusVector,idbg) << std::endl;
      }
      
      
#endif
      
      // Compute delta radii vectorized
      auto deltaRs = isBottomCandidate 
	? hn::Sub(rMs, radiusVector)   // rM - otherSP radius
	: hn::Sub(radiusVector, rMs);  // otherSP radius - rM
      
#ifdef _DEBUG_
      std::cout << "  Delta Radii Vector: " << std::endl;
      for (size_t idbg = 0; idbg < lanes; ++idbg) {
	std::cout << "    Lane " << idbg << ": " << hn::ExtractLane(deltaRs, idbg) << std::endl;
      }
      
#endif
      
      // For bottom candidates I need to check that deltaRMinSP < deltaR < deltaRMaxSP
      // For top candidates I need to check    that deltaRMinSP < deltaR < deltaRMaxSP  
      
      // Vectorized comparison for early termination
      auto minThresholdVector = hn::Set(d, deltaRMinSP);
      auto maxThresholdVector = hn::Set(d, deltaRMaxSP);
      
#ifdef _DEBUG_
      std::cout << "  Threshold Vectors (min/max):" << std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(minThresholdVector,idbg) << std::endl;
      }
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(maxThresholdVector,idbg) << std::endl;
      }
      
#endif
      
      
      //PF:: TODO TEST:M MaskedGt(M m, V a, V b): returns a[i] > b[i] or false if m[i] is false.
      //Which might be more efficient after calling Lt with the first mask
      auto ltMask = hn::Lt(deltaRs, maxThresholdVector); //true if deltaR < deltaRMaxSP
      auto gtMask = hn::Gt(deltaRs, minThresholdVector); //true if deltaR > deltaRMinSP
      auto deltaRMask = hn::And(ltMask,gtMask);
      
      //Alternative:
      //auto deltaRMask = hn::MaskedGt(ltMask, deltaRs, minThreasholdVector);
      
#ifdef _DEBUG_
      
      std::cout << "  Termination Mask:" << std::endl;


      // up to 8 lanes
      uint8_t mask_dbg = 0;
      hn::StoreMaskBits(d, deltaRMask, &mask_dbg);
            
      std::cout << "Mask values:" << std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	bool active = (mask_dbg >> idbg) & 1;
	std::cout << "Lane " << idbg << ": " << active << std::endl;
      }
      std::cout<<"Now storing in the output mask"<<std::endl;
      
#endif
      
      // calc of deltaZ
      //pit zm in a vector
      auto zMs     = hn::Set(d, zM);
      
      //TODO:::Here I could use a masked load?
      auto zVector = hn::Load(d, zBuffer);
      
      // compute deltaZ
      auto deltaZs = isBottomCandidate
	? hn::Sub(zMs, zVector)
	: hn::Sub(zVector, zMs);
      
      // compute zOriginTimesDeltaR
      // Will use masked aritmetic, where the mask is the deltaR mask from above
      // V MaskedNegMulAdd(M m, V a, V b, V c): returns -a[i] * b[i] + c[i] or 0 if m[i] is false.
      // V MaskedMul(M m, V a, V b): returns a[i] * b[i] or 0 if m[i] is false.
      
      //zM*deltaR - rM*deltaZ = 
      //-rM*deltaZ + zM*deltaR
      
      // Not available in the tag, only in master
      //auto zOriginTimesDeltaRs = hn::MaskedNegMulAdd(deltaRMask, rMs, deltaZs,
      //						   hn::MaskedMul(deltaRMask,zMs,deltaRs));
      //HWY_API V MaskedNegMulAdd(M m, V mul, V x, V add) {
      //return IfThenElseZero(m, NegMulAdd(mul, x, add));
      //}
      //template <class V, class M>
      //HWY_API V MaskedMul(M m, V a, V b) {
      //  return IfThenElseZero(m, Mul(a, b));
      //}
      

      // Easy computations: just do them.
      auto zMdeltaRs = hn::Mul(zMs,deltaRs);
      
      auto zOriginTimesDeltaRs = hn::NegMulAdd(rMs, deltaZs, zMdeltaRs);
      
      // Prepare collision region checks
      auto collRegMin = hn::Set(d, collisionRegionMin);
      
      //V MaskedMul(M m, V a, V b): returns a[i] * b[i] or 0 if m[i] is false.
      //auto collRegMinDeltaR = hn::MaskedMul(deltaRMask,collRegionMin,deltaRs);
      auto collRegMinDeltaR = hn::Mul(collRegMin,deltaRs);
      
      auto collRegMax = hn::Set(d, collisionRegionMax);
      //auto collRegMaxDeltaR = hn::MaskedMul(deltaRMask,collRegionMin,deltaRs);
      auto collRegMaxDeltaR = hn::Mul(collRegMax,deltaRs);
      
      
      // TODO:: maybe combine the masks differently?
      // M MaskedLt(M m, V a, V b): returns a[i] < b[i] or false if m[i] is false.
      
      //MaskedGt return And(m, Gt(a, b));
      //true if zOriginTimesDeltaRs > collRegMinDeltaR
      auto collMaskMin = hn::And(deltaRMask,hn::Gt(zOriginTimesDeltaRs, collRegMinDeltaR));
      
      // collMaskMax combines, deltaRmask, collMaskMin and collMaskMax
      //HWY_API auto MaskedLt(M m, V a, V b) -> decltype(a == b) {
      //return And(m, Lt(a, b));
      
      //true if zOriginTimesDeltaRs < collRegMaxDeltaR
      auto combinedMask = hn::And(collMaskMin,hn::Lt(zOriginTimesDeltaRs, collRegMaxDeltaR)); // true if 
      
      
      
#ifdef _DEBUG_
      
      std::cout<<"zMs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(zMs,idbg) << std::endl;
      }
      
      std::cout<<"zVector"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(zVector,idbg) << std::endl;
      }
      
      std::cout<<" deltaZs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(deltaZs,idbg) << std::endl;
      }
      
      std::cout<<" zMdeltaRs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(zMdeltaRs,idbg) << std::endl;
      }
      
      //std::cout<<" zMdeltaRs"<<std::endl;
      //for (size_t idbg = 0; idbg < lanes; idbg++) {
      //  std::cout << "Lane " << idbg << ": " << hn::ExtractLane(zMdeltaRs,idbg) << std::endl;
      //}
      
      std::cout<<" zOriginTimesDeltaRs:"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(zOriginTimesDeltaRs,idbg) << std::endl;
      }
      
      std::cout<<"collRegMinDeltaR"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(collRegMinDeltaR,idbg) << std::endl;
      }
      
      std::cout<<"collRegMaxDeltaR"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(collRegMaxDeltaR,idbg) << std::endl;
      }
      
      std::cout<<"The combined mask"<<std::endl;
      
      uint8_t combinedMask_dbg = 0;
      hn::StoreMaskBits(d, combinedMask, &combinedMask_dbg);
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	bool active = (combinedMask_dbg >> idbg) & 1;
      std::cout << "Lane " << idbg << ": " << active << std::endl;
      }
      
#endif
      
      // The interactionPointCut check not implemented for the moment
      
      auto xMs = hn::Set(d, xM);
      auto yMs = hn::Set(d, yM);
      auto cosPhiMs = hn::Set(d,cosPhiM);
      auto sinPhiMs = hn::Set(d,sinPhiM);

      auto xs = hn::Load(d,xBuffer);
      auto ys = hn::Load(d,yBuffer);
      
      auto deltaXs = hn::Sub(xs, xMs);
      auto deltaYs = hn::Sub(ys, yMs);
      
      // masked aritmetics
      //auto deltaXcosPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaXs,cosPhiMs));
      //auto deltaYsinPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaYs,sinPhiMs));
      //auto deltaYcosPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaYs,cosPhiMs));
      //auto deltaXsinPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaXs,sinPhiMs);)

      auto deltaXcosPhiMs = hn::Mul(deltaXs,cosPhiMs);
      auto deltaYsinPhiMs = hn::Mul(deltaYs,sinPhiMs);
      auto deltaYcosPhiMs = hn::Mul(deltaYs,cosPhiMs);
      auto deltaXsinPhiMs = hn::Mul(deltaXs,sinPhiMs);
            
      auto xNewFrames = hn::Add(deltaXcosPhiMs,deltaYsinPhiMs);
      auto yNewFrames = hn::Sub(deltaYcosPhiMs,deltaXsinPhiMs);
      
      auto deltaX2s = hn::Mul(deltaXs,deltaXs);
      auto deltaY2s = hn::Mul(deltaYs,deltaYs);
      auto deltaR2s = hn::Add(deltaX2s,deltaY2s);
      auto ones = hn::Set(d,1.);
      
      auto iDeltaR2s = hn::IfThenElseZero(combinedMask,hn::Div(ones,deltaR2s));
      auto uTs = hn::Mul(xNewFrames,iDeltaR2s);
      auto vTs = hn::Mul(yNewFrames,iDeltaR2s);

      auto impactMaxs = hn::Set(d, impactMax);
      auto impactMaxTimesxNewFrames = hn::IfThenElseZero(combinedMask,hn::Mul(impactMaxs,xNewFrames));
      auto rMTimesyNewFrames = hn::Abs(hn::IfThenElseZero(combinedMask, hn::Mul(rMs,yNewFrames)));

      // I make sure I combine it with the combined Mask to avoid miscomputations in masked lines
      auto frameMask = hn::And(combinedMask,hn::Lt(rMTimesyNewFrames,impactMaxTimesxNewFrames));
      
      
      auto cotThetaMaxs = hn::Set(d,cotThetaMax);
      auto deltaRCotThetaMaxs = hn::IfThenElseZero(frameMask,hn::Mul(deltaRs,cotThetaMaxs));


      // -cotThetaMax*deltaR< deltaZ < cotThetaMax*deltaR
      auto cotThetaMask = hn::And(hn::Gt(deltaZs,hn::Neg(deltaRCotThetaMaxs)),hn::Lt(deltaZs,deltaRCotThetaMaxs));
      
      //auto deltaZMaxs = hn::Set(d,deltaZMax);

      // -deltaZMax<deltaZ<deltaZMax
      //auto deltaZMask = hn::And(hn::Gt(deltaZs,hn::Neg(deltaZMaxs)),hn::Lt(deltaZs,deltaZMaxs));

      // combine the deltaZ masks
      //auto combinedZmasks = hn::And(cotThetaMask,deltaZMask);

      // combine with the previously combined mask
      
      auto mask = hn::And(frameMask, cotThetaMask);

      // Can be optimized
      auto iDeltaRs  = hn::IfThenElseZero(mask, hn::Sqrt(iDeltaR2s));
      auto cotThetas = hn::IfThenElseZero(mask, hn::Mul(deltaZs,iDeltaRs)); 

      // For the moment do not implement the bottom candidate cut
      
#ifdef _DEBUG_
      
      std::cout<<"iDeltaRs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(iDeltaRs,idbg) << std::endl;
      }

      std::cout<<"cotThetas"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(cotThetas,idbg) << std::endl;
      }
      
      
#endif

      
      
  } //loop on other SPs


  return valid_mask;
  } // optimizeSpacePointSearch
}




