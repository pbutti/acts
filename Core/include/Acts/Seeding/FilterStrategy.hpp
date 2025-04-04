#pragma once

#include <vector>
#include <algorithm>


#include <iostream>
#include <hwy/highway.h>
//#include "hwy/print.h"
#include <hwy/aligned_allocator.h>
#include <boost/align/aligned_allocator.hpp>

//#define _DEBUG_
//#define _DEBUG2_

namespace hn = hwy::HWY_NAMESPACE;

namespace Acts {
  
  /*
template <typename T, std::size_t Alignment>
struct AlignedAllocator {
    using value_type = T;

    T* allocate(std::size_t n) {
        if (auto ptr = static_cast<T*>(std::aligned_alloc(Alignment, n * sizeof(T)))) {
            return ptr;
        }
        throw std::bad_alloc();
    }

    void deallocate(T* ptr, std::size_t) noexcept {
        std::free(ptr);
    }
};

  */
  
  // Use this allocator with std::vector
  //using AlignedVector = std::vector<float, AlignedAllocator<float, 64>>;

  // Define an aligned vector using Highway's allocator
  
  template <typename T>
  using hwy_aligned_vector = std::vector<T, hwy::AlignedAllocator<T>>;

  // Define an alignment vector using Boost's allocator
  
  template <typename T>
  using boost_aligned_vector = std::vector<T, boost::alignment::aligned_allocator<T, 64>>;
  
  // Create SOA layout for better vectorization
  template <typename external_spacepoint_t>
  struct SpacePointsSOA {
    boost_aligned_vector<float> radius;
    boost_aligned_vector<float> x;
    boost_aligned_vector<float> y;
    boost_aligned_vector<float> z;
    boost_aligned_vector<float> varianceR;
    boost_aligned_vector<float> varianceZ;
    	
    // Constructor to preload from pointer-based structure
    explicit SpacePointsSOA(const std::vector<const external_spacepoint_t*>& points) {
      size_t size = points.size();
      radius.reserve(size);
      x.reserve(size);
      y.reserve(size);
      z.reserve(size);
      varianceR.reserve(size);
      varianceZ.reserve(size);
      for (size_t ipoint=0; ipoint<size; ipoint++) {
	radius.push_back(points[ipoint]->radius());
	x.push_back(points[ipoint]->x());
	y.push_back(points[ipoint]->y());
	z.push_back(points[ipoint]->z());
	varianceR.push_back(points[ipoint]->varianceR());
	varianceZ.push_back(points[ipoint]->varianceZ());
	
      }
    }
  };
  
  template<typename external_spacepoint_t, bool isBottomCandidate>
  bool optimizeSpacePointSearchSOA(
	  const SpacePointsSOA<external_spacepoint_t>& pointsSOA,
	  float rM, float xM, float yM, float zM, float uIP2, float cosPhiM, float sinPhiM,
	  float varianceRM, float varianceZM, float deltaRMaxSP, float deltaRMinSP,
	  float impactMax, float cotThetaMax, float deltaZMax,
	  float collisionRegionMin, float collisionRegionMax) {

   
    // Highway SIMD vector configuration for float
    const hn::ScalableTag<float> d;
    const auto lanes = hn::Lanes(d);
    
    constexpr unsigned int align = 64;
    // ref radii
    auto rMs = hn::Set(d, rM);
    auto zMs = hn::Set(d, zM);
    auto xMs = hn::Set(d, xM);
    auto yMs = hn::Set(d, yM);

    
    // Vectorized comparison for early termination
    auto minThresholdVector = hn::Set(d, deltaRMinSP);
    auto maxThresholdVector = hn::Set(d, deltaRMaxSP);
    // Prepare collision region checks
    auto collRegMin = hn::Set(d, collisionRegionMin);
    auto collRegMax = hn::Set(d, collisionRegionMax);

        
    
    auto cosPhiMs = hn::Set(d,cosPhiM);
    auto sinPhiMs = hn::Set(d,sinPhiM);
    auto impactMaxs = hn::Set(d, impactMax);
    auto cotThetaMaxs = hn::Set(d,cotThetaMax);
    //auto deltaZMaxs = hn::Set(d,deltaZMax);
    auto vZms = hn::Set(d,varianceZM);
    auto vRms = hn::Set(d,varianceRM);


    auto ones = hn::Set(d,1.f);

    
    
    // Process in SIMD chunks
    for (size_t i = 0; i < pointsSOA.radius.size(); i += lanes) {
      // Directly load from SOA vectors - no scatter-gather overhead
      auto radiusVector = hn::LoadN(d, pointsSOA.radius.data() + i, 
				    std::min(lanes, pointsSOA.radius.size() - i));
      auto zVector = hn::LoadN(d, pointsSOA.z.data() + i,
      			       std::min(lanes, pointsSOA.z.size() - i));
      auto xVector = hn::LoadN(d, pointsSOA.x.data() + i,
      			       std::min(lanes, pointsSOA.x.size() - i));
      auto yVector = hn::LoadN(d, pointsSOA.y.data() + i, 
			       std::min(lanes, pointsSOA.y.size() - i));
      auto vZs = hn::LoadN(d, pointsSOA.varianceZ.data() + i,
				std::min(lanes, pointsSOA.varianceZ.size() - i));
      auto vRs = hn::LoadN(d, pointsSOA.varianceR.data() + i,
				std::min(lanes, pointsSOA.varianceR.size() - i));


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
      
      //V MaskedMul(M m, V a, V b): returns a[i] * b[i] or 0 if m[i] is false.
      //auto collRegMinDeltaR = hn::MaskedMul(deltaRMask,collRegionMin,deltaRs);
      auto collRegMinDeltaR = hn::Mul(collRegMin,deltaRs);
      
      
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
      
      auto deltaXs = hn::Sub(xVector, xMs);
      auto deltaYs = hn::Sub(yVector, yMs);
      
      // masked aritmetics
      //auto deltaYsinPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaYs,sinPhiMs));
      //auto deltaYcosPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaYs,cosPhiMs));
      //auto deltaXsinPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaXs,sinPhiMs);)

      auto xNewFrames = hn::Add(hn::Mul(deltaXs,cosPhiMs),hn::Mul(deltaYs,sinPhiMs));
      auto yNewFrames = hn::Sub(hn::Mul(deltaYs,cosPhiMs),hn::Mul(deltaXs,sinPhiMs));
      
      auto deltaX2s = hn::Mul(deltaXs,deltaXs);
      auto deltaY2s = hn::Mul(deltaYs,deltaYs);
      auto deltaR2s = hn::Add(deltaX2s,deltaY2s);
      
      
      auto iDeltaR2s = hn::IfThenElseZero(combinedMask,hn::Div(ones,deltaR2s));
      auto uTs = hn::Mul(xNewFrames,iDeltaR2s);
      auto vTs = hn::Mul(yNewFrames,iDeltaR2s);

      
      auto impactMaxTimesxNewFrames = hn::Mul(impactMaxs,xNewFrames);
      auto rMTimesyNewFrames = hn::Abs(hn::Mul(rMs,yNewFrames));

      // I make sure I combine it with the combined Mask to avoid miscomputations in masked lines
      auto frameMask = hn::And(combinedMask,hn::Lt(rMTimesyNewFrames,impactMaxTimesxNewFrames));

      auto negFrameMask = hn::Not(frameMask);
      
      
      auto deltaRCotThetaMaxs = hn::Mul(deltaRs,cotThetaMaxs);

      
      // -cotThetaMax*deltaR< deltaZ < cotThetaMax*deltaR
      auto cotThetaMask = hn::And(hn::Gt(deltaZs,hn::Neg(deltaRCotThetaMaxs)),hn::Lt(deltaZs,deltaRCotThetaMaxs));
      
      // -deltaZMax<deltaZ<deltaZMax
      //auto deltaZMask = hn::And(hn::Gt(deltaZs,hn::Neg(deltaZMaxs)),hn::Lt(deltaZs,deltaZMaxs));

      // combine the deltaZ masks
      //auto combinedZmasks = hn::And(cotThetaMask,deltaZMask);

      // combine with the previously combined mask
      
      auto mask = hn::And(frameMask, cotThetaMask);

      // Can be optimized
      auto iDeltaRs  = hn::IfThenElseZero(mask, hn::Sqrt(iDeltaR2s));
      auto cotThetas = hn::Mul(deltaZs,iDeltaRs); 
      
      // Prepare the errors
      
      auto Ers = hn::IfThenElseZero(mask,hn::Add(hn::Add(vZms,vZs),hn::Mul(hn::Mul(hn::Mul(cotThetas,cotThetas), hn::Add(vRms,vRs)),iDeltaR2s)));
      
      // For the moment do not implement the bottom candidate cut
      
#ifdef _DEBUG2_

      std::cout<<"cotThetas"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(cotThetas,idbg) << std::endl;
      }
            
      std::cout<<"iDeltaRs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(iDeltaRs,idbg) << std::endl;
      }

      std::cout<<"Ers"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(Ers,idbg) << std::endl;
      }

      std::cout<<"uTs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(uTs,idbg) << std::endl;
      }

      std::cout<<"vTs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(vTs,idbg) << std::endl;
      }

      std::cout<<"xNewFrames"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(xNewFrames,idbg) << std::endl;
      }

      std::cout<<"yNewFrames"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(yNewFrames,idbg) << std::endl;
      }
#endif
      
      
    } // loop on other SPs
    

    return true;
  }

  

  // This returns a mask of the SPs that pass the condition.
  template<typename external_spacepoint_t, bool isBottomCandidate>
  bool optimizeSpacePointSearch(
	  const std::vector<const external_spacepoint_t*>& otherSPs,
	  float rM, float xM, float yM, float zM, float uIP2, float cosPhiM, float sinPhiM,
	  float varianceRM, float varianceZM, float deltaRMaxSP, float deltaRMinSP,
	  float impactMax, float cotThetaMax, float deltaZMax,
	  float collisionRegionMin, float collisionRegionMax) {

    
    // Initialize:
    
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

    constexpr unsigned int align = 64;
    
    // Output mask. Use uint8_t for better memory layout.
    // Alternatively consider <bool>
    //std::vector<uint8_t> valid_mask;
    //valid_mask.reserve(otherSPs.size());
    size_t paddingSize = 0;

    // Compare constants only once
    
    // ref radii
    auto rMs = hn::Set(d, rM);
    auto zMs = hn::Set(d, zM);
    auto xMs = hn::Set(d, xM);
    auto yMs = hn::Set(d, yM);
    
    // Vectorized comparison for early termination
    auto minThresholdVector = hn::Set(d, deltaRMinSP);
    auto maxThresholdVector = hn::Set(d, deltaRMaxSP);
    // Prepare collision region checks
    auto collRegMin = hn::Set(d, collisionRegionMin);
    auto collRegMax = hn::Set(d, collisionRegionMax);
    
    auto cosPhiMs = hn::Set(d,cosPhiM);
    auto sinPhiMs = hn::Set(d,sinPhiM);
    auto impactMaxs = hn::Set(d, impactMax);
    auto cotThetaMaxs = hn::Set(d,cotThetaMax);
    //auto deltaZMaxs = hn::Set(d,deltaZMax);
    auto vZms = hn::Set(d,varianceZM);
    auto vRms = hn::Set(d,varianceRM);
    auto ones = hn::Set(d,1.f);

    // Prepare buffers for vectorized processing
    alignas(align) float radiusBuffer[lanes];
    alignas(align) float zBuffer[lanes];
    alignas(align) float xBuffer[lanes];
    alignas(align) float yBuffer[lanes];
    alignas(align) float vZBuffer[lanes];
    alignas(align) float vRBuffer[lanes];
      
      
    // Prepare vectorized processing
    for (size_t i = 0; i < otherSPs.size(); i += lanes) {    
      
      // Load radii into buffer
      for (size_t j = 0; j < lanes; ++j) {
	if (i + j < otherSPs.size()) {
	  radiusBuffer[j] = otherSPs[i + j]->radius();
	  
	  // I load them all in the first pass
	  // TODO consider masked loads
	  // Vec<D> MaskedLoad(M mask, D d, const T* p): equivalent to MaskedLoadOr(Zero(d), mask, d, p), but potentially slightly more efficient.
	  //zBuffer[j] = otherSPs[i + j]->z();
	  //xBuffer[j] = otherSPs[i + j]->x();
	  //yBuffer[j] = otherSPs[i + j]->y();
	  //vZBuffer[j] = otherSPs[i + j]->varianceZ();
	  //vRBuffer[j] = otherSPs[i + j]->varianceR();
	  
	} else {
	  // Pad with sentinel values to avoid processing invalid elements
	  radiusBuffer[j] = isBottomCandidate 
	    ? std::numeric_limits<float>::max() 
	    : std::numeric_limits<float>::lowest();
	  
	  
	  //zBuffer[j] = 0.f;
	  //xBuffer[j] = 0.f;
	  //yBuffer[j] = 0.f;
	  //vZBuffer[j] = 0.f;
	  //vRBuffer[j] = 0.f;
	  
	  
	  // Last iteration
	  paddingSize = otherSPs.size() - i;
	}
      } // buffer preparation
      
      // Load the buffer into the scalable tag
      auto radiusVector = hn::Load(d, radiusBuffer);
      
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

      // Skip computations of all deltaR fail
      if (hn::AllFalse(d, deltaRMask))
	continue;
            
      
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


      // Load radii into buffer
      for (size_t j = 0; j < lanes; ++j) {
	if (i + j < otherSPs.size()) {
	  zBuffer[j] = otherSPs[i + j]->z();
	} else {
	  zBuffer[j] = 0.f;
	}
      } // buffer preparation
           
      
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
      
      //V MaskedMul(M m, V a, V b): returns a[i] * b[i] or 0 if m[i] is false.
      //auto collRegMinDeltaR = hn::MaskedMul(deltaRMask,collRegionMin,deltaRs);
      auto collRegMinDeltaR = hn::Mul(collRegMin,deltaRs);
      
      
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

      if (hn::AllFalse(d, combinedMask))
	continue;
            
      
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


      // Load radii into buffer
      for (size_t j = 0; j < lanes; ++j) {
	if (i + j < otherSPs.size()) {
	  xBuffer[j] = otherSPs[i + j]->x();
	  yBuffer[j] = otherSPs[i + j]->y();
	} else {
	  xBuffer[j] = 0.f;
	  yBuffer[j] = 0.f;
	}
      } // buffer preparation

      
      auto xs = hn::Load(d,xBuffer);
      auto ys = hn::Load(d,yBuffer);
      
      auto deltaXs = hn::Sub(xs, xMs);
      auto deltaYs = hn::Sub(ys, yMs);
      
      // masked aritmetics
      //auto deltaXcosPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaXs,cosPhiMs));
      //auto deltaYsinPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaYs,sinPhiMs));
      //auto deltaYcosPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaYs,cosPhiMs));
      //auto deltaXsinPhiMs = hn::IfThenElseZero(combinedMask, hn::Mul(deltaXs,sinPhiMs);)

      auto xNewFrames = hn::Add(hn::Mul(deltaXs,cosPhiMs),hn::Mul(deltaYs,sinPhiMs));
      auto yNewFrames = hn::Sub(hn::Mul(deltaYs,cosPhiMs),hn::Mul(deltaXs,sinPhiMs));
      
      auto deltaR2s = hn::Add(hn::Mul(deltaXs,deltaXs),hn::Mul(deltaYs,deltaYs));
            
      auto iDeltaR2s = hn::IfThenElseZero(combinedMask,hn::Div(ones,deltaR2s));
      auto uTs = hn::Mul(xNewFrames,iDeltaR2s);
      auto vTs = hn::Mul(yNewFrames,iDeltaR2s);

      
      auto impactMaxTimesxNewFrames = hn::Mul(impactMaxs,xNewFrames);
      auto rMTimesyNewFrames = hn::Abs(hn::Mul(rMs,yNewFrames));

      // I make sure I combine it with the combined Mask to avoid miscomputations in masked lines
      auto frameMask = hn::And(combinedMask,hn::Lt(rMTimesyNewFrames,impactMaxTimesxNewFrames));

      auto negFrameMask = hn::Not(frameMask);
      
      
      auto deltaRCotThetaMaxs = hn::Mul(deltaRs,cotThetaMaxs);

      
      // -cotThetaMax*deltaR< deltaZ < cotThetaMax*deltaR
      auto cotThetaMask = hn::And(hn::Gt(deltaZs,hn::Neg(deltaRCotThetaMaxs)),hn::Lt(deltaZs,deltaRCotThetaMaxs));
      
      // -deltaZMax<deltaZ<deltaZMax
      //auto deltaZMask = hn::And(hn::Gt(deltaZs,hn::Neg(deltaZMaxs)),hn::Lt(deltaZs,deltaZMaxs));

      // combine the deltaZ masks
      //auto combinedZmasks = hn::And(cotThetaMask,deltaZMask);

      // combine with the previously combined mask
      
      auto mask = hn::And(frameMask, cotThetaMask);
      
      if (hn::AllFalse(d, mask))
	continue;
      
      
      // Can be optimized
      auto iDeltaRs  = hn::IfThenElseZero(mask, hn::Sqrt(iDeltaR2s));
      auto cotThetas = hn::Mul(deltaZs,iDeltaRs); 
      
      // Prepare the errors

      
      for (size_t j = 0; j < lanes; ++j) {
	if (i + j < otherSPs.size()) {
	  vZBuffer[j] = otherSPs[i + j]->varianceR();
	  vRBuffer[j] = otherSPs[i + j]->varianceZ();
	} else {
	  vZBuffer[j] = 0.f;
	  vRBuffer[j] = 0.f;
	}
      } // buffer preparation
      
      
      auto vZs = hn::Load(d,vZBuffer);
      auto vRs = hn::Load(d,vRBuffer);
      
      auto Ers = hn::IfThenElseZero(mask,hn::Add(hn::Add(vZms,vZs),hn::Mul(hn::Mul(hn::Mul(cotThetas,cotThetas), hn::Add(vRms,vRs)),iDeltaR2s)));
      
      // For the moment do not implement the bottom candidate cut
      
#ifdef _DEBUG2_

      std::cout<<"cotThetas"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(cotThetas,idbg) << std::endl;
      }
            
      std::cout<<"iDeltaRs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(iDeltaRs,idbg) << std::endl;
      }

      std::cout<<"Ers"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(Ers,idbg) << std::endl;
      }

      std::cout<<"uTs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(uTs,idbg) << std::endl;
      }

      std::cout<<"vTs"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(vTs,idbg) << std::endl;
      }

      std::cout<<"xNewFrames"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(xNewFrames,idbg) << std::endl;
      }

      std::cout<<"yNewFrames"<<std::endl;
      for (size_t idbg = 0; idbg < lanes; idbg++) {
	std::cout << "Lane " << idbg << ": " << hn::ExtractLane(yNewFrames,idbg) << std::endl;
      }
            
#endif
            
  } //loop on other SPs


    return true;
  } // optimizeSpacePointSearch
}




