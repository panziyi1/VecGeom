#ifndef VECCORE_BACKEND_UMESIMD_H
#define VECCORE_BACKEND_UMESIMD_H

#ifdef VECCORE_ENABLE_UMESIMD

#include <umesimd/UMESimd.h>

namespace vecCore {

namespace backend {

// a UME::SIMD backend with automatically chosen SIMD sizes
class UMESimd {
public:
  using Real_v   = UME::SIMD::SIMDVec<Real_s, SIMDWidth<Real_s>()>;
  using Float_v  = UME::SIMD::SIMDVec<Float_s, SIMDWidth<Float_s>()>;
  using Double_v = UME::SIMD::SIMDVec<Double_s, SIMDWidth<Double_s>()>;

  using Int_v   = UME::SIMD::SIMDVec<Int_s, SIMDWidth<Int_s>()>;
  using Int16_v = UME::SIMD::SIMDVec<Int16_s, SIMDWidth<Int16_s>()>;
  using Int32_v = UME::SIMD::SIMDVec<Int32_s, SIMDWidth<Int32_s>()>;
  using Int64_v = UME::SIMD::SIMDVec<Int64_s, SIMDWidth<Int64_s>()>;

  using UInt_v   = UME::SIMD::SIMDVec<UInt_s, SIMDWidth<UInt_s>()>;
  using UInt16_v = UME::SIMD::SIMDVec<UInt16_s, SIMDWidth<UInt16_s>()>;
  using UInt32_v = UME::SIMD::SIMDVec<UInt32_s, SIMDWidth<UInt32_s>()>;
  using UInt64_v = UME::SIMD::SIMDVec<UInt64_s, SIMDWidth<UInt64_s>()>;
};


#define UMESIMD_REAL_FUNC(f, name)                                                                 \
  VECCORE_FORCE_INLINE UMESimd::Real_v f(const UMESimd::Real_v &x) \
  {                                                                                                \
    return x.name();                                                                               \
  }

UMESIMD_REAL_FUNC(Abs, abs)
UMESIMD_REAL_FUNC(Exp, exp)
UMESIMD_REAL_FUNC(Log, log)
UMESIMD_REAL_FUNC(Sin, sin)
UMESIMD_REAL_FUNC(Cos, cos)
UMESIMD_REAL_FUNC(Tan, tan)
UMESIMD_REAL_FUNC(Atan, atan)
UMESIMD_REAL_FUNC(Sqrt, sqrt)
UMESIMD_REAL_FUNC(Round, round)
UMESIMD_REAL_FUNC(Floor, floor)
UMESIMD_REAL_FUNC(Ceil, ceil)

#undef UMESIMD_REAL_FUNC

} // namespace backend

} // namespace vecCore

// bring in common implementation
#include "UMESimdCommon.h"

#endif
#endif
