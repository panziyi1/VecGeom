#ifndef VECCORE_BACKEND_BASIC_H
#define VECCORE_BACKEND_BASIC_H

#include <cassert>
#include <cstdint>
#include <type_traits>

#include "ScalarWrapper.h"

namespace vecCore {

template <typename T>
struct ScalarType<ScalarWrapper<T>> { using Type = T; };

namespace backend {

class Basic {
public:
  using Real_v   = ScalarWrapper<Real_t>;
  using Float_v  = ScalarWrapper<Float_t>;
  using Double_v = ScalarWrapper<Double_t>;

  using Int_v    = ScalarWrapper<Int_t>;
  using Int16_v  = ScalarWrapper<Int16_t>;
  using Int32_v  = ScalarWrapper<Int32_t>;
  using Int64_v  = ScalarWrapper<Int64_t>;

  using UInt_v   = ScalarWrapper<UInt_t>;
  using UInt16_v = ScalarWrapper<UInt16_t>;
  using UInt32_v = ScalarWrapper<UInt32_t>;
  using UInt64_v = ScalarWrapper<UInt64_t>;

  template <typename T> using Mask_v  = BoolWrapper;
  template <typename T> using Index_v = ScalarWrapper<Int_t>;

  template <typename T>
  static constexpr UInt64_t VectorSize(const T& x)
  {
    return 1UL;
  }
};

} // namespace backend
} // namespace vecCore

#endif
