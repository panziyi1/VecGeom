#ifndef VECCORE_BACKEND_SIMD_H
#define VECCORE_BACKEND_SIMD_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>

#define __noinline__     __attribute__((noinline))
#define __force_inline__ inline __attribute__((always_inline))

namespace simd {
  template<typename T, size_t N>
    class VECCORE_ALIGNED array {
    public:
      typedef T                 value_type;
      typedef value_type*       pointer;
      typedef value_type*       iterator;
      typedef value_type&       reference;
      typedef const value_type* const_pointer;
      typedef const value_type* const_iterator;
      typedef const value_type& const_reference;
      typedef std::size_t       size_type;
      typedef std::ptrdiff_t    difference_type;

      value_type fData[N];

      constexpr array() { }

      template<typename _T>
      array(const _T& s)
      { fill(s); }

      template<typename _T>
      array(const _T* __restrict__ ptr)
      { for (size_t i = 0; i < N; i++) (*this)[i] = static_cast<T>(ptr[i]); }

      template<typename _T>
      array(const array<_T, N>& a)
      { for (size_t i = 0; i < N; i++) (*this)[i] = T(a[i]); }

      void
      load(const value_type* __restrict__ ptr) noexcept
      { for (size_t i = 0; i < N; i++) (*this)[i] = ptr[i]; }

      void
      store(value_type* __restrict__ ptr) const noexcept
      { for (size_t i = 0; i < N; i++) ptr[i] = (*this)[i]; }

      void
      fill(const value_type& __u)
      { std::fill_n(begin(), size(), __u); }

      void
      swap(array& __other)
      noexcept(noexcept(swap(std::declval<T&>(), std::declval<T&>())))
      { std::swap_ranges(begin(), end(), __other.begin()); }

      iterator
      begin() noexcept
      { return iterator(data()); }

      const_iterator
      begin() const noexcept
      { return const_iterator(data()); }

      iterator
      end() noexcept
      { return iterator(data() + N); }

      const_iterator
      end() const noexcept
      { return const_iterator(data() + N); }

      const_iterator
      cbegin() const noexcept
      { return const_iterator(data()); }

      const_iterator
      cend() const noexcept
      { return const_iterator(data() + N); }

      constexpr size_type
      size() const noexcept { return N; }

      constexpr size_type
      max_size() const noexcept { return N; }

      constexpr bool
      empty() const noexcept { return size() == 0; }

      reference
      operator[](size_type __n) noexcept
      { return const_cast<T&>(fData[__n]); }

      constexpr const_reference
      operator[](size_type __n) const noexcept
      { return fData[__n]; }

      reference
      at(size_type __n)
      {
        static_assert(__n >= N, "index out of range");
        return const_cast<T&>(fData[__n]);
      }

      constexpr const_reference
      at(size_type __n) const
      {
        static_assert(__n >= N, "index out of range");
        return fData[__n];
      }

      reference
      front() noexcept
      { return *begin(); }

      constexpr const_reference
      front() const noexcept
      { return *cbegin(); }

      reference
      back() noexcept
      { return *(end() - 1); }

      constexpr const_reference
      back() const noexcept
      { return *(cend() - 1); }

      pointer
      data() noexcept
      { return const_cast<T*>(fData); }

      const_pointer
      data() const noexcept
      { return fData; }

      bool
      all_false() const noexcept
      {
        const_iterator it = begin();
        for (; it != end(); ++it)
          if (*it) break;
        return it == end();
      }

      bool
      all_true() const noexcept
      {
        const_iterator it = begin();
        for (; it != end(); ++it)
          if (!(*it)) break;
        return it == end();
      }

      inline array&
      assign(const array& mask, const array& val)
      {
        for (size_t i = 0; i < N; ++i)
          if(mask[i]) (*this)[i] = val[i];
        return *this;
      }

      static __force_inline__
      array blend(const array& mask, const array& a, const array& b)
      {
        array result;
        for (size_t i = 0; i < N; ++i)
          result[i] = mask[i] ? a[i] : b[i];
        return result;
      }
    };

  #define UNARY_OPERATOR(OP)                                  \
  template<typename T, size_t N> __force_inline__             \
  array<T, N> operator OP(const array<T, N>& a)               \
  {                                                           \
    array<T, N> result;                                       \
    for (size_t i = 0; i < N; i++)                            \
      result[i] = OP a[i];                                    \
    return result;                                            \
  }

  #define BINARY_OPERATOR(OP)                                 \
  template<typename T, size_t N> __force_inline__             \
  array<T, N> operator OP(const array<T, N>& a,               \
                          const array<T, N>& b)               \
  {                                                           \
    array<T, N> result;                                       \
    for (size_t i = 0; i < N; i++)                            \
      result[i] = a[i] OP b[i];                               \
    return result;                                            \
  }                                                           \
                                                              \
  template<typename T, size_t N> __force_inline__             \
  array<T, N> operator OP(const T& a, const array<T, N>& b)   \
  { return array<T, N>(array<T, N>(a) OP b); }                \
                                                              \
  template<typename T, size_t N> __force_inline__             \
  array<T, N> operator OP(const array<T, N>& a, const T& b)   \
  { return array<T, N>(a OP array<T, N>(b)); }                \
                                                              \
  template<typename T1, typename T2, typename R = T2, size_t N>    \
  __force_inline__                                            \
  array<R, N> operator OP(const array<T1, N>& a,              \
                          const array<T2, N>& b)              \
  {                                                           \
    array<R, N> result;                                       \
    for (size_t i = 0; i < N; i++)                            \
      result[i] = a[i] OP b[i];                               \
    return result;                                            \
  }                                                           \
                                                              \
  template<typename T1, typename T2, typename R = T2, size_t N>    \
  __force_inline__                                            \
  array<R, N> operator OP(const T1& a, const array<T2, N>& b) \
  { return array<R, N>(array<T2, N>(a) OP b); }               \
                                                              \
  template<typename T1, typename T2, typename R = T2, size_t N>    \
  __force_inline__                                            \
  array<R, N> operator OP(const array<T1, N>& a, const T2& b) \
  { return array<R, N>(a OP array<T1, N>(b)); }               \

  // arithmetic
  UNARY_OPERATOR(+)
  UNARY_OPERATOR(-)
  BINARY_OPERATOR(+)
  BINARY_OPERATOR(-)
  BINARY_OPERATOR(*)
  BINARY_OPERATOR(/)
  BINARY_OPERATOR(%)

  // bitwise
  UNARY_OPERATOR(~)
  BINARY_OPERATOR(&)
  BINARY_OPERATOR(|)
  BINARY_OPERATOR(^)
  BINARY_OPERATOR(>>)
  BINARY_OPERATOR(<<)

  // comparison
  BINARY_OPERATOR(<)
  BINARY_OPERATOR(>)
  BINARY_OPERATOR(==)
  BINARY_OPERATOR(!=)
  BINARY_OPERATOR(<=)
  BINARY_OPERATOR(>=)

  // logical
  UNARY_OPERATOR(!)
  BINARY_OPERATOR(&&)
  BINARY_OPERATOR(||)

  #undef UNARY_OPERATOR
  #undef BINARY_OPERATOR
} // namespace simd

namespace vecCore {

template<typename T, size_t N>
struct TypeTraits<simd::array<T, N>> {
  using ScalarType = T;
  using MaskType   = simd::array<uint32_t, N>;
  using IndexType  = simd::array<uint32_t, N>;
};

namespace backend {
  template<size_t N>
  class SIMD {
  public:
    using Real_v   = simd::array<Real_s, N>;
    using Float_v  = simd::array<float, N>;
    using Double_v = simd::array<double, N>;

    using Int_v    = simd::array<int, N>;
    using Int16_v  = simd::array<int16_t, N>;
    using Int32_v  = simd::array<int32_t, N>;
    using Int64_v  = simd::array<int64_t, N>;

    using UInt_v   = simd::array<unsigned int, N>;
    using UInt16_v = simd::array<uint16_t, N>;
    using UInt32_v = simd::array<uint32_t, N>;
    using UInt64_v = simd::array<uint64_t, N>;
  };
} // namespace backend

template<typename T, size_t N>
VECCORE_FORCE_INLINE
VECCORE_CUDA_HOST_DEVICE
bool MaskEmpty(const simd::array<T, N>& mask)
{ return mask.all_false(); }

template<typename T, size_t N>
VECCORE_FORCE_INLINE
VECCORE_CUDA_HOST_DEVICE
bool MaskFull(const simd::array<T, N>& mask)
{ return mask.all_true(); }

template <typename T, size_t N>
VECCORE_FORCE_INLINE
VECCORE_CUDA_HOST_DEVICE
simd::array<T, N> Blend(const simd::array<uint32_t, N>& mask,
                        const simd::array<T, N>& tval,
                        const simd::array<T, N>& fval)
{ return simd::array<T, N>::blend(mask, tval, fval); }

template <typename T, size_t N>
VECCORE_FORCE_INLINE
VECCORE_CUDA_HOST_DEVICE
void MaskedAssign(simd::array<T, N>& dest,
                  const simd::array<uint32_t, N>& mask,
                  const simd::array<T, N>& src)
{ dest.assign(mask, src); }

namespace math {

#define SIMD_FUNCTION(function, f)      \
  template<typename T, size_t N>        \
  simd::array<T, N>                     \
  function(const simd::array<T, N>& x)  \
  {                                     \
    simd::array<T, N> result;           \
    for (size_t i = 0; i < N; ++i)      \
      result[i] = f(x[i]);              \
    return result;                      \
  }

#define SIMD_FUNCTION_2(function, f)    \
  template<typename T, size_t N>        \
  simd::array<T, N>                     \
  function (const simd::array<T, N>& x, \
            const simd::array<T, N>& y) \
  {                                     \
    simd::array<T, N> result;           \
    for (size_t i = 0; i < N; i++)      \
      result[i] = f(x[i], y[i]);        \
    return result;                      \
  }

SIMD_FUNCTION(Abs, std::abs)
SIMD_FUNCTION(Floor, std::floor)
SIMD_FUNCTION(Ceil, std::ceil)

SIMD_FUNCTION(Sin, std::sin)
SIMD_FUNCTION(Asin, std::asin)

SIMD_FUNCTION(Cos, std::cos)
SIMD_FUNCTION(Acos, std::acos)

SIMD_FUNCTION(Tan, std::tan)
SIMD_FUNCTION(Atan, std::atan)
SIMD_FUNCTION_2(Atan2, std::atan2)

SIMD_FUNCTION(Sinh, std::sinh)
SIMD_FUNCTION(Asinh, std::asinh)

SIMD_FUNCTION(Cosh, std::cosh)
SIMD_FUNCTION(Acosh, std::acosh)

SIMD_FUNCTION(Exp, std::exp)
SIMD_FUNCTION(Exp2, std::exp2)
SIMD_FUNCTION(Expm1, std::expm1)

SIMD_FUNCTION(Log, std::log)
SIMD_FUNCTION(Log2, std::log2)
SIMD_FUNCTION(Log10, std::log10)
SIMD_FUNCTION(Sqrt, std::sqrt)

SIMD_FUNCTION_2(Pow, std::pow)

template<typename T, size_t N>
void SinCos(const simd::array<T, N>& x,
            simd::array<T, N>* sin,
            simd::array<T, N>* cos)
{
  for (size_t i = 0; i < N; i++)
    sincos(x[i], &((*sin)[i]), &((*cos)[i]));
}

}

} // namespace vecCore

#endif
