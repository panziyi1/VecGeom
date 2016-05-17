#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>

#include <malloc.h>

#include <VecCore/VecCore>
#include "VecCore/Timer.h"

using namespace vecCore;

#ifdef VECCORE_DISABLE_SIMD
using Backend = backend::Scalar;
#elif defined(VECCORE_ENABLE_VC)
using Backend = backend::VcVector;
#elif defined(VECCORE_ENABLE_UMESIMD)
using Backend = backend::UMESimd;
#else
using Backend = backend::ScalarWrapper;
#endif

static const Int_s N = (8 * 1024 * 1024);

// solve ax2 + bx + c = 0, return number of roots found

template <typename T>
void QuadSolve(const T& a, const T& b, const T& c, T &x1, T &x2, int& roots) {
  T a_inv = T(1.0) / a;
  T delta = b * b - T(4.0) * a * c;
  T s = (b >= 0) ? T(1.0) : T(-1.0);

  roots = delta > NumericLimits<T>::Epsilon() ? 2 : delta < T(0.0) ? 0 : 1;

  switch (roots) {
    case 2:
      x1 = T(0.5) * (-b + s * math::Sqrt(delta));
      x2 = c / x1;
      x1 *= a_inv;
      return;

    case 0: return;

    case 1:
      x1 = x2 = T(-0.5) * b * a_inv;
      return;

    default: return;
  }
}

// SIMD version

// solve ax2 + bx + c = 0, return number of roots found
// solves many equation simultaneously, depending on value of Float_v::Size

template <class Backend>
void QuadSolveSIMD(typename Backend::Float_v const &a,
                   typename Backend::Float_v const &b,
                   typename Backend::Float_v const &c,
                   typename Backend::Float_v &x1,
                   typename Backend::Float_v &x2,
                   typename Backend::Int32_v &roots)
{
  using Float_v = typename Backend::Float_v;
  using Int32_v = typename Backend::Int32_v;
  using FMask = Mask_v<Float_v>;
  using IMask = Mask_v<Int32_v>;

  roots = 0;
  Float_v delta = b * b - Float_v(4.0) * a * c;

  FMask mask0(delta < Float_v(0.0));
  FMask mask2(delta >= NumericLimits<Float_v>::Epsilon());

  MaskedAssign(roots, (IMask)mask2, Int32_v(2));

  FMask mask = (b >= Float_v(0.0));

  Float_v sign = Blend(mask, Float_v(-1.0), Float_v(1.0));
  Float_v root = (-b + sign * math::Sqrt(delta)) / (Float_v(2.0) * a);

  MaskedAssign(x1, mask2 && mask,  root);
  MaskedAssign(x2, mask2 && !mask, root);

  root = c / (a * Blend(mask, x1, x2));

  MaskedAssign(x1, mask2 && !mask, root);
  MaskedAssign(x2, mask2 && mask,  root);

  FMask mask1 = !(mask0 || mask2);

  if (MaskEmpty(mask1))
    return;

  MaskedAssign(roots, (IMask)mask1, Int32_v(1));

  root = -Float_v(0.5) * b / a;

  MaskedAssign(x1, mask1, root);
  MaskedAssign(x2, mask1, root);
}

int main(int argc, char *argv[]) {
  Float_s *a = (Float_s *)memalign(64, N * sizeof(Float_s));
  Float_s *b = (Float_s *)memalign(64, N * sizeof(Float_s));
  Float_s *c = (Float_s *)memalign(64, N * sizeof(Float_s));

  Int_s *roots = (Int_s *)memalign(64, N * sizeof(Int_s));
  Float_s *x1 = (Float_s *)memalign(64, N * sizeof(Float_s));
  Float_s *x2 = (Float_s *)memalign(64, N * sizeof(Float_s));

  srand48(time(NULL));
  Int_s index = (Int_s)((N - 100) * drand48());

  for (Int_s i = 0; i < N; i++) {
    a[i] = 10.0 * (drand48() - 0.5);
    b[i] = 10.0 * (drand48() - 0.5);
    c[i] = 50.0 * (drand48() - 0.5);
    x1[i] = 0.0;
    x2[i] = 0.0;
    roots[i] = 0;
  }

  Timer<milliseconds> timer;

  for (Int_s i = 0; i < N; i++) {
    QuadSolve<float>(a[i], b[i], c[i], x1[i], x2[i], roots[i]);
  }

  Double_s t = timer.Elapsed();

  // print random result to ensure scalar and vector backends give same result

  for (Int_s i = index; i < index + 10; i++) {
    printf("%d: a = % 8.3f, b = % 8.3f, c = % 8.3f,"
           " roots = %d, x1 = % 8.3f, x2 = % 8.3f\n",
            i, a[i], b[i], c[i], roots[i],
            roots[i] > 0 ? x1[i] : 0,
            roots[i] > 1 ? x2[i] : 0);
  }

  printf("\nelapsed time = %.3lfms (scalar code)\n", t);

  timer.Start();

  for (Int_s i = 0; i < N; i += Int_s(VectorSize<Backend::Float_v>())) {
    QuadSolveSIMD<Backend>(
      (Backend::Float_v &)(a[i]),
      (Backend::Float_v &)(b[i]),
      (Backend::Float_v &)(c[i]),
      (Backend::Float_v &)(x1[i]),
      (Backend::Float_v &)(x2[i]),
      (Backend::Int32_v &)(roots[i])
    );
  }

  t = timer.Elapsed();

  for (Int_s i = index; i < index + 10; i++) {
    printf("%d: a = % 8.3f, b = % 8.3f, c = % 8.3f,"
           " roots = %d, x1 = % 8.3f, x2 = % 8.3f\n",
            i, a[i], b[i], c[i], roots[i],
            roots[i] > 0 ? x1[i] : 0,
            roots[i] > 1 ? x2[i] : 0);
  }

  printf("\nelapsed time = %.3lfms (vector backend)\n", t);

  return 0;
}
