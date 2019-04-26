/// @file TorusImplementation.h

#ifndef VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_

#include "base/Vector3D.h"
#include "volumes/TorusStruct.h"
#include "volumes/kernel/GenericKernels.h"
#include "VecCore/VecCore"

#include "volumes/kernel/TubeImplementation.h"

#include <cstdio>
#include <VecCore/VecCore>

namespace vecgeom {

// VECGEOM_DEVICE_FORWARD_DECLARE(struct TorusImplementation;);
// VECGEOM_DEVICE_DECLARE_CONV(struct, TorusImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

//_____________________________________________________________________________
template <typename T>
VECCORE_ATT_HOST_DEVICE
unsigned int SolveCubic(T a, T b, T c, T *x)
{
  // Find real solutions of the cubic equation : x^3 + a*x^2 + b*x + c = 0
  // Input: a,b,c
  // Output: x[3] real solutions
  // Returns number of real solutions (1 or 3)
  const T ott        = 1. / 3.;
  const T sq3        = Sqrt(3.);
  const T inv6sq3    = 1. / (6. * sq3);
  unsigned int ireal = 1;
  T p                = b - a * a * ott;
  T q                = c - a * b * ott + 2. * a * a * a * ott * ott * ott;
  T delta            = 4 * p * p * p + 27. * q * q;
  T t, u;

  if (vecCore::MaskFull(delta >= 0)) {
    delta = Sqrt(delta);
    t     = (-3 * q * sq3 + delta) * inv6sq3;
    u     = (3 * q * sq3 + delta) * inv6sq3;
    x[0]  = CopySign(T(1.), t) * Cbrt(Abs(t)) - CopySign(T(1.), u) * Cbrt(Abs(u)) - a * ott;
  } else {
    delta = Sqrt(-delta);
    t     = -0.5 * q;
    u     = delta * inv6sq3;
    x[0]  = 2. * Pow(t * t + u * u, 0.5 * ott) * Cos(ott * ATan2(u, t)) - a * ott;
  }

  t     = x[0] * x[0] + a * x[0] + b;
  u     = a + x[0];
  delta = u * u - 4. * t;
  if (!vecCore::MaskEmpty(delta >= 0)) {
    ireal = 3;
    delta = Sqrt(Max(delta, T(0.0)));
    x[1]  = 0.5 * (-u - delta);
    x[2]  = 0.5 * (-u + delta);
  }

  return ireal;
}

template <typename T, unsigned int i, unsigned int j>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void CmpAndSwap(T *array)
{
  if (vecCore::MaskFull(array[i] > array[j])) {
    T c      = array[j];
    array[j] = array[i];
    array[i] = c;
  }
}

// a special function to sort a 4 element array
// sorting is done inplace and in increasing order
// implementation comes from a sorting network
template <typename T>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Sort4(T *array)
{
  CmpAndSwap<T, 0, 2>(array);
  CmpAndSwap<T, 1, 3>(array);
  CmpAndSwap<T, 0, 1>(array);
  CmpAndSwap<T, 2, 3>(array);
  CmpAndSwap<T, 1, 2>(array);
}

// solve quartic taken from ROOT/TGeo and adapted
//_____________________________________________________________________________

template <typename T>
VECCORE_ATT_HOST_DEVICE
int SolveQuartic(T a, T b, T c, T d, T *x)
{
  // Find real solutions of the quartic equation : x^4 + a*x^3 + b*x^2 + c*x + d = 0
  // Input: a,b,c,d
  // Output: x[4] - real solutions
  // Returns number of real solutions (0 to 3)
  T zero(0.0);
  T e     = b - T(3.) * a * a / T(8.);
  T f     = c + a * a * a / T(8.) - T(0.5) * a * b;
  T g     = d - T(3.) * a * a * a * a / T(256.) + a * a * b / T(16.) - a * c / T(4.);
  T xx[4] = {vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength};
  T delta;
  T h(0.0);
  unsigned int ireal = 0;

  // special case when f is zero
  using Bool_v = vecCore::Mask_v<T>;
  Bool_v done(false);
  Bool_v small = (Abs(f) < T(1.0e-6));
  Bool_v posDelta;
  if (!vecCore::MaskEmpty(small)) {
    delta    = e * e - 4. * g;
    posDelta = (delta >= zero);

    delta        = Sqrt(Max(delta, zero));
    h            = 0.5 * (-e - delta);
    Bool_v valid = !done && posDelta && (h > zero);

    if (!vecCore::MaskEmpty(valid)) {
      h = Sqrt(Max(h, zero));
      vecCore__MaskedAssignFunc(x[ireal++], valid, -h - 0.25 * a);
      vecCore__MaskedAssignFunc(x[ireal++], valid, h - 0.25 * a);
    }
    done = done | valid;

    h     = 0.5 * (-e + delta);
    valid = !done && posDelta && (h > zero);
    h     = Sqrt(Max(h, zero));
    if (!vecCore::MaskEmpty(valid)) {
      vecCore__MaskedAssignFunc(x[ireal++], valid, -h - 0.25 * a);
      vecCore__MaskedAssignFunc(x[ireal++], valid, h - 0.25 * a);
    }
    done = done | valid;

    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) {
      Sort4(x);
      return ireal;
    }
  } // end if(Abs(f) < 1.e-6)

  small = !done && (Abs(g) < T(1.e-6));
  if (!vecCore::MaskEmpty(small)) {
    vecCore__MaskedAssignFunc(x[ireal++], small, T(-0.25) * a);
    // this actually wants to solve a second order equation
    // we should specialize if it happens often
    unsigned int ncubicroots = SolveCubic<T>(0, e, f, xx);
    // this loop is not nice
    for (unsigned int i = 0; i < ncubicroots; ++i) {
      vecCore__MaskedAssignFunc(x[ireal++], small, xx[i] - T(0.25) * a);
    }
    done = done | small;

    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) {
      Sort4(x); // could be Sort3
      return ireal;
    }
  }

  ireal = SolveCubic<T>(T(2.) * e, e * e - T(4.) * g, -f * f, xx);
  h     = zero;
  if (ireal == 1) {
    done = done | (xx[0] <= h); // <= 0.
    vecCore::MaskedAssign(h, !done, Sqrt(Max(xx[0], h)));
  } else {
    // 3 real solutions of the cubic
    h = zero;
    for (unsigned int i = 0; i < 3; i++) {
      vecCore__MaskedAssignFunc(h, !done, xx[i]);
      if (vecCore::MaskFull(h >= zero)) break;
    }
    if (vecCore::MaskFull(done | (h <= zero))) {
      return 0;
    }
    vecCore__MaskedAssignFunc(h, !done, Sqrt(Max(h, zero)));
  }
  T j   = T(0.5) * (e + h * h - f / h);
  ireal = 0;
  delta = h * h - 4. * j;
  posDelta = !done && (delta >= zero);

  delta = Sqrt(Max(delta, zero));
  vecCore__MaskedAssignFunc(x[ireal++], posDelta, T(0.5) * (-h - delta) + T(-0.25) * a);
  vecCore__MaskedAssignFunc(x[ireal++], posDelta, T(0.5) * (-h + delta) + T(-0.25) * a);

  delta    = h * h - 4. * g / j;
  posDelta = !done && (delta >= zero);
  delta    = Sqrt(Max(delta, zero));
  vecCore__MaskedAssignFunc(x[ireal++], posDelta, T(0.5) * (h - delta) + T(-0.25) * a);
  vecCore__MaskedAssignFunc(x[ireal++], posDelta, T(0.5) * (h + delta) + T(-0.25) * a);

  Sort4(x);
  return ireal;
}

// template <typename T>
// class SIMDPlacedTorus;
// class SIMDUnplacedTorus;

class PlacedTorus;
class UnplacedTorus;


struct TorusImplementation {
  // using UnplacedStruct_t = ::vecgeom::TorusStruct<double>;
  // using UnplacedVolume_t = SIMDUnplacedTorus;
  // using PlacedShape_t    = SIMDPlacedTorus<UnplacedVolume_t>;
  using PlacedShape_t    = PlacedTorus;
  using UnplacedStruct_t = TorusStruct<double>;
  using UnplacedVolume_t = UnplacedTorus;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType() { printf("SpecializedTorus"); }

  template <typename Stream>
  static void PrintType(Stream &s, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    s << "SpecializedTorus<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &s)
  {
    s << "TorusImplementation";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s)
  {
    s << "UnplacedTorus";
  }

  template <typename Real_v, bool ForInside, bool doPhiTreatment>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void GenericKernelForContainsAndInside(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point,
                                                typename vecCore::Mask_v<Real_v> &completelyinside,
                                                typename vecCore::Mask_v<Real_v> &completelyoutside)

  {
    // using vecgeom::GenericKernels;
    // here we are explicitely unrolling the loop since a for statement will likely be a penalty
    // check if second call to Abs is compiled away
    // and it can anyway not be vectorized
    /* rmax */
    using Bool_v            = vecCore::Mask_v<Real_v>;
    constexpr Precision tol = 100. * vecgeom::kTolerance;

    Real_v rxy   = Sqrt(point[0] * point[0] + point[1] * point[1]);
    Real_v radsq = (rxy - torus.rtor()) * (rxy - torus.rtor()) + point[2] * point[2];

    if (ForInside) {
      completelyoutside = radsq > (tol * torus.rmax() + torus.rmax2()); // rmax
      completelyinside  = radsq < (-tol * torus.rmax() + torus.rmax2());
    } else {
      completelyoutside = radsq > torus.rmax2();
    }

    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(completelyoutside)) {
        return;
      }
    }
    /* rmin */
    if (ForInside) {
      completelyoutside |= radsq < (-tol * torus.rmin() + torus.rmin2()); // rmin
      completelyinside &= radsq > (tol * torus.rmin() + torus.rmin2());
    } else {
      completelyoutside |= radsq < torus.rmin2();
    }

    // NOT YET NEEDED WHEN NOT PHI TREATMENT
    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(completelyoutside)) {
        return;
      }
    }

    /* phi */
    if ((torus.dphi() < kTwoPi) && (doPhiTreatment)) {
      Bool_v completelyoutsidephi;
      Bool_v completelyinsidephi;
      torus.GetWedge().GenericKernelForContainsAndInside<Real_v, ForInside>(point, completelyinsidephi,
                                                                            completelyoutsidephi);

      completelyoutside |= completelyoutsidephi;
      if (ForInside) completelyinside &= completelyinsidephi;
    }
  }

  template <typename Real_v, typename Bool_v, bool doPhiTreatment>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void ContainsKernel(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    Bool_v unused;
    Bool_v outside;
    TorusImplementation::GenericKernelForContainsAndInside<Real_v, false, doPhiTreatment>(torus, point, unused, outside);
    inside = !outside;
  }

  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void InsideKernel(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Inside_t &inside)
  {

    using Bool_v = vecCore::Mask_v<Real_v>;
    //
    Bool_v completelyinside, completelyoutside;
    TorusImplementation::GenericKernelForContainsAndInside<Real_v, true, true>(torus, point, completelyinside,
                                                                               completelyoutside);
    inside = Inside_t(EInside::kSurface);
    vecCore::MaskedAssign(inside, completelyoutside, Inside_t(EInside::kOutside));
    vecCore::MaskedAssign(inside, completelyinside, Inside_t(EInside::kInside));
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void UnplacedContainsDisk(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    ContainsKernel<Real_v, Bool_v, false>(torus, point, inside);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point,
                       typename vecCore::Mask_v<Real_v> &contains)
  {
    using Bool_v = vecCore::Mask_v<Real_v>;
    Bool_v unused, outside;
    TorusImplementation::GenericKernelForContainsAndInside<Real_v, true, false>(torus, point, unused, outside);
    contains = !outside;
  }

  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Inside_t &inside)
  {
    TorusImplementation::InsideKernel<Real_v, Inside_t>(torus, point, inside);
  }

  template <class Real_v>
  VECCORE_ATT_HOST_DEVICE
  static Real_v DistSqrToTorusR(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point,
                                Vector3D<Real_v> const &dir, Real_v dist)
  {
    // Computes the squared distance to "axis" or "defining ring" of the torus from point point + t*dir;
    Vector3D<Real_v> p = point + dir * dist;
    Real_v rxy         = p.Perp();
    return (rxy - torus.rtor()) * (rxy - torus.rtor()) + p.z() * p.z();
  }

  template <typename Real_v, bool ForRmin>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static Real_v ToBoundary(UnplacedStruct_t const &torus, Vector3D<Real_v> const &pt, Vector3D<Real_v> const &dir,
                           Real_v radius, bool out)
  {
    // to be taken from ROOT
    // Returns distance to the surface or the torus from a point, along
    // a direction. Point is close enough to the boundary so that the distance
    // to the torus is decreasing while moving along the given direction.

    // Compute coeficients of the quartic
    Real_v s         = vecgeom::kInfLength;
    const Real_v tol = 100. * vecgeom::kTolerance;
    Real_v r0sq      = pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2];
    Real_v rdotn     = pt[0] * dir[0] + pt[1] * dir[1] + pt[2] * dir[2];
    Real_v rsumsq    = torus.rtor2() + radius * radius;
    Real_v a         = 4. * rdotn;
    Real_v b         = 2. * (r0sq + 2. * rdotn * rdotn - rsumsq + 2. * torus.rtor2() * dir[2] * dir[2]);
    Real_v c         = 4. * (r0sq * rdotn - rsumsq * rdotn + 2. * torus.rtor2() * pt[2] * dir[2]);
    Real_v d         = r0sq * r0sq - 2. * r0sq * rsumsq + 4. * torus.rtor2() * pt[2] * pt[2] +
               (torus.rtor2() - radius * radius) * (torus.rtor2() - radius * radius);

    Real_v x[4] = {vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength};
    int nsol    = 0;

    // special condition
    if (vecCore::MaskFull(Abs(dir[2]) < 1E-3 && Abs(pt[2]) < 0.1 * radius)) {
      Real_v r0        = torus.rtor() - Sqrt((radius - pt[2]) * (radius + pt[2]));
      Real_v invdirxy2 = 1. / (1 - dir.z() * dir.z());
      Real_v b0        = (pt[0] * dir[0] + pt[1] * dir[1]) * invdirxy2;
      Real_v c0        = (pt[0] * pt[0] + (pt[1] - r0) * (pt[1] + r0)) * invdirxy2;
      Real_v delta     = b0 * b0 - c0;
      if (vecCore::MaskFull(delta > 0)) {
        x[nsol] = -b0 - Sqrt(delta);
        if (vecCore::MaskFull(x[nsol] > -tol)) nsol++;
        x[nsol] = -b0 + Sqrt(delta);
        if (vecCore::MaskFull(x[nsol] > -tol)) nsol++;
      }
      r0    = torus.rtor() + Sqrt((radius - pt[2]) * (radius + pt[2]));
      c0    = (pt[0] * pt[0] + (pt[1] - r0) * (pt[1] + r0)) * invdirxy2;
      delta = b0 * b0 - c0;
      if (vecCore::MaskFull(delta > 0)) {
        x[nsol] = -b0 - Sqrt(delta);
        if (vecCore::MaskFull(x[nsol] > -tol)) nsol++;
        x[nsol] = -b0 + Sqrt(delta);
        if (vecCore::MaskFull(x[nsol] > -tol)) nsol++;
      }
      if (nsol) {
        Sort4(x);
      }
    } else { // generic case
      nsol = SolveQuartic(a, b, c, d, x);
    }
    if (!nsol) {
      return vecgeom::kInfLength;
    }

    // look for first positive solution
    Real_v ndotd;
    bool inner = vecCore::MaskFull(Abs(radius - torus.rmin()) < vecgeom::kTolerance);
    for (int i = 0; i < nsol; i++) {
      if (vecCore::MaskFull(x[i] < -10)) continue;

      Vector3D<Real_v> r0   = pt + x[i] * dir;
      Vector3D<Real_v> norm = r0;
      r0.z()                = 0.;
      r0.Normalize();
      r0 *= torus.rtor();
      norm -= r0;
      // norm = pt
      // for (unsigned int ipt = 0; ipt < 3; ipt++)
      //   norm[ipt] = pt[ipt] + x[i] * dir[ipt] - r0[ipt];
      // ndotd = norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2];
      ndotd = norm.Dot(dir);
      if (inner ^ out) {
        if (vecCore::MaskFull(ndotd < 0)) continue; // discard this solution
      } else {
        if (vecCore::MaskFull(ndotd > 0)) continue; // discard this solution
      }

      // The crossing point should be in the phi wedge
      if (torus.dphi() < vecgeom::kTwoPi) {
        if (!vecCore::MaskFull(torus.GetWedge().ContainsWithBoundary<Real_v>(r0))) continue;
      }

      s = x[i];
      // refine solution with Newton iterations
      Real_v eps   = vecgeom::kInfLength;
      Real_v delta = s * s * s * s + a * s * s * s + b * s * s + c * s + d;
      Real_v eps0  = -delta / (4. * s * s * s + 3. * a * s * s + 2. * b * s + c);
      int ntry     = 0;
      while (vecCore::MaskFull(Abs(eps) > vecgeom::kTolerance)) {
        if (vecCore::MaskFull(Abs(eps0) > 100)) break;
        s += eps0;
        if (vecCore::MaskFull(Abs(s + eps0) < vecgeom::kTolerance)) break;
        delta = s * s * s * s + a * s * s * s + b * s * s + c * s + d;
        eps   = -delta / (4. * s * s * s + 3. * a * s * s + 2. * b * s + c);
        if (vecCore::MaskFull(Abs(eps) >= Abs(eps0))) break;
        ntry++;
        // Avoid infinite recursion
        if (ntry > 100) break;
        eps0 = eps;
      }
      // discard this solution
      if (vecCore::MaskFull(s < -tol)) continue;
      return Max(Real_v(0.), s);
    }
    return vecgeom::kInfLength;
  }


  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point,
                           Vector3D<Real_v> const &dir, Real_v const &stepMax, Real_v &distance)
  {
    // Vector3D<Float_t> localPoint     = transformation.Transform<transCodeT, rotCodeT>(point);
    // Vector3D<Float_t> localDirection = transformation.TransformDirection<rotCodeT>(dir);

    using Bool_v   = vecCore::Mask_v<Real_v>;
    using Inside_v = vecCore::Index_v<Real_v>;

    ////////First naive implementation
    distance = kInfLength;

    // Check Bounding Cylinder first
    Bool_v inBounds;
    Bool_v done         = Bool_v(false);
    Inside_v inside     = Inside_v(EInside::kOutside);
    Real_v tubeDistance = kInfLength;

#ifndef VECGEOM_NO_SPECIALIZATION
    // call the tube functionality -- first of all we check whether we are inside
    // bounding volume
    TubeImplementation<TubeTypes::HollowTube>::Contains(torus.GetBoundingTube().GetStruct(), point, inBounds);

    // only need to do this check if all particles (in vector) are outside ( otherwise useless )
    TubeImplementation<TubeTypes::HollowTube>::DistanceToIn(torus.GetBoundingTube().GetStruct(), point,
                                                            dir, stepMax, tubeDistance);
#else
    // call the tube functionality -- first of all we check whether we are inside the bounding volume
    TubeImplementation<TubeTypes::UniversalTube>::Contains(torus.GetBoundingTube().GetStruct(), point, inBounds);

    // only need to do this check if all particles (in vector) are outside ( otherwise useless )
    // vecCore::Mask_v<Real_v> notInBounds { !inBounds };
    if (!vecCore::MaskFull(inBounds)) {
      TubeImplementation<TubeTypes::UniversalTube>::DistanceToIn(torus.GetBoundingTube().GetStruct(), point,
                                                                 dir, stepMax, tubeDistance);
    } else {
      tubeDistance = 0.;
    }
#endif // VECGEOM_NO_SPECIALIZATION

    if (!vecCore::MaskEmpty(inBounds)) {
      // Check points on the wrong side (inside torus)
      TorusImplementation::InsideKernel<Real_v, Inside_v>(torus, point, inside);
      if (vecCore::MaskFull(inside == Inside_v(EInside::kInside))) {
        done     = Bool_v(true);
        distance = Real_v(-1.);
      }
    } else {
      done = Bool_v(vecCore::MaskFull(tubeDistance == kInfLength));
    }

    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(done)) {
        return;
      }
    }

    // Propagate the point to the bounding tube, as this will reduce the
    // coefficients of the quartic and improve precision of the solutions
    Vector3D<Real_v> ptAtSurf = point + tubeDistance * dir;
    Bool_v hasphi = Bool_v(torus.dphi() < vecgeom::kTwoPi);
    if (vecCore::MaskFull(hasphi)) {
      Real_v d1, d2;

      auto wedge = torus.GetWedge();
      // checking distance to phi wedges
      // NOTE: if the tube told me its hitting surface, this would be unnecessary
      wedge.DistanceToIn<Real_v>(ptAtSurf, dir, d1, d2);

      // check phi intersections if bounding tube intersection is due to phi in which case we are done
      if (vecCore::MaskFull(d1 != kInfLength)) {
        Real_v daxis = DistSqrToTorusR(torus, ptAtSurf, dir, d1);
        if (vecCore::MaskFull(daxis >= torus.rmin2() && daxis < torus.rmax2())) {
          distance = d1;
          // check if tube intersections is due to phi in which case we are done
          if (vecCore::MaskFull(Abs(distance) < kTolerance)) {
            distance += tubeDistance;
            return;
          }
        }
      }

      if (vecCore::MaskFull(d2 != kInfLength)) {
        Real_v daxis = DistSqrToTorusR(torus, ptAtSurf, dir, d2);
        if (vecCore::MaskFull(daxis >= torus.rmin2() && daxis < torus.rmax2())) {
          distance = Min(d2, distance);
          // check if tube intersections is due to phi in which case we are done
          if (vecCore::MaskFull(Abs(distance) < kTolerance)) {
            distance += tubeDistance;
            return;
          }
        }
      }
      distance = kInfLength;
    }

    Real_v dd = ToBoundary<Real_v, false>(torus, ptAtSurf, dir, torus.rmax(), false);

    // in case of a phi opening we also need to check the Rmin surface
    if (torus.rmin() > 0.) {
      Real_v ddrmin = ToBoundary<Real_v, true>(torus, ptAtSurf, dir, torus.rmin(), false);
      dd            = Min(dd, ddrmin);
    }
    distance = Min(distance, dd);
    distance += tubeDistance;
    // This has to be added because distance can become > kInfLength due to
    // missing early returns in CUDA. This makes comparisons to kInfLength fail.
    if (vecCore::MaskFull(Abs(distance) > kInfLength)) distance = kInfLength;

    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                            Real_v const & /*stepMax*/, Real_v &distance)
  {
    using Inside_v = vecCore::Index_v<Real_v>;
    using Bool_v   = vecCore::Mask_v<Real_v>;

    bool hasphi  = (torus.dphi() < kTwoPi);
    bool hasrmin = (torus.rmin() > 0);

    //=== First step: for points outside --> return infinity
    Bool_v done = Bool_v(false);
    distance    = kInfLength;

    // // Check points on the wrong side (inside torus)
    // Inside_v inside;s
    // TorusImplementation::InsideKernel<Real_v, Inside_v>(torus, point, inside);
    // if (vecCore::MaskFull(inside == Inside_v(EInside::kOutside))) {
    //   distance = -1.;
    //   return;
    // }

    if (vecCore::EarlyReturnAllowed()) {
      // very simple calculations -- only if can save some time
      Real_v distz = Abs(point.z()) - torus.rmax();
      done |= distz > kHalfTolerance;

      // outside of bounding tube?
      Real_v rsq = point.x() * point.x() + point.y() * point.y();
      // Real_v rdotv = point.x()*dir.x() + point.y()*dir.y();
      Precision outerExclRadius = torus.rtor() + torus.rmax() + kHalfTolerance;
      done |= rsq > outerExclRadius * outerExclRadius;
      Precision innerExclRadius = torus.rtor() - torus.rmax() - kHalfTolerance;
      done |= rsq < innerExclRadius * innerExclRadius;
      vecCore__MaskedAssignFunc(distance, done, Real_v(-1.0));
    }

    //=== Use InsideKernel() for a quick check, and if outside --> return -1
    // Bool_t inside=false, outside=false;
    // GenericKernelForContainsAndInside<Backend,true,true>(torus, point, inside, outside);
    // MaskedAssign( inside, -1.0, &distance );
    // done |= inside;
    Inside_v locus;
    TorusImplementation::InsideKernel<Real_v, Inside_v>(torus, point, locus);
    vecCore__MaskedAssignFunc(distance, locus == EInside::kOutside, Real_v(-1.0));
    done |= locus == EInside::kOutside;
    vecCore__MaskedAssignFunc(distance, done, Real_v(-1.0));
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    Real_v dout = ToBoundary<Real_v, false>(torus, point, dir, torus.rmax(), true);
    // ToBoundary<Backend, false, true>(torus, point, dir, torus.rmax());
    Real_v din(kInfLength);
    if (hasrmin) {
      din = ToBoundary<Real_v, true>(torus, point, dir, torus.rmin(), true);
      // ToBoundary<Backend, true, true>(torus, point, dir, torus.rmin());
    }
    distance = Min(dout, din);
    // std::cout << "dout, din: " << dout << ", " << din << '\n';
    // std::cout << "distance = Min(dout, din): " << distance << '\n';

    if (hasphi) {
      Real_v distPhi1;
      Real_v distPhi2;
      // torus.GetWedge().DistanceToOut<Backend>(point, dir, distPhi1, distPhi2);
      torus.GetWedge().DistanceToOut<Real_v>(point, dir, distPhi1, distPhi2);
      Bool_v smallerphi = distPhi1 < distance;
      if (!vecCore::MaskEmpty(smallerphi)) {
        Vector3D<Real_v> intersectionPoint = point + dir * distPhi1;
        Bool_v insideDisk;
        // UnplacedContainsDisk<Backend>(torus, intersectionPoint, insideDisk);
        UnplacedContainsDisk<Real_v, Bool_v>(torus, intersectionPoint, insideDisk);

        if (!vecCore::MaskEmpty(insideDisk)) // Inside Disk
        {
          Real_v diri = intersectionPoint.x() * torus.GetWedge().GetAlong1().x() +
                        intersectionPoint.y() * torus.GetWedge().GetAlong1().y();
          Bool_v rightside = (diri >= 0);

          vecCore__MaskedAssignFunc(distance, rightside && smallerphi && insideDisk, distPhi1);
        }
      }
      smallerphi = distPhi2 < distance;
      if (!vecCore::MaskEmpty(smallerphi)) {

        Vector3D<Real_v> intersectionPoint = point + dir * distPhi2;
        Bool_v insideDisk;
        // UnplacedContainsDisk<Backend>(torus, intersectionPoint, insideDisk);
        UnplacedContainsDisk<Real_v, Bool_v>(torus, intersectionPoint, insideDisk);
        if (!vecCore::MaskEmpty(insideDisk)) // Inside Disk
        {
          Real_v diri2 = intersectionPoint.x() * torus.GetWedge().GetAlong2().x() +
                         intersectionPoint.y() * torus.GetWedge().GetAlong2().y();
          Bool_v rightside = (diri2 >= 0);
          vecCore__MaskedAssignFunc(distance, rightside && (distPhi2 < distance) && smallerphi && insideDisk, distPhi2);
        }
      }
    }

    vecCore__MaskedAssignFunc(distance, distance >= kInfLength, Real_v(-1.));
    vecCore__MaskedAssignFunc(distance, vecCore::math::Abs(distance) < vecgeom::kTolerance, Real_v(0.));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Real_v &safety)
  {
    Real_v rxy = Sqrt(point[0] * point[0] + point[1] * point[1]);
    Real_v rad = Sqrt((rxy - torus.rtor()) * (rxy - torus.rtor()) + point[2] * point[2]);
    safety     = torus.rmax() - rad;
    if (torus.rmin()) {
      safety = Min(rad - torus.rmin(), torus.rmax() - rad);
    }

    // TODO: extend implementation for phi sector case
    bool hasphi = (torus.dphi() < kTwoPi);
    if (hasphi) {
      Real_v safetyPhi = torus.GetWedge().SafetyToOut<Real_v>(point);
      safety           = Min(safetyPhi, safety);
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(UnplacedStruct_t const &torus, Vector3D<Real_v> const &point, Real_v &safety)
  {
    // implementation taken from TGeoTorus
    Real_v rxy = Sqrt(point[0] * point[0] + point[1] * point[1]);
    Real_v rad = Sqrt((rxy - torus.rtor()) * (rxy - torus.rtor()) + point[2] * point[2]);
    safety     = rad - torus.rmax();
    if (torus.rmin()) {
      safety = Max(torus.rmin() - rad, rad - torus.rmax());
    }

    bool hasphi = (torus.dphi() < kTwoPi);
    if (hasphi && vecCore::MaskFull(rxy != 0.)) {
      Real_v safetyPhi = torus.GetWedge().SafetyToIn<Real_v>(point);
      safety           = Max(safetyPhi, safety);
    }
  }

}; // end struct
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_
