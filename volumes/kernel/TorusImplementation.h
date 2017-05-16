/// @file TorusImplementation.h

#ifndef VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_

//#include "base/Global.h"
#include "base/Transformation3D.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/kernel/TubeImplementation.h"
#include "volumes/UnplacedTorus.h"

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(struct, TorusImplementation, TranslationCode, translation::kGeneric,
                                        RotationCode, rotation::kGeneric);

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
  if (delta >= 0) {
    delta = Sqrt(delta);
    t     = (-3 * q * sq3 + delta) * inv6sq3;
    u     = (3 * q * sq3 + delta) * inv6sq3;
    x[0]  = CopySign(1., t) * Cbrt(Abs(t)) - CopySign(1., u) * Cbrt(Abs(u)) - a * ott;
  } else {
    delta = Sqrt(-delta);
    t     = -0.5 * q;
    u     = delta * inv6sq3;
    x[0]  = 2. * Pow(t * t + u * u, 0.5 * ott) * cos(ott * ATan2(u, t));
    x[0] -= a * ott;
  }

  t     = x[0] * x[0] + a * x[0] + b;
  u     = a + x[0];
  delta = u * u - 4. * t;
  if (delta >= 0) {
    ireal = 3;
    delta = Sqrt(delta);
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
  if (array[i] > array[j]) {
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
  T e     = b - 3. * a * a / 8.;
  T f     = c + a * a * a / 8. - 0.5 * a * b;
  T g     = d - 3. * a * a * a * a / 256. + a * a * b / 16. - a * c / 4.;
  T xx[4] = {vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength};
  T delta;
  T h                = 0.;
  unsigned int ireal = 0;

  // special case when f is zero
  if (Abs(f) < 1E-6) {
    delta = e * e - 4. * g;
    if (delta < 0) return 0;
    delta = Sqrt(delta);
    h     = 0.5 * (-e - delta);
    if (h >= 0) {
      h          = Sqrt(h);
      x[ireal++] = -h - 0.25 * a;
      x[ireal++] = h - 0.25 * a;
    }
    h = 0.5 * (-e + delta);
    if (h >= 0) {
      h          = Sqrt(h);
      x[ireal++] = -h - 0.25 * a;
      x[ireal++] = h - 0.25 * a;
    }
    Sort4(x);
    return ireal;
  }

  if (Abs(g) < 1E-6) {
    x[ireal++] = -0.25 * a;
    // this actually wants to solve a second order equation
    // we should specialize if it happens often
    unsigned int ncubicroots = SolveCubic<T>(0, e, f, xx);
    // this loop is not nice
    for (unsigned int i = 0; i < ncubicroots; i++)
      x[ireal++]        = xx[i] - 0.25 * a;
    Sort4(x); // could be Sort3
    return ireal;
  }

  ireal = SolveCubic<T>(2. * e, e * e - 4. * g, -f * f, xx);
  if (ireal == 1) {
    if (xx[0] <= 0) return 0;
    h = Sqrt(xx[0]);
  } else {
    // 3 real solutions of the cubic
    for (unsigned int i = 0; i < 3; i++) {
      h = xx[i];
      if (h >= 0) break;
    }
    if (h <= 0) return 0;
    h = Sqrt(h);
  }
  T j   = 0.5 * (e + h * h - f / h);
  ireal = 0;
  delta = h * h - 4. * j;
  if (delta >= 0) {
    delta      = Sqrt(delta);
    x[ireal++] = 0.5 * (-h - delta) - 0.25 * a;
    x[ireal++] = 0.5 * (-h + delta) - 0.25 * a;
  }
  delta = h * h - 4. * g / j;
  if (delta >= 0) {
    delta      = Sqrt(delta);
    x[ireal++] = 0.5 * (h - delta) - 0.25 * a;
    x[ireal++] = 0.5 * (h + delta) - 0.25 * a;
  }
  Sort4(x);
  return ireal;
}

class PlacedTorus;

template <TranslationCode transCodeT, RotationCode rotCodeT>
struct TorusImplementation {

  static const int transC = transCodeT;
  static const int rotC   = rotCodeT;

  using PlacedShape_t   = PlacedTorus;
  using UnplacedShape_t = UnplacedTorus;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType() { printf("SpecializedTorus<%i, %i>", transCodeT, rotCodeT); }

  template <typename Stream>
  static void PrintType(Stream &s)
  {
    s << "SpecializedTorus<" << transCodeT << "," << rotCodeT << ","
      << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &s)
  {
    s << "TorusImplemenation<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &s)
  {
    s << "UnplacedTorus";
  }

  /// GenericKernel Contains/Inside implementation
  template <typename Backend, bool ForInside, bool doPhiTreatment>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void GenericKernelForContainsAndInside(UnplacedTorus const &torus,
                                                Vector3D<typename Backend::precision_v> const &point,
                                                typename Backend::bool_v &completelyinside,
                                                typename Backend::bool_v &completelyoutside)

  {
    // using vecgeom::GenericKernels;
    // here we are explicitely unrolling the loop since  a for statement will likely be a penality
    // check if second call to Abs is compiled away
    // and it can anyway not be vectorized
    /* rmax */
    typedef typename Backend::precision_v Float_t;
    typedef typename Backend::bool_v Bool_t;
    constexpr Float_t tol = 100. * vecgeom::kTolerance;

    // very fast check on z-height
    completelyoutside = point[2] > MakePlusTolerant<ForInside>(torus.rmax());
    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(completelyoutside)) {
        return;
      }
    }

    Float_t rxy   = Sqrt(point[0] * point[0] + point[1] * point[1]);
    Float_t radsq = (rxy - torus.rtor()) * (rxy - torus.rtor()) + point[2] * point[2];

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

    /* phi */
    if ((torus.dphi() < kTwoPi) && (doPhiTreatment)) {
      if (vecCore::EarlyReturnAllowed()) {
        if (vecCore::MaskFull(completelyoutside)) {
          return;
        }
      }

      Bool_t completelyoutsidephi;
      Bool_t completelyinsidephi;
      torus.GetWedge().GenericKernelForContainsAndInside<Backend, ForInside>(point, completelyinsidephi,
                                                                             completelyoutsidephi);

      completelyoutside |= completelyoutsidephi;
      if (ForInside) completelyinside &= completelyinsidephi;
    }
  }

  template <class Backend, bool doPhiTreatment>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void ContainsKernel(UnplacedTorus const &torus, Vector3D<typename Backend::precision_v> const &point,
                             typename Backend::bool_v &inside)
  {
    typedef typename Backend::bool_v Bool_t;

    Bool_t unused = false, outside = false;
    GenericKernelForContainsAndInside<Backend, false, doPhiTreatment>(torus, point, unused, outside);
    inside = !outside;
  }

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  static void InsideKernel(UnplacedTorus const &torus, Vector3D<typename Backend::precision_v> const &point,
                           typename Backend::inside_v &inside)
  {

    typedef typename Backend::bool_v Bool_t;
    //
    Bool_t completelyinside = false, completelyoutside = false;
    GenericKernelForContainsAndInside<Backend, true, true>(torus, point, completelyinside, completelyoutside);
    inside = EInside::kSurface;
    vecCore::MaskedAssign(inside, completelyoutside, EInside::kOutside);
    vecCore::MaskedAssign(inside, completelyinside, EInside::kInside);
  }

  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void UnplacedContains(UnplacedTorus const &torus, Vector3D<typename Backend::precision_v> const &point,
                               typename Backend::bool_v &inside)
  {
    ContainsKernel<Backend, true>(torus, point, inside);
  }
  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void UnplacedContainsDisk(UnplacedTorus const &torus, Vector3D<typename Backend::precision_v> const &point,
                                   typename Backend::bool_v &inside)
  {
    ContainsKernel<Backend, false>(torus, point, inside);
  }
  template <typename Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(UnplacedTorus const &unplaced, Transformation3D const &transformation,
                       Vector3D<typename Backend::precision_v> const &masterPoint,
                       Vector3D<typename Backend::precision_v> &localPoint, typename Backend::bool_v &inside)
  {
    localPoint = transformation.Transform<transCodeT, rotCodeT>(masterPoint);
    UnplacedContains<Backend>(unplaced, localPoint, inside);
  }
  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedTorus const &torus, Transformation3D const &transformation,
                     Vector3D<typename Backend::precision_v> const &point, typename Backend::inside_v &inside)
  {
    InsideKernel<Backend>(torus, transformation.Transform<transCodeT, rotCodeT>(point), inside);
  }

  /////End GenericKernel Contains/Inside implementation

  template <class T>
  VECCORE_ATT_HOST_DEVICE
  static T CheckZero(T b, T c, T d, T e, T x)
  {
    T x2 = x * x;
    return x2 * x2 + b * x2 * x + c * x2 + d * x + e;
  }

  template <class T>
  VECCORE_ATT_HOST_DEVICE
  static T NewtonIter(T b, T c, T d, T /*e*/, T x, T fold)
  {
    T x2     = x * x;
    T fprime = 4 * x2 * x + 3 * b * x2 + 2 * c * x + d;
    return x - fold / fprime;
  }

  //_____________________________________________________________________________
  template <class T>
  VECCORE_ATT_HOST_DEVICE
  static T DistSqrToTorusR(UnplacedTorus const &torus, Vector3D<T> const &point, Vector3D<T> const &dir, T dist)
  {
    // Computes the squared distance to "axis" or "defining ring" of the torus from point point + t*dir;
    Vector3D<T> p = point + dir * dist;
    T rxy         = p.Perp();
    return (rxy - torus.rtor()) * (rxy - torus.rtor()) + p.z() * p.z();
  }

  template <class Backend, bool ForRmin, bool distToOut>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static typename Backend::precision_v ToBoundary(UnplacedTorus const &torus,
                                                  Vector3D<typename Backend::precision_v> const &pt,
                                                  Vector3D<typename Backend::precision_v> const &dir, Precision radius)
  {
    // to be taken from ROOT
    // Returns distance to the surface or the torus from a point, along
    // a direction. Point is close enough to the boundary so that the distance
    // to the torus is decreasing while moving along the given direction.

    typedef typename Backend::precision_v Real_v;
    typedef typename vecCore::Mask<Real_v> Bool_t;

    // Compute coeficients of the quartic
    Real_v s             = vecgeom::kInfLength;
    Real_v result        = s;
    constexpr Real_v tol = 100. * vecgeom::kTolerance;
    Real_v r0sq          = pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2];
    Real_v rdotv         = pt[0] * dir[0] + pt[1] * dir[1] + pt[2] * dir[2];
    Real_v rsumsq        = Real_v(torus.rtor2() + radius * radius);
    Real_v rdiffsq       = Real_v(torus.rtor2() - radius * radius);
    Real_v a             = 4. * rdotv;
    Real_v b             = 2. * (r0sq + 2. * rdotv * rdotv - rsumsq + 2. * torus.rtor2() * dir[2] * dir[2]);
    Real_v c             = 4. * (r0sq * rdotv - rsumsq * rdotv + 2. * torus.rtor2() * pt[2] * dir[2]);
    Real_v d             = r0sq * r0sq - 2. * r0sq * rsumsq + 4. * torus.rtor2() * pt[2] * pt[2] + rdiffsq * rdiffsq;

    Real_v x[4] = {vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength, vecgeom::kInfLength};
    int nsol    = 0;

    // special condition
    if (Abs(dir[2]) < 1E-3 && Abs(pt[2]) < 0.1 * radius) {
      Real_v r0        = torus.rtor() - Sqrt((radius - pt[2]) * (radius + pt[2]));
      Real_v invdirxy2 = 1. / (1 - dir.z() * dir.z());
      Real_v b0        = (pt[0] * dir[0] + pt[1] * dir[1]) * invdirxy2;
      Real_v c0        = (pt[0] * pt[0] + (pt[1] - r0) * (pt[1] + r0)) * invdirxy2;
      Real_v delta     = b0 * b0 - c0;
      if (delta > 0) {
        x[nsol] = -b0 - Sqrt(delta);
        if (x[nsol] > -tol) nsol++;
        x[nsol] = -b0 + Sqrt(delta);
        if (x[nsol] > -tol) nsol++;
      }
      r0    = torus.rtor() + Sqrt((radius - pt[2]) * (radius + pt[2]));
      c0    = (pt[0] * pt[0] + (pt[1] - r0) * (pt[1] + r0)) * invdirxy2;
      delta = b0 * b0 - c0;
      if (delta > 0) {
        x[nsol] = -b0 - Sqrt(delta);
        if (x[nsol] > -tol) nsol++;
        x[nsol] = -b0 + Sqrt(delta);
        if (x[nsol] > -tol) nsol++;
      }
      if (nsol) {
        Sort4(x);
      }
    } else {
      // generic case
      nsol = SolveQuartic(a, b, c, d, x);
    }
    if (!nsol) return vecgeom::kInfLength;

    // look for first positive solution
    Real_v ndotd;
    bool inner  = Abs(radius - torus.rmin()) < vecgeom::kTolerance;
    Bool_t done = Backend::kFalse;
    for (int i = 0; i < nsol; i++) {
      if (x[i] < -10) continue;

      Vector3D<Precision> r0   = pt + x[i] * dir;
      Vector3D<Precision> norm = r0;
      r0.z()                   = 0.;
      r0.Normalize();
      r0 *= torus.rtor();
      norm -= r0;
      // norm = pt
      // for (unsigned int ipt = 0; ipt < 3; ipt++)
      //   norm[ipt] = pt[ipt] + x[i] * dir[ipt] - r0[ipt];
      // ndotd = norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2];
      ndotd       = norm.Dot(dir);
      Bool_t skip = done;
      // if (inner ^ distToOut) {
      //   if (ndotd < 0) continue; // discard this solution
      // } else {
      //   if (ndotd > 0) continue; // discard this solution
      // }
      skip = skip || ((inner ^ distToOut) && (ndotd < 0.0));
      skip = skip || (!(inner ^ distToOut) && (ndotd > 0.0));

      // The crossing point could be in the phi wedge
      if (torus.dphi() < vecgeom::kTwoPi) {
        skip = skip || (!torus.GetWedge().ContainsWithBoundary<Backend>(r0));
      }

      vecCore__MaskedAssignFunc(s, !skip, x[i]);
      // refine solution with Newton iterations
      Real_v eps0 = vecgeom::kInfLength;
      Real_v delta, eps;
      // Real_v delta = s * s * s * s + a * s * s * s + b * s * s + c * s + d;
      // Real_v eps0  = -delta / (4. * s * s * s + 3. * a * s * s + 2. * b * s + c);
      Bool_t doneWhile = done || skip;
      int ntry         = 0;
      while (doneWhile && !vecCore::MaskEmpty(Abs(eps) > 2. * vecgeom::kTolerance)) {
        // delta = s * s * s * s + a * s * s * s + b * s * s + c * s + d;
        // eps   = -delta / (4. * s * s * s + 3. * a * s * s + 2. * b * s + c);
        vecCore__MaskedAssignFunc(delta, !done, (((s + a) * s + b) * s + c) * s + d);
        vecCore__MaskedAssignFunc(eps, !done, -delta / (((4. * s + 3. * a) * s + 2. * b) * s + c));

        doneWhile |= Abs(eps) > 100;
        doneWhile |= Abs(eps) >= Abs(eps0);
        vecCore__MaskedAssignFunc(s, !doneWhile, s + eps0);
        // std::cout<<" i="<< i <<", ntry="<< ntry <<", s="<< s <<", eps="<< eps <<", eps0="<< eps0 <<"\n";

        doneWhile |= delta < 0.1 * kTolerance;
        eps0 = eps;
        ++ntry;
      }
      // ignore skipped lanes
      const Bool_t converged = doneWhile && !skip && (s > -tol);
      vecCore__MaskedAssignFunc(result, converged, s);
      done = done || converged;

      if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return result;
    }
    return result;
  }

  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(UnplacedTorus const &torus, Transformation3D const &transformation,
                           Vector3D<typename Backend::precision_v> const &masterPoint,
                           Vector3D<typename Backend::precision_v> const &masterDirection,
                           typename Backend::precision_v const &stepMax, typename Backend::precision_v &distance)
  {
    typedef typename Backend::precision_v Float_t;
    typedef typename Backend::bool_v Bool_t;
    // std::cout<<" torus D2I(): masterPoint="<< masterPoint <<", masterDir="<< masterDirection <<"\n";

    Vector3D<Float_t> point = transformation.Transform<transCodeT, rotCodeT>(masterPoint);
    Vector3D<Float_t> dir   = transformation.TransformDirection<rotCodeT>(masterDirection);
    // std::cout<<" torus D2I(): point="<< point <<", dir="<< dir <<"\n";

    //=== First, for points outside and moving away --> return infinity
    Bool_t done = Backend::kFalse;
    distance    = kInfLength;

    // outside of Z range and going away?
    Float_t distz = Abs(point.z()) - torus.rmax();
    done |= distz > kHalfTolerance && point.z() * dir.z() >= 0;

    // outside of bounding tube and going away?
    Float_t rsq           = point.x() * point.x() + point.y() * point.y();
    Float_t rdotv         = point.x() * dir.x() + point.y() * dir.y();
    Precision outerRadius = torus.rtor() + torus.rmax() + kHalfTolerance;
    done |= rsq > outerRadius * outerRadius && rdotv >= 0;
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    //=== Next, check all dimensions of the tube, whether points are inside --> return -1
    // Bool_t inside=false, outside=false;
    // GenericKernelForContainsAndInside<Backend,true,true>(torus, point, inside, outside);
    // MaskedAssign( inside, -1.0, &distance );
    // done |= inside;
    typedef typename Backend::precision_v Real_v;
    typename Backend::inside_v locus;
    InsideKernel<Backend>(torus, point, locus);
    vecCore__MaskedAssignFunc(distance, locus == EInside::kInside, Real_v(-1.0));
    done |= locus == EInside::kInside;

    // and for points on surface and going away --> ignore distance to surface

    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(done)) return;
    }

    // if on surface and going out --> return infinity
    // if( !(inside || outside) )

    // Check Bounding Cylinder first
    Bool_t inBounds;
    typename Backend::inside_v inside = EInside::kOutside;
    Float_t tubeDistance              = kInfLength;

#ifndef VECGEOM_NO_SPECIALIZATION
    // call the tube functionality -- first of all we check whether we are inside bounding volume
    TubeImplementation<TubeTypes::HollowTube>::Contains(torus.GetBoundingTube().GetStruct(), point, inBounds);

    // only need to do this check if all particles (in vector) are outside ( otherwise useless )
    TubeImplementation<TubeTypes::HollowTube>::DistanceToIn(torus.GetBoundingTube().GetStruct(), point, dir, stepMax,
                                                            tubeDistance);
#else
    // call the tube functionality -- first of all we check whether we are inside bounding volume
    TubeImplementation<TubeTypes::UniversalTube>::Contains(torus.GetBoundingTube().GetStruct(), point, inBounds);

    // only need to do this check if all particles (in vector) are outside ( otherwise useless )
    if (!inBounds)
      TubeImplementation<TubeTypes::UniversalTube>::DistanceToIn(torus.GetBoundingTube().GetStruct(), point, dir,
                                                                 stepMax, tubeDistance);
    else
      tubeDistance = -1.;
#endif // VECGEOM_NO_SPECIALIZATION

    if (inBounds) {
      // Check points on the wrong side (inside torus)
      TorusImplementation::InsideKernel<Backend>(torus, point, inside);
      if (inside == EInside::kInside) {
        done     = true;
        distance = -1.;
      }
    } else {
      done = (tubeDistance == kInfLength);
    }

    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(done)) {
        return;
      }
    }

    // Propagate the point to the bounding tube, as this will reduce the
    // coefficients of the quartic and improve precision of the solutions
    point += tubeDistance * dir;
    Bool_t hasphi = (torus.dphi() < vecgeom::kTwoPi);
    if (hasphi) {
      Float_t d1, d2;

      auto wedge = torus.GetWedge();
      // checking distance to phi wedges
      // NOTE: if the tube told me its hitting surface, this would be unnecessary
      wedge.DistanceToIn<Backend>(point, dir, d1, d2);

      // check phi intersections if bounding tube intersection is due to phi in which case we are done
      if (d1 != kInfLength) {
        Float_t daxis    = DistSqrToTorusR(torus, point, dir, d1);
        Bool_t contained = (daxis >= torus.rmin2() && daxis < torus.rmax2());
        vecCore::MaskedAssign(distance, contained, d1);
        // // check if tube intersections is due to phi in which case we are done
        // if (Abs(distance) < kTolerance) {
        //   distance += tubeDistance;
        //   return;
        // }
      }

      if (d2 != kInfLength) {
        Float_t daxis    = DistSqrToTorusR(torus, point, dir, d2);
        Bool_t contained = (daxis >= torus.rmin2() && daxis < torus.rmax2());
        vecCore::MaskedAssign(distance, contained && d2 < distance, d2);
        // // check if tube intersections is due to phi in which case we are done
        // if (Abs(distance) < kTolerance) {
        //   distance += tubeDistance;
        //   return;
        // }
      }
      // distance = kInfLength;
    }

    Float_t dd = ToBoundary<Backend, false, false>(torus, point, dir, torus.rmax());

    // in case of a phi opening we also need to check the Rmin surface
    if (torus.rmin() > 0.) {
      Float_t ddrmin = ToBoundary<Backend, true, false>(torus, point, dir, torus.rmin());
      dd             = Min(dd, ddrmin);
    }
    distance = Min(distance, dd);
    distance += tubeDistance;
    // This has to be added because distance can become > kInfLength due to
    // missing early returns in CUDA. This makes comparisons to kInfLength fail.
    vecCore::MaskedAssign(distance, distance > kInfLength, kInfLength);
    return;
  }

  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(UnplacedTorus const &torus, Vector3D<typename Backend::precision_v> const &point,
                            Vector3D<typename Backend::precision_v> const &dir,
                            typename Backend::precision_v const & /*stepMax*/, typename Backend::precision_v &distance)
  {

    typedef typename Backend::precision_v Float_t;
    typedef typename Backend::bool_v Bool_t;
    distance = kInfLength;

    bool hasphi  = (torus.dphi() < kTwoPi);
    bool hasrmin = (torus.rmin() > 0);

    // Check points on the wrong side (inside torus)
    typename Backend::inside_v inside;
    TorusImplementation::InsideKernel<Backend>(torus, point, inside);
    if (inside == EInside::kOutside) {
      distance = -1.;
      return;
    }

    Float_t dout = ToBoundary<Backend, false, true>(torus, point, dir, torus.rmax());
    Float_t din(kInfLength);
    if (hasrmin) {
      din = ToBoundary<Backend, true, true>(torus, point, dir, torus.rmin());
    }
    distance = Min(dout, din);

    if (hasphi) {
      Float_t distPhi1;
      Float_t distPhi2;
      torus.GetWedge().DistanceToOut<Backend>(point, dir, distPhi1, distPhi2);
      Bool_t smallerphi = distPhi1 < distance;
      if (!vecCore::MaskEmpty(smallerphi)) {
        Vector3D<Float_t> intersectionPoint = point + dir * distPhi1;
        Bool_t insideDisk;
        UnplacedContainsDisk<Backend>(torus, intersectionPoint, insideDisk);

        if (!vecCore::MaskEmpty(insideDisk)) // Inside Disk
        {
          Float_t diri = intersectionPoint.x() * torus.GetWedge().GetAlong1().x() +
                         intersectionPoint.y() * torus.GetWedge().GetAlong1().y();
          Bool_t rightside = (diri >= 0);

          vecCore::MaskedAssign(distance, rightside && smallerphi && insideDisk, distPhi1);
        }
      }
      smallerphi = distPhi2 < distance;
      if (!vecCore::MaskEmpty(smallerphi)) {

        Vector3D<Float_t> intersectionPoint = point + dir * distPhi2;
        Bool_t insideDisk;
        UnplacedContainsDisk<Backend>(torus, intersectionPoint, insideDisk);
        if (!vecCore::MaskEmpty(insideDisk)) // Inside Disk
        {
          Float_t diri2 = intersectionPoint.x() * torus.GetWedge().GetAlong2().x() +
                          intersectionPoint.y() * torus.GetWedge().GetAlong2().y();
          Bool_t rightside = (diri2 >= 0);
          vecCore::MaskedAssign(distance, rightside && (distPhi2 < distance) && smallerphi && insideDisk, distPhi2);
        }
      }
    }
    if (distance >= kInfLength) distance                             = -1.;
    if (vecCore::math::Abs(distance) < vecgeom::kTolerance) distance = 0.;
  }

  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(UnplacedTorus const &torus, Transformation3D const &transformation,
                         Vector3D<typename Backend::precision_v> const &masterPoint,
                         typename Backend::precision_v &safety)
  {

    typedef typename Backend::precision_v Float_t;
    Vector3D<Float_t> point = transformation.Transform<transCodeT, rotCodeT>(masterPoint);

    // implementation taken from TGeoTorus
    Float_t rxy = Sqrt(point[0] * point[0] + point[1] * point[1]);
    Float_t rad = Sqrt((rxy - torus.rtor()) * (rxy - torus.rtor()) + point[2] * point[2]);
    safety      = rad - torus.rmax();
    if (torus.rmin()) {
      safety = Max(torus.rmin() - rad, rad - torus.rmax());
    }

    bool hasphi = (torus.dphi() < kTwoPi);
    if (hasphi && (rxy != 0.)) {
      Float_t safetyPhi = torus.GetWedge().SafetyToIn<Backend>(point);
      safety            = Max(safetyPhi, safety);
    }
  }

  template <class Backend>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(UnplacedTorus const &torus, Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety)
  {

    typedef typename Backend::precision_v Float_t;
    Float_t rxy = Sqrt(point[0] * point[0] + point[1] * point[1]);
    Float_t rad = Sqrt((rxy - torus.rtor()) * (rxy - torus.rtor()) + point[2] * point[2]);
    safety      = torus.rmax() - rad;
    if (torus.rmin()) {
      safety = Min(rad - torus.rmin(), torus.rmax() - rad);
    }

    // TODO: extend implementation for phi sector case
    bool hasphi = (torus.dphi() < kTwoPi);
    if (hasphi) {
      Float_t safetyPhi = torus.GetWedge().SafetyToOut<Backend>(point);
      safety            = Min(safetyPhi, safety);
    }
  }
}; // end struct
}
} // end namespace

#endif // VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_
