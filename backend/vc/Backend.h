/// \file vc/backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BACKEND_VCBACKEND_H_
#define VECGEOM_BACKEND_VCBACKEND_H_

#include "base/Global.h"

#include "backend/scalar/Backend.h"

#include <Vc/Vc>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

struct kVc {
  typedef Vc::int_v int_v;
  typedef Vc::Vector<Precision> precision_v;
  typedef Vc::Vector<Precision>::Mask bool_v;
  typedef Vc::Vector<int> inside_v;
  const static precision_v kOne;
  const static precision_v kZero;
  const static bool_v kTrue;
  const static bool_v kFalse;
  // alternative typedefs ( might supercede above typedefs )
  typedef Vc::int_v Int_t;
  typedef Vc::Vector<Precision> Double_t;
  typedef Vc::Vector<Precision>::Mask Bool_t;
  typedef Vc::Vector<Precision> Index_t;
};

#ifdef VECGEOM_VC
constexpr int kVectorSize = kVc::precision_v::Size;
#define VECGEOM_BACKEND_TYPE vecgeom::kVc
#define VECGEOM_BACKEND_PRECISION_FROM_PTR(P) vecgeom::VcPrecision(P)
#define VECGEOM_BACKEND_PRECISION_TYPE vecgeom::VcPrecision
#define VECGEOM_BACKEND_PRECISION_TYPE_SIZE vecgeom::VcPrecision::Size
#define VECGEOM_BACKEND_PRECISION_NOT_SCALAR
#define VECGEOM_BACKEND_BOOL vecgeom::VcBool
#define VECGEOM_BACKEND_INSIDE vecgeom::kVc::inside_v
#endif

typedef kVc::int_v VcInt;
typedef kVc::precision_v VcPrecision;
typedef kVc::bool_v VcBool;
typedef kVc::inside_v VcInside;

// stores a vector type into a memory position ( normally an array ) toaddr
// toaddr has to be properly aligned
// this function is an abstraction for the Vc API "store"
template <typename Type>
VECGEOM_FORCE_INLINE
void StoreTo(typename Vc::Vector<Type> const &what, Type *toaddr)
{
  what.store(toaddr);
}

VECGEOM_FORCE_INLINE
void StoreTo(VcBool const &what, bool *toaddr)
{
  what.store(toaddr);
}

} // End inline namespace

} // End global namespace

#endif // VECGEOM_BACKEND_VCBACKEND_H_
