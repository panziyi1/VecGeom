#pragma once

#include "VecGeom/base/Global.h"
#include "VecGeom/base/SOA3D.h"

#include <algorithm>

#ifdef VECGEOM_DISTANCE_DEBUG
#include "VecGeom/volumes/utilities/ResultComparator.h"
#endif

namespace vecgeom {

// putting a forward declaration by hand
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(class, SpecializedVolImplHelper, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

template <class Specialization>
class SpecializedVolImplHelper : public Specialization::PlacedShape_t {

  using PlacedShape_t    = typename Specialization::PlacedShape_t;
  using UnplacedVolume_t = typename Specialization::UnplacedVolume_t;

public:
  using PlacedShape_t::Contains;
  using PlacedShape_t::DistanceToIn;
  using PlacedShape_t::DistanceToOut;
  using PlacedShape_t::Inside;
  using PlacedShape_t::PlacedShape_t;
  using PlacedShape_t::SafetyToIn;
  using PlacedShape_t::SafetyToOut;
  using PlacedShape_t::UnplacedContains;

}; // End class SpecializedVolImplHelper

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom
