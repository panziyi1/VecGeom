#ifndef VECGEOM_VOLUMES_PLACEDTORUS_H_
#define VECGEOM_VOLUMES_PLACEDTORUS_H_

#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/TorusImplementation.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedTorus.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedTorus;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedTorus);

inline namespace VECGEOM_IMPL_NAMESPACE {

// the base class of all placed torus
// exists for stronger typing reasons and to be able
// to do runtime type inference on placed volumes
class PlacedTorus : public PlacedVolumeImplHelper<UnplacedTorus, VPlacedVolume> {
  // some common functionality for all placed tori
  // like constructors
  using Base = PlacedVolumeImplHelper<UnplacedTorus, VPlacedVolume>;

public:
  // inherit base-class constructors
  using Base::Base;

#ifndef VECCORE_CUDA
  PlacedTorus(char const *const label, LogicalVolume const *const logical_volume,
	      Transformation3D const *const transformation, ::vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logical_volume, transformation, boundingBox)
  { }

  PlacedTorus(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
               ::vecgeom::PlacedBox const *const boundingBox)
      : PlacedTorus("", logical_volume, transformation, boundingBox)
  { }

#else
  VECCORE_ATT_DEVICE
  PlacedTorus(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
               PlacedBox const *const boundingBox, const int id)
      : Base(logical_volume, transformation, boundingBox, id)
  { }
#endif

  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedTorus() {}

#ifndef VECCORE_CUDA
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;

#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif

#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif
};

class UnplacedTorus;
// a placed tube knowing abouts its volume/structural specialization
template <typename UnplacedTorus_t>
class SIMDPlacedTorus : public PlacedVolumeImplHelper<UnplacedTorus_t, PlacedTorus> {
  using Base = PlacedVolumeImplHelper<UnplacedTorus_t, PlacedTorus>;

public:
  typedef UnplacedTorus UnplacedShape_t;
  using Base::Base;
};

////  Maybe this is needed for TorusTypes specializations
// // a placed tube knowing abouts its volume/structural specialization
// template <typename UnplacedTorus_t>
// class SIMDPlacedTorus : public PlacedVolumeImplHelper<UnplacedTorus_t, PlacedTorus> {
//   using Base = PlacedVolumeImplHelper<UnplacedTorus_t, PlacedTorus>;

// public:
//   typedef UnplacedTorus UnplacedShape_t;
//   using Base::Base;
// };

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_PLACEDTORUS_H_
