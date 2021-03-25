#ifndef VECGEOM_VOLUMES_PLACEDTUBE_H_
#define VECGEOM_VOLUMES_PLACEDTUBE_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/kernel/TubeImplementation.h"
#include "VecGeom/volumes/PlacedVolImplHelper.h"
#include "VecGeom/volumes/UnplacedTube.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedTube;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedTube);

inline namespace VECGEOM_IMPL_NAMESPACE {

// the base class of all placed tubes
// exists for stronger typing reasons and to be able
// to do runtime type inference on placed volumes

class PlacedTube : public PlacedVolumeImplHelper<UnplacedTube> {
  using Base = PlacedVolumeImplHelper<UnplacedTube, VPlacedVolume>;

public:
  using Base::Base;

#ifndef VECCORE_CUDA
  PlacedTube(char const *const label, LogicalVolume const *const logical_volume,
             Transformation3D const *const transformation)
      : Base(label, logical_volume, transformation)
  {
  }

  PlacedTube(LogicalVolume const *const logical_volume, Transformation3D const *const transformation)
      : PlacedTube("", logical_volume, transformation)
  {
  }
#else
  VECCORE_ATT_DEVICE PlacedTube(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
                                const int id, const int copy_no, const int child_id)
      : Base(logical_volume, transformation, id, copy_no, child_id)
  {
  }
#endif
  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedTube() {}

  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif

#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_PLACEDTUBE_H_
