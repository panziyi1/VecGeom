/// \file UnplacedRootVolume.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "VecGeom/volumes/UnplacedRootVolume.h"

#include "VecGeom/volumes/PlacedRootVolume.h"

#include <stdio.h>

namespace vecgeom {

void UnplacedRootVolume::Print() const
{
  printf("UnplacedRootVolume");
}

void UnplacedRootVolume::Print(std::ostream &os) const
{
  os << "UnplacedRootVolume";
}

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedRootVolume::PlaceVolume(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                               const int id, const int copy_no, const int child_id,
#endif
                                               VPlacedVolume *const placement) const
{
  if (placement) {
    return new (placement) PlacedRootVolume(logical_volume, transformation
#ifdef VECCORE_CUDA
                                            , id, copy_no, child_id
#endif
    );
    return placement;
  }
  return new PlacedRootVolume(logical_volume, transformation
#ifdef VECCORE_CUDA
                              , id, copy_no, child_id
#endif
  );
}

#ifdef VECGEOM_CUDA_INTERFACE
DevicePtr<cuda::VUnplacedVolume> UnplacedRootVolume::CopyToGpu() const
{
  assert(0 && "Attempted to copy unsupported ROOT volume to GPU.");
  return DevicePtr<cuda::VUnplacedVolume>(nullptr);
}
DevicePtr<cuda::VUnplacedVolume> UnplacedRootVolume::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> /* gpu_ptr */) const
{
  assert(0 && "Attempted to copy unsupported ROOT volume to GPU.");
  return DevicePtr<cuda::VUnplacedVolume>(nullptr);
}
#endif

} // End namespace vecgeom
