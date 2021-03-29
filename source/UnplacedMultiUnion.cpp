#include "VecGeom/volumes/PlacedMultiUnion.h"
#include "VecGeom/volumes/UnplacedMultiUnion.h"
#include "VecGeom/base/RNG.h"
#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

Precision UnplacedMultiUnion::Capacity() const
{
  if (fMultiUnion.fCapacity >= 0.) {
    return fMultiUnion.fCapacity;
  }
  // Sample in the solid extent and estimate capacity by counting how many points are
  // sampled inside
  const size_t nsamples = 100000;
  fMultiUnion.fCapacity = EstimateCapacity(nsamples);
  return fMultiUnion.fCapacity;
}

Precision UnplacedMultiUnion::SurfaceArea() const
{
  if (fMultiUnion.fSurfaceArea >= 0.) {
    return fMultiUnion.fSurfaceArea;
  }
  // Sample points on components, count how many are on the solid surface to correct
  // the surface area of each component
  const size_t nsamples    = 10000;
  fMultiUnion.fSurfaceArea = 0.;
  for (size_t i = 0; i < GetNumberOfSolids(); ++i) {
    size_t nsurf = 0;
    for (size_t ip = 0; ip < nsamples; ++ip) {
      Vector3D<Precision> point =
          GetNode(i)->GetTransformation()->InverseTransform(GetNode(i)->GetUnplacedVolume()->SamplePointOnSurface());
      if (Inside(point) == vecgeom::kSurface) nsurf++;
    }
    fMultiUnion.fSurfaceArea += GetNode(i)->SurfaceArea() * nsurf / nsamples;
  }
  return fMultiUnion.fSurfaceArea;
}

Vector3D<Precision> UnplacedMultiUnion::SamplePointOnSurface() const
{
  // Select a random component
  auto counter                = 0;
  VPlacedVolume const *volume = nullptr;
  Vector3D<Precision> point;
  size_t id = 0;
  do {
    if (counter == 0) {
      id     = (size_t)RNG::Instance().uniform(0., GetNumberOfSolids());
      volume = GetNode(id);
    }
    point   = volume->GetTransformation()->InverseTransform(volume->GetUnplacedVolume()->SamplePointOnSurface());
    counter = (counter + 1) % 1000;
  } while (Inside(point) != vecgeom::kSurface);
  return point;
}

bool UnplacedMultiUnion::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal to solid in a point
  bool valid = false;
  normal     = MultiUnionImplementation::NormalKernel<Precision>(fMultiUnion, point, valid);
  return valid;
}

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedMultiUnion::PlaceVolume(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                               const int id, const int copy_no, const int child_id,
#endif
                                               VPlacedVolume *const placement) const
{
  if (placement) {
    return new (placement) PlacedMultiUnion(logical_volume, transformation
#ifdef VECCORE_CUDA
                                            , id, copy_no, child_id
#endif
    );
    return placement;
  }
  return new PlacedMultiUnion(logical_volume, transformation
#ifdef VECCORE_CUDA
                              , id, copy_no, child_id
#endif
  );
}

#if defined(VECGEOM_CUDA_INTERFACE) && defined(VECGEOM_CUDA_HYBRID2)

DevicePtr<cuda::VUnplacedVolume> UnplacedMultiUnion::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<UnplacedMultiUnion>(in_gpu_ptr);
}

DevicePtr<cuda::VUnplacedVolume> UnplacedMultiUnion::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedMultiUnion>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedMultiUnion>::SizeOf();
template void DevicePtr<cuda::UnplacedMultiUnion>::Construct() const;

} // namespace cxx

#elif defined(VECGEOM_CUDA_INTERFACE) && !defined(VECGEOM_CUDA_HYBRID2)

namespace cuda {
// class UnplacedMultiUnion {};
}
inline namespace cxx {

template <>
size_t DevicePtr<cuda::PlacedMultiUnion>::SizeOf()
{
  return 0;
}
// template size_t DevicePtr<cuda::PlacedMultiUnion>::SizeOf();

template <>
template <>
void DevicePtr<cuda::PlacedMultiUnion>::Construct(DevicePtr<vecgeom::cuda::LogicalVolume>,
                                                  DevicePtr<vecgeom::cuda::Transformation3D>, unsigned int, int, int) const
{
  return;
}
// template void DevicePtr<cuda::PlacedMultiUnion>>::Construct() const;

} // namespace cxx
#endif

} // namespace vecgeom
