/// \file PlacedScaledShape.h
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_PLACEDSCALEDSHAPE_H_
#define VECGEOM_VOLUMES_PLACEDSCALEDSHAPE_H_

#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedScaledShape.h"
#include "volumes/kernel/ScaledShapeImplementation.h"
#include "volumes/PlacedVolImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedScaledShape;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedScaledShape);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedScaledShape : public PlacedVolumeImplHelper<UnplacedScaledShape, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedScaledShape, VPlacedVolume>;

public:
#ifndef VECGEOM_NVCC

  using Base::Base;
  PlacedScaledShape(char const *const label, LogicalVolume const *const logicalVolume,
                    Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedScaledShape(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                    vecgeom::PlacedBox const *const boundingBox)
      : PlacedScaledShape("", logicalVolume, transformation, boundingBox)
  {
  }

#else

  __device__ PlacedScaledShape(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                               PlacedBox const *const boundingBox, const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }

#endif
  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedScaledShape() {}

  // Accessors

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedScaledShape const *GetUnplacedVolume() const
  {
    return static_cast<UnplacedScaledShape const *>(GetLogicalVolume()->GetUnplacedVolume());
  }

#if !defined(VECGEOM_NVCC)
  virtual Precision Capacity() override { return GetUnplacedVolume()->Volume(); }

  virtual void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    GetUnplacedVolume()->Extent(aMin, aMax);
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override
  {
    return GetUnplacedVolume()->Normal(point, normal);
  }

  virtual Vector3D<Precision> GetPointOnSurface() const override { return GetUnplacedVolume()->GetPointOnSurface(); }

  virtual double SurfaceArea() override { return GetUnplacedVolume()->SurfaceArea(); }

#if defined(VECGEOM_USOLIDS)
  virtual std::string GetEntityType() const override { return GetUnplacedVolume()->GetEntityType(); }
#endif
#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual void PrintType() const override;

  /** @brief Print type name */
  void PrintType(std::ostream &os) const override;

  // CUDA specific

  virtual int memory_size() const override { return sizeof(*this); }

// Comparison specific

#ifndef VECGEOM_NVCC
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
  virtual ::VUSolid const *ConvertToUSolids() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECGEOM_NVCC
};
}
} // End global namespace

#endif // VECGEOM_VOLUMES_PLACEDSCALEDSHAPE_H_
