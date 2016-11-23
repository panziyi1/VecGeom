/*
 * PlacedCutTube.h
 *
 *  Created on: 03.11.2016
 *      Author: mgheata
 */
#ifndef VECGEOM_VOLUMES_PLACEDCUTTUBE_H_
#define VECGEOM_VOLUMES_PLACEDCUTTUBE_H_

#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/CutTubeImplementation.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedCutTube.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedCutTube;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedCutTube);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedCutTube : public PlacedVolumeImplHelper<UnplacedCutTube, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedCutTube, VPlacedVolume>;

public:
#ifndef VECGEOM_NVCC
  // constructor inheritance;
  using Base::Base;
  PlacedCutTube(char const *const label, LogicalVolume const *const logicalVolume,
                Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedCutTube(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                vecgeom::PlacedBox const *const boundingBox)
      : PlacedCutTube("", logicalVolume, transformation, boundingBox)
  {
  }
#else
  __device__ PlacedCutTube(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                           PlacedBox const *const boundingBox, const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }
#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedCutTube() {}

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedCutTube const *GetUnplacedVolume() const
  {
    return static_cast<UnplacedCutTube const *>(GetLogicalVolume()->GetUnplacedVolume());
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision rmin() const { return GetUnplacedVolume()->rmin(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision rmax() const { return GetUnplacedVolume()->rmax(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision z() const { return GetUnplacedVolume()->z(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision sphi() const { return GetUnplacedVolume()->sphi(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision dphi() const { return GetUnplacedVolume()->dphi(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Vector3D<Precision> BottomNormal() const { return GetUnplacedVolume()->BottomNormal(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Vector3D<Precision> TopNormal() const { return GetUnplacedVolume()->TopNormal(); }

#if !defined(VECGEOM_NVCC)
  /** @brief Interface method for computing capacity */
  virtual Precision Capacity() override { return GetUnplacedVolume()->volume(); }
  /** @brief Computes the extent on X/Y/Z of the trapezoid */
  virtual void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    GetUnplacedVolume()->Extent(aMin, aMax);
  }

  /** @brief Shortcut for computing the normal */
  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override
  {
    return GetUnplacedVolume()->Normal(point, normal);
  }

  /** @brief Generates randomly a point on the surface of the trapezoid */
  virtual Vector3D<Precision> GetPointOnSurface() const override { return GetUnplacedVolume()->GetPointOnSurface(); }

  /** @brief Implementation of surface area computation */
  virtual double SurfaceArea() override { return GetUnplacedVolume()->SurfaceArea(); }

#if defined(VECGEOM_USOLIDS)
  /** @brief Get type name */
  virtual std::string GetEntityType() const override { return GetUnplacedVolume()->GetEntityType(); }
#endif
#endif

  // CUDA specific
  /** @brief Memory size in bytes */
  VECGEOM_FORCE_INLINE
  virtual int memory_size() const override { return sizeof(*this); }

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

} // end inline namespace
} // End global namespace

#endif
