// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Declaration of the unplaced volume interfaces.
/// \file volumes/UnplacedVolume.h
/// \author created by Sandro Wenzel

#ifndef VECGEOM_VOLUMES_UNPLACEDVOLUME_H_
#define VECGEOM_VOLUMES_UNPLACEDVOLUME_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Transformation3D.h"
#include "VecGeom/volumes/VolumeTypes.h"
#include <string>
#include <ostream>

#ifndef VECCORE_CUDA
#include "VecGeom/volumes/SolidMesh.h"
#endif

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class VUnplacedVolume;);
VECGEOM_DEVICE_DECLARE_CONV(class, VUnplacedVolume);

inline namespace VECGEOM_IMPL_NAMESPACE {

class LogicalVolume;
class VPlacedVolume;

/**
 * The abstract interface class for unplaced volumes.
 *
 * An unplaced volume represents a geometry shape (primitive) and offers
 * interfaces to query distance, location, containment, etc. in its "natural"
 * system of coordinates.
 */
class VUnplacedVolume {

private:
  friend class CudaManager;

protected:
  VolumeTypes fType = VolumeTypes::kUnknown;
  bool fGlobalConvexity;
  bool fIsAssembly = false; // indicates if this volume is an assembly

public:
  VECCORE_ATT_HOST_DEVICE
  virtual ~VUnplacedVolume() {}

  // ---------------- Contains --------------------------------------------------------------------

  /*!
   * Returns whether a space point pos is contained or not in the shape.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual bool Contains(Vector3D<Precision> const &pos) const;

  /*!
   * Returns whether a space point pos is inside, on the surface or outside
   * the shape. The surface is defined by a thickness constant.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual EnumInside Inside(Vector3D<Precision> const &pos) const;

  // ---------------- DistanceToOut functions -----------------------------------------------------

  /*!
   * Returns the distance from an internal or surface space point pos to the surface
   * of the shape along the normalized direction dir.
   * Does not have to look for surfaces beyond an optional distance of step_max.
   * Calling it with an outside point might result in undefined behaviour.
   *
   * TODO: Clarify return value in case step_max is non-default.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToOut(Vector3D<Precision> const &pos, Vector3D<Precision> const &dir,
                                  Precision step_max = kInfLength) const;

  // ---------------- SafetyToOut functions -----------------------------------------------------

  /*!
   * Returns the estimated minimum distance from an internal or surface space point pos to the
   * boundary of the shape. The estimate will be strictly smaller or equal to the true value.
   * Calling it with an outside point might result in undefined behaviour.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToOut(Vector3D<Precision> const &pos) const;

  // ---------------- DistanceToIn functions -----------------------------------------------------

  /*!
   * Returns the distance from an outside space point pos to the surface
   * of the shape along the normalized direction dir.
   * Does not have to look for surfaces beyond an optional distance of step_max.
   * Calling it with an inside point might result in undefined behaviour.
   *
   * TODO: Clarify return value in case step_max is non-default.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToIn(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                                 const Precision step_max = kInfLength) const;

  // ---------------- SafetyToIn functions -------------------------------------------------------

  /*!
   * Returns the estimated minimum distance from an outside or surface space point pos to the
   * boundary of the shape. The estimate will be strictly smaller or equal to the true value.
   * Calling it with an inside point is undefined behaviour.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToIn(Vector3D<Precision> const &pos) const;

  // ---------------- Normal ---------------------------------------------------------------------

  /*!
   * Calculates the surface normal unit vector for a space point pos, assuming
   * that pos is on the surface (i.e. Inside(pos) == kSurface).
   * The behaviour for a point not on the surface is undefined.
   * TODO: Clarify whether normal always points outwards.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual bool Normal(Vector3D<Precision> const &pos, Vector3D<Precision> &normal) const /* = 0 */;

  // ---------------- SamplePointOnSurface ----------------------------------------------------------
  /*!
   * Generates random point pos on the surface of the shape.
   * The returned point satisfies Inside(pos)==kSurface.
   */
  virtual Vector3D<Precision> SamplePointOnSurface() const /* = 0 */;

  // ----------------- Extent --------------------------------------------------------------------

  /*!
   * Returns the extent of the shape as corner points of the enclosing
   * bounding box.
   * @param aMin point of bounding box corner with minimum coordinates
   * @param aMax point of bounding box corner with maximum coordinates
   */
  VECCORE_ATT_HOST_DEVICE
  virtual void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const /* = 0 */;

  /*!
   *  Returns whether the shape is (globally) convex or not.
   *  If not known, returns false.
   */
  VECCORE_ATT_HOST_DEVICE
  bool IsConvex() const { return fGlobalConvexity; }

  /*!
   *  Returns whether the shape is an assembly
   */
  VECCORE_ATT_HOST_DEVICE
  bool IsAssembly() const { return fIsAssembly; }

  // ----------------- Capacity --------------------------------------------------------------------
  /*!
   *  Returns the (exact or estimated) cubic volume/capacity of the shape.
   */
  virtual Precision Capacity() const = 0;

  /*!
   *  Calculates an estimate of the cubic volume of the shape via a sampling technique.
   *  @param nStat number of sample points to be used
   */
  Precision EstimateCapacity(int nStat = 100000) const;

  // ----------------- Surface Area ----------------------------------------------------------------
  /*!
   *  Returns the (exact or estimated) surface area of the shape.
   */
  virtual Precision SurfaceArea() const = 0;

  /*!
   *  Calculates an estimate of the surface area of the shape via a sampling technique.
   *  @param nStat number of sample points to be used
   */
  Precision EstimateSurfaceArea(int nStat = 100000) const;

  /*!
   * Standard output operator for a textual representation.
   * (Uses the virtual method print(std::ostream &ps))
   */
  friend std::ostream &operator<<(std::ostream &os, VUnplacedVolume const &vol);

  /*!
   * Return the size of the deriving class in bytes. Necessary for
   * copying to the GPU.
   */
  virtual int MemorySize() const = 0;

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const = 0;

  /*!
   * Constructs the deriving class on the GPU and returns a pointer to GPU
   * memory where the object has been instantiated.
   */
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const                                               = 0;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const = 0;

  template <typename Derived, typename... ArgsTypes>
  DevicePtr<cuda::VUnplacedVolume> CopyToGpuImpl(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr,
                                                 ArgsTypes... params) const
  {
    DevicePtr<CudaType_t<Derived>> gpu_ptr(in_gpu_ptr);
    gpu_ptr.Construct(params...);
    CudaAssertError();
    // Need to go via the void* because the regular c++ compilation
    // does not actually see the declaration for the cuda version
    // (and thus can not determine the inheritance).
    return DevicePtr<cuda::VUnplacedVolume>((void *)gpu_ptr);
  }
  template <typename Derived>
  DevicePtr<cuda::VUnplacedVolume> CopyToGpuImpl() const
  {
    DevicePtr<CudaType_t<Derived>> gpu_ptr;
    gpu_ptr.Allocate();
    return this->CopyToGpu(DevicePtr<cuda::VUnplacedVolume>((void *)gpu_ptr));
  }

#endif

  /*!
   * Print a textual representation of the shape to a given outstream os.
   * This should typically tell the parameters, class, etc. of the shape.
   */
  virtual void Print(std::ostream &os) const = 0;

  /**
   * C-style printing for CUDA purposes.
   * TODO: clarify relation to other Print.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const = 0;

/// Generates mesh representation of the solid
#ifndef VECCORE_CUDA
  virtual SolidMesh *CreateMesh3D(Transformation3D const & /*trans*/, const size_t /*nSegments*/) const
  {
    return nullptr;
  };
#endif

#ifndef VECCORE_CUDA
  VPlacedVolume *PlaceVolume(char const *const label, LogicalVolume const *const volume,
                             Transformation3D const *const transformation, VPlacedVolume *const placement = NULL) const;
#endif

  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *PlaceVolume(LogicalVolume const *const volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                     const int id, const int copy_no, const int child_id,
#endif
                                     VPlacedVolume *const placement = NULL) const = 0;

};

std::ostream &operator<<(std::ostream &os, VUnplacedVolume const &vol);

} // namespace VECGEOM_IMPL_NAMESPACE

} // namespace vecgeom

#endif // VECGEOM_VOLUMES_UNPLACEDVOLUME_H_
