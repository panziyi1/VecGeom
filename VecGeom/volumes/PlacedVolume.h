// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \file PlacedVolume.h
/// \author created by Johannes de Fine Licht, Sandro Wenzel (CERN)

#ifndef VECGEOM_VOLUMES_PLACEDVOLUME_H_
#define VECGEOM_VOLUMES_PLACEDVOLUME_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/VolumeTypes.h"
#include <string>

#ifdef VECGEOM_GEANT4
#include <G4VSolid.hh>
#endif

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class VPlacedVolume;);
VECGEOM_DEVICE_DECLARE_CONV(class, VPlacedVolume);
#ifndef VECCORE_CUDA
template <>
struct kCudaType<const cxx::VPlacedVolume *> {
  using type_t = const cuda::VPlacedVolume *;
};
#endif

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedBox;
class GeoManager;

/*!
 * \brief A placed volume is a positioned logical volume.
 *
 * Special features of VecGeom are that:
 *
 *   1. Placed volumes offer
 *      similar geometry interfaces (distance, etc) as unplaced volumes for user convenience.
 *   2. We can have sub-types of placed volumes, specialized according
 *      to underlying unplaced volume type as well as per categories
 *      of placing transformations.
 *
 * Given a logical volume and a transformation, factory methods
 * will generate the most suitable sub-type of a placed volume so that
 * geometry APIs, such as VPlacedVolume::DistanceToIn(position, direction)
 * are optimized as much as possible. Some CMake flags allow to select how
 * far this is done or not.
 */
class VPlacedVolume {
  friend class GeoManager;

private:
  unsigned int id_;  ///< Integer id
  int copy_no_ = 0;  ///< Copy number for the physical volume, used by transport
  int ichild_  = -1; ///< Index in the mother volume list;

  // Use a pointer so the string won't be constructed on the GPU
  std::string *label_;            ///< Label/name of placed volume
  static unsigned int g_id_count; ///< Static instance counter

protected:
  VolumeTypes type = VolumeTypes::kUnknown;
  LogicalVolume const *logical_volume_; ///< Pointer to positioned logical volume
#ifdef VECGEOM_INPLACE_TRANSFORMATIONS
  Transformation3D fTransformation; ///< The positioning transformation
#else
  Transformation3D const *fTransformation; ///< The positioning transformation
#endif

#ifndef VECCORE_CUDA

  /** Constructor
   * \param label Name of logical volume.
   * \param logical_vol The logical volume to be positioned.
   * \param transform The positioning transformation.
   */
  VPlacedVolume(char const *const label, LogicalVolume const *const logical_vol,
                Transformation3D const *const transform);

  /** Constructor
   * \param logical_vol The logical volume to be positioned.
   * \param transform The positioning transformation.
   */
  VPlacedVolume(LogicalVolume const *const logical_vol, Transformation3D const *const transform)
      : VPlacedVolume("", logical_vol, transform)
  {
  }

#else
  /// CUDA version of constructor
  VECCORE_ATT_DEVICE VPlacedVolume(LogicalVolume const *const logical_vol, Transformation3D const *const transformation,
                                   unsigned int id, int copy_no, int ichild)
#ifdef VECGEOM_INPLACE_TRANSFORMATIONS
      : logical_volume_(logical_vol), fTransformation(*transformation), id_(id), copy_no_(copy_no), ichild_(ichild),
        label_(NULL)
  {
  }
#else
      : logical_volume_(logical_vol), fTransformation(transformation), id_(id), copy_no_(copy_no), ichild_(ichild),
        label_(NULL)
  {
  }
#endif
#endif

  VECGEOM_FORCE_INLINE
  void SetChildId(int index) { ichild_ = index; }

public:
  VECCORE_ATT_HOST_DEVICE
  VPlacedVolume(VPlacedVolume const &);
  VECCORE_ATT_HOST_DEVICE
  VPlacedVolume *operator=(VPlacedVolume const &);

  VECCORE_ATT_HOST_DEVICE
  virtual ~VPlacedVolume();

  /// Returns integer index associated to this volume.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  unsigned int id() const { return id_; }

  /// Returns copy number.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  int GetChildId() const { return ichild_; }

  /// LogicalVolume::PlaceDaughters is a friend that can set the child index
  friend void LogicalVolume::PlaceDaughter(VPlacedVolume *const placed);

  /// Returns copy number.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  int GetCopyNo() const { return copy_no_; }

  /// Returns value of static instance counter
  static unsigned int GetIdCount() { return g_id_count; }

  /// Returns name/label.
  std::string const &GetLabel() const { return *label_; }

  /// Returns underlying logical volume.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  LogicalVolume const *GetLogicalVolume() const { return logical_volume_; }

  /// Returns daughter container of logical volume.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Vector<Daughter> const &GetDaughters() const { return logical_volume_->GetDaughters(); }

  /// Finds the index of a given daughter having its pointer (linear complexity)
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  int IndexOf(Daughter daughter) const
  {
    int id = 0;
    for (auto d : logical_volume_->GetDaughters()) {
      if (d == daughter) return id;
      id++;
    }
    return -1;
  }

  /// Returns name/label.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  const char *GetName() const { return (*label_).c_str(); }

  /// Returns unplaced volume encapsulated in the logical volume.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  VUnplacedVolume const *GetUnplacedVolume() const { return logical_volume_->GetUnplacedVolume(); }

  /// Returns if underlying unplaced volume is an assembly.
  VECCORE_ATT_HOST_DEVICE
  bool IsAssembly() const { return GetUnplacedVolume()->IsAssembly(); }

  /// Returns underlying transformation.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Transformation3D const *GetTransformation() const
  {
#ifdef VECGEOM_INPLACE_TRANSFORMATIONS
    return &fTransformation;
#else
    return fTransformation;
#endif
  }
#ifndef VECCORE_CUDA
  SolidMesh *CreateMesh3D(size_t nFaces) const
  {
    return GetUnplacedVolume()->CreateMesh3D(*this->GetTransformation(), nFaces);
  }
#endif
  /// Sets logical volume.
  VECCORE_ATT_HOST_DEVICE
  void SetLogicalVolume(LogicalVolume const *const logical_vol) { logical_volume_ = logical_vol; }

  /// Sets transformation.
  VECCORE_ATT_HOST_DEVICE
  void SetTransformation(Transformation3D const *const transform)
  {
#ifdef VECGEOM_INPLACE_TRANSFORMATIONS
    fTransformation = *transform;
#else
    fTransformation = transform;
#endif
  }

  /// Sets name/label.
  void set_label(char const *label)
  {
    if (label_) delete label_;
    label_ = new std::string(label);
  }

  /// Sets copy number.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetCopyNo(int copy_no) { copy_no_ = copy_no; }

  friend std::ostream &operator<<(std::ostream &os, VPlacedVolume const &vol);

  /// Returns in-memory size in bytes of deriving objects (used to copy to GPU).
  virtual int MemorySize() const = 0;

  /// Print info about placed volume.
  VECCORE_ATT_HOST_DEVICE
  virtual void Print(const int indent = 0) const;

  /// Print info about placed volume.
  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const = 0;

  // some functions allowing for some very basic "introspection"

  /// Print the actual volume type to an outstream
  virtual void PrintType(std::ostream &os) const = 0;

  /// Recursively prints contained volumes to standard output.
  VECCORE_ATT_HOST_DEVICE
  void PrintContent(const int depth = 0) const;

  // Geometry functionality like in unplaced volume but taking the placement
  // into account.

  /*!
   * Returns whether a space point is contained or not in the placed volume.
   * This is similar to the functionality in VUnplacedVolume but taking into account
   * the positioning of the shape due to the placement.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual bool Contains(Vector3D<Precision> const &point) const;

  /*!
   * Returns whether a space point is contained or not in the placed volume.
   * Also returns the transformed position.
   *
   * \param point A given space point.
   * \param localPoint The point in the natural reference frame of the shape.
   *
   * This is similar to the functionality in VUnplacedVolume but taking into account
   * the positioning of the shape due to the placement.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual bool Contains(Vector3D<Precision> const &point, Vector3D<Precision> &localPoint) const;

  /// Direct dispatch to Contains of underlying unplaced volume without coordinate/placement transformation.
  VECCORE_ATT_HOST_DEVICE
  virtual bool UnplacedContains(Vector3D<Precision> const &localPoint) const;

  /**
   * Like similar function in VUnplacedVolume but taking into account
   * the positioning of the shape due to the placement.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual EnumInside Inside(Vector3D<Precision> const &point) const;

  /**
   * Like similar function in VUnplacedVolume but taking into account
   * the positioning of the shape due to the placement.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToIn(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                                 const Precision step_max = kInfLength) const;

  /**
   * Like similar function in VUnplacedVolume. Here position and direction are supposed to be
   * in the frame of the placed volume!
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToOut(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                                  Precision const step_max = kInfLength) const;

  /** A "placed" version of the DistanceToOut function; here
   * the point and direction are first of all transformed into the reference frame of the
   * shape. So given a position and direction in the reference frame in which the placed volume
   * is positioned, we transform everything into the coordinate system of the placed volume and
   * calculate DistanceToOut from there.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision PlacedDistanceToOut(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                                        Precision const step_max = kInfLength) const;

  /**
   * Like similar function in VUnplacedVolume but taking into account
   * the positioning of the shape due to the placement.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToIn(Vector3D<Precision> const &position) const;

  /**
   * Like similar function in VUnplacedVolume. Here position is supposed to be
   * in the frame of the placed volume.
   */
  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToOut(Vector3D<Precision> const &position) const;

  /// Simple forward to capacity on VUnplacedVolume
  virtual Precision Capacity();

  /// Get Extent of VUnplacedVolume, then apply transformation and recalculate
  VECCORE_ATT_HOST_DEVICE
  virtual void Extent(Vector3D<Precision> & /* min */, Vector3D<Precision> & /* max */) const;

  /// Get Normal of VUnplacedVolume, then apply transformation
  VECCORE_ATT_HOST_DEVICE
  virtual bool Normal(Vector3D<Precision> const & /*point*/, Vector3D<Precision> & /*normal*/) const;

  /// Like SurfaceArea on VUnplacedVolume
  virtual Precision SurfaceArea() const = 0;

public:
#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const                                                                       = 0;
  virtual DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                                   DevicePtr<cuda::Transformation3D> const transform,
                                                   DevicePtr<cuda::VPlacedVolume> const gpu_ptr) const      = 0;
  virtual DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                                   DevicePtr<cuda::Transformation3D> const transform) const = 0;

  template <typename Derived>
  DevicePtr<cuda::VPlacedVolume> CopyToGpuImpl(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                               DevicePtr<cuda::Transformation3D> const transform,
                                               DevicePtr<cuda::VPlacedVolume> const in_gpu_ptr) const
  {
    DevicePtr<CudaType_t<Derived>> gpu_ptr(in_gpu_ptr);
    gpu_ptr.Construct(logical_volume, transform, nullptr, this->id(), this->GetCopyNo(), this->GetChildId());
    CudaAssertError();
    // Need to go via the void* because the regular c++ compilation
    // does not actually see the declaration for the cuda version
    // (and thus can not determine the inheritance).
    return DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr);
  }
  template <typename Derived>
  DevicePtr<cuda::VPlacedVolume> CopyToGpuImpl(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                               DevicePtr<cuda::Transformation3D> const transform) const
  {
    DevicePtr<CudaType_t<Derived>> gpu_ptr;
    gpu_ptr.Allocate();
    return this->CopyToGpuImpl<Derived>(logical_volume, transform, DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr));
  }

#endif

#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
  /// A conversion function to a TGeoShape (when TGeo support is available)
  virtual TGeoShape const *ConvertToRoot() const = 0;
#endif
#ifdef VECGEOM_GEANT4
  /// A conversion function to a Geant4 G4VSolid (when G4 support is available)
  virtual G4VSolid const *ConvertToGeant4() const;
#endif
#endif // VECCORE_CUDA
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#ifdef VECCORE_CUDA

#define VECGEOM_DEVICE_INST_PLACED_VOLUME(PlacedVol)                                                       \
  namespace cxx {                                                                                          \
  template size_t DevicePtr<cuda::PlacedVol>::SizeOf();                                                    \
  template void DevicePtr<cuda::PlacedVol>::Construct(DevicePtr<cuda::LogicalVolume> const logical_volume, \
                                                      DevicePtr<cuda::Transformation3D> const transform,   \
                                                      const unsigned int id, const int copy_no,            \
                                                      const int child_id) const;                           \
  }

#endif

#endif // VECGEOM_VOLUMES_PLACEDVOLUME_H_
