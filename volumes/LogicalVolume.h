/**
 * @file logical_volume.h
 * @author Johannes de Fine Licht (johannes.definelicht@cern.ch)
 */

#ifndef VECGEOM_VOLUMES_LOGICALVOLUME_H_
#define VECGEOM_VOLUMES_LOGICALVOLUME_H_

#include "base/Global.h"
#include "base/Vector.h"
#include "volumes/UnplacedVolume.h"
//#include "volumes/PlacedVolume.h"

#include <string>
#include <list>
#include <set>
#include <typeindex>
#include <typeinfo>

class TGeoShape;
// class VUSolid;

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class LogicalVolume;);
VECGEOM_DEVICE_FORWARD_DECLARE(class VPlacedVolume;);
VECGEOM_DEVICE_FORWARD_DECLARE(class VLevelLocator;);
VECGEOM_DEVICE_FORWARD_DECLARE(class VSafetyEstimator;);
VECGEOM_DEVICE_FORWARD_DECLARE(class VNavigator;);

VECGEOM_DEVICE_DECLARE_CONV(class, LogicalVolume);

inline namespace VECGEOM_IMPL_NAMESPACE {

class VLevelLocator;
class VSafetyEstimator;
class VNavigator;
typedef VPlacedVolume const *Daughter;
class GeoManager;
/**
 * @brief Class responsible for storing the unplaced volume, material and
 *        daughter volumes of a mother volume.
 */
class LogicalVolume {
  friend class GeoManager;

private:
  // pointer to concrete unplaced volume/shape
  VUnplacedVolume const *fUnplacedVolume;

  int fId; // global id of logical volume object

  std::string *fLabel; // name of logical volume

  static int gIdCount; // static class counter

  /** a pointer member to register arbitrary objects with logical volume;
        included for the moment to model UserExtension like in TGeoVolume
  */
  void *fUserExtensionPtr;
  /** some specific pointers used by Geant-V
   *
   */
  void *fTrackingMediumPtr;
  void *fBasketManagerPtr;

  void const *fRegion = nullptr;        // pointer to a region object
  static std::type_index gRegionTypeId; // stores a type-index of the concrete type for fRegion
                                        // in order to disallow accidential Setting and Getting with
                                        // different types in one user application
  static bool gRegionTypeInitialized;

  VLevelLocator const *fLevelLocator;       // a locator class for this logical volume
  VSafetyEstimator const *fSafetyEstimator; // a safety estimator class for this logical volume
  VNavigator const *fNavigator;             // the attached navigator

  // the container of daughter (placed) volumes which are placed inside this logical
  // Volume
  Vector<Daughter> *fDaughters;

  using CudaDaughter_t = cuda::VPlacedVolume const *;
  friend class CudaManager;
  //  friend class GeoManager;

  // possibility to change pointer of daughter volumes ( can be used by GeoManager )
  //  void SetDaughter(unsigned int i, VPlacedVolume const *pvol);

public:
#ifndef VECGEOM_NVCC
  // Standard constructor when constructing geometries. Will initiate an empty
  // daughter list which can be populated by placing daughters.
  // \sa PlaceDaughter()
  LogicalVolume(char const *const label, VUnplacedVolume const *const unplaced_vol);

  LogicalVolume(VUnplacedVolume const *const unplaced_vol) : LogicalVolume("", unplaced_vol) {}

  //
  // copy operator since we have pointer data members
  //
  LogicalVolume(LogicalVolume const &other) = delete;
  LogicalVolume *operator=(LogicalVolume const &other) = delete;

#else
  VECGEOM_CUDA_HEADER_DEVICE
  LogicalVolume(VUnplacedVolume const *const unplaced_vol, Vector<Daughter> *GetDaughter);
#endif

  ~LogicalVolume();

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  VUnplacedVolume const *GetUnplacedVolume() const { return fUnplacedVolume; }

  // will be deprecated in favour of better encapsulation of internal storage
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Vector<Daughter> const &GetDaughters() const { return *fDaughters; }

  // will be deprecated in favour of better encapsulation of internal storage
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Vector<Daughter> const *GetDaughtersp() const { return fDaughters; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Vector<Daughter> *GetDaughtersp() { return fDaughters; }

  //  VECGEOM_CUDA_HEADER_BOTH
  //  VECGEOM_FORCE_INLINE
  //  VPlacedVolume const* GetDaughter(unsigned int i) const { return daughters_->operator[](i); }
  //
  //  VECGEOM_CUDA_HEADER_BOTH
  //  VECGEOM_FORCE_INLINE
  //  unsigned int GetNDaughters() const { return daughters_->size(); }

  // returns the total number of placed volume contained in this logical volume AND below
  size_t GetNTotal() const;

  VECGEOM_FORCE_INLINE
  void *GetUserExtensionPtr() const { return fUserExtensionPtr; }
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  void *GetTrackingMediumPtr() const { return fTrackingMediumPtr; }
  VECGEOM_FORCE_INLINE
  void *GetBasketManagerPtr() const { return fBasketManagerPtr; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  VLevelLocator const *GetLevelLocator() const { return fLevelLocator; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  void SetLevelLocator(VLevelLocator const *locator) { fLevelLocator = locator; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  VSafetyEstimator const *GetSafetyEstimator() const { return fSafetyEstimator; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  void SetSafetyEstimator(VSafetyEstimator const *est) { fSafetyEstimator = est; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  VNavigator const *GetNavigator() const { return fNavigator; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  void SetNavigator(VNavigator const *n) { fNavigator = n; }

  int id() const { return fId; }

  const char *GetName() const { return fLabel->c_str(); }
  std::string GetLabel() const { return *fLabel; }

  void SetLabel(char const *const label)
  {
    if (fLabel) delete fLabel;
    fLabel = new std::string(label);
  }

  VECGEOM_FORCE_INLINE
  void SetUserExtensionPtr(void *userpointer) { fUserExtensionPtr = userpointer; }

  VECGEOM_FORCE_INLINE
  void SetTrackingMediumPtr(void *tmediumpointer) { fTrackingMediumPtr = tmediumpointer; }

  VECGEOM_FORCE_INLINE
  void SetBasketManagerPtr(void *basketpointer) { fBasketManagerPtr = basketpointer; }

  VECGEOM_CUDA_HEADER_BOTH
  void Print(const int indent = 0) const;

  VECGEOM_CUDA_HEADER_BOTH
  void PrintContent(const int depth = 0) const;

  VPlacedVolume *Place(char const *const label, Transformation3D const *const transformation) const;

  VPlacedVolume *Place(Transformation3D const *const transformation) const;

  VPlacedVolume *Place(char const *const label) const;

  VPlacedVolume *Place() const;

  VPlacedVolume const *PlaceDaughter(char const *const label, LogicalVolume const *const volume,
                                     Transformation3D const *const transformation);

  VPlacedVolume const *PlaceDaughter(LogicalVolume const *const volume, Transformation3D const *const transformation);

  void PlaceDaughter(VPlacedVolume const *const placed);

  // returns true of at least one of the daughters is an assembly
  bool ContainsAssembly() const;

  friend std::ostream &operator<<(std::ostream &os, LogicalVolume const &vol);

  template <typename Region>
  Region const *GetRegion() const
  {
    // note that the use of a simple assert still does not prevent getting a different region class
    // than we wrote (in case of a release build)
    assert(std::type_index(typeid(Region)) == gRegionTypeId);
    return reinterpret_cast<Region const *>(fRegion);
  }

  // Region is a any object which will be associated
  // to this logical volume
  // *and* if pushdown=true to all logical volumes below (in the geom hierarchy)
  // set region will override any existing region
  // the user has to make sure to call in the right order
  template <typename Region>
  void SetRegion(Region const *region, bool pushdown = true);

#ifdef VECGEOM_CUDA_INTERFACE
  DevicePtr<cuda::LogicalVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const unplaced_vol,
                                           DevicePtr<cuda::Vector<CudaDaughter_t>> GetDaughter) const;
  DevicePtr<cuda::LogicalVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const unplaced_vol,
                                           DevicePtr<cuda::Vector<CudaDaughter_t>> GetDaughter,
                                           DevicePtr<cuda::LogicalVolume> const gpu_ptr) const;
#endif

private:
  std::set<LogicalVolume *> GetSetOfDaughterLogicalVolumes() const;

}; // End class

template <typename Region>
void LogicalVolume::SetRegion(Region const *region, bool pushdown)
{
  // make sure we only store Regions of the same type-id
  if (!gRegionTypeInitialized) {
    // first use of this function
    gRegionTypeId          = std::type_index(typeid(Region));
    gRegionTypeInitialized = true;
  }
  // throw error if trying to push different types
  if (std::type_index(typeid(Region)) != gRegionTypeId) {
    std::runtime_error("SetRegion called with unexpected type");
  }

  fRegion = reinterpret_cast<void const *>(region);
  if (pushdown) {
    for (auto &lv : GetSetOfDaughterLogicalVolumes()) {
      lv->SetRegion<Region>(region, true);
    }
  }
}

} // End inline namespace
} // End global namespace

#endif // VECGEOM_VOLUMES_LOGICALVOLUME_H_
