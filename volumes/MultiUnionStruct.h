#ifndef VECGEOM_VOLUMES_MULTIUNIONSTRUCT_H_
#define VECGEOM_VOLUMES_MULTIUNIONSTRUCT_H_

#include "base/Global.h"
#include "management/HybridManager2.h"
#include "navigation/HybridNavigator2.h"
#include "management/ABBoxManager.h"
#ifndef VECCORE_CUDA
#include <atomic>
#endif

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 @brief Struct containing placed volume components representing a multiple union solid
 @author mihaela.gheata@cern.ch
*/

struct MultiUnionStruct {
  // template <typename U>
  // using vector_t     = vecgeom::Vector<U>;
  using BVHStructure = HybridManager2::HybridBoxAccelerationStructure;

  vector_t<VPlacedVolume const *> fVolumes; ///< Component placed volumes
  BVHStructure *fNavHelper = nullptr;       ///< Navigation helper using bounding boxes

  Vector3D<double> fMinExtent;      ///< Minimum extent
  Vector3D<double> fMaxExtent;      ///< Maximum extent
  mutable double fCapacity    = -1; ///< Capacity of the multiple union
  mutable double fSurfaceArea = -1; ///< Surface area of the multiple union
#ifndef VECCORE_CUDA
  mutable std::atomic<size_t> fLast; ///< Last located component for opportunistic relocation
#endif
  unsigned int fNboxes = 0;
  size_t **fNeighbours = nullptr; ///< Array of lists of neigbours
  size_t *fNneighbours = nullptr; ///< Number of neighbours for each component

  size_t *fBuffer = nullptr; ///< Scratch space for storing neighbours

  VECCORE_ATT_HOST_DEVICE
  MultiUnionStruct()
  {
    fMinExtent.Set(kInfLength);
    fMaxExtent.Set(-kInfLength);
#ifndef VECCORE_CUDA
    fLast.store(0);
#endif
  }

  VECCORE_ATT_HOST_DEVICE
  ~MultiUnionStruct()
  {
    delete[] fNeighbours;
    delete[] fNneighbours;
    delete[] fBuffer;
  }

  VECCORE_ATT_HOST_DEVICE
  void AddNode(VPlacedVolume const *volume)
  {
    using vecCore::math::Min;
    using vecCore::math::Max;
    Vector3D<double> amin, amax;
    ABBoxManager::ComputeABBox(volume, &amin, &amax);
    fMinExtent.Set(Min(fMinExtent.x(), amin.x()), Min(fMinExtent.y(), amin.y()), Min(fMinExtent.z(), amin.z()));
    fMaxExtent.Set(Max(fMaxExtent.x(), amax.x()), Max(fMaxExtent.y(), amax.y()), Max(fMaxExtent.z(), amax.z()));
    fVolumes.push_back(volume);
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool ABBoxOverlap(Vector3D<double> const &amin1, Vector3D<double> const &amax1, Vector3D<double> const &amin2,
                    Vector3D<double> const &amax2)
  {
    // Check if two aligned boxes overlap
    if ((amax1 - amin2).Min() < -kTolerance || (amax2 - amin1).Min() < -kTolerance) return false;
    return true;
  }

  VECCORE_ATT_HOST_DEVICE
  void Close()
  {
    // This method prepares the navigation structure
    using Boxes_t           = ABBoxManager::ABBoxContainer_t;
    using BoxCorner_t       = ABBoxManager::ABBox_s;
    fNboxes                 = fVolumes.size();
    BoxCorner_t *boxcorners = new BoxCorner_t[2 * fNboxes];
    Vector3D<double> amin, amax;
    for (size_t i = 0; i < fNboxes; ++i)
      ABBoxManager::ComputeABBox(fVolumes[i], &boxcorners[2 * i], &boxcorners[2 * i + 1]);
    Boxes_t boxes = &boxcorners[0];
    fNavHelper    = HybridManager2::Instance().BuildStructure(boxes, fNboxes);
    // Compute the lists of possibly overlapping neighbours
    fNeighbours  = new size_t *[fNboxes];
    fNneighbours = new size_t[fNboxes];
    memset(fNneighbours, 0, fNboxes * sizeof(size_t));
    fBuffer = new size_t[fNboxes * fNboxes];
    for (size_t i = 0; i < fNboxes; ++i) {
      fNeighbours[i] = fBuffer + i * fNboxes;
    }
    size_t newsize = 0;
    for (size_t i = 0; i < fNboxes - 1; ++i) {
      for (size_t j = i + 1; j < fNboxes; ++j) {
        if (ABBoxOverlap(boxcorners[2 * i], boxcorners[2 * i + 1], boxcorners[2 * j], boxcorners[2 * j + 1])) {
          fNeighbours[i][fNneighbours[i]++] = j;
          fNeighbours[j][fNneighbours[j]++] = i;
          newsize += 2;
        }
      }
    }
    // Compacting buffer of neighbours
    size_t *buffer  = new size_t[newsize];
    size_t *nextloc = buffer;
    for (size_t i = 0; i < fNboxes; ++i) {
      memcpy(nextloc, fNeighbours[i], fNneighbours[i] * sizeof(size_t));
      fNeighbours[i] = nextloc;
      nextloc += fNneighbours[i];
    }
    delete[] fBuffer;
    fBuffer = buffer;
  }

}; // End struct

} // End impl namespace
} // End global namespace

#endif /* VECGEOM_VOLUMES_MULTIUNIONSTRUCT_H_ */
