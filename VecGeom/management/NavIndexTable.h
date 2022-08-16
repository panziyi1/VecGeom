// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Class storing the global navigation index lookup table.
/// \file management/NavIndexTable.h
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#ifndef VECGEOM_MANAGEMENT_NAVINDEXTABLE_H_
#define VECGEOM_MANAGEMENT_NAVINDEXTABLE_H_

#include "VecGeom/management/GeoVisitor.h"
#include "VecGeom/navigation/NavigationState.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

class NavIndexTable;

class BuildNavIndexVisitor : public GeoVisitorNavIndex {
private:
  size_t fTableSize = sizeof(
      NavIndex_t); ///< The iteration needs to be done in two passes, the first one is just to compute the array size
  NavIndex_t *fNavInd = nullptr; ///< Array storing the navigation info related to a given state
  int fLimitDepth     = 0;       ///< limit depth to scache the transformations, 0 means unlimited
  // NavIndex_t fWorld   = 1;       ///< index for the world volume
  NavIndex_t fCurrent = 1;     ///< Current navigation index being filled.
  bool fDoCount       = true;  ///< First pass to compute the table size
  bool fValidate      = false; ///< If set, this flag will force a thorough validation upon visiting

public:
  BuildNavIndexVisitor(int depth_limit, bool do_count)
      : GeoVisitorNavIndex(), fLimitDepth(depth_limit), fDoCount(do_count)
  {
  }

  size_t GetTableSize() const { return fTableSize; }
  void SetTable(NavIndex_t *table) { fNavInd = table; }
  void SetDoCount(bool flag) { fDoCount = flag; }
  void SetValidate(bool flag) { fValidate = flag; }

  NavIndex_t apply(NavStatePath *state, int level, NavIndex_t mother, int dind, NavIndex_t &id);
};

class NavIndexTable {
private:
  NavIndex_t *fNavInd       = nullptr; ///< address of the table
  VPlacedVolume *fVolBuffer = nullptr; ///< placed volume buffer
  NavIndex_t fWorld         = 1;       ///< index for the world volume
  size_t fTableSize         = 0;       ///< table size in bytes
  size_t fDepthLimit        = 0;       ///< depth limnit to which transformations will be cached

  NavIndexTable(NavIndex_t *table, size_t table_size) : fNavInd(table), fTableSize(table_size) {}

public:
  ~NavIndexTable() { CleanTable(); }
  static NavIndexTable *Instance(NavIndex_t *nav_table = nullptr, size_t table_size = 0)
  {
    static NavIndexTable instance(nav_table, table_size);
    return &instance;
  }

  void CleanTable()
  {
    delete[] fNavInd;
    fNavInd = nullptr;
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetTableSize() const { return fTableSize; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  NavIndex_t *GetTable() const { return fNavInd; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetVolumeBuffer(VPlacedVolume *buffer) { fVolBuffer = buffer; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  NavIndex_t GetWorld() const { return fWorld; }

  bool AllocateTable(size_t bytes)
  {
    bool success       = true;
    size_t nelem       = bytes / sizeof(NavIndex_t) + 1;
    NavIndex_t *buffer = nullptr;
    try {
      buffer = new NavIndex_t[nelem];
    } catch (std::bad_alloc &) {
      success = false;
      std::cout << "=== EEE === AlocateTable bad_alloc intercepted while trying to allocate " << bytes << " bytes\n";
    }
    if (success) {
      fNavInd    = buffer;
      fNavInd[0] = 0;
      fTableSize = bytes;
    }
    return success;
  }

  static size_t ComputeTableSize(VPlacedVolume const *top, int maxdepth, int depth_limit)
  {
    NavStatePath *state = NavStatePath::MakeInstance(maxdepth);
    BuildNavIndexVisitor visitor(depth_limit, true); // just count table size
    NavIndex_t id = 1;
    visitAllPlacedVolumesNavIndex(top, &visitor, state, id);
    size_t table_size = visitor.GetTableSize();
    NavStatePath::ReleaseInstance(state);
    return table_size;
  }

  bool CreateTable(VPlacedVolume const *top, int maxdepth, int depth_limit)
  {
    fDepthLimit         = depth_limit;
    NavStatePath *state = NavStatePath::MakeInstance(maxdepth);
    BuildNavIndexVisitor visitor(depth_limit, true); // just count table size

    NavIndex_t id = 1;
    visitAllPlacedVolumesNavIndex(top, &visitor, state, id);
    bool hasTable = AllocateTable(visitor.GetTableSize());
    if (!hasTable) return false;

    visitor.SetTable(fNavInd);
    visitor.SetDoCount(false);

    state->Clear();

    id = 1;
    visitAllPlacedVolumesNavIndex(top, &visitor, state, id);
    NavStatePath::ReleaseInstance(state);

    return true;
  }

  bool Validate(VPlacedVolume const *top, int maxdepth) const
  {
    NavStatePath *state = NavStatePath::MakeInstance(maxdepth);
    state->Clear();
    auto visitor = new BuildNavIndexVisitor(0, false);
    visitor->SetValidate(true);
    NavIndex_t id = 1;
    visitAllPlacedVolumesNavIndex(top, visitor, state, id);
    NavStatePath::ReleaseInstance(state);
    return true;
  }

  NavIndex_t ValidateState(NavStatePath *state);

  // vecgeom::cuda::NavIndexTable *CopyToGPU() const

  /// Traverses the geometry tree keeping track of the state context (volume path or navigation state)
  /// and applies the injected Visitor
  template <typename Visitor>
  static void visitAllPlacedVolumesNavIndex(VPlacedVolume const *currentvolume, Visitor *visitor, NavStatePath *state,
                                            NavIndex_t &id, int level = 0, NavIndex_t mother = 0, int dind = 0)
  {
    if (currentvolume != NULL) {
      state->Push(currentvolume);
      NavIndex_t new_mother = visitor->apply(state, level, mother, dind, id);
      int size              = currentvolume->GetDaughters().size();
      for (int i = 0; i < size; ++i) {
        visitAllPlacedVolumesNavIndex(currentvolume->GetDaughters().operator[](i), visitor, state, id, level + 1,
                                      new_mother, i);
      }
      state->Pop();
    }
  }
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
