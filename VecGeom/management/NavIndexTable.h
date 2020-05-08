// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Class storing the global navigation index lookup table.
/// \file management/NavIndexTable.h
/// \author Andrei Gheata (andrei.gheata@cern.ch)

#ifndef VECGEOM_MANAGEMENT_NAVINDEXTABLE_H_
#define VECGEOM_MANAGEMENT_NAVINDEXTABLE_H_

#include "VecGeom/navigation/NavigationState.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

class BuildNavIndexVisitor : public GeoVisitorNavIndex {
private:
  size_t fTableSize =
      0; ///< The iteration needs to be done in two passes, the first one is just to compute the array size
  NavIndex_t *fNavInd = nullptr; ///< Array storing the navigation info related to a given state
  int fLimitDepth     = 0;       ///< limit depth to scache the transformations, 0 means unlimited
  NavIndex_t fCurrent = 0;       ///< Current navigation index being filled.
  bool fDoCount       = true;    ///< First pass to compute the table size

public:
  BuildNavIndexVisitor(int depth_limit, bool do_count)
      : GeoVisitorNavIndex(), fLimitDepth(depth_limit), fDoCount(do_count)
  {
  }

  size_t GetTableSize() const { return fTableSize; }
  void SetTable(NavIndex_t *table) { fNavInd = table; }
  void SetDoCount(bool flag) { fDoCount = flag; }

  NavIndex_t apply(NavigationState *state, int level, NavIndex_t mother, int dind)
  {

    bool cacheTrans = true;
    size_t nd       = 0;
    // Check if matrix has to be cached for this node
    if (fLimitDepth > 0) {
      auto lv = state->Top()->GetLogicalVolume();
      nd      = lv->GetDaughters().size();
      if (level > fLimitDepth && !lv->IsReqCaching()) cacheTrans = false;
    }
    // Size in bytes of the current node data
    size_t current_size =
        3 * sizeof(unsigned int) + int(cacheTrans) * 12 * sizeof(double) + (nd + nd % 2) * sizeof(unsigned short);
    if (fDoCount) {
      fTableSize += current_size;
      return 0;
    }

    // Add data for the current element.

    // Fill the mother index for the current node
    fNavInd[fCurrent] = mother;

    // Fill the node index in the mother list of daughters
    if (mother > 0) {
      auto content_dind_mother  = (unsigned short *)(&fNavInd[mother + 3]);
      content_dind_mother[dind] = fCurrent;
    }

    // Physical volume index
    fNavInd[fCurrent + 1] = (level > 0) ? state->ValueAt(level - 1) : 0;

    // Write current level in next byte
    auto content_ddt = (unsigned char *)(&fNavInd[fCurrent + 2]);
    assert(level < std::numeric_limits<unsigned char>::max() && "fatal: geometry deph more than 255 not supported");
    *content_ddt = (unsigned char)level;

    // Write number of daughters in next 2 bytes
    auto content_nd = (unsigned short *)(content_ddt + 1);
    assert(nd < std::numeric_limits<unsigned short>::max() && "fatal: not supporting more than 65535 daughters");
    *content_nd = (unsigned short)nd;

    // Write the flag if matrix is stored in the next byte
    auto content_hasm = (bool *)(content_ddt + 3);
    *content_hasm     = cacheTrans;

    // Prepare the space for the daughter indices
    auto content_dind = (unsigned short *)(&fNavInd[fCurrent + 3]);
    for (auto i = 0; i < nd; ++i)
      content_dind[i] = 0;

    fCurrent += 3 + ((nd + nd % 2) >> 1);

    if (!cacheTrans) return fCurrent;

    // Write the transformation elements
    Transformation3D mat;
    state->TopMatrix(mat);

    auto content_mat = (double *)(&fNavInd[fCurrent]);
    for (auto i = 0; i < 3; ++i)
      content_mat[i] = mat.Translation(i);
    for (auto i = 0; i < 9; ++i)
      content_mat[i + 3] = mat.Rotation(i);

    // Set new value for fCurrent
    fCurrent += 24;
    return fCurrent;
  }
};

class NavIndexTable {
private:
  NavIndex_t *fNavInd = nullptr; ///< address of the table
  size_t fTableSize   = 0;       ///< table size in bytes
  size_t fDepthLimit  = 0;       ///< depth limnit to which transformations will be cached

  NavIndexTable(NavIndex_t *table, size_t table_size) : fNavInd(table), fTableSize(table_size) {}

public:
  ~NavIndexTable() { delete[](char *) fNavInd; }
  static NavIndexTable *Instance(NavIndex_t *nav_table = nullptr, size_t table_size = 0)
  {
    static NavIndexTable instance(nav_table, table_size);
    return &instance;
  }

  bool AllocateTable(size_t bytes)
  {
    bool success = true;
    char *buffer = nullptr;
    try {
      buffer = new char[bytes];
    } catch (std::bad_alloc &) {
      success = false;
    }
    if (success) {
      fNavInd    = (NavIndex_t *)buffer;
      fTableSize = bytes;
    }
    return success;
  }

  bool CreateTable(VPlacedVolume const *top, int maxdepth, int depth_limit)
  {
    fDepthLimit            = depth_limit;
    NavigationState *state = NavigationState::MakeInstance(maxdepth);
    state->Clear();
    auto visitor = new BuildNavIndexVisitor(depth_limit, true); // just count table size

    visitAllPlacedVolumesNavIndex(top, visitor, state);
    bool hasTable = AllocateTable(visitor->GetTableSize());
    if (!hasTable) return false;

    std::cout << "Navigation index table size is " << fTableSize << "bytes\n";

    visitor->SetTable(fNavInd);
    visitor->SetDoCount(false);

    state->Clear();
    visitAllPlacedVolumesNavIndex(top, visitor, state);
    return true;
  }

  // vecgeom::cuda::NavIndexTable *CopyToGPU() const

  /// Traverses the geometry tree keeping track of the state context (volume path or navigation state)
  /// and applies the injected Visitor
  template <typename Visitor>
  void visitAllPlacedVolumesNavIndex(VPlacedVolume const *currentvolume, Visitor *visitor, NavigationState *state,
                                     int level = 0, NavIndex_t mother = 0, int dind = 0) const
  {
    if (currentvolume != NULL) {
      state->Push(currentvolume);
      NavIndex_t new_mother = visitor->apply(state, level, mother, dind);
      int size              = currentvolume->GetDaughters().size();
      for (int i = 0; i < size; ++i) {
        visitAllPlacedVolumesNavIndex(currentvolume->GetDaughters().operator[](i), visitor, state, level + 1,
                                      new_mother, i);
      }
      state->Pop();
    }
  }
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
