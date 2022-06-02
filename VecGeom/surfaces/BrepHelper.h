#ifndef VECGEOM_SURFACE_BREPHELPER_H_
#define VECGEOM_SURFACE_BREPHELPER_H_

#include <cassert>
#include <functional>
#include <VecGeom/surfaces/Model.h>

#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/management/GeoManager.h"


namespace vgbrep {

template <typename Real_t>
constexpr Real_t Tolerance() { return 0; }

template <>
constexpr double Tolerance() { return 1.e-9; }

template <>
constexpr float Tolerance() { return 1.e-4; }

template <typename Real_t>
bool ApproxEqual(Real_t t1, Real_t t2)
{
  return std::abs(t1 - t2) <= Tolerance<Real_t>();
}

template <typename Real_t>
bool ApproxEqual(vecgeom::Vector3D<Real_t> const &v1, vecgeom::Vector3D<Real_t> const &v2)
{
  return ApproxEqual(v1[0], v2[0]) && ApproxEqual(v1[1], v2[1]) && ApproxEqual(v1[2], v2[2]);
}

// Placeholder (on host) for all surfaces belonging to a volume. An array of those will be indexed
// by the logical volume id. Also an intermediate helper for building portals.
// Note: the local surfaces defined by solids will have local references that will be changed by
// the flattening process, depending on the scene on which the parent volume will be flattened
struct VolumeShell {
  std::vector<int>         fSurfaces; ///< Local surface id's for this volume
};

// We would need masks, imprints, portals created contiguously, by some factory
// Host helper for filling the SurfData store
template <typename Real_t>
class BrepHelper {
  using SurfData_t     = SurfData<Real_t>;
  using CylData_t      = CylData<Real_t>;
  using ConeData_t     = ConeData<Real_t>;
  using SphData_t      = SphData<Real_t>;
  using Transformation = vecgeom::Transformation3D;
  using Vector         = vecgeom::Vector3D<Real_t>;
 
private:
  SurfData_t *fSurfData{nullptr};            ///< Surface data
  SurfData_t *fSurfDataGPU{nullptr};         ///< Surface data on device

  std::vector<RangeMask<Real_t>> fRanges;    ///< List of frame ranges
  std::vector<Frame>     fFrames;        ///< vector of masks
  std::vector<CylData_t>   fCylSphData;      ///< data for cyl surfaces
  std::vector<ConeData_t>  fConeData;        ///< data for conical surfaces
  std::vector<Transformation> fLocalTrans;   ///< local transformations
  std::vector<Transformation> fGlobalTrans;  ///< global transformations for surfaces in the scene
  std::vector<PlacedSurface> fLocalSurfaces; ///< local surfaces
  std::vector<PlacedSurface> fGlobalSurfaces; ///< global surfaces
  std::vector<CommonSurface> fCommonSurfaces; ///< common surfaces
  std::vector<VolumeShell>   fShells;        ///< vector of local volume surfaces
  std::vector<std::vector<int>> fCandidates; ///< candidate lists for each state

  BrepHelper() : fSurfData(new SurfData_t()) {}

public:
  
  /// Returns the singleton instance (CPU only)
  static BrepHelper &Instance()
  {
    static BrepHelper instance;
    return instance;
  }

  void ClearData()
  {
    // Dispose of surface data and shrink the container
    fRanges.clear(); std::vector<RangeMask<Real_t>>().swap(fRanges);
    fFrames.clear(); std::vector<Frame>().swap(fFrames);
    fCylSphData.clear(); std::vector<CylData_t>().swap(fCylSphData);
    fConeData.clear(); std::vector<ConeData_t>().swap(fConeData);
    fLocalTrans.clear(); std::vector<Transformation>().swap(fLocalTrans);
    fGlobalTrans.clear(); std::vector<Transformation>().swap(fGlobalTrans);
    fLocalSurfaces.clear(); std::vector<PlacedSurface>().swap(fLocalSurfaces);
    fGlobalSurfaces.clear(); std::vector<PlacedSurface>().swap(fGlobalSurfaces);
    fCommonSurfaces.clear(); std::vector<CommonSurface>().swap(fCommonSurfaces);
    for (size_t i = 0; i < fShells.size(); ++i) {
      fShells[i].fSurfaces.clear(); std::vector<int>().swap(fShells[i].fSurfaces);
    }
    fShells.clear(); std::vector<VolumeShell>().swap(fShells);

    delete fSurfData;
    fSurfData = nullptr;
    // cudaDelete(fSurfDataGPU);
    fSurfDataGPU = nullptr;
  }

  ~BrepHelper()
  {
    delete fSurfData;
  }

  void SortSides(int shared_id)
  {
    // lambda to remove a surface from a side
    auto removeSurface = [&](Side &side, int ind)
    {
      for (int i = ind + 1; i < side.fNsurf; ++i)
        side.fSurfaces[i-1] = side.fSurfaces[i];
      if (ind < side.fNsurf) side.fNsurf--;
    };

    // lambda to detect surfaces on a side that have identical frame
    auto sortAndRemoveCommonFrames = [&](Side &side)
    {
      if (!side.fNsurf) return;
      std::sort(side.fSurfaces, side.fSurfaces + side.fNsurf,
                [&](int i, int j) {return fGlobalSurfaces[i] < fGlobalSurfaces[j];});
      for (int i = 0; i < side.fNsurf - 1; ++i) {
        for (int j = side.fNsurf - 1; j > i; --j) {
          if (EqualFrames(side, i, j)) removeSurface(side, j);
        }
      }
    };

    sortAndRemoveCommonFrames(fCommonSurfaces[shared_id].fLeftSide);
    sortAndRemoveCommonFrames(fCommonSurfaces[shared_id].fRightSide);
  }

  void ComputeDefaultStates(int shared_id)
  {
    using vecgeom::NavStateIndex;
    // Computes the default states for each side of a common surface
    Side &left  = fCommonSurfaces[shared_id].fLeftSide;
    Side &right = fCommonSurfaces[shared_id].fRightSide;
    assert(left.fNsurf > 0 || right.fNsurf > 0);

    NavIndex_t default_ind;

    // A lambda that finds the deepest common ancestor between 2 states
    auto getCommonState = [&](NavIndex_t const &s1, NavIndex_t const &s2)
    {
      NavIndex_t a1(s1), a2(s2);
      // Bring both states at the same level
      while (NavStateIndex::GetLevelImpl(a1) > NavStateIndex::GetLevelImpl(a2))
        a1 = NavStateIndex::PopImpl(a1);
      while (NavStateIndex::GetLevelImpl(a2) > NavStateIndex::GetLevelImpl(a1))
        a2 = NavStateIndex::PopImpl(a2);

      // Pop until we reach the same state
      while (a1 != a2) {
        a1 = NavStateIndex::PopImpl(a1);
        a2 = NavStateIndex::PopImpl(a2);
      }
      return a1;
    };

    int minlevel = 10000; // this is a big-enough number as level
    for (int isurf = 0; isurf < left.fNsurf; ++isurf) {
      auto navind = fGlobalSurfaces[left.fSurfaces[isurf]].fState;
      minlevel = std::min(minlevel, (int)NavStateIndex::GetLevelImpl(navind));
    }
    for (int isurf = 0; isurf < right.fNsurf; ++isurf) {
      auto navind = fGlobalSurfaces[right.fSurfaces[isurf]].fState;
      minlevel = std::min(minlevel, (int)NavStateIndex::GetLevelImpl(navind));
    }

    // initialize the default state
    if (left.fNsurf > 0)
      default_ind = fGlobalSurfaces[left.fSurfaces[0]].fState;
    else if (right.fNsurf > 0)
      default_ind = fGlobalSurfaces[right.fSurfaces[0]].fState;

    for (int isurf = 0; isurf < left.fNsurf; ++isurf)
      default_ind = getCommonState(default_ind, fGlobalSurfaces[left.fSurfaces[isurf]].fState);
    for (int isurf = 0; isurf < right.fNsurf; ++isurf)
      default_ind = getCommonState(default_ind, fGlobalSurfaces[right.fSurfaces[isurf]].fState);
    
    if (NavStateIndex::GetLevelImpl(default_ind) == minlevel)
      default_ind = NavStateIndex::PopImpl(default_ind);
    fCommonSurfaces[shared_id].fDefaultState = default_ind;
  }

  void PrintCommonSurface(int id)
  {
    auto const &surf = fCommonSurfaces[id];
    printf("== common surface %d: default state: ", id);
    vecgeom::NavStateIndex default_state(surf.fDefaultState);
    default_state.Print();
    printf("   side1: %d surfaces\n", surf.fLeftSide.fNsurf);
    for (int i = 0; i < surf.fLeftSide.fNsurf; ++i) {
      int idglob = surf.fLeftSide.fSurfaces[i];
      auto const &placed = fGlobalSurfaces[idglob];
      printf("    surf %d: trans: ", idglob);
      fGlobalTrans[placed.fTrans].Print();
      printf(", ");
      vecgeom::NavStateIndex state(placed.fState);
      state.Print();
    }
    printf("   side2: %d surfaces\n", surf.fRightSide.fNsurf);
    for (int i = 0; i < surf.fRightSide.fNsurf; ++i) {
      int idglob = surf.fRightSide.fSurfaces[i];
      auto const &placed = fGlobalSurfaces[idglob];
      printf("    surf %d: trans: ", idglob);
      fGlobalTrans[placed.fTrans].Print();
      vecgeom::NavStateIndex state(placed.fState);
      state.Print();
    }
  }

  void SetNvolumes(int nvolumes)
  {
    if (fShells.size() > 0) {
      std::cout << "BrepHelper::SetNvolumes already called for this instance.\n";
      return; 
    }
    fShells.resize(nvolumes);
  }

  bool CreateLocalSurfaces()
  {
    // Iterate logical volumes and create local surfaces
    std::vector<vecgeom::LogicalVolume*> volumes;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(volumes);
    SetNvolumes(volumes.size());
    // creat a placeholder for surface data
    for (auto volume : volumes) {
      vecgeom::VUnplacedVolume const *solid = volume->GetUnplacedVolume();
      vecgeom::UnplacedBox const *box = dynamic_cast<vecgeom::UnplacedBox const*>(solid);
      if (!box) {
        std::cout << "testEm3: solid type not supported for volume: " << volume->GetName() << "\n";
        return false;
      }
      CreateBoxSurfaces(*box, volume->id());
    }
    return true;
  }

  bool CreateCommonSurfacesFlatTop()
  {
    // Iterate the geometry tree and flatten surfaces at top level
    int nphysical = 0;
    vecgeom::NavStateIndex state;
    
    // recursive geometry visitor lambda creating the common surfaces for the current placed volume
    typedef std::function<void(vecgeom::VPlacedVolume const *)> func_t;
    func_t createCommonSurfaces = [&](vecgeom::VPlacedVolume const *pvol) {
      state.Push(pvol);
      const auto vol   = pvol->GetLogicalVolume();
      auto daughters   = vol->GetDaughters();
      int nd           = daughters.size();
      nphysical++;
      Transformation trans;
      state.TopMatrix(trans);
      VolumeShell const &shell = fShells[vol->id()];
      for (int lsurf_id : shell.fSurfaces) {
        PlacedSurface const &lsurf = fLocalSurfaces[lsurf_id];
        Transformation global(trans);
        global.MultiplyFromRight(fLocalTrans[lsurf.fTrans]);
        int trans_id = fGlobalTrans.size();
        fGlobalTrans.push_back(global);
        // Create the global surface
        int id_glob = fGlobalSurfaces.size();
        fGlobalSurfaces.push_back({lsurf.fSurface, lsurf.fFrame, trans_id, state.GetNavIndex()});
        CreateCommonSurface(id_glob);
      }

      // Now do the daughters
      for (int id = 0; id < nd; ++id) {
        createCommonSurfaces(daughters[id]);
      }
      state.Pop();
    };

    // add a dummy common surface since index 0 is not allowed for correctly handling sides
    fCommonSurfaces.push_back({});

    createCommonSurfaces(vecgeom::GeoManager::Instance().GetWorld());

    for (size_t isurf = 1; isurf < fCommonSurfaces.size(); ++isurf) {
      ComputeDefaultStates(isurf);
      SortSides(isurf);
      PrintCommonSurface(isurf);
    }

    std::cout << "Visited " << nphysical << " physical volumes, created " << fCommonSurfaces.size() - 1 << " common surfaces\n";

    CreateCandidateLists();
    return true;
  }

  ///< This method creates helper lists of candidate surfaces for each navigation state
  void CreateCandidateLists()
  {
    int numNodes = vecgeom::GeoManager::Instance().GetTotalNodeCount() + 1; // count also outside state
    fCandidates.reserve(numNodes);

    // Lambda adding the surface id as candidate to all states from a side
    auto addSurfToSideStates = [&](int isurf, int iside)
    {
      Side const &side = (iside > 0) ? fCommonSurfaces[isurf].fLeftSide
                                     : fCommonSurfaces[isurf].fRightSide;
      for (int i = 0; i < side.fNsurf; ++i) {
        int idglob = side.fSurfaces[i];
        auto const &placed = fGlobalSurfaces[idglob];
        vecgeom::NavStateIndex state(placed.fState);
        int state_id = state.GetId();
        fCandidates[state_id].push_back(isurf * iside);
      }
    };

    // prepare all lists
    for (int i = 0; i < numNodes; ++i) fCandidates.push_back({});
  
    // loop over all common surfaces and add their index in the appropriate list
    for (size_t isurf = 1; isurf < fCommonSurfaces.size(); ++isurf) {
      auto const &surf = fCommonSurfaces[isurf];
      // Add to default surface state
      vecgeom::NavStateIndex state(surf.fDefaultState);
      fCandidates[state.GetId()].push_back(-isurf);
      // Add to side states
      addSurfToSideStates(isurf, 1);
      addSurfToSideStates(isurf, -1);
    }

    vecgeom::NavStateIndex state;

    // recursive geometry visitor lambda printing the candidates lists
    // We have no direct access from a state (contiguous) id to the actual state index
    typedef std::function<void(vecgeom::VPlacedVolume const *)> func_t;
    func_t printCandidates = [&](vecgeom::VPlacedVolume const *pvol) {
      state.Push(pvol);
      const auto vol   = pvol->GetLogicalVolume();
      auto daughters   = vol->GetDaughters();
      int nd           = daughters.size();
      state.Print();
      printf(" %lu candidates: ", fCandidates[state.GetId()].size());
      for (auto isurf : fCandidates[state.GetId()]) printf("%d ", isurf);
      printf("\n");

      // do daughters
      for (int id = 0; id < nd; ++id) {
        printCandidates(daughters[id]);
      }
      state.Pop();
    };

    printf("\nCandidate surfaces per state:");
    state.Print();
    printf(" %lu candidates: ", fCandidates[state.GetId()].size());
    for (auto isurf : fCandidates[state.GetId()]) printf("%d ", isurf);
    printf("\n");

    printCandidates(vecgeom::GeoManager::Instance().GetWorld());
  }

private:
  int AddSurfaceToShell(int logical_id, int isurf)
  {
    if (fShells.size() == 0) {
      std::cout << "BrepHelper::AddSurfaceToShell: need to call SetNvolumes first\n";
      return -1;
    }
    assert(logical_id < (int)fShells.size() && "surface shell id exceeding number of volumes");
    int id = fShells[logical_id].fSurfaces.size();
    fShells[logical_id].fSurfaces.push_back(isurf);
    return id;
  }

  UnplacedSurface CreateUnplacedSurface(SurfaceType type, Real_t *data = nullptr)
  {
    switch (type) {
      case kPlanar:
        return UnplacedSurface(type, -1);
      case kCylindrical:
      case kSpherical:
        fCylSphData.push_back({data[0]});
        return UnplacedSurface(type, fCylSphData.size() - 1);
      case kConical:
        fConeData.push_back({data[0], data[1]});
        return UnplacedSurface(type, fConeData.size() - 1);
      case kTorus:
      case kGenSecondOrder:
        std::cout << "kTorus, kGenSecondOrder unhandled\n";
        return UnplacedSurface(type, -1);
    };
    return UnplacedSurface(type, -1);
  }

  Frame CreateFrame(FrameType type, RangeMask<Real_t> mask)
  {
    int id = fRanges.size();
    fRanges.push_back(mask);
    return Frame(type, id);
  }

  int CreateLocalTransformation(Transformation const &trans)
  {
    int id = fLocalTrans.size();
    fLocalTrans.push_back(trans);
    return id;
  }

  int CreateLocalSurface(UnplacedSurface const &unplaced, Frame const &frame, int trans)
  {
    int id = fLocalSurfaces.size();
    fLocalSurfaces.push_back({unplaced, frame, trans});
    return id;
  }

  int CreateCommonSurface(int idglob)
  {
    bool flip;

    auto approxEqual = [&](int idglob1, int idglob2)
    {
      flip = false;
      PlacedSurface const &s1 = fGlobalSurfaces[idglob1];
      PlacedSurface const &s2 = fGlobalSurfaces[idglob2];
      // Surfaces may be in future "compatible" even if they are not the same, for now enforce equality
      if (s1.fSurface.type != s2.fSurface.type)
        return false;

      // Check if the 2 surfaces are parallel
      Transformation const &t1 = fGlobalTrans[s1.fTrans];
      Transformation const &t2 = fGlobalTrans[s2.fTrans];
      // Calculate normalized connection vector between the two transformations
      // Use double precision explicitly
      vecgeom::Vector3D<double> tdiff = t1.Translation() - t2.Translation();
      bool same_tr = ApproxEqual(tdiff, {0, 0, 0});
      vecgeom::Vector3D<double> ldir;
      switch (s1.fSurface.type) {
        case kPlanar:
          if (same_tr) break;
          // For planes to match, the connecting vector must be along the planes
          tdiff.Normalize();
          t1.TransformDirection(tdiff, ldir);
          if (std::abs(ldir[2]) > vecgeom::kTolerance)
            return false;
          break;
        case kCylindrical:
          if (std::abs( fCylSphData[s1.fSurface.id].radius - fCylSphData[s1.fSurface.id].radius) > vecgeom::kTolerance)
            return false;
          if (same_tr) break;
          tdiff.Normalize();
          t1.TransformDirection(tdiff, ldir);
          // For connected cylinders, the connecting vector must be along the Z axis
          if (!ApproxEqual(ldir, {0, 0, ldir[2]}))
            return false;
          break;
        case kConical:
        case kSpherical:
        case kTorus:
        case kGenSecondOrder:
        default:
          printf("CreateCommonSurface: case not implemented\n");
          return false;
      };

      // Now check if the rotations are matching. The z axis transformed
      // with the two rotations should end up in aligned vectors. This is
      // true for planes (Z is the normal) but also for tubes/cones where
      // Z is the axis of symmetry
      vecgeom::Vector3D<double> zaxis(0, 0, 1);
      auto z1 = t1.InverseTransformDirection(zaxis);
      auto z2 = t2.InverseTransformDirection(zaxis);
      if (!ApproxEqual(z1.Cross(z2), {0, 0, 0}))
        return false;
      
      flip = z1.Dot(z2) < 0;
      return true;
    };
    
    // this may be slow
    auto it = std::find_if(std::begin(fCommonSurfaces), std::end(fCommonSurfaces), [&](const CommonSurface &t)
                          {return (t.fLeftSide.fNsurf > 0) ? approxEqual(t.fLeftSide.fSurfaces[0], idglob) : false;});
    int id = -1;
    if (it != std::end(fCommonSurfaces)) {
      id = int(it - std::begin(fCommonSurfaces));
      // Add the global surface to the appropriate side
      if (flip)
        (*it).fRightSide.AddSurface(idglob);
      else
        (*it).fLeftSide.AddSurface(idglob);

    } else {
      // Construct a new common surface from the current placed global surface
      id = fCommonSurfaces.size();
      fCommonSurfaces.push_back({fGlobalSurfaces[idglob].fSurface.type, idglob});
    }
    return id;
  }

  // The code for creating solid-specific surfaces should sit in the specific solid struct type
  void CreateBoxSurfaces(vecgeom::UnplacedBox const &box, int logical_id)
  {
    int isurf;
    // surface at -dx:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateFrame(kWindow, {box.z(), box.y()}),
                                CreateLocalTransformation({-box.x(), 0, 0, 90, -90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at +dx:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateFrame(kWindow, {box.z(), box.y()}),
                                CreateLocalTransformation({box.x(), 0, 0, 90, 90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at -dy:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateFrame(kWindow, {box.x(), box.z()}),
                                CreateLocalTransformation({0, -box.y(), 0, 0, 90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at +dy:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateFrame(kWindow, {box.x(), box.z()}),
                                CreateLocalTransformation({0, box.y(), 0, 0, -90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at -dz:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateFrame(kWindow, {box.x(), box.y()}),
                                CreateLocalTransformation({0, 0, -box.z(), 0, 180, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at +dz:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateFrame(kWindow, {box.x(), box.y()}),
                                CreateLocalTransformation({0, 0, box.z(), 0, 0, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
  }

  ///< This function evaluates if the frames of two placed surfaces on the same side
  ///< of a common surface are matching
  bool EqualFrames(Side const &side, int i1, int i2)
  {
    PlacedSurface const &s1 = fGlobalSurfaces[side.fSurfaces[i1]];
    PlacedSurface const &s2 = fGlobalSurfaces[side.fSurfaces[i2]];
    if (s1.fFrame.type != s2.fFrame.type)
      return false;
    // Get displacement vector between the 2 frame centers and check if it has null length
    Transformation const &t1 = fGlobalTrans[s1.fTrans];
    Transformation const &t2 = fGlobalTrans[s2.fTrans];
    vecgeom::Vector3D<double> tdiff = t1.Translation() - t2.Translation();
    if (!ApproxEqual(tdiff, {0, 0, 0}))
      return false;
    
    auto frameData1 = fRanges[s1.fFrame.id];
    auto frameData2 = fRanges[s2.fFrame.id];
    vecgeom::Vector3D<double> v1(frameData1.range[0], frameData1.range[1], 0);
    vecgeom::Vector3D<double> v1_transposed(frameData1.range[1], frameData1.range[0], 0);
    vecgeom::Vector3D<double> v2(frameData2.range[0], frameData2.range[1], 0);
    
    // Different treatment of different frame types
    switch (s1.fFrame.type) {
      case kRangeZ:
      case kRangeCyl:
      case kRangeSph:
        if (ApproxEqual(v1, v2))
          return true;
        return false;
      case kWindow:
        // Here we can have matching also if the window is rotated
        if (!ApproxEqual(v1, v2) && !ApproxEqual(v1_transposed, v2))
          return false;
        {
          // The vectors connecting the local window center to the closest edge
          // must be aligned in the transformed frames
          vecgeom::Vector3D<double> vshort;
          if (v1[0] > v1[1]) vshort.Set(0, v1[1], 0);
          else               vshort.Set(v1[0], 0, 0);       
          auto v1glob = t1.InverseTransform(vshort);
          if (v2[0] > v2[1]) vshort.Set(0, v2[1], 0);
          else               vshort.Set(v2[0], 0, 0);       
          auto v2glob = t2.InverseTransform(vshort);
          if (ApproxEqual(v1glob.Cross(v2glob), {0, 0, 0}))
            return true;
        } 
        break;
      case kTriangle:
        // to be implemented
        break;
    };
    return false; 
  }

  ///< The method updates the SurfData storage
  void UpdateSurfData()
  {
    // Create and copy surface data
    fSurfData->fCylSphData = new CylData_t[fCylSphData.size()];
    for (size_t i = 0; i < fCylSphData.size(); ++i)
      fSurfData->fCylSphData[i] = fCylSphData[i];

    fSurfData->fConeData = new ConeData_t[fConeData.size()];
    for (size_t i = 0; i < fConeData.size(); ++i)
      fSurfData->fConeData[i] = fConeData[i];

    fSurfData->fConeData = new ConeData_t[fConeData.size()];
    for (size_t i = 0; i < fConeData.size(); ++i)
      fSurfData->fConeData[i] = fConeData[i];
  }
  
};

} // namespace vgbrep
#endif

