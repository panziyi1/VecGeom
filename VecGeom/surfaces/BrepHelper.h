#ifndef VECGEOM_SURFACE_BREPHELPER_H_
#define VECGEOM_SURFACE_BREPHELPER_H_

#include <cassert>
#include <functional>
#include <VecGeom/surfaces/Model.h>

#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/management/GeoManager.h"


namespace vgbrep {

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

  bool ApproxEqual(Real_t t1, Real_t t2)
  {
    return std::abs(t1 - t2) <= vecgeom::kTolerance;
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
      auto tdiff = t1.Translation() - t2.Translation();
      tdiff.Normalize();
      decltype(tdiff) ldir;
      t1.TransformDirection(tdiff, ldir);
      switch (s1.fSurface.type) {
        case kPlanar:
          // For planes to match, the connecting vector must be along the planes
          if (std::abs(ldir[2]) > vecgeom::kTolerance)
            return false;
          break;
        case kCylindrical:
          if (std::abs( fCylSphData[s1.fSurface.id].radius - fCylSphData[s1.fSurface.id].radius) > vecgeom::kTolerance)
            return false;
          // For connected cylinders, the connecting vector must be along the Z axis
          if (!ApproxEqual(ldir[0], 0) || !ApproxEqual(ldir[1], 0))
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

      // Now check if the rotations are matching
      // Check if translations are equal (but this is just a sub-case)
      //for (int i = 0; i < 3; ++i)
      //  if (!ApproxEqual(t1.Translation(i), t2.Translation(i)))
      //    return false;

      // Check rotation. Two rotations are considered equal if we can reach one from
      // the other by composing with a flip operation:
      //   R2 = R1 * F    =>    R1.Inverse * R2 = F
      //   R1 = R2 * F    =>    R2.Inverse * R1 = F
      // so: R1.Inverse * R2 = R2.Inverse * R1
      // which has a trivial solution (F = I) for R1 = R2 and a non-trivial one.
      bool equal = true;
      for (int i = 0; i < 9; ++i) {
        if (!ApproxEqual(t1.Rotation(i), t2.Rotation(i))) {
          equal = false;
          break;
        }
      }
      if (equal) return true;
      // We need to check the flipped solution
      Transformation flip1, flip2;
      t1.Inverse(flip1);
      t2.Inverse(flip2);
      flip1.MultiplyFromRight(t2);
      flip2.MultiplyFromRight(t1);
      
      for (int i = 0; i < 9; ++i)
        if (!ApproxEqual(flip1.Rotation(i), flip2.Rotation(i)))
          return false;
      
      flip = true;
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
  
};

} // namespace vgbrep
#endif

