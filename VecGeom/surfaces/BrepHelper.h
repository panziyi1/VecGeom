#ifndef VECGEOM_SURFACE_BREPHELPER_H_
#define VECGEOM_SURFACE_BREPHELPER_H_

#include <cassert>
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

  std::vector<RangeMask<Real_t>> fRanges;    ///< List of outline ranges
  std::vector<Outline>     fOutlines;        ///< vector of masks
  std::vector<CylData_t>   fCylSphData;      ///< data for cyl surfaces
  std::vector<ConeData_t>  fConeData;        ///< data for conical surfaces
  std::vector<Transformation> fLocalTrans;   ///< local transformations
  std::vector<Transformation> fGlobalTrans;  ///< global transformations for surfaces in the scene
  std::vector<PlacedSurface> fLocalSurfaces; ///< local surfaces
  std::vector<PlacedSurface> fGlobalSurfaces; ///< global surfaces
  std::vector<VolumeShell>   fShells;        ///< vector of local volume surfaces

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
    fOutlines.clear(); std::vector<Outline>().swap(fOutlines);
    fCylSphData.clear(); std::vector<CylData_t>().swap(fCylSphData);
    fConeData.clear(); std::vector<ConeData_t>().swap(fConeData);
    fLocalTrans.clear(); std::vector<Transformation>().swap(fLocalTrans);
    fGlobalTrans.clear(); std::vector<Transformation>().swap(fGlobalTrans);
    fLocalSurfaces.clear(); std::vector<PlacedSurface>().swap(fLocalSurfaces);
    fGlobalSurfaces.clear(); std::vector<PlacedSurface>().swap(fGlobalSurfaces);
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

  bool CreateGlobalSurfacesFlatTop()
  {
    // Iterate the geometry tree and flatten surfaces at top level
    int nphysical = 0;
    int nglobalsurf = 0;
    vecgeom::NavStateIndex state;
    
    // recursive geometry visitor lambda creating the global surfaces for the current placed volume
    typedef std::function<void(vecgeom::VPlacedVolume const *)> func_t;
    func_t createGlobalSurfaces = [&](vecgeom::VPlacedVolume const *pvol) {
      state.Push(pvol);
      const auto vol   = pvol->GetLogicalVolume();
      auto daughters   = vol->GetDaughters();
      int nd           = daughters.size();
      nphysical++;
      //state.Print();
      Transformation trans;
      state.TopMatrix(trans);
      VolumeShell const &shell = fShells[vol->id()];
      for (int lsurf_id : shell.fSurfaces) {
        PlacedSurface const &lsurf = fLocalSurfaces[lsurf_id];
        Transformation global(trans);
        global.MultiplyFromRight(fLocalTrans[lsurf.fTrans]);
        state.Print();
        global.Print();
        std::cout << "\n";
        int global_trans_id = CreateGlobalTransformation(global);
        CreateGlobalSurface(lsurf.fSurface, lsurf.fOutline, global_trans_id, state.GetNavIndex());
        nglobalsurf++;
      }

      // Now do the daughters
      for (int id = 0; id < nd; ++id) {
        createGlobalSurfaces(daughters[id]);
      }
      state.Pop();
    };

    createGlobalSurfaces(vecgeom::GeoManager::Instance().GetWorld());

    std::cout << "Visited " << nphysical << " physical volumes, created " << nglobalsurf << " global surfaces\n";
    return true;
  }

private:
  int AddSurfaceToShell(int logical_id, int isurf)
  {
    if (fShells.size() == 0) {
      std::cout << "BrepHelper::AddSurfaceToShell: need to call SetNvolumes first\n";
      return -1;
    }
    assert(logical_id < fShells.size() && "surface shell id exceeding number of volumes");
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

  Outline CreateOutline(OutlineType type, RangeMask<Real_t> mask)
  {
    int id = fRanges.size();
    fRanges.push_back(mask);
    return Outline(type, id);
  }

  int CreateLocalTransformation(Transformation const &trans)
  {
    int id = fLocalTrans.size();
    fLocalTrans.push_back(trans);
    return id;
  }

  int CreateGlobalTransformation(Transformation const &trans)
  {
    int id = fGlobalTrans.size();
    fGlobalTrans.push_back(trans);
    return id;
  }

  int CreateLocalSurface(UnplacedSurface const &unplaced, Outline const &outline, int trans)
  {
    int id = fLocalSurfaces.size();
    fLocalSurfaces.push_back({unplaced, outline, trans});
    return id;
  }

  int CreateGlobalSurface(UnplacedSurface const &unplaced, Outline const &outline, int trans, NavIndex_t navindex)
  {
    int id = fGlobalSurfaces.size();
    fGlobalSurfaces.push_back({unplaced, outline, trans, navindex});
    return id;
  }

  // The code for creating solid-specific surfaces should sit in the specific solid struct type
  void CreateBoxSurfaces(vecgeom::UnplacedBox const &box, int logical_id)
  {
    int isurf;
    // surface at -dx:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateOutline(kWindow, {box.z(), box.y()}),
                                CreateLocalTransformation({-box.x(), 0, 0, 0, 90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at +dx:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateOutline(kWindow, {box.z(), box.y()}),
                                CreateLocalTransformation({box.x(), 0, 0, 0, -90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at -dy:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateOutline(kWindow, {box.x(), box.z()}),
                                CreateLocalTransformation({0, -box.y(), 0, 90, -90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at +dy:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateOutline(kWindow, {box.x(), box.z()}),
                                CreateLocalTransformation({0, box.y(), 0, 90, 90, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at -dz:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateOutline(kWindow, {box.x(), box.y()}),
                                CreateLocalTransformation({0, 0, -box.z(), 0, 180, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
    // surface at +dz:
    isurf = CreateLocalSurface( CreateUnplacedSurface(kPlanar),
                                CreateOutline(kWindow, {box.x(), box.y()}),
                                CreateLocalTransformation({0, 0, box.z(), 0, 0, 0})
                              );
    AddSurfaceToShell(logical_id, isurf);
  }
  
};

} // namespace vgbrep
#endif

