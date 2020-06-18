#include "VecGeom/volumes/Box.h"
#include "VecGeom/base/Transformation3D.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/management/GeoManager.h"
#include "VecGeom/navigation/VLevelLocator.h"

#ifdef VECGEOM_ROOT
#include "TGeoManager.h"
#endif

namespace vecgeom {

// Purpose: parse a .gdml file into VecGeom.
//   If VECGEOM_GEANT4 is defined, loads also a Geant4 geometry.
//   If VECGEOM_ROOT is defined, loads also a Root geometry.
//
// input: a .gdml filename char*
//
bool SetupGDMLGeometry(const char* fname) {

  bool loaded = false;
  std::string gdmlName(fname);
  if ( gdmlName.find("gdml") >= gdmlName.size() ) {
    gdmlName += ".gdml";
  }

#ifdef VECGEOM_GDML
  std::cout<<"*** Using native GDML parser...\n";
  loaded = vgdml::Frontend::Load(gdmlName);
#endif

#ifdef VECGEOM_ROOT
  if (loaded) {
    std::cout<<" *** Exporting VecGeom geometry to .root format...\n";
    RootGeoManager::Instance().ExportToROOTGeometry(GeoManager::Instance().GetWorld(), "navBench.root");
    std::cout<<" *** Cleaning up RootGeoManager...\n";
    RootGeoManager::Instance().Clear();

    // read it back to fill Root geometry
    std::cout<<" *** Reading back .root file to fill TGeoManager geometry...\n";
    RootGeoManager::Instance().set_verbose(1);
    RootGeoManager::Instance().LoadRootGeometry("navBench.root");
  }
  else {
    std::cout<<"*** Using Root's GDML parser...\n";
    loaded = TGeoManager::Import(gdmlName.c_str());
    RootGeoManager::Instance().LoadRootGeometry();
  }

  #ifdef VECGEOM_GEANT4
  // load Geant4 geometry from ROOT
  G4GeoManager::Instance().GetG4GeometryFromROOT();
  #endif

#else
#ifdef VECGEOM_GEANT4xb
  std::cout<<"*** Using G4 GDML parser...\n";
  G4GeoManager::Instance().LoadG4Geometry(gdmlName.c_str(), true);
#endif

#endif

  const VPlacedVolume *w = RootGeoManager::Instance().world();
  GeoManager::Instance().SetWorld(w);
  GeoManager::Instance().CloseGeometry();

  return loaded;
}

/*
bool SetupGDMLGeometry()
{
#ifdef VERBOSE
  //.. print geometry details
  for (auto &element : GeoManager::Instance().GetLogicalVolumesMap()) {
    auto *lvol = element.second;
    // lvol->SetLevelLocator(nullptr); // force-disable locators to test default GlobalLocator implementations
    std::cerr << "SetupBoxGeom(): logVol=" << lvol << ", name=" << lvol->GetName()
	      << ", locator=" << (lvol->GetLevelLocator() ? lvol->GetLevelLocator()->GetName() : "NULL") << "\n";
  }

  std::vector<VPlacedVolume *> v1;
  GeoManager::Instance().getAllPlacedVolumes(v1);
  for (auto &plvol : v1) {
    std::cerr << "placedVol=" << plvol << ", name=" << plvol->GetName() << ", world=" << world << ", <"
	      << world->GetName() << ", " << GeoManager::Instance().GetWorld() << ">\n";
  }

  //.. more details
  world->PrintContent();
  std::cerr << "\n";
#endif

#ifdef VECGEOM_ROOT
  // Exporting to ROOT file
  RootGeoManager::Instance().ExportToROOTGeometry(GeoManager::Instance().GetWorld(), "gdml2root.root");
  RootGeoManager::Instance().Clear();
#endif

  return true;
}
*/

} // namespace vecgeom
