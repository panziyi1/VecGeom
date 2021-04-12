#include <stdio.h>

#include "VecGeom/base/Config.h"
#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/management/CudaManager.h"
#endif
#include "VecGeom/navigation/GlobalLocator.h"
#include "VecGeom/navigation/NavigationState.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/Box.h"
#include "TGeoManager.h"
#ifdef VECGEOM_GDML
  #include "Frontend.h"
#endif

using namespace vecgeom;

int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cerr << std::endl;
    std::cerr << "Need to give gdmlfile and volumename" << std::endl;
    std::cerr << "USAGE : test/check-parser [gdmlfile] [VolumeName] [--parser=ROOT(default:VecGeom)] "
              << std::endl;
    std::cerr << "  ex) test/check-parser cms2018.gdml MUON --vgdml (default)" << std::endl;
    std::cerr << "      test/check-parser cms2018.gdml MUON --root" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  // flag to choose GDML parser
  bool parseRoot = false;
  for (auto i = 3; i < argc; ++i)
  {
    if (!strcmp(argv[i], "--vgdml")) parseRoot = false;
    if (!strcmp(argv[i], "--root")) parseRoot = true;
  }

  if (parseRoot) {
#ifdef VECGEOM_ROOT
    // Use Root parser
    std::cout<<" *** Using ROOT GDML parser...\n";
    TGeoManager::Import(argv[1]);
    // convert current gGeoManager to a VecGeom geometry
    RootGeoManager::Instance().LoadRootGeometry();
#else
    std::cerr<<" *** ERROR: Root GDML parser requested in a ROOT=OFF build configuration.\n";
    return 1;
#endif
  }
  else {
#ifdef VECGEOM_GDML
    // Use VecGeom parser
    std::cout<<" *** Using VGDML parser...\n";
    constexpr bool validate_xml_schema = false;
    vgdml::Frontend::Load(argv[1], validate_xml_schema, 1.0);
#else
    std::cerr<<" *** ERROR: VGDML parser requested in a GDML=OFF build configuration.\n";
    return 1;
#endif
  }

  std::string testvolume(argv[2]);

  Vector3D<Precision> point(-6954, 10, 8562.5);
  std::cout<<" --- Locating point "<< point << std::endl;

  NavigationState *path = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  GlobalLocator::LocateGlobalPoint(GeoManager::Instance().GetWorld(), point, *path, true);
  path->Print();

  //GeoManager::Instance().FindLogicalVolume("YE2p_b0x7f4a8eed33c0");
  //GeoManager::Instance().FindPlacedVolume("YE2p_b0x7f4a8eed33c0");
  for (unsigned int i = 0; i <= path->GetLevel(); ++i) {
    auto pvol = path->At(i);
    pvol->Print(3*i);
  }

  NavigationState::ReleaseInstance(path);

  return 0;
}
