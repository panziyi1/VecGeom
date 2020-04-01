/*
 *  File: NavigationBenchmark.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: swenzel, lima
 */

#undef VERBOSE
#define VERBOSE

//#include "VecGeom/base/Config.h"
#ifdef VECGEOM_GDML
#include "persistency/gdml/source/include/Frontend.h"
#include "persistency/gdml/source/include/Backend.h"
#endif

#ifdef VECGEOM_ROOT
#include "TGeoManager.h"
#include "VecGeom/management/RootGeoManager.h"
#include "utilities/Visualizer.h"
#endif

#ifdef VECGEOM_GEANT4
#include "VecGeom/management/G4GeoManager.h"
#include "G4ThreeVector.hh"
// #include "G4TouchableHistoryHandle.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#endif

#include "VecGeom/benchmarking/NavigationBenchmarker.h"
#include "test/benchmark/ArgParser.h"
#include "VecGeom/volumes/utilities/VolumeUtilities.h"

#include "VecGeom/management/GeoManager.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/volumes/Orb.h"
#include "VecGeom/volumes/Trapezoid.h"

using namespace VECGEOM_NAMESPACE;

VPlacedVolume *SetupGeometry()
{
  std::cout<<"\nBuilding default navBench geometry...\n";
  UnplacedBox *worldUnplaced      = new UnplacedBox(10, 10, 10);
  UnplacedTrapezoid *trapUnplaced = new UnplacedTrapezoid(4, 0, 0, 4, 4, 4, 0, 4, 4, 4, 0);
  UnplacedBox *boxUnplaced        = new UnplacedBox(2, 2, 2);
  UnplacedOrb *orbUnplaced        = new UnplacedOrb(2.8);

  LogicalVolume *world = new LogicalVolume("world", worldUnplaced);
  LogicalVolume *trap  = new LogicalVolume("trap", trapUnplaced);
  LogicalVolume *box   = new LogicalVolume("box", boxUnplaced);
  LogicalVolume *orb   = new LogicalVolume("orb", orbUnplaced);

  Transformation3D *ident = new Transformation3D(0, 0, 0, 0, 0, 0);
  orb->PlaceDaughter("box1", box, ident);
  trap->PlaceDaughter("orb1", orb, ident);

  Transformation3D *placement1 = new Transformation3D(5, 5, 5, 0, 0, 0);
  Transformation3D *placement2 = new Transformation3D(-5, 5, 5, 0, 0, 0);   // 45,  0,  0);
  Transformation3D *placement3 = new Transformation3D(5, -5, 5, 0, 0, 0);   // 0, 45,  0);
  Transformation3D *placement4 = new Transformation3D(5, 5, -5, 0, 0, 0);   // 0,  0, 45);
  Transformation3D *placement5 = new Transformation3D(-5, -5, 5, 0, 0, 0);  // 45, 45,  0);
  Transformation3D *placement6 = new Transformation3D(-5, 5, -5, 0, 0, 0);  // 45,  0, 45);
  Transformation3D *placement7 = new Transformation3D(5, -5, -5, 0, 0, 0);  // 0, 45, 45);
  Transformation3D *placement8 = new Transformation3D(-5, -5, -5, 0, 0, 0); // 45, 45, 45);

  world->PlaceDaughter("trap_0", trap, placement1);
  world->PlaceDaughter("trap_1", trap, placement2);
  world->PlaceDaughter("trap_2", trap, placement3);
  world->PlaceDaughter("trap_3", trap, placement4);
  world->PlaceDaughter("trap_4", trap, placement5);
  world->PlaceDaughter("trap_5", trap, placement6);
  world->PlaceDaughter("trap_6", trap, placement7);
  world->PlaceDaughter("trap_7", trap, placement8);

  VPlacedVolume *w = world->Place();
  GeoManager::Instance().SetWorld(w);
  GeoManager::Instance().CloseGeometry();

  // cleanup
  delete ident;
  delete placement1;
  delete placement2;
  delete placement3;
  delete placement4;
  delete placement5;
  delete placement6;
  delete placement7;
  delete placement8;
  return w;
}

void VerifyVecGeomGeometry() {
  size_t numPlaced = GeoManager::Instance().GetPlacedVolumesCount();
  size_t numLogical = GeoManager::Instance().GetRegisteredVolumesCount();
  size_t numNodes = GeoManager::Instance().GetTotalNodeCount();
  std::cout<<" *** VecGeom volumes: logical="<< numLogical <<", placed="<< numPlaced <<", #nodes="<< numNodes <<"\n";
}

#ifdef VECGEOM_ROOT
void VerifyRootGeometry() {
  size_t rootNumShapes = gGeoManager->GetListOfShapes()->GetEntries();
  size_t rootNumVols   = gGeoManager->GetListOfVolumes()->GetEntries();
  size_t rootNumNodes  = gGeoManager->GetListOfNodes()->GetEntries();
  size_t rootNumPhys   = gGeoManager->GetListOfPhysicalNodes()->GetEntries();
  std::cout<<" *** Root geometry: #shapes="<< rootNumShapes <<", #vols="<< rootNumVols
	   <<", #phys="<< rootNumPhys <<", #nodes="<< rootNumNodes <<"\n\n";
}
#endif

#ifdef VECGEOM_GEANT4
void VerifyGeant4Geometry() {
  G4VPhysicalVolume* g4vol = G4GeoManager::Instance().GetNavigator()->GetWorldVolume();
  size_t g4NumPhysVol = G4PhysicalVolumeStore::GetInstance()->size();
  size_t g4NumLogVol = G4LogicalVolumeStore::GetInstance()->size();
  std::cout<<" *** Geant4 volumes: World="<< g4vol->GetName() <<
    " -- #logVols="<< g4NumLogVol <<", #physVols="<< g4NumPhysVol <<"\n";
}
#endif

// Purpose: parse a .gdml file into VecGeom.
//   If VECGEOM_GEANT4 is defined, loads also a Geant4 geometry.
//   If VECGEOM_ROOT is defined, loads also a Root geometry.
//
// input: a .gdml filename char*
//
bool LoadFromGDML(const char* fname) {

  bool loaded = false;
  std::string gdmlName(fname);
  if ( gdmlName.find("gdml") >= gdmlName.size() ) {
    gdmlName += ".gdml";
  }
  std::cout<<" Using Geant4 GDML parser on file "<< gdmlName <<"\n";

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
  }

  #ifdef VECGEOM_GEANT4
  // load Geant4 geometry from ROOT
  G4GeoManager::Instance().GetG4GeometryFromROOT();
  #endif

#else
  #ifdef VECGEOM_GEANT4
    std::cout<<"*** Using Geant4 GDML parser...\n";
    G4GeoManager::Instance().LoadG4Geometry(gdmlName.c_str(), true);
  #endif
#endif

  return loaded;
}

// Saves a VecGeom geometry from memory to a GDML file
void SaveGDML(const char* filename) {
  auto aBackend       = vgdml::Backend();
  auto const loaded   = aBackend.Load(filename);
  aBackend.Save(loaded, filename);
}

#ifdef VECGEOM_ROOT
bool LoadFromRoot(const char* fname) {
  // loads both Root and VecGeom geometries
  RootGeoManager::Instance().LoadRootGeometry( fname );

#ifdef VECGEOM_GEANT4
  // load Geant4 geometry directly from Root
  G4VPhysicalVolume *world = G4GeoManager::Instance().GetG4GeometryFromROOT();
  G4GeoManager::Instance().LoadG4Geometry(world);
#endif
  return true;
}
#endif

//====================

int main(int argc, char *argv[])
{
  OPTION_INT(ntracks, 1024);
  OPTION_INT(nreps, 3);
  OPTION_STRING(geometry, "navBench");
  OPTION_STRING(logvol, "world");
  OPTION_DOUBLE(bias, 0.8f);
#ifdef VECGEOM_ROOT
  OPTION_BOOL(vis, false);
#endif

  // default values used above are always printed.  If help true, stop now, so user will know which options
  // are available, and what the default values are.
  OPTION_BOOL(help, false);
  if (help) return 0;

#ifdef VECGEOM_ENABLE_CUDA
  // If CUDA enabled, then GPU hardware is required!
  int nDevice;
  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
  }
  else {
    std::cout << "\n ***** No Cuda Capable Device!!! *****\n" << std::endl;
    return 0;
  }
#endif

  std::cout<<" Setting up geometry: "<< geometry <<"\n";
  bool loaded = false;
  if (geometry.find("gdml") < geometry.size()) {
    // read GDML file -- load VG, G4 and Root geometries
    loaded = LoadFromGDML( geometry.c_str() );
  }
  else if (geometry.find("root") < geometry.size()) {
    loaded = LoadFromRoot( geometry.c_str() );
  }
  else {

    if (geometry.compare("navBench") == 0) {
      SetupGeometry();

      //SaveGDML("navBench.out.gdml");

#ifdef VERBOSE
      auto *world = GeoManager::Instance().GetWorld();

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
      RootGeoManager::Instance().ExportToROOTGeometry(GeoManager::Instance().GetWorld(), "navBench.root");
      RootGeoManager::Instance().Clear();
#endif

      //#ifdef VECGEOM_ROOT
      // Now try to read back in.  This is needed to make comparisons to VecGeom easily,
      // since it builds VecGeom geometry based on the ROOT geometry and its TGeoNodes.
      // auto rootgeom = geometry + ".root";
      // RootGeoManager::Instance().set_verbose(0);
      // RootGeoManager::Instance().LoadRootGeometry(rootgeom.c_str());
      //#endif
    } // endif("navBench")

  } // else -- no GDML or ROOT file provided

  // verify volumes
  VerifyVecGeomGeometry();
#ifdef VECGEOM_GEANT4
  VerifyGeant4Geometry();
#endif
#ifdef VECGEOM_ROOT
  VerifyRootGeometry();
#endif

  // Visualization
#ifdef VECGEOM_ROOT
  if (vis) { // note that visualization block returns, excluding the rest of benchmark
    Visualizer visualizer;
    visualizer.SetVerbosity(1);
    const VPlacedVolume *world = GeoManager::Instance().GetWorld();
    visualizer.AddVolume(*world);

    Vector<Daughter> const *daughters = world->GetLogicalVolume()->GetDaughtersp();
    for (size_t i = 0; i < daughters->size(); ++i) {
      VPlacedVolume const *daughter = (*daughters)[i];
      Transformation3D const &trf1  = *(daughter->GetTransformation());
      visualizer.AddVolume(*daughter, trf1);

      Vector<Daughter> const* daughters2 = daughter->GetLogicalVolume()->GetDaughtersp();
      for (int ii=0; ii < (int)daughters2->size(); ++ii) {
	VPlacedVolume const* daughter2 = (*daughters2)[ii];
	Transformation3D const& trf2 = *(daughter2->GetTransformation());
	Transformation3D comb = trf1;
	comb.MultiplyFromRight(trf2);
	visualizer.AddVolume(*daughter2, comb);
      }
    }

    std::vector<VPlacedVolume *> v1;
    GeoManager::Instance().getAllPlacedVolumes(v1);
    for (auto &plvol : v1) {
      std::cerr << "placedVol=" << plvol << ", name=" << plvol->GetName() << ", world=" << world << ", <"
                << world->GetName() << ", " << GeoManager::Instance().GetWorld() << ">\n";
    }

    visualizer.Show();
    return 0;
  }
#endif

  std::cout << "\n*** Validating VecGeom navigation...\n";

  const LogicalVolume *startVolume = GeoManager::Instance().GetWorld()->GetLogicalVolume();
  if (logvol.compare("world") != 0) {
    startVolume = GeoManager::Instance().FindLogicalVolume(logvol.c_str());
  }

#ifdef VERBOSE
  std::cout << "NavigationBenchmark: logvol=<" << logvol << ">, startVolume=<"
            << (startVolume ? startVolume->GetLabel() : "NULL") << ">\n";
  if (startVolume) std::cout << *startVolume << "\n";
#endif

  // no more than about 1000 points used for validation
  int np = Min(ntracks, 1024);

  // prepare tracks to be used for benchmarking
  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);

  Vector3D<Precision> samplingVolume(10, 10, 10);
  vecgeom::volumeUtilities::FillRandomPoints(samplingVolume, points);
  vecgeom::volumeUtilities::FillRandomDirections(dirs);

  // run validation on subsample of np tracks
  Precision *maxSteps = (Precision *)vecCore::AlignedAlloc(32, sizeof(Precision) * np);
  for (int i = 0; i < np; ++i)
    maxSteps[i] = 10. * RNG::Instance().uniform();

  // Must be validated before being benchmarked
  bool ok = validateVecGeomNavigation(np, points, dirs, maxSteps);
  if (!ok) {
    std::cout << "VecGeom validation failed." << std::endl;
    //return 1;
  }
  else {
    std::cout << "VecGeom validation passed." << std::endl;
  }

  // on mic.fnal.gov CPUs, loop execution takes ~70sec for ntracks=10M
  while (np <= ntracks) {
    std::cout << "\n*** Running navigation benchmarks with ntracks=" << np << " and nreps=" << nreps << ".\n";
    runNavigationBenchmarks(startVolume, np, nreps, maxSteps, bias);
    np *= 8;
  }

  // cleanup
  vecCore::AlignedFree(maxSteps);
#ifdef VECGEOM_ROOT
  RootGeoManager::Instance().Clear();
#endif
  return 0;
}
