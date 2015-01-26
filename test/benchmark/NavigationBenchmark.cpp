/*
 * testGPUNavigation.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: swenzel, lima
 */

#include "benchmarking/NavigationBenchmarker.h"
#include "ArgParser.h"

#ifdef VECGEOM_ROOT
#include "management/RootGeoManager.h"
#include "utilities/Visualizer.h"
#endif

#include "management/GeoManager.h"
#include "volumes/Box.h"

using namespace VECGEOM_NAMESPACE;

VPlacedVolume* SetupBoxGeometry() {
  UnplacedBox *worldUnplaced = new UnplacedBox(10, 10, 10);
  UnplacedBox *boxUnplaced = new UnplacedBox(2.5, 2.5, 2.5);
  Transformation3D *placement1 = new Transformation3D( 5,  5,  5,  0,  0,  0);
  Transformation3D *placement2 = new Transformation3D(-5,  5,  5, 45,  0,  0);
  Transformation3D *placement3 = new Transformation3D( 5, -5,  5,  0, 45,  0);
  Transformation3D *placement4 = new Transformation3D( 5,  5, -5,  0,  0, 45);
  Transformation3D *placement5 = new Transformation3D(-5, -5,  5, 45, 45,  0);
  Transformation3D *placement6 = new Transformation3D(-5,  5, -5, 45,  0, 45);
  Transformation3D *placement7 = new Transformation3D( 5, -5, -5,  0, 45, 45);
  Transformation3D *placement8 = new Transformation3D(-5, -5, -5, 45, 45, 45);
  LogicalVolume *world = new LogicalVolume("world",worldUnplaced);
  LogicalVolume *box = new LogicalVolume("box",boxUnplaced);
  world->PlaceDaughter("box0",box, placement1);
  world->PlaceDaughter("box1",box, placement2);
  world->PlaceDaughter("box2",box, placement3);
  world->PlaceDaughter("box3",box, placement4);
  world->PlaceDaughter("box4",box, placement5);
  world->PlaceDaughter("box5",box, placement6);
  world->PlaceDaughter("box6",box, placement7);
  world->PlaceDaughter("box7",box, placement8);
  VPlacedVolume  * w = world->Place();
  GeoManager::Instance().SetWorld(w);
  GeoManager::Instance().CloseGeometry();

  return w;
}


int main(int argc, char* argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_INT(nreps, 3);
#ifdef VECGEOM_ROOT
  OPTION_BOOL(vis, false);
#endif

  VPlacedVolume *world = SetupBoxGeometry();

#ifdef VECGEOM_ROOT
  if(vis) {  // note that visualization block returns, excluding the rest of benchmark
    Visualizer visualizer;
    visualizer.AddVolume(*world);

    Vector<Daughter> const* daughters = world->logical_volume()->daughtersp();
    for(int i=0; i<daughters->size(); ++i) {
      VPlacedVolume const* daughter = (*daughters)[i];
      visualizer.AddVolume(*daughter, *(daughter->transformation()));
    }

    visualizer.Show();
    return 0;
  }
#endif

  testVectorSafety(world);

#ifdef VECGEOM_ROOT
  // Exporting to ROOT file
  RootGeoManager::Instance().ExportToROOTGeometry( world, "geom1.root" );
  RootGeoManager::Instance().Clear();

  // Now try to read back in.  This is needed to make comparisons to VecGeom easily,
  // since it builds VecGeom geometry based on the ROOT geometry and its TGeoNodes.
  RootGeoManager::Instance().set_verbose(0);
  RootGeoManager::Instance().LoadRootGeometry("geom1.root");
#endif

  std::cout<<"\n*** Validating VecGeom navigation..."<< std::endl;
  bool ok = validateVecGeomNavigation(np, points, dirs);

  // Must be validated before being benchmarked
  if(!ok) {
    std::cout<<"VecGeom validation failed."<< std::endl;
    return 1;
  }

  std::cout<<"VecGeom validation passed."<< std::endl;

  // on FNAL GPUs, loop execution takes ~70sec for npoints=10M
  while(npoints<=1000000) {
    std::cout<<"\n*** Running navigation benchmarks with npoints="<<npoints<<" and nreps="<< nreps <<".\n";
    runNavigationBenchmarks(GeoManager::Instance().GetWorld(), npoints, nreps);
    npoints*=10;
  }

/*
// GPU part
  int nDevice;
  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
  }
  else {
    std::cout << "No Cuda Capable Device ... " << std::endl;
    return 0;
  }
*/

}
