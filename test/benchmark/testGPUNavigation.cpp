/*
 * testGPUNavigation.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: swenzel, lima
 */

#include "benchmarking/NavigationBenchmarker.h"

#include "volumes/utilities/VolumeUtilities.h"
#include "volumes/Box.h"
#include "volumes/PlacedVolume.h"
#include "base/Transformation3D.h"
#include "base/SOA3D.h"
#include "management/RootGeoManager.h"
#include "navigation/NavigationState.h"
#include "navigation/SimpleNavigator.h"
#include "management/GeoManager.h"
#include "base/Global.h"
#include "ArgParser.h"

#include "navigationgpu.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "utilities/Visualizer.h"

using namespace VECGEOM_NAMESPACE;

// Prepare VecGeom vs. root mapping for comparison purposes
static std::map<VPlacedVolume const*, TGeoNode const*> vg2rootMap;

VPlacedVolume* SetupBoxGeometry() {
  UnplacedBox *worldUnplaced = new UnplacedBox(10, 10, 10);
  UnplacedBox *boxUnplaced = new UnplacedBox(3, 3, 3);
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

  // Visualizer visualizer;
  // visualizer.AddVolume(*w);

  // Vector<Daughter> const* daughters = w->logical_volume()->daughtersp();
  // for(int i=0; i<daughters->size(); ++i) {
  //   VPlacedVolume const* daughter = (*daughters)[i];
  //   visualizer.AddVolume(*daughter, *(daughter->transformation()));
  // }

  // visualizer.Show();

  return w;
}

// function to test safety
void testVectorSafety( VPlacedVolume const* world ){
   SOA3D<Precision> points(1024);
   SOA3D<Precision> workspace(1024);
   Precision * safeties = (Precision *) _mm_malloc(sizeof(Precision)*1024,32);
   vecgeom::volumeUtilities::FillUncontainedPoints( *world, points );

   // now setup all the navigation states
   NavigationState ** states = new NavigationState*[1024];

   vecgeom::SimpleNavigator nav;

   for (int i=0;i<1024;++i){
     states[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
     nav.LocatePoint( world, points[i], *states[i], true);
   }

    // calculate safeties with vector interface
    nav.GetSafeties(points, states, workspace, safeties );

    // verify against serial interface
    for (int i=0;i<1024;++i){
        vecgeom::Assert( safeties[i] == nav.GetSafety( points[i], *states[i] ), ""
                " Problem in VectorSafety (in SimpleNavigator)" );
    }
    std::cout << "Safety test passed\n";
   _mm_free(safeties);
}

// function to test vector navigator
void runNavigationTests( VPlacedVolume const* world, int np, int nreps) {

  std::cerr << "\n\n =============== Running navigation tests...\n";

  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);

  Precision * maxSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
  Precision * refSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
  Precision * rootSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);

  for (int i=0;i<np;++i) maxSteps[i] = kInfinity;

  vecgeom::volumeUtilities::FillUncontainedPoints( *world, points );
  vecgeom::volumeUtilities::FillRandomDirections( dirs );

  // now setup all the navigation states
  NavigationState ** states = new NavigationState*[np];
  NavigationState ** vgVectorStates = new NavigationState*[np];
  NavigationState ** vgSerialStates = new NavigationState*[np];

  // one loop at a time, for better data locality

  for (int i=0;i<np;++i) states[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
  for (int i=0;i<np;++i) vgSerialStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
  for (int i=0;i<np;++i) vgVectorStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  memset(refSteps, 0, sizeof(Precision)*np);
  memset(rootSteps, 0, sizeof(Precision)*np);

  vecgeom::SimpleNavigator nav;
  Precision cputime;

  cputime = benchmarkLocatePoint(world,np,nreps,points);
  printf("CPU elapsed time (locating and setting steps) %f ms\n", 1000.*cputime);

  cputime = benchmarkSerialNavigation(world,np,nreps,points, dirs);
  printf("CPU elapsed time (serialized navigation) %f ms\n", 1000.*cputime);

  cputime = benchmarkVectorNavigation(world,np,nreps,points, dirs);
  printf("CPU elapsed time (vectorized navigation) %f ms\n", 1000.*cputime);

  cputime = benchmarkROOTNavigation(world,np,nreps,points, dirs);
  printf("CPU elapsed time (ROOT navigation) %f ms\n", 1000.*cputime);


//=== Comparing with
return;
#ifdef VECGEOM_CUDA
  // Precision * gpuSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   printf("Start GPU\n");
   RunNavigationCuda(world, np,
                     points.x(),  points.y(), points.z(),
                     dirs.x(), dirs.y(), dirs.z(), maxSteps, gpuSteps );
#endif
   // unsigned mismatches=0;
   // for (int i=0;i<np;++i) {
     // if( abs(vecSteps[i]-gpuSteps[i]) > kTolerance ) {
     //   ++mismatches;
     //   std::cerr << "i " << i << " steps " << vecSteps[i] << " CUDA steps " << gpuSteps[i]
     //             << " - diff="<< (vecSteps[i]-gpuSteps[i]) <<"\n";
     //  vecgeom::Assert( vecSteps[i] == gpuSteps[i], " Problem in CUDA Navigator " );
     // }
   // }

   // if(mismatches>0) std::cout << "Navigation test failed with "<< mismatches <<" mismatches\n";
   // else std::cout<<"Navigation test passed.\n";
   for(int i=0; i<np; ++i) NavigationState::ReleaseInstance( vgSerialStates[i] );
}

bool validateVecGeomNavigation( VPlacedVolume const* top, int npoints) {

  bool result = true;
  Precision maxStep = kInfinity;
  int np = Min( npoints, 4096 );

  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);
  vecgeom::volumeUtilities::FillUncontainedPoints( *top, points );
  vecgeom::volumeUtilities::FillRandomDirections( dirs );

  TGeoNavigator * rootnav = ::gGeoManager->GetCurrentNavigator();
  vecgeom::SimpleNavigator nav;

  Precision vgstep;
  NavigationState* curState = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
  NavigationState* newState = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
  for (int i=0;i<np;++i) {
    Vector3D<Precision> const& pos = points[i];
    Vector3D<Precision> const& dir = dirs[i];

    curState->Clear();
    newState->Clear();
    nav.LocatePoint( top, pos, *curState, true);
    nav.FindNextBoundaryAndStep( pos, dir, *curState, *newState, maxStep, vgstep);

    // vecgeom::Assert( vecSteps[i] == refSteps[i] ,
    //                  " Problem in VectorNavigation (steps) (in SimpleNavigator)" );
    // vecgeom::Assert( vgSerialStates[i]->Top() == vgVectorStates[i]->Top() ,
    //                  " Problem in VectorNavigation (states) (in SimpleNavigator)" );
    // vecgeom::Assert( vgSerialStates[i]->IsOnBoundary() == vgVectorStates[i]->IsOnBoundary(),
    //                  " Problem in VectorNavigation (boundary) (in SimpleNavigator)" );
    // vecgeom::Assert( safeties[i] == nav.GetSafety( points[i], *states[i] ),
    //                  " Problem with safety" );

    rootnav->ResetState();
    rootnav->FindNode( pos.x(), pos.y(), pos.z() );
    rootnav->SetCurrentPoint( pos.x(), pos.y(), pos.z() );
    rootnav->SetCurrentDirection( dir.x(), dir.y(), dir.z() );
    rootnav->FindNextBoundaryAndStep( maxStep );

    std::cout<<"testGPUNav: rootNode="<< rootnav->GetCurrentNode()->GetName() << std::endl;
    if(newState->Top())
       std::cout<<"           vecGnode="<< RootGeoManager::Instance().tgeonode(newState->Top())->GetName() << std::endl;
     else
       std::cout<<"           vecGnode=NULL"<< std::endl;

     if( newState->Top() == NULL ) {
       if (! rootnav->IsOutside() ) {
         result = false;
         std::cerr << " OUTSIDEERROR \n";
       }
     }
     else if( Abs(vgstep - rootnav->GetStep()) > kTolerance
              || rootnav->GetCurrentNode() != RootGeoManager::Instance().tgeonode(newState->Top()) )
     {
       result = false;
       std::cerr << "INSIDE ERROR ON ITERATION " << i <<":"
                 <<" ROOT node="<< rootnav->GetCurrentNode()->GetName()
                 <<" outside="<< rootnav->IsOutside()
                 <<" step="<< rootnav->GetStep()
                 <<" <==> VecGeom node="<< (newState->Top() ? newState->Top()->GetLabel() : "NULL")
                 <<" step="<< vgstep
                 << std::endl;
       std::cerr<< rootnav->GetCurrentNode() <<' '<< RootGeoManager::Instance().tgeonode(newState->Top()) << std::endl;

       std::cerr << "\n=======> Summary: ITERATION " << i;
       std::cerr <<" pos = " << pos;
       std::cerr <<" dir = " << dir << "\n";
       std::cerr << "ROOT GOES HERE: " << rootnav->GetCurrentNode()->GetName();
       std::cerr << " with step: "<< rootnav->GetStep() << "\n";
       if(newState->Top()==NULL) {
         std::cerr<<"VECGEOM GOES TO <NULL> with step: "<< vgstep <<"\n";
       }
       else {
         std::cerr<<"VECGEOM GOES TO <" << newState->Top()->GetLabel() << "> with step: "<< vgstep <<"\n";
         nav.InspectEnvironmentForPointAndDirection( pos, dir, *curState );
       }
     }
  }
  return result;
}

int main(int argc, char* argv[])
{
  OPTION_INT(npoints, 1000);
  OPTION_INT(nreps, 1000);

  VPlacedVolume *world = SetupBoxGeometry();

  testVectorSafety(world);

  // Exporting to ROOT file
  RootGeoManager::Instance().ExportToROOTGeometry( world, "geom1.root" );
  RootGeoManager::Instance().Clear();

  // Now try to read back in.  This is needed to make comparisons to VecGeom easily,
  // since it builds VecGeom geometry based on the ROOT geometry and its TGeoNodes.
  RootGeoManager::Instance().set_verbose(1);
  RootGeoManager::Instance().LoadRootGeometry("geom1.root");

//  validateVecGeomNavigation(world, 1024);

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

  runNavigationTests(GeoManager::Instance().GetWorld(), npoints, nreps);
}
