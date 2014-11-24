/*
 * testGPUNavigation.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: swenzel, lima
 */

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
#include "test/benchmark/ArgParser.h"
#include "base/Stopwatch.h"

#include "navigationgpu.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "utilities/Visualizer.h"

using namespace VECGEOM_NAMESPACE;

// Prepare VecGeom vs. root mapping for comparison purposes
static std::map<VPlacedVolume const*, TGeoNode const*> vg2rootMap;

/*
Precision benchVecGeomSerialNavigation(int nreps) {
  StopWatch timer;
  timer.Start();
  for(int n=0; n<nreps; ++n) {
    for (int i=0;i<np;++i) {
      vgSerialStates->Clear();
      nav.FindNextBoundaryAndStep( points[i], dirs[i], *(states[i]),
                                   *(vgSerialStates[i]), pSteps[i], serialStep );
    }
  }
  return timer.Stop();
}


Precision benchVecGeomVectorNavigation(int nreps) {
  StopWatch timer;
  timer.Start();

  for(int n=0; n<nreps; ++n) {
    nav.FindNextBoundaryAndStep( points, dirs, workspace1, workspace2,
                                 states, vgVectorStates, pSteps, safeties, steps, intworkspace );
  }
  return timer.Stop();
}

*/

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

void RecursiveBuildMap(VPlacedVolume const* mother, std::map<VPlacedVolume const*, TGeoNode const*>& aMap) {

  // insert mother entry
  const TGeoNode* node = RootGeoManager::Instance().tgeonode(mother);
  std::cout<<"RecursiveMapping: adding to map: "<< mother <<" <"<< mother->GetLabel() <<">"
           <<" and "<< node <<" <name="<< node->GetName()<<">"<< std::endl;
  aMap.insert( std::pair<VPlacedVolume const*, TGeoNode const*>(mother, RootGeoManager::Instance().tgeonode(mother)) );

  // recursive call for each daughter
  Vector<Daughter> const* daughters = mother->logical_volume()->daughtersp();
  for(int i=0; i<daughters->size(); ++i) {
    VPlacedVolume const* daughter = (*daughters)[i];
    RecursiveBuildMap(daughter, aMap);
  }
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
void testVectorNavigator( VPlacedVolume const* world, int np, int nreps ){
   SOA3D<Precision> points(np);
   SOA3D<Precision> dirs(np);

   Precision * steps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   Precision * pSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   Precision * GPUSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   Precision * safeties = (Precision *) _mm_malloc(sizeof(Precision)*np,32);


   vecgeom::volumeUtilities::FillUncontainedPoints( *world, points );
   vecgeom::volumeUtilities::FillRandomDirections( dirs );

   vecgeom::SimpleNavigator nav;
   Stopwatch timer;

   // now setup all the navigation states
   NavigationState ** states = new NavigationState*[np];
   NavigationState ** vgSerialStates = new NavigationState*[np];
   NavigationState ** vgVectorStates = new NavigationState*[np];

   // one loop at a time, for better data locality
   for (int i=0;i<np;++i) pSteps[i] = kInfinity;
   for (int i=0;i<np;++i) states[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
   for (int i=0;i<np;++i) vgVectorStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
   for (int i=0;i<np;++i) vgSerialStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

   //========= timing 1
   timer.Start();
   for(int n=0; n<nreps; ++n) {
     for (int i=0;i<np;++i) {
       nav.LocatePoint( world, points[i], *states[i], true);
       steps[i] = (Precision)i;
     }
   }
   Precision elapsedCPU = timer.Stop();
   printf("CPU elapsed time (locating and setting steps) %f ms\n", 1000.*elapsedCPU);

   //========= timing 2
   timer.Start();
   for(int n=0; n<nreps; ++n) {
     for (int i=0;i<np;++i) {
       nav.LocatePoint( world, points[i], *states[i], true);
     }
   }
   elapsedCPU = timer.Stop();
   printf("CPU elapsed time (only locating points) %f ms\n", 1000.*elapsedCPU);

   //========= timing 3
   timer.Start();
   for(int n=0; n<nreps; ++n) {
     for (int i=0;i<np;++i) {
       steps[i] = (Precision)i;
     }
   }
   elapsedCPU = timer.Stop();
   printf("CPU elapsed time (only setting steps) %f ms\n", 1000.*elapsedCPU);

   //========== timing 4

   SOA3D<Precision> workspace1(np);
   SOA3D<Precision> workspace2(np);
   int * intworkspace = (int *) _mm_malloc(sizeof(int)*np,32);

   timer.Start();

   // calculate steps with vector interface
   for(int n=0; n<nreps; ++n) {
     printf("vecNavigation loop: n=%d\n",n);
     nav.FindNextBoundaryAndStep( points, dirs, workspace1, workspace2,
                                  states, vgVectorStates, pSteps, safeties, steps, intworkspace );
   }
   elapsedCPU = timer.Stop();
   printf("CPU elapsed time (vectorized navigation): %f ms\n", 1000.*elapsedCPU);
   //===========================================


   //======== timing: serial navigation on CPU
   timer.Start();
   Precision serialStep;
   for(int n=0; n<nreps; ++n) {
     for (int i=0;i<np;++i) {
       vgSerialStates[i]->Clear();
       nav.FindNextBoundaryAndStep( points[i], dirs[i], *(states[i]),
                                    *(vgSerialStates[i]), pSteps[i], serialStep );
     }
   }
   elapsedCPU = timer.Stop();
   printf("CPU elapsed time (serialized navigation): %f ms\n", 1000.*elapsedCPU);
   //===========================================

//=== ROOT Navigation

   TGeoNavigator * rootnav = ::gGeoManager->GetCurrentNavigator();
   TGeoNode ** rootNodes = new TGeoNode*[np];

   //====== timing: vectorized navigation on CPU
   timer.Start();
   for(int n=0; n<nreps; ++n) {
     for (int i=0;i<np;++i) {
       Vector3D<Precision> const& pos = points[i];
       Vector3D<Precision> const& dir = dirs[i];

       rootnav->ResetState();
       rootnav->FindNode( pos.x(), pos.y(), pos.z() );

       rootnav->SetCurrentPoint( pos.x(), pos.y(), pos.z() );
       rootnav->SetCurrentDirection( dir.x(), dir.y(), dir.z() );
//     path->UpdateNavigator( rootnav );
       rootnav->FindNextBoundaryAndStep( pSteps[i] );

       // store node at end of step, to compare with VecGeom results
       rootNodes[i] = rootnav->GetCurrentNode();
     }
   }
   Precision elapsedROOT = timer.Stop();
   printf("ROOT elapsed time (scalar): %f ms\n", 1000.*elapsedROOT);


   // Validation is done at the end, so that any mismatches will more benchmark
   // results off screen (note it doesn't help much if results are on a file)

   // VecGeom validation:  vectorized vs. serial navigation
   NavigationState* cmp = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
   for (int i=0;i<np;++i) {
     cmp->Clear();
     nav.FindNextBoundaryAndStep( points[i], dirs[i], *(states[i]),
               *cmp, pSteps[i], serialStep );

     vecgeom::Assert( steps[i] == serialStep ,
                      " Problem in VectorNavigation (steps) (in SimpleNavigator)" );
     vecgeom::Assert( cmp->Top() == vgVectorStates[i]->Top() ,
                      " Problem in VectorNavigation (states) (in SimpleNavigator)" );
     vecgeom::Assert( cmp->IsOnBoundary() == vgVectorStates[i]->IsOnBoundary(),
                      " Problem in VectorNavigation (boundary) (in SimpleNavigator)" );
     vecgeom::Assert( safeties[i] == nav.GetSafety( points[i], *states[i] ),
                      " Problem with safety " );
   }


   // compare ROOT and vectorized VecGeom results
   for (int i=0;i<np;++i) {
     if( cmp->Top() == NULL ) {
        if (! rootnav->IsOutside() )
            std::cerr << " OUTSIDEERROR \n";
     }
     else if( rootnav->GetCurrentNode() != RootGeoManager::Instance().tgeonode(cmp->Top()) ){

         std::cerr << "INSIDE ERROR ON ITERATION " << i <<":"
                   <<" ROOT node="<< rootnav->GetCurrentNode()->GetName()
                   <<" ROOT outside"<< rootnav->IsOutside()
                   <<" VecGeom node="<< (cmp->Top() ? cmp->Top()->GetLabel() : "NULL")
                   << std::endl;
         std::cerr<< rootnav->GetCurrentNode() <<' '<< RootGeoManager::Instance().tgeonode(cmp->Top()) << std::endl;

         Vector3D<Precision> const& pos = points[i];
         Vector3D<Precision> const& dir = dirs[i];
         std::cerr << "\n=======> Summary: ITERATION " << i;
         std::cerr <<" pos = " << pos;
         std::cerr <<" dir = " << dir << "\n";
         std::cerr << "ROOT GOES HERE: " << rootnav->GetCurrentNode()->GetName();
         std::cerr << " with step: "<< rootnav->GetStep() << "\n";
         if( !cmp->Top() ) {
           std::cerr<<"VECGEOM GOES TO <NULL> with step: "<< serialStep <<"\n";
         }
         else {
           std::cerr << "VECGEOM GOES TO <" << cmp->Top()->GetLabel() << "> with step: "<< serialStep <<"\n";
           nav.InspectEnvironmentForPointAndDirection( pos, dir, *states[i] );
         }
     }
   }

//=== Comparing with 
return;
#ifdef VECGEOM_CUDA
   printf("Start GPU\n");
   RunNavigationCuda(world, np,
                     points.x(),  points.y(), points.z(),
                     dirs.x(), dirs.y(), dirs.z(), pSteps, GPUSteps );
#endif
   unsigned mismatches=0;
   for (int i=0;i<np;++i) {
     if( abs(steps[i]-GPUSteps[i]) > kTolerance ) {
       ++mismatches;
       std::cerr << "i " << i << " steps " << steps[i] << " CUDA steps " << GPUSteps[i]
                 << " - diff="<< (steps[i]-GPUSteps[i]) <<"\n";
     //  vecgeom::Assert( steps[i] == GPUSteps[i], " Problem in CUDA Navigator " );
     }
   }

   if(mismatches>0) std::cout << "Navigation test failed with "<< mismatches <<" mismatches\n";
   else std::cout<<"Navigation test passed.\n";
   _mm_free(steps);
   _mm_free(intworkspace);
   _mm_free(pSteps);
   _mm_free(safeties);
}


int main(int argc, char* argv[])
{
  OPTION_INT(npoints, 100);
  OPTION_INT(nreps, 1024);
  VPlacedVolume *world = SetupBoxGeometry();

  testVectorSafety(world);

  // exporting to ROOT file
  RootGeoManager::Instance().ExportToROOTGeometry( world, "geom1.root" );
  RootGeoManager::Instance().Clear();

  // now try to read back in.  This is needed to make comparisons to VecGeom easier
  RootGeoManager::Instance().set_verbose(1);
  RootGeoManager::Instance().LoadRootGeometry("geom1.root");

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

  std::cerr << "\n\n =============== testing vector navigation...\n";
  testVectorNavigator(GeoManager::Instance().GetWorld(), npoints, nreps);
}
