/// \file NavigationBenchmarker.h
/// \author Guilherme Lima (lima at fnal dot gov)
//
// 2014-11-26 G.Lima - created, by adapting Johannes' Benchmarker for navigation

#include "benchmarking/NavigationBenchmarker.h"

#include "base/SOA3D.h"
#include "base/Stopwatch.h"
#include "volumes/utilities/VolumeUtilities.h"

#include "volumes/PlacedVolume.h"
#include "navigation/NavigationState.h"
#include "navigation/SimpleNavigator.h"

#ifdef VECGEOM_ROOT
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#endif

namespace vecgeom {

//==================================

double benchmarkLocatePoint( VPlacedVolume const* world, int nPoints, int nReps,
                             SOA3D<Precision> const& points) {

  NavigationState * state = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
  vecgeom::SimpleNavigator nav;

  Stopwatch timer;
  timer.Start();
  for(int n=0; n<nReps; ++n) {
    for( int i=0; i<nPoints; ++i ) {
      state->Clear();
      nav.LocatePoint( world, points[i], *state, true);
    }
  }
  Precision elapsed = timer.Stop();

  NavigationState::ReleaseInstance( state );
  return elapsed;
}

//==================================

double benchmarkSerialNavigation( VPlacedVolume const* world, int nPoints, int nReps,
                                  SOA3D<Precision> const& points,
                                  SOA3D<Precision> const& dirs ) {

  NavigationState ** curStates = new NavigationState*[nPoints];
  for( int i=0; i<nPoints; ++i) curStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  NavigationState ** newStates = new NavigationState*[nPoints];
  for( int i=0; i<nPoints; ++i) newStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  vecgeom::SimpleNavigator nav;

  for( int i=0; i<nPoints; ++i ) {
    curStates[i]->Clear();
    nav.LocatePoint( world, points[i], *curStates[i], true);
  }

  Precision maxStep = kInfinity, step=0.0;
  Stopwatch timer;
  timer.Start();
  for(int n=0; n<nReps; ++n) {
    for( int i=0; i<nPoints; ++i ) {
         nav.FindNextBoundaryAndStep( points[i], dirs[i], *(curStates[i]),
                                      *(newStates[i]), maxStep, step );
    }
  }
  Precision elapsed = timer.Stop();

  // cleanup
  for( int i=0; i<nPoints; ++i) {
    NavigationState::ReleaseInstance( curStates[i] );
    NavigationState::ReleaseInstance( newStates[i] );
  }
  delete[] curStates;
  delete[] newStates;

  return elapsed;
}

//==================================

double benchmarkVectorNavigation( VPlacedVolume const* world, int nPoints, int nReps,
                                  SOA3D<Precision> const& points,
                                  SOA3D<Precision> const& dirs ) {

  NavigationState ** curStates = new NavigationState*[nPoints];
  for( int i=0; i<nPoints; ++i) curStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  NavigationState ** newStates = new NavigationState*[nPoints];
  for( int i=0; i<nPoints; ++i) newStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  SOA3D<Precision> workspace1(nPoints);
  SOA3D<Precision> workspace2(nPoints);

  int * intworkspace = (int *) _mm_malloc(sizeof(int)*nPoints,32);
  Precision * maxSteps = (Precision *) _mm_malloc(sizeof(Precision)*nPoints,32);
  Precision * vecSteps = (Precision *) _mm_malloc(sizeof(Precision)*nPoints,32);
  Precision * safeties = (Precision *) _mm_malloc(sizeof(Precision)*nPoints,32);

  for (int i=0;i<nPoints;++i) maxSteps[i] = kInfinity;
  memset(vecSteps, 0, sizeof(Precision)*nPoints);
  memset(safeties, 0, sizeof(Precision)*nPoints);

  vecgeom::SimpleNavigator nav;

  for( int i=0; i<nPoints; ++i ) {
    curStates[i]->Clear();
    nav.LocatePoint( world, points[i], *curStates[i], true);
  }

  Stopwatch timer;
  timer.Start();
  for(int n=0; n<nReps; ++n) {
    nav.FindNextBoundaryAndStep( points, dirs, workspace1, workspace2, curStates, newStates,
                                 maxSteps, safeties, vecSteps, intworkspace );
  }
  Precision elapsed = timer.Stop();

  // cleanup
  for( int i=0; i<nPoints; ++i) {
    NavigationState::ReleaseInstance( curStates[i] );
    NavigationState::ReleaseInstance( newStates[i] );
  }
  delete[] curStates;
  delete[] newStates;

  _mm_free(intworkspace);
  _mm_free(maxSteps);
  _mm_free(vecSteps);
  _mm_free(safeties);

  return elapsed;
}

//==================================

#ifdef VECGEOM_ROOT
double benchmarkROOTNavigation( VPlacedVolume const* world, int nPoints, int nReps,
                                SOA3D<Precision> const& points,
                                SOA3D<Precision> const& dirs ) {

  TGeoNavigator * rootnav = ::gGeoManager->GetCurrentNavigator();
  TGeoNode ** rootNodes = new TGeoNode*[nPoints];

  Precision * maxSteps = (Precision *) _mm_malloc(sizeof(Precision)*nPoints,32);
  for (int i=0;i<nPoints;++i) maxSteps[i] = kInfinity;

  Stopwatch timer;
  timer.Start();
  for(int n=0; n<nReps; ++n) {
    for (int i=0;i<nPoints;++i) {
      Vector3D<Precision> const& pos = points[i];
      Vector3D<Precision> const& dir = dirs[i];

      rootnav->ResetState();
      rootNodes[i] = rootnav->FindNode( pos.x(), pos.y(), pos.z() );

      rootnav->SetCurrentPoint( pos.x(), pos.y(), pos.z() );
      rootnav->SetCurrentDirection( dir.x(), dir.y(), dir.z() );
      rootnav->FindNextBoundaryAndStep( maxSteps[i] );
    }
  }

  // cleanup
  delete[] rootNodes;
  _mm_free(maxSteps);

  return timer.Stop();
}
#endif

//==================================
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

   for (int i=0;i<1024;++i) NavigationState::ReleaseInstance( states[i] );
   delete states;
   _mm_free(safeties);
}

//=======================================
/// Function to run navigation benchmarks
void runNavigationBenchmarks( VPlacedVolume const* world, int np, int nreps) {

  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);

  vecgeom::volumeUtilities::FillUncontainedPoints( *world, points );
  vecgeom::volumeUtilities::FillRandomDirections( dirs );

  Precision cputime;

  cputime = benchmarkLocatePoint(world,np,nreps,points);
  printf("CPU elapsed time (locating and setting steps) %f ms\n", 1000.*cputime);

  cputime = benchmarkSerialNavigation(world,np,nreps,points, dirs);
  printf("CPU elapsed time (serialized navigation) %f ms\n", 1000.*cputime);

  cputime = benchmarkVectorNavigation(world,np,nreps,points, dirs);
  printf("CPU elapsed time (vectorized navigation) %f ms\n", 1000.*cputime);

#ifdef VECGEOM_ROOT
  cputime = benchmarkROOTNavigation(world,np,nreps,points, dirs);
  printf("CPU elapsed time (ROOT navigation) %f ms\n", 1000.*cputime);
#endif

  return;
}

//=======================================

#ifdef VECGEOM_ROOT
bool validateNavigationStepAgainstRoot( Vector3D<Precision> const& pos, Vector3D<Precision> const& dir,
                                        NavigationState const& testState, Precision maxStep, Precision testStep)
{
  bool result = true;

  TGeoNavigator* rootnav = ::gGeoManager->GetCurrentNavigator();
  rootnav->ResetState();
  rootnav->FindNode( pos.x(), pos.y(), pos.z() );
  rootnav->SetCurrentPoint( pos.x(), pos.y(), pos.z() );
  rootnav->SetCurrentDirection( dir.x(), dir.y(), dir.z() );
  rootnav->FindNextBoundaryAndStep( maxStep );

  if( testState.Top() == NULL ) {
    if (! rootnav->IsOutside() ) {
      result = false;
      std::cerr << " OUTSIDEERROR \n";
    }
  }
  else if( Abs(testStep - rootnav->GetStep()) > kTolerance
           || rootnav->GetCurrentNode() != RootGeoManager::Instance().tgeonode(testState.Top()) )
  {
    result = false;
    std::cerr << "*** ERROR on validateAgainstROOT: "
              <<" ROOT node="<< rootnav->GetCurrentNode()->GetName()
              <<" outside="<< rootnav->IsOutside()
              <<" step="<< rootnav->GetStep()
              <<" <==> VecGeom node="<< (testState.Top() ? testState.Top()->GetLabel() : "NULL")
              <<" step="<< testStep
              << std::endl;

    std::cerr<< rootnav->GetCurrentNode() <<' '<< RootGeoManager::Instance().tgeonode(testState.Top()) << std::endl;
  }

  return result;
}
#endif // VECGEOM_ROOT

//=======================================

bool validateVecGeomNavigation( VPlacedVolume const* top, int npoints) {

  bool result = true;
  int errorCount = 0;
  int np = Min( npoints, 4096 );  // no more than 4096 points used for validation

  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);
  vecgeom::volumeUtilities::FillUncontainedPoints( *top, points );
  vecgeom::volumeUtilities::FillRandomDirections( dirs );

  // now setup all the navigation states - one loop at a time for better data locality
  NavigationState ** origStates = new NavigationState*[np];
  NavigationState ** vgSerialStates = new NavigationState*[np];

  for (int i=0;i<np;++i) origStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
  for (int i=0;i<np;++i) vgSerialStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  vecgeom::SimpleNavigator nav;
  Precision * maxSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
  Precision * refSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);

  for (int i=0;i<np;++i) maxSteps[i] = kInfinity;
  memset(refSteps, 0, sizeof(Precision)*np);

  // navigation using the serial interface
  for (int i=0; i<np; ++i) {
    Vector3D<Precision> const& pos = points[i];
    Vector3D<Precision> const& dir = dirs[i];

    nav.LocatePoint( top, pos, *origStates[i], true);
    nav.FindNextBoundaryAndStep( pos, dir, *origStates[i], *vgSerialStates[i], maxSteps[i], refSteps[i]);

#ifdef VECGEOM_ROOT
    // validate serial interface agains ROOT, if available
    bool ok = validateNavigationStepAgainstRoot(pos, dir, *vgSerialStates[i], maxSteps[i], refSteps[i]);
    result &= ok;
    if( !ok ) {
      ++errorCount;
      TGeoNavigator * rootnav = ::gGeoManager->GetCurrentNavigator();
      std::cout << "\n=======> Summary: ITERATION " << i
                <<" pos = " << pos <<" dir = " << dir << "\n";
      std::cout << "ROOT GOES HERE: " << rootnav->GetCurrentNode()->GetName()
                << " with step: "<< rootnav->GetStep() << "\n";
      if(vgSerialStates[i]->Top()==NULL) {
        std::cout<<"VECGEOM GOES TO <NULL> with step: "<< refSteps[i] <<"\n";
      }
      else {
        std::cout<<"VECGEOM GOES TO <" << vgSerialStates[i]->Top()->GetLabel() << "> with step: "<< refSteps[i] <<"\n";
        // nav.InspectEnvironmentForPointAndDirection( pos, dir, *origState );
      }
    }
#endif // VECGEOM_ROOT
  }

  std::cout<<"VecGeom navigation - serial interface: #errors = "<< errorCount <<" / "<< np << std::endl;
  //=== Vector interface

  NavigationState ** vgVectorStates = new NavigationState*[np];
  for (int i=0;i<np;++i) vgVectorStates[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

  SOA3D<Precision> workspace1(np);
  SOA3D<Precision> workspace2(np);
  int * intworkspace = (int *) _mm_malloc(sizeof(int)*np,32);

  Precision * vecSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
  Precision * safeties = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
  memset(vecSteps, 0, sizeof(Precision)*np);
  memset(safeties, 0, sizeof(Precision)*np);

  nav.FindNextBoundaryAndStep( points, dirs, workspace1, workspace2, origStates, vgVectorStates,
                               maxSteps, safeties, vecSteps, intworkspace );

  // compare with serial results
  //  ok = compareNavigationResults(refSteps, vgSerialStates, vecSteps, vgVectorStates);
  errorCount = 0;
  for(int i=0; i<np; ++i) {
    bool mismatch = false;
    void* void1 = (void*)vgSerialStates[i]->Top();
    void* void2 = (void*)vgVectorStates[i]->Top();
    if( Abs( vecSteps[i] - refSteps[i] ) > kTolerance )                   mismatch = true;
    if( void1 && void2 && void1 != void2 )                                mismatch = true;
    if( vgSerialStates[i]->IsOnBoundary() != vgVectorStates[i]->IsOnBoundary()) mismatch = true;
    if( safeties[i] != nav.GetSafety( points[i], *origStates[i] ))   mismatch = true;
    if(mismatch) {
      ++errorCount;
      std::cout<<"Vector navigation problems: iteration "<< i
               <<" steps: "<< refSteps[i] <<" / "<< vecSteps
               <<" navStates: "<< ( void1 ? vgSerialStates[i]->Top()->GetLabel() : "NULL")
               << (vgSerialStates[i]->IsOnBoundary() ? "*" : "")
               <<" / "<< (void2 ? vgVectorStates[i]->Top()->GetLabel() : "NULL")
               << (vgVectorStates[i]->IsOnBoundary() ? "*" : "")
               << std::endl;
    }
  }

  std::cout<<"VecGeom navigation - vector interface: #errors = "<< errorCount <<" / "<< np << std::endl;

#ifdef VECGEOM_NVCC
  Precision * gpuSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
  printf("Start GPU navigation...\n");
  runNavigationCuda(GeoManager::Instance().GetWorld(), np,
                    points.x(),  points.y(), points.z(),
                    dirs.x(), dirs.y(), dirs.z(), maxSteps, gpuSteps );

  //=== Comparing results from GPU with serialized navigation
  errorCount = 0;
  for(int i=0; i<np; ++i) {
    bool mismatch = false;
    if( Abs( gpuSteps[i] - refSteps[i] ) > kTolerance )                         mismatch = true;
    // if( vgSerialStates[i]->Top() != vgVectorStates[i]->Top() )                  mismatch = true;
    // if( vgSerialStates[i]->IsOnBoundary() != vgVectorStates[i]->IsOnBoundary()) mismatch = true;
    // if( safeties[i] != nav.GetSafety( points[i], *origStates[i] ))              mismatch = true;
    if(mismatch) {
      ++errorCount;
      std::cout<<"GPU navigation problems: iteration "<< i
               <<" steps: "<< refSteps[i] <<" / "<< gpuSteps
               // <<" navStates: "<< vgSerialStates[i]->Top()->GetLabel()
               // << (vgSerialStates[i]->IsOnBoundary() ? "*" : "")
               // <<" / "<< vgVectorStates[i]->Top()->GetLabel()
               // << (vgVectorStates[i]->IsOnBoundary() ? "*" : "")
               << std::endl;
    }
  }

  std::cout<<"VecGeom navigation on the GPUs: #errors = "<< errorCount <<" / "<< np << std::endl;
#endif // VECGEOM_NVCC

  // if(mismatches>0) std::cout << "Navigation test failed with "<< mismatches <<" mismatches\n";
  // else std::cout<<"Navigation test passed.\n";


  //=== cleanup
  for(int i=0; i<np; ++i) NavigationState::ReleaseInstance( origStates[i] );
  for(int i=0; i<np; ++i) NavigationState::ReleaseInstance( vgSerialStates[i] );
  for(int i=0; i<np; ++i) NavigationState::ReleaseInstance( vgVectorStates[i] );
  delete[] origStates;
  delete[] vgSerialStates;
  delete[] vgVectorStates;
  _mm_free(intworkspace);
  _mm_free(maxSteps);
  _mm_free(refSteps);
  _mm_free(vecSteps);
  _mm_free(safeties);
#ifdef VECGEOM_NVCC
  _mm_free(gpuSteps);
#endif
  return result;
}

} // End namespace vecgeom
