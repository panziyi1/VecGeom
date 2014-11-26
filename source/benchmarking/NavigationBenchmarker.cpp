/// \file NavigationBenchmarker.h
/// \author Guilherme Lima (lima at fnal dot gov)
//
// 2014-11-26 G.Lima - created, by adapting Johannes' Benchmarker for navigation

#ifndef VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_
#define VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_

#include "benchmarking/NavigationBenchmarker.h"

// #include "base/Global.h"
#include "base/SOA3D.h"
#include "base/Stopwatch.h"

#include "volumes/PlacedVolume.h"
#include "navigation/NavigationState.h"
#include "navigation/SimpleNavigator.h"
//#include "management/RootGeoManager.h"

#ifdef VECGEOM_ROOT
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#endif

namespace vecgeom {

double benchmarkLocatePoint( VPlacedVolume const* world, int nPoints, int nReps,
                             SOA3D<Precision> const& points) {

   NavigationState ** states = new NavigationState*[nPoints];
   for( int i=0; i<nPoints; ++i)
     states[i] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );

   vecgeom::SimpleNavigator nav;

   Stopwatch timer;
   timer.Start();
   for(int n=0; n<nReps; ++n) {
     for( int i=0; i<nPoints; ++i ) {
       states[i]->Clear();
       nav.LocatePoint( world, points[i], *states[i], true);
     }
   }
   Precision elapsed = timer.Stop();

   for( int i=0; i<nPoints; ++i)
     NavigationState::ReleaseInstance( states[i] );

   return elapsed;
}

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

  for( int i=0; i<nPoints; ++i) {
    NavigationState::ReleaseInstance( curStates[i] );
    NavigationState::ReleaseInstance( newStates[i] );
  }

  return elapsed;
}

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

  _mm_free(intworkspace);
  _mm_free(maxSteps);
  _mm_free(vecSteps);
  _mm_free(safeties);

  return elapsed;
}


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
  return timer.Stop();
}

} // End namespace vecgeom

#endif // VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_
