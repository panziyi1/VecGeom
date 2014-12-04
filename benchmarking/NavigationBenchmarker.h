/// \file NavigationBenchmarker.h
/// \author Guilherme Lima (lima at fnal dot gov)
//
// 2014-11-26 G.Lima - created with help from Sandro

#ifndef VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_
#define VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_

#include "base/SOA3D.h"

namespace vecgeom {

  class NavigationState;
  class VPlacedVolume;

  double benchmarkLocatePoint(
    VPlacedVolume const* top,
    int nPoints,
    int nReps,
    SOA3D<Precision> const& points
    );

  double benchmarkSerialNavigation(
    VPlacedVolume const* top,
    int nPoints,
    int nReps,
    SOA3D<Precision> const& points,
    SOA3D<Precision> const& dirs
    );

  double benchmarkVectorNavigation(
    VPlacedVolume const* world,
    int nPoints,
    int nReps,
    SOA3D<Precision> const& points,
    SOA3D<Precision> const& dirs
    );

  double benchmarkROOTNavigation(
    VPlacedVolume const* world,
    int nPoints,
    int nReps,
    SOA3D<Precision> const& points,
    SOA3D<Precision> const& dirs );

  void testVectorSafety( VPlacedVolume const* top );

  void runNavigationBenchmarks( VPlacedVolume const* top, int np, int nreps);

  bool validateNavigationStepAgainstRoot(
    Vector3D<Precision> const& pos,
    Vector3D<Precision> const& dir,
    NavigationState const& testState,
    Precision maxStep,
    Precision testStep
    );

  bool validateVecGeomNavigation( VPlacedVolume const* top, int npoints);

} // End namespace vecgeom

#endif // VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_
