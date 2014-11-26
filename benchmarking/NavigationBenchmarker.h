/// \file NavigationBenchmarker.h
/// \author Guilherme Lima (lima at fnal dot gov)
//
// 2014-11-26 G.Lima - created, by adapting Johannes' Benchmarker for navigation

#ifndef VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_
#define VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_

#include "volumes/PlacedVolume.h"
// #include "base/SOA3D.h"

namespace vecgeom {

  double benchmarkLocatePoint( VPlacedVolume const* world, int nPoints, int nReps,
                               SOA3D<Precision> const& points);

  double benchmarkSerialNavigation( VPlacedVolume const* world, int nPoints, int nReps,
                                    SOA3D<Precision> const& points,
                                    SOA3D<Precision> const& dirs );

  double benchmarkVectorNavigation( VPlacedVolume const* world, int nPoints, int nReps,
                                    SOA3D<Precision> const& points,
                                    SOA3D<Precision> const& dirs );

  double benchmarkROOTNavigation( VPlacedVolume const* world, int nPoints, int nReps,
                                  SOA3D<Precision> const& points,
                                  SOA3D<Precision> const& dirs );

} // End namespace vecgeom

#endif // VECGEOM_BENCHMARKING_NAVIGATIONBENCHMARKER_H_
