#ifndef VECGEOM_ENABLE_CUDA

#include "volumes/LogicalVolume.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"
#include "base/Stopwatch.h"
#include <iostream>
#include "volumes/Extruded.h"
#include "volumes/SExtru.h"

using namespace vecgeom;

#endif

int main(int argc, char *argv[])
{
#ifndef VECGEOM_ENABLE_CUDA
  OPTION_INT(npoints, 1024);
  OPTION_INT(nrep, 4);
  OPTION_INT(nvert, 8);
  OPTION_INT(nsect, 2);
  OPTION_BOOL(convex, 't');
  OPTION_BOOL(tsl, 't');

  double rmin = 10.;
  double rmax = 20.;

  vecgeom::XtruVertex2 *vertices = new vecgeom::XtruVertex2[nvert];
  vecgeom::XtruSection *sections = new vecgeom::XtruSection[nsect];
  double *x                      = new double[nvert];
  double *y                      = new double[nvert];

  double phi = 2. * kPi / nvert;
  double r;
  for (int i = 0; i < nvert; ++i) {
    r                           = rmax;
    if (i % 2 > 0 && !convex) r = rmin;
    vertices[i].x               = r * vecCore::math::Cos(i * phi);
    vertices[i].y               = r * vecCore::math::Sin(i * phi);
    x[i]                        = vertices[i].x;
    y[i]                        = vertices[i].y;
  }
  for (int i = 0; i < nsect; ++i) {
    sections[i].fOrigin.Set(0, 0, -20. + i * 40. / (nsect - 1));
    sections[i].fScale = 1;
  }

  Transformation3D placement(0, 0, 0);
  UnplacedBox worldUnplaced = UnplacedBox(20., 20., 20.);
  LogicalVolume world("world", &worldUnplaced);

  UnplacedSExtruVolume sxtruv(nvert, x, y, -20, 20);
  LogicalVolume sxtruVol("xtru", &sxtruv);
  if (!tsl) {
    world.PlaceDaughter("extruded", &sxtruVol, &placement);
    std::cout << "Benchmarking simple extruded polygon (SExtru) having " << nvert << " vertices and " << nsect
              << " sections\n";
  }

  UnplacedExtruded xtru(nvert, vertices, nsect, sections);
  LogicalVolume xtruVol("xtru", &xtru);

  if (tsl) {
    world.PlaceDaughter("extruded", &xtruVol, &placement);
    std::cout << "Benchmarking extruded polygon having " << nvert << " vertices and " << nsect << " sections\n";
  }

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetVerbosity(2);
  tester.SetRepetitions(nrep);
  tester.SetPointCount(npoints);
  tester.SetToInBias(0.8);
  tester.SetPoolMultiplier(1);
  delete[] x;
  delete[] y;

  //  tester.RunToInBenchmark();
  //  tester.RunToOutBenchmark();
  return tester.RunBenchmark();
#else
  return 0;
#endif
}