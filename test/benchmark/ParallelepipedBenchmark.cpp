#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "volumes/Parallelepiped.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"

using namespace vecgeom;

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 16384);
  OPTION_INT(nrep, 1);
/*
  OPTION_DOUBLE(dx, 20.);
  OPTION_DOUBLE(dy, 30.);
  OPTION_DOUBLE(dz, 40.);
  OPTION_DOUBLE(alpha, 30. / 180. * kPi);
  OPTION_DOUBLE(theta, 15. / 180. * kPi);
  OPTION_DOUBLE(phi, 30. / 180. * kPi);
*/

  OPTION_DOUBLE(dx, 7.5);
  OPTION_DOUBLE(dy, 10.);
  OPTION_DOUBLE(dz, 15.);
  OPTION_DOUBLE(alpha, 10. / 180. * kPi);
  OPTION_DOUBLE(theta, 20. / 180. * kPi);
  OPTION_DOUBLE(phi, 5. / 180. * kPi);

/*
  const double &deg = kDegToRad;

  OPTION_DOUBLE(dx, 7.5);
  OPTION_DOUBLE(dy, 10.);
  OPTION_DOUBLE(dz, 15.);
  OPTION_DOUBLE(alpha, 10.*deg);
  OPTION_DOUBLE(theta, 20.*deg);
  OPTION_DOUBLE(phi, 5.*deg);
*/



  UnplacedBox worldUnplaced = UnplacedBox(dx * 4, dy * 4, dz * 4);
  UnplacedParallelepiped paraUnplaced(dx, dy, dz, alpha, theta, phi);
  LogicalVolume world("w0rld", &worldUnplaced);
  LogicalVolume para("p4r4", &paraUnplaced);
  Transformation3D placement = Transformation3D(5, 5, 5);
  world.PlaceDaughter(&para, &placement);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetVerbosity(3);
  tester.SetPoolMultiplier(1);
  tester.SetRepetitions(nrep);
  tester.SetPointCount(npoints);
  return tester.RunBenchmark();
}
