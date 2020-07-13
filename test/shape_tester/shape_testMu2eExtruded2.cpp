#ifndef VECGEOM_ENABLE_CUDA

#include "test/benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "VecGeom/volumes/Extruded.h"
//#include "test/core/TessellatedOrb.h"

using namespace vecgeom;
#endif

int main(int argc, char *argv[])
{
#ifndef VECGEOM_ENABLE_CUDA
  OPTION_INT(npoints, 1000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);

  ///.. Test an extruded volume from Mu2e experiment
  constexpr int nvertices = 10;
  XtruVertex2 *vertices = new XtruVertex2[nvertices];
  vertices[0]= { -5461.0, 8128.0 };
  vertices[1]= { -5003.8, 8128.0 };
  vertices[2]= { -5003.8, 6858.0 };
  vertices[3]= { 1468.12, 6858.0 };
  vertices[4]= { 3251.20, 7289.8 };
  vertices[5]= { 3251.20, 8674.1 };
  vertices[6]= { 3708.40, 8674.1 };
  vertices[7]= { 3708.40, 6771.64 };
  vertices[8]= { 1544.32, 6248.4 };
  vertices[9]= { -5461.0, 6248.4 };

  constexpr int nsections = 2;
  XtruSection* sections = new XtruSection[nsections];
  sections[0].fOrigin.Set(0, 0, -2209.8);
  sections[0].fScale = 1;
  sections[1].fOrigin.Set(0, 0,  2209.8);
  sections[1].fScale = 1;

  UnplacedExtruded mu2eXtru( nvertices, vertices, nsections, sections);
  LogicalVolume vol( &mu2eXtru );
  VPlacedVolume *placed = vol.Place();

  ShapeTester<vecgeom::VPlacedVolume> tester;
  tester.setDebug(debug);
  tester.setStat(stat);
  tester.SetMaxPoints(npoints);
  tester.SetSolidTolerance(1.e-9);
  tester.SetTestBoundaryErrors(false);
  int errCode = tester.Run(placed);

  std::cout << "Final Error count for Shape *** " << placed->GetName() << "*** = " << errCode << "\n";
  std::cout << "=========================================================" << std::endl;

#endif
  return 0;
}
