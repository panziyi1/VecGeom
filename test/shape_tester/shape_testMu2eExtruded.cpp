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
  constexpr int nsections = 2;
  constexpr int nvertices = 6;
  XtruVertex2 *vertices = new XtruVertex2[nvertices];
  vertices[0]= { 13343.2, 15903.4 };
  vertices[1]= { 15614.2, 16453.6 };
  vertices[2]= { 15758.0, 15861.7 };
  vertices[3]= { 17585.5, 16306.1 };
  vertices[4]= { 17082.1, 18377.4 };
  vertices[5]= { 12984.2, 17384.5 };

  XtruSection* sections = new XtruSection[nsections];
  sections[0].fOrigin.Set(0, 0, -304.8);
  sections[0].fScale = 1;
  sections[1].fOrigin.Set(0, 0,  304.8);
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
