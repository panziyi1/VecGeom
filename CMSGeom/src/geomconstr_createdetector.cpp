#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/LogicalVolume.h"
#include "base/Transformation3D.h"
#include "management/GeoManager.h"
#include "base/Stopwatch.h"
using namespace vecgeom;extern void GenerateTransformations_part0();
extern void GenerateTransformations_part1();
extern void GenerateTransformations_part2();
extern void GenerateTransformations_part3();
extern void CreateLogicalVolumes();
extern void GeneratePlacedVolumes_part0();
extern void GeneratePlacedVolumes_part1();
extern void GeneratePlacedVolumes_part2();
extern void GeneratePlacedVolumes_part3();
extern void GeneratePlacedVolumes_part4();
extern LogicalVolume * lvol0;
extern Transformation3D * idtrans;
VPlacedVolume const * generateDetector() {
  GenerateTransformations_part0();
  GenerateTransformations_part1();
  GenerateTransformations_part2();
  GenerateTransformations_part3();
CreateLogicalVolumes();
  GeneratePlacedVolumes_part0();
  GeneratePlacedVolumes_part1();
  GeneratePlacedVolumes_part2();
  GeneratePlacedVolumes_part3();
  GeneratePlacedVolumes_part4();
VPlacedVolume const * world = lvol0->Place( idtrans ); 
return world;
}
int main(){
// function could be used like this 
 GeoManager & geom = GeoManager::Instance();
 Stopwatch timer;
 timer.Start();
 geom.SetWorld( generateDetector() );
 geom.CloseGeometry();
 timer.Stop();
 std::cerr << "loading took  " << timer.Elapsed() << " s " << std::endl;
 std::cerr << "loaded geometry has " << geom.getMaxDepth() << " levels " << std::endl;
 return 0;}
