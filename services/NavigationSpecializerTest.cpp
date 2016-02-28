/*
 * NavigationSpecializerTest.cpp
 *
 *  Created on: 11.09.2015
 *      Author: swenzel
 */

#include "services/NavigationSpecializer.h"
#include "management/RootGeoManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace vecgeom;

int main(int argc, char * argv[]){
  if (argc < 3) {
    std::cerr << "usage: " << argv[0] << " geometryfile.root volumename [--basenav BasicNavigator] [--loopunroll]\n";
    return 1;
  }
  vecgeom::NavigationSpecializer specializer;
  for (auto i = 3; i < argc; i++) {
     if (strcmp(argv[i], "--basenav") == 0) {
       std::cerr << "setting a basenav\n";
       specializer.SetBaseNavigator(argv[i + 1]);
     }
   }

  vecgeom::RootGeoManager::Instance().LoadRootGeometry(argv[1]);

  for (auto i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--loopunroll") == 0)
      specializer.EnableLoopUnrolling();
  }



  std::ofstream outputfile;
  std::stringstream outputname;
  outputname << argv[2] << "Navigator.h";
  outputfile.open(outputname.str());
  specializer.ProduceSpecializedNavigator(vecgeom::GeoManager::Instance().FindLogicalVolume(argv[2]), outputfile);
  outputfile.close();
  return 0;
}
