// this file is part of VecGeom
// a service that generates a cpp function/file which instantiates a given list of
// specialized navigators with the purpose to link them into a shared library
// the service also generates a CMakeLists.txt file to facilitate this process
// started 27.2.2016; sandro.wenzel@cern.ch

#include "SpecNavServices.h"

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "usage : " << argv[0] << " LVolumeName [LVolumeName ...] \n";
    return 1;
  }

  std::vector<std::string> volumenames;

  for (int i = 1; i < argc; ++i) {
    volumenames.push_back(argv[i]);
  }

  bool success = vecgeom::SpecNavServices::GenerateNavigatorLibrary(volumenames);

  return (success) ? 0 : 1;
}
