/// @file SpecNavServices.h
/// @author andrei.gheata@cern.ch
///
/// Services used for generating a specialized navigator library
/// a service that generates a cpp function/file which instantiates a given list of
/// specialized navigators with the purpose to link them into a shared library
/// the service also generates a CMakeLists.txt file to facilitate this process
/// started 27.2.2016; sandro.wenzel@cern.ch

#ifndef VECGEOM_SERVICES_LIBRARYGENERATOR_H_
#define VECGEOM_SERVICES_LIBRARYGENERATOR_H_

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

namespace vecgeom {
namespace SpecNavServices {

void GenerateHeaderIncludes(std::ostream &ss, std::vector<std::string> const &navigatornames);
void GenerateNavigatorInstantiationFunction(std::ostream &ss, std::vector<std::string> const &volumenames,
                                            std::vector<std::string> const &navigatornames);

// to generate the CMakeFile in order to compile and link this
void GenerateCMakeFile(std::ostream &ss);

// to generate the specialized navigators library
bool GenerateNavigatorLibrary(std::vector<std::string> const &volumenames);

} // namespace SpecNavServices
} // namespace vecgeom

#endif // VECGEOM_SERVICES_LIBRARYGENERATOR_H_
