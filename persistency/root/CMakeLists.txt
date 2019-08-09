# CMAKE for root_persistency

cmake_minimum_required (VERSION 3.3)
project (root_persistency)

message(STATUS "root_persistency CMAKE_GENERATOR is ${CMAKE_GENERATOR}")
message(STATUS "root_persistency CMAKE_MAKE_PROGRAM is ${CMAKE_MAKE_PROGRAM}")

option(ROOT_PERSISTENCY_TESTING "enable testing of ROOT persitency" ON)

set(VECGEOM_LIBRARIES vecgeom)
get_target_property(VECGEOM_INCLUDE_DIR vecgeom INCLUDE_DIRECTORIES)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) #see https://root.cern.ch/how/integrate-root-my-project-cmake
include(${ROOT_USE_FILE})

set(VECGEOM_DEFINITIONS_LIST "${VECGEOM_DEFINITIONS} --debug")
separate_arguments(VECGEOM_DEFINITIONS_LIST)

#generating the dictionary
REFLEX_GENERATE_DICTIONARY(RootPClasses include/Classes.h SELECTION xml/Selection.xml OPTIONS ${VECGEOM_DEFINITIONS_LIST} )
add_library(RootPClassesDict SHARED RootPClasses.cxx)
set_target_properties(RootPClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})
target_link_libraries(RootPClassesDict ${VECGEOM_LIBRARIES} ${ROOT_LIBRARIES})

#copying the dictionary to installation_dir/lib
install(TARGETS RootPClassesDict DESTINATION lib)
install(FILES ${PROJECT_BINARY_DIR}/RootPClasses_rdict.pcm DESTINATION lib)

# adding test directory
if(ROOT_PERSISTENCY_TESTING)
    add_subdirectory(test)
endif()