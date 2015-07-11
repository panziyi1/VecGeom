#----------------------------------------------------------------------------
# Setup the project
#

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(USolids)

option(USOLIDS_USE_UVECTOR3
       "Use the UVector3 implementation instead of Vector3D." OFF)

if (NOT USOLIDS_DIR)
  if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
  endif()
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.icc)

if (NOT USOLIDS_DIR)
  include_directories(${PROJECT_SOURCE_DIR}/include)
  find_package(VecGeom REQUIRED)
  include_directories(${VECGEOM_INCLUDE_DIR})
endif()

# avoid "file has no symbols" warning
if(NOT USOLIDS_USE_UVECTOR3)
  list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/UVector3.cc)
endif()

source_group("Header Files" FILES ${headers})
source_group("Source Files" FILES ${sources})

if (NOT USOLIDS_DIR)
  set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -std=c++11)
endif()

add_library(${PROJECT_NAME} STATIC ${sources} ${headers})
set_target_properties( ${PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

install(FILES ${headers} DESTINATION include)
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
