#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(CMSGeom)

if (NOT CMSGEOM_DIR)
  if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
  endif()

  # Fix bug on APPLE, this is the default everywhere
  if(APPLE)
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
  endif()

#----------------------------------------------------------------------------
# Additional options
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

#----------------------------------------------------------------------------
# Find VecGeom geometry headers library 
#
  find_package(VecGeom REQUIRED)
  if( VECGEOM_FOUND )
    find_package(Vc)
    if ( Vc_FOUND )
        include_directories( AFTER ${Vc_INCLUDE_DIR})
    endif()
    include_directories(AFTER ${VECGEOM_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "Required VecGeom library not found")
  endif()

  include_directories("${VECGEOM_INSTALL_DIR}/include/VecCore")

  message(STATUS "VECGEOM_INSTALL_DIR is -${VECGEOM_INSTALL_DIR}-")
  message(STATUS "VECGEOM_INCLUDE_DIR is -${VECGEOM_INCLUDE_DIR}-")

#----------------------------------------------------------------------------
# Add rules to compile lib
#
# We don't want to have trailing CMAKE_CXX_FLAGS_"BUILD_TYPE" options on the command line
# but rather have it as first set of argument.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}} ${VECGEOM_COMPILATION_FLAGS} ${VECGEOM_DEFINITION_FLAGS}")
set(CMAKE_CXX_FLAGS_${_build_type} "" CACHE STRING "" FORCE)

message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

if (NOT CMSGEOM_DIR)
set_source_files_properties(${SOURCES} PROPERTIES COMPILE_FLAGS -std=c++11)
endif()

add_library(${PROJECT_NAME} SHARED ${sources})
target_link_libraries(${PROJECT_NAME} vecgeom)

if (Vc_FOUND)
  if (MIC)
    target_link_libraries(${PROJECT_NAME} ${Vc_MIC_LIBRARIES})
  else()
    target_link_libraries(${PROJECT_NAME} ${Vc_LIBRARIES})
  endif()
endif()

set_target_properties( ${PROJECT_NAME}
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)


#set(TESTEXEC ${CMAKE_SOURCE_DIR}/test/TestDynamicLoad.cpp)
#set_source_files_properties(${TESTEXEC} PROPERTIES COMPILE_FLAGS -std=c++11)

#add_executable(TestDynamicLoad ${TESTEXEC})
#if (MIC)
#  target_link_libraries(TestDynamicLoad ${VECGEOM_LIBRARIES} ${Vc_MIC_LIBRARIES})
#else()
#  target_link_libraries(TestDynamicLoad ${VECGEOM_LIBRARIES} ${Vc_LIBRARIES})
#endif()
