#----------------------------------------------------------------------------
# Setup the project
#

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(USolids)

option(USOLIDS_USE_UVECTOR3
       "Use the UVector3 implementation instead of Vector3D." OFF)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

if (NOT USOLIDS_DIR)
  include_directories(${PROJECT_SOURCE_DIR}/include)
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.icc)

# avoid "file has no symbols" warning
if(NOT USOLIDS_USE_UVECTOR3)
  list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/UVector3.cc)
endif()

source_group("Header Files" FILES ${headers})
source_group("Source Files" FILES ${sources})

add_library(${PROJECT_NAME} STATIC ${sources} ${headers})
set_target_properties( ${PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

install(FILES ${headers} DESTINATION include)
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
