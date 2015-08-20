cmake_minimum_required(VERSION 2.8.8)
project(VecGeom)

enable_language(CXX)

set(VECGEOM ON)

################################################################################
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)

# Setting the C++ compiler to use
set("${CMAKE_CXX_COMPILER_ID}" TRUE CACHE STRING "C++ Compiler")
if (NOT GNU AND NOT Clang AND NOT Intel)
  message(WARNING "Unsupported compiler. Build will likely fail.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

include(MacroUtilities)
################################################################################

# Configuration options

option(CUDA "Enable compilation for CUDA." OFF)

option(CUDA_VOLUME_SPECIALIZATION "Use specialized volumes for CUDA." ON)

option(ROOT "Include ROOT." OFF)

option(USOLIDS "Include USolids." OFF)

option(GEANT4 "Include Geant4." OFF)

option(BENCHMARK "Enable performance comparisons." OFF)

option(NO_SPECIALIZATION "Disable specialization of volumes." ON)

option(COVERAGE_TESTING "Enable coverage testing flags." OFF)

option(CTEST "Enable CTest when building." ON)

option(VALIDATION "Enable validation tests from CMS geometry." OFF)

option(PLANESHELL "Enable the use of PlaneShell class for the trapezoid." ON)

option(DISTANCE_DEBUG "Enable comparison of calculated distances againt ROOT/Geant4 behind the scenes" OFF)

option(DATA_DOWNLOAD "Enable downloading of data for tests" OFF)

set(_vector native sse2 sse3 ssse3 sse4.1 sse4.2 avx avx2 mic avx512)

enum_option(VECGEOM_VECTOR
      DOC "Vector instruction set to be used"
      TYPE STRING
      VALUES ${_vector}
      CASE_INSENSITIVE
    )

set(_backend default Scalar Vc Cilk micvec)

enum_option(BACKEND
            DOC "Backend infrastructure to be used"
            TYPE STRING
            VALUES ${_backend}
            CASE_INSENSITIVE
           )

option(VC_ACCELERATION "Internal vectorization using Vc." OFF)

option(VTUNE "Enable the use of profiler Vtune" OFF)

option(OFFLOAD_MODE "Enable offload mode for Intel Xeon Phi." OFF)

if (USOLIDS)
  set(USolids ON)
endif()

if (GEANT4)
#  if (BENCHMARK)
#    set(Geant4 ON)
#  else()
#    message(WARNING
#      "Enabling Geant4 without benchmarks has no effect and will be disabled.")
#    set(GEANT4 OFF)
#    set(Geant4 OFF)
#  endif()
   set(Geant4 ON)
endif()

if (VC_ACCELERATION)
  set (Vc ON)
endif()

if (CTEST)
  set (DATA_DOWNLOAD ON)
endif()

################################################################################

# Set default C++ compiler flags for each build type
# C++ standard set to c++11 by default
# some compilers do not default to C++11 standard yet
if (NOT CUDA)
  set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -std=c++11")
endif()

if (GNU)
  set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -Wall -fPIC")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -Weffc++ -ggdb -O0")

  else()
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -finline-limit=10000000 -ffast-math -ftree-vectorize")
  endif()
  if (Vc)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -fabi-version=6")
  endif()
elseif(Intel)
  set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -Wall -fPIC")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -g -O0")
  else()
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -fno-alias")
  endif()
else(Clang)
  set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -Wall -fPIC")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -Weffc++ -ggdb -O0")
  else()
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -ffast-math -ftree-vectorize")
  endif()
endif()

if (CUDA)
if (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  set(VECGEOM_NVCC_COMPILATION_FLAGS -g -G ${VECGEOM_NVCC_COMPILATION_FLAGS})
endif()
if (CMAKE_BUILD_TYPE MATCHES MinSizeRel)
  set(VECGEOM_NVCC_COMPILATION_FLAGS -use_fast_math ${VECGEOM_NVCC_COMPILATION_FLAGS})
endif()
endif()

string(TOLOWER ${VECGEOM_VECTOR} _arch_lo)
string(TOUPPER ${VECGEOM_VECTOR} _arch_up)

# Set C++ SIMD compiler flags.

if (${_arch_lo} MATCHES mic)
  if (Intel)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -m${_arch_lo}")
    set(MIC TRUE)
  else()
    message(FATAL_ERROR "MIC native or offload compilation requires Intel Compiler.")
  endif()
elseif (${_arch_lo} MATCHES native)
  if (Intel)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -xHost")
  else()
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -march=${_arch_lo}")
  endif()
  set(VECGEOM_NVCC_COMPILATION_FLAGS -g -G ${VECGEOM_NVCC_COMPILATION_FLAGS})
else()
  if (Intel)
    if (${_arch_lo} MATCHES avx2)
      set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -xCORE-${_arch_up}")
    else()
      set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -x${_arch_up}")
    endif()
  else()
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -m${_arch_lo}")
  endif()
endif()
message(STATUS "Compiling for ${_arch_up} SIMD achitecture")

# Set C++ compiler definition flags. Will also be forwarded by CUDA when compiling C++.

if ((NOT PLANESHELL))
  message(STATUS "Disabling PlaneShell class")
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_PLANESHELL_DISABLE")
endif()


if (CUDA)
  set(CUDA_ARCH 20 CACHE STRING "CUDA device architecture.")
  set(CUDA_ARCH "-arch=sm_${CUDA_ARCH}")
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_CUDA")
  if (CUDA_VOLUME_SPECIALIZATION)
    set(VECGEOM_DEFINITION_FLAGS
        "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_CUDA_VOLUME_SPECIALIZATION")
  else()
    set(VECGEOM_DEFINITION_FLAGS
        "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_CUDA_NO_VOLUME_SPECIALIZATION")
  endif()
  set(VECGEOM_NVCC_COMPILATION_FLAGS ${VECGEOM_NVCC_COMPILATION_FLAGS} 
      -std=c++11
      -Xcudafe "--diag_suppress=code_is_unreachable"
      -Xcudafe "--diag_suppress=initialization_not_reachable")
  set(VECGEOM_NVCC_CONFIGS
        DEBUG -g -G
        RELEASE  -O3 -use_fast_math
        RELWITHDEBINFO -g -G -O3 -use_fast_math
        MINSIZEREL -use_fast_math )
endif()

################################################################################

# Determine backend

string(TOUPPER ${BACKEND} _backend_up)

if (_backend_up MATCHES VC)
  set(Vc ON)
  message(STATUS "Configuring with \"Vc\" vector backend")
elseif (_backend_up MATCHES SCALAR)
  set(Scalar ON)
  message(STATUS "Configuring with \"Scalar\" vector backend")
elseif (_backend_up MATCHES CILK)
  if (NOT Intel)
    message(FATAL_ERROR "Must use Intel C++ compiler (icpc) for Cilk+ backend.")
  endif()
  set(Cilk ON)
  message(STATUS "Configuring with \"Cilk+\" vector backend")
elseif (_backend_up MATCHES MICVEC)
  if (NOT Intel OR NOT MIC)
    message(FATAL_ERROR "Must use Intel C++ compiler (icpc) with MIC SIMD for MICVEC backend.")
  endif()
  set(MICVEC ON)
  message(STATUS "Configuring with \"Micvec\" vector backend")
elseif (_backend_up MATCHES DEFAULT)
  set(Vc ON)
  message(WARNING "Backend not specified. Defaulting to Vc.")
else()
  message(FATAL_ERROR "Value '${BACKEND}' for variable BACKEND is not allowed. ")
endif()

# Backend specific
if (Scalar)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_SCALAR")
endif()
if (Vc)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_VC")
endif()
if (Cilk)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_CILK")
endif()
if (MICVEC)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_MICVEC")
endif()

# Enable pretty diagnostics coloring if gcc 4.9+ is being used
if (GNU)
  # gcc needs -ldl to be added explicitly
  set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} dl)
  set(VECGEOM_ERROR_LIMIT 20 CACHE STRING
      "Limit number of errors output by diagnostics.")
  set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -fmax-errors=${VECGEOM_ERROR_LIMIT}")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  if (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} -fdiagnostics-color=auto")
  endif()
endif()

# Coverage testing support
if (COVERAGE_TESTING)
  if (GNU)
    set(VECGEOM_CXX_FLAGS
        "${VECGEOM_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} gcov)
  endif()
  if (Clang)
    set(VECGEOM_CXX_FLAGS "${VECGEOM_CXX_FLAGS} --coverage")
    set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} profile_rt)
  endif()
  if (Intel)
    message(FATAL_ERROR "Coverage testing not supported for icc.")
  endif()
endif()

# Forward CMake options
if (VC_ACCELERATION)
  set(VECGEOM_DEFINITION_FLAGS
      "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_VC_ACCELERATION")
endif()
if (BENCHMARK)
  set(VECGEOM_DEFINITION_FLAGS
      "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_BENCHMARK")
endif()
if (NO_SPECIALIZATION)
  set(VECGEOM_DEFINITION_FLAGS
      "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_NO_SPECIALIZATION")
endif()
if (ROOT)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_ROOT")
endif()
if (USolids)
  set(VECGEOM_DEFINITION_FLAGS
      "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_USOLIDS")
endif()
if (Geant4)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_GEANT4")
endif()
if (DISTANCE_DEBUG)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_DISTANCE_DEBUG")
endif()
if (OFFLOAD_MODE)
  if (NOT Intel)
    message(FATAL_ERROR "Must use Intel C++ compiler (icc) for Xeon Phi Offload Mode.")
  endif()
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DOFFLOAD_MODE")
  message(STATUS "Preparing for a offload compiler mode.")
  set(MICVEC ON)
  set(VTUNE OFF)
endif()
if (VTUNE)
  set(VECGEOM_DEFINITION_FLAGS "${VECGEOM_DEFINITION_FLAGS} -DVECGEOM_VTUNE")
endif()
#################################################################################

include_directories(AFTER ${CMAKE_SOURCE_DIR})

# Backend configuration

if (VC_ACCELERATION OR Vc)

  find_package(Vc REQUIRED)
  if (MIC)
    set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} ${Vc_MIC_LIBRARIES})
  else()
    set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} ${Vc_LIBRARIES})
  endif()
  include_directories(AFTER ${Vc_INCLUDE_DIR})

  if (Vc)
    set(SRC_EXTERNAL ${SRC_EXTERNAL} source/backend/vc/Backend.cpp)
  endif()

endif()

if (MICVEC)
  set(SRC_EXTERNAL ${SRC_EXTERNAL} source/backend/mic/Backend.cpp)
endif()

if (Cilk)
  set(SRC_EXTERNAL ${SRC_EXTERNAL} source/backend/cilk/Backend.cpp)
endif()

if (CUDA)

  find_package(CUDA REQUIRED)
  include(FindCUDA)

  set(CUDA_SEPARABLE_COMPILATION ON)

  set(SRC_EXTERNAL ${SRC_EXTERNAL} source/CudaManager.cpp source/backend/cuda/Interface.cpp)
  set(SRC_CUDA ${SRC_CUDA}
      source/CudaManager.cu)

  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${VECGEOM_NVCC_COMPILATION_FLAGS})
endif()

### CXXFLAGS needs to be set early, so that VECGEOM_GEANT4, VECGEOM_ROOT,
### and VECGEOM_USOLIDS are already defined when compiling USolids, or it breaks

# Pass flags to compilers
# We may have addition in compiler flags from the above included packages
# We don't want to have trailing CMAKE_CXX_FLAGS_"BUILD_TYPE" options on the command line
# but rather have it as first set of argument.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}} ${VECGEOM_CXX_FLAGS} ${VECGEOM_DEFINITION_FLAGS}")
set(VECGEOM_COMPILATION_FLAGS ${VECGEOM_CXX_FLAGS})
#set(VECGEOM_CXX_FLAGS ${CMAKE_CXX_FLAGS_${_build_type}})
#unset(CMAKE_CXX_FLAGS_${_build_type} CACHE)

message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
if (CUDA)
  message(STATUS "Compiling with NVCC flags: ${CUDA_NVCC_FLAGS}")
endif()

# Uncomment line below for building/installing VecCore
# add_subdirectory(VecCore)

if (ROOT)

  find_package(ROOT REQUIRED)

  include_directories(AFTER ${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR})
  set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL}
      ${ROOT_LIBRARIES} Geom)
  set(SRC_EXTERNAL ${SRC_EXTERNAL}
      source/PlacedRootVolume.cpp
      source/RootGeoManager.cpp
      source/UnplacedRootVolume.cpp
      source/Visualizer.cpp)

endif()

if (USolids)

  find_package(USolids REQUIRED)
  include_directories(AFTER ${USOLIDS_INCLUDE_DIR})

  set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} ${USOLIDS_LIBRARIES})

endif()

if (Geant4)

  find_package(Geant4 REQUIRED)
  include_directories(AFTER ${Geant4_INCLUDE_DIRS})
  set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL}
      ${Geant4_LIBRARIES})
  message(STATUS "Geant4 found in: ${Geant4_INCLUDE_DIRS}")

endif()

if (VTUNE)

  find_package(VTUNE REQUIRED)
  include_directories(AFTER ${VTUNE_INCLUDE_DIR})
  set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} ${VTUNE_LIBRARIES} -lpthread -ldl)

endif()


#include_directories(AFTER ${CMAKE_SOURCE_DIR}/VecCore/inc )
################################################################################


set(SRC_CPP_RELATIVE
  source/RNG.cpp
  source/AOS3D.cpp
  source/LogicalVolume.cpp
  source/PlacedCone.cpp
  source/PlacedBox.cpp
  source/PlacedTube.cpp
  source/PlacedTorus.cpp
  source/PlacedTorus2.cpp
  source/PlacedTrd.cpp
  source/PlacedParallelepiped.cpp
  source/PlacedParaboloid.cpp
  source/PlacedPolyhedron.cpp
  source/PlacedTrapezoid.cpp
  source/PlacedVolume.cpp
  source/Planes.cpp
  source/Quadrilaterals.cpp
  source/Rectangles.cpp
  source/SOA3D.cpp
  source/Transformation3D.cpp
  source/UnplacedCone.cpp
  source/UnplacedBox.cpp
  source/UnplacedTube.cpp
  source/UnplacedTorus.cpp
  source/UnplacedTorus2.cpp
  source/UnplacedTrd.cpp
  source/UnplacedParaboloid.cpp
  source/UnplacedParallelepiped.cpp
  source/UnplacedPolyhedron.cpp
  source/UnplacedTrapezoid.cpp
  source/UnplacedVolume.cpp
  source/Vector.cpp
  source/NavigationState.cpp
  source/SimpleNavigator.cpp

  source/UnplacedOrb.cpp
  source/PlacedOrb.cpp

  source/UnplacedSphere.cpp
  source/PlacedSphere.cpp

  source/UnplacedBooleanVolume.cpp
  source/PlacedBooleanVolume.cpp

  source/UnplacedPolycone.cpp
  source/PlacedPolycone.cpp
  source/Wedge.cpp
  source/ABBoxManager.cpp

  source/USolidsInterfaceHelper.cpp
)

# avoid "file has no symbols" warning
if(NOT CUDA)
  list(REMOVE_ITEM SRC_CPP_RELATIVE source/RNG.cpp)
  list(REMOVE_ITEM SRC_CPP_RELATIVE source/AOS3D.cpp)
  list(REMOVE_ITEM SRC_CPP_RELATIVE source/SOA3D.cpp)
  list(REMOVE_ITEM SRC_CPP_RELATIVE source/Vector.cpp)
endif()

foreach(SRC ${SRC_CPP_RELATIVE})
  set(SRC_CPP ${SRC_CPP} ${CMAKE_SOURCE_DIR}/${SRC})
endforeach()

# Copy all source files to .cu-files in order for NVCC to compile them as CUDA
# code and not regular C++ files.

if (CUDA)

  foreach(SRC_FILE ${SRC_CPP})

    get_filename_component(SRC_FILENAME ${SRC_FILE} NAME_WE)

    if (NOT "${SRC_FILENAME}" STREQUAL "ABBoxManager")
      ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC_FILE}
            ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
        DEPENDS ${SRC_FILE}
      )

      set(
        SRC_CUDA ${SRC_CUDA}
        ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      )
    endif()
  endforeach()

endif()

# Add files that won't be compiled for CUDA

if (USolids)
  set(SRC_CPP ${SRC_CPP}
    test/shape_tester/ShapeTester.cpp
  )
endif()
if (BENCHMARK)
  set(SRC_CPP ${SRC_CPP}
    source/benchmarking/BenchmarkResult.cpp
    source/benchmarking/Benchmarker.cpp
    source/benchmarking/VolumePointers.cpp
  )
  set(SRC_CUDA ${SRC_CUDA}
    source/benchmarking/Benchmarker.cu
    source/benchmarking/NavigationBenchmarker.cu
  )
endif()
set(SRC_CPP ${SRC_CPP} ${SRC_EXTERNAL}
  source/GeoManager.cpp
  source/ShapeDebugger.cpp
  source/benchmarking/NavigationBenchmarker.cpp
  source/CppExporter.cpp
)

# some CUDA kernel code in "userspace"
if (CUDA)
    set(SRC_USER_CUDA_CPP
        userexamples/src/TestNavigationStatePool.cu
    )
endif()

# Compile cpp-files for C++11

if (CUDA)
  set_source_files_properties(${SRC_CPP} PROPERTIES COMPILE_FLAGS -std=c++11)
endif()

################################################################################
# Build libraries
################################################################################

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

if (USolids)
  if (USOLIDS_EMBEDDED)
    add_subdirectory(${USOLIDS_DIR})
  endif()
endif()

add_library(vecgeom STATIC ${SRC_CPP})
target_link_libraries(vecgeom ${VECGEOM_LIBRARIES}
                      ${VECGEOM_LIBRARIES_EXTERNAL})

set(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} vecgeom)
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(VECGEOM_LIBRARIES_EXTERNAL ${VECGEOM_LIBRARIES_EXTERNAL} rt)
endif()

if (CUDA)
  set_source_files_properties(
    ${SRC_EXECUTABLES}
    PROPERTIES COMPILE_FLAGS -std=c++11
  )
endif()

# build the CUDA version of the library
if (CUDA)
  cuda_add_library(
    vecgeomcuda
    ${SRC_CUDA}
    OPTIONS ${CUDA_ARCH}
    ${VECGEOM_NVCC_CONFIGS}
  )
  target_link_libraries(vecgeomcuda ${VECGEOM_LIBRARIES}
                        ${VECGEOM_LIBRARIES_EXTERNAL})
  set(VECGEOM_LIBRARIES ${VECGEOM_LIBRARIES} vecgeomcuda)
endif()


# build the CUDA user-code library
if (CUDA)
  cuda_add_library(
    cudauserlib
    ${SRC_USER_CUDA_CPP}
    OPTIONS ${CUDA_ARCH}
  )
  target_link_libraries(cudauserlib ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL})
  set(USERKERNELLIB cudauserlib)
endif()
################################################################################

# Define executables

set(TEST_EXECUTABLES_CORE
  ${CMAKE_SOURCE_DIR}/test/core/ContainerTest.cpp
  ${CMAKE_SOURCE_DIR}/test/core/create_geometry.cpp
  ${CMAKE_SOURCE_DIR}/test/core/testVectorSafety.cpp
  ${CMAKE_SOURCE_DIR}/test/core/PlanesTest.cpp
  ${CMAKE_SOURCE_DIR}/test/core/QuadrilateralTest.cpp
#  ${CMAKE_SOURCE_DIR}/test/core/SOATest.cpp
  ${CMAKE_SOURCE_DIR}/test/core/Transformation3DTest.cpp
#  ${CMAKE_SOURCE_DIR}/test/core/boolminustest.cpp
#  ${CMAKE_SOURCE_DIR}/test/core/boolminustest2.cpp
  ${CMAKE_SOURCE_DIR}/test/core/PhiWedgeTest.cpp
  ${CMAKE_SOURCE_DIR}/test/core/ThetaConeTest.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestVecGeomPolycone.cpp
)


set(TEST_EXECUTABLES_ROOT
  ${CMAKE_SOURCE_DIR}/test/core/TestNavigationStatePool.cpp
  ${CMAKE_SOURCE_DIR}/test/root/root_geometry.cpp
  ${CMAKE_SOURCE_DIR}/test/root/complex_test1.cpp
  ${CMAKE_SOURCE_DIR}/test/root/E03Test.cpp
  ${CMAKE_SOURCE_DIR}/test/root/ImportFromRootFileTest.cpp
  ${CMAKE_SOURCE_DIR}/test/root/ImportTGeoPgon.cpp
  ${CMAKE_SOURCE_DIR}/test/root/TestExportToROOT.cpp
  ${CMAKE_SOURCE_DIR}/test/root/ExitingOrEntering.cpp
)

# Tests below fail to build when both ROOT and USolids are enabled
if(NOT USolids)
  set(TEST_EXECUTABLES_CORE
    ${TEST_EXECUTABLES_CORE}
    ${CMAKE_SOURCE_DIR}/test/core/BitSetTest.cpp
  )
  # higher level benchmarks
  set(TEST_EXECUTABLES_ROOT
    ${TEST_EXECUTABLES_ROOT}
    ${CMAKE_SOURCE_DIR}/test/globalbenchmarks/LocatePointsBenchmark.cpp
    ${CMAKE_SOURCE_DIR}/test/globalbenchmarks/NavigationBenchmark.cpp
    ${CMAKE_SOURCE_DIR}/test/globalbenchmarks/TraceTrack.cpp
    ${CMAKE_SOURCE_DIR}/test/globalbenchmarks/XRayBenchmarkFromROOTFile.cpp
  )
endif()

# separate list for ROOT UNIT tests
set(TEST_UNITTESTEXECUTABLES_ROOT
  ${CMAKE_SOURCE_DIR}/test/root/complex_test1.cpp
  ${CMAKE_SOURCE_DIR}/test/root/E03Test.cpp
  ${CMAKE_SOURCE_DIR}/test/root/TestExportToROOT.cpp
  ${CMAKE_SOURCE_DIR}/test/root/ImportTGeoPgon.cpp
  ${CMAKE_SOURCE_DIR}/test/core/TestNavigationStatePool.cpp
)

if(BENCHMARK)
  set(TEST_EXECUTABLES_ROOT
  ${TEST_EXECUTABLES_ROOT}
  ${CMAKE_SOURCE_DIR}/test/root/BenchmarkShapeFromROOTFile.cpp
  ${CMAKE_SOURCE_DIR}/test/root/BenchmarkShapeFromROOTFile_WithVisualization.cpp
  ${CMAKE_SOURCE_DIR}/test/root/CompareDistances.cpp
)
endif()

set(TEST_EXECUTABLES_VISUALIZATION
  ${CMAKE_SOURCE_DIR}/test/visualization/VisualizeTrap.cpp
  ${CMAKE_SOURCE_DIR}/test/visualization/VisualizeTube.cpp
  ${CMAKE_SOURCE_DIR}/test/visualization/DebugTube.cpp
  ${CMAKE_SOURCE_DIR}/test/visualization/DebugPolyhedron.cpp
)

set(TEST_EXECUTABLES_USOLIDS
  ${CMAKE_SOURCE_DIR}/test/usolids/CompatibilityTest.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestBox.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestCons.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestTube.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestTrd.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestTrap.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestPolycone.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestPolyhedra.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestTet.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestOrb.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestSphere.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/TestBoxlikeTrapezoid.cpp
  ${CMAKE_SOURCE_DIR}/test/unit_tests/usolidTest.cpp
)

set(TEST_EXECUTABLES_BENCHMARK
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParaboloidBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParaboloidScriptBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParallelepipedBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/PolyhedronBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TubeBenchmark.cpp
 # ${CMAKE_SOURCE_DIR}/test/benchmark/BoxMinusHollowTubeBenchmark.cpp
 # ${CMAKE_SOURCE_DIR}/test/benchmark/BoxMinusTubeSegBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxMinusHollowTubeBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxUnionHollowTubeBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxIntersectHollowTubeBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxUnionBoxBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/3LevelSubtractionBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TorusBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrapezoidBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrapezoidBenchmarkScript.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/OrbBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/SphereBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrdBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ConeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/PolyconeBenchmark.cpp
)
set(TEST_EXECUTABLES_CUDA
  # Files go here
)
set(TEST_EXECUTABLES_SHAPETESTER
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testBox.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testOrb.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testSphere.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testCone.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testTube.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testTrd.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testTrapezoid.cpp
  # ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testParaboloid.cpp
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testPolycone.cpp 
)

if(ROOT) # add shape testers that depend on ROOT
  set(TEST_EXECUTABLES_SHAPETESTER
  ${TEST_EXECUTABLES_SHAPETESTER}
  ${CMAKE_SOURCE_DIR}/test/shape_tester/shape_testFromROOTFile.cpp )
endif()

set(TEST_EXECUTABLES
  ${TEST_EXECUTABLES_CORE}
  ${TEST_EXECUTABLES_ROOT}
  ${TEST_EXECUTABLES_VISUALIZATION}
  ${TEST_EXECUTABLES_USOLIDS}
  ${TEST_EXECUTABLES_BENCHMARK}
  ${TEST_EXECUTABLES_CUDA}
  ${TEST_EXECUTABLES_SHAPETESTER}
)

if (CUDA)
  set_source_files_properties(
    ${TEST_EXECUTABLES}
    PROPERTIES COMPILE_FLAGS -std=c++11
  )
endif()

# Build executables

if (CTEST)
  message(STATUS "Testing with CTest enabled.")
  enable_testing()
endif()

macro(build_executables EXECUTABLES)
  foreach(_exec ${EXECUTABLES})
    get_filename_component(_target ${_exec} NAME_WE)
    add_executable(${_target} ${_exec})
    set_target_properties(${_target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    target_link_libraries(${_target} ${VECGEOM_LIBRARIES}
                          ${VECGEOM_LIBRARIES_EXTERNAL}
                          ${USERKERNELLIB})
  endforeach()
endmacro()

macro(add_to_ctest EXECUTABLES)
  foreach(_exec ${EXECUTABLES})
    get_filename_component(_target ${_exec} NAME_WE)
    add_test(NAME ${_target} COMMAND ${_target})
  endforeach()
endmacro()

build_executables("${TEST_EXECUTABLES_CORE}")
add_to_ctest("${TEST_EXECUTABLES_CORE}")

if (BENCHMARK)
  build_executables("${TEST_EXECUTABLES_BENCHMARK}")
endif()

if (ROOT)
  build_executables("${TEST_EXECUTABLES_ROOT}")
  build_executables("${TEST_EXECUTABLES_VISUALIZATION}")
  add_to_ctest("${TEST_UNITTESTEXECUTABLES_ROOT}")
endif()

if (USolids)
  build_executables("${TEST_EXECUTABLES_USOLIDS}")
  build_executables("${TEST_EXECUTABLES_SHAPETESTER}")
  #add_to_ctest("${TEST_EXECUTABLES_USOLIDS}")
  #add_to_ctest("${TEST_EXECUTABLES_SHAPETESTER}")

  # add unit tests
  add_test(NAME USolidBox COMMAND TestBox --usolids)
  add_test(NAME VecGeomBox COMMAND TestBox --vecgeom)
  add_test(NAME USolidTube COMMAND TestTube --usolids)
  add_test(NAME VecTube COMMAND TestTube --vecgeom)

endif()

if (CUDA)
  build_executables("${TEST_EXECUTABLES_CUDA}")
#  add_to_ctest("${TEST_EXECUTABLES_CUDA}")
endif()

if (CTEST AND VALIDATION AND ROOT)
  macro(add_cmsshapevalidation_test TESTNAME SHAPEFILE)
    add_test(NAME ${TESTNAME} COMMAND bash -c "${CMAKE_SOURCE_DIR}/test/scripts/RunRandomValidation ${CMAKE_BINARY_DIR}/BenchmarkShapeFromROOTFile  ${CMAKE_SOURCE_DIR}/test/cmstestdata/cms2015.root  ${CMAKE_SOURCE_DIR}/test/cmstestdata/${SHAPEFILE}" )
  endmacro()

  # Adding various shapes tests for nightlies
  add_cmsshapevalidation_test( tubevalidation cmstubes.txt )
  add_cmsshapevalidation_test( trapvalidation cmstraps.txt )
  add_cmsshapevalidation_test( polyconevalidation cmspolycones.txt )
  add_cmsshapevalidation_test( polyhedravalidation cmspolyhedra.txt )
  add_cmsshapevalidation_test( conevalidation cmscones.txt )
  add_cmsshapevalidation_test( boxvalidation cmsboxes.txt )
  #taken out due to problems: add_cmsshapevalidation_test( booleanvalidation cmsbooleans.txt )
  #taken out due to Issue-133: add_cmsshapevalidation_test( torusvalidation cmstori.txt )
endif()



#---------------------------------------------------------------------------
# Include creation and submission to dashboard cdash.cern.ch

include(CTest)

#---------------------------------------------------------------------------
# Set data to be dowloaded

if(DATA_DOWNLOAD)
  message(STATUS "Downloading data files")

  # define a function downloading from a URL into a local file LOCALFILE
  function(FILE_DOWNLOAD FILE_URL LOCALFILE ) 
   if(APPLE)
       execute_process(COMMAND curl -o ${LOCALFILE} ${FILE_URL})
   else()  
      execute_process(COMMAND wget -q ${FILE_URL} -O ${LOCALFILE})
   endif()
  endfunction(FILE_DOWNLOAD)
  # end of function FILE DOWNLOAD

  # define a function checking md5 hashes
  # result is stored in MD5MATCHES ( 1 == true, 0 == false )
  function(CHECKMD5 FILETOCHECK EXPECTEDMD5HASH MD5MATCHES)
      if(APPLE)
          execute_process(COMMAND md5 ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
          string(LENGTH ${MD5SUM} MD5LENGTH)
          MATH(EXPR START "${MD5LENGTH} - 33")
          string(SUBSTRING ${MD5SUM} ${START} 32 MD5SUM)
      else()
          execute_process(COMMAND md5sum ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
          string(SUBSTRING ${MD5SUM} 0 32 MD5SUM)
      endif()
      if(MD5SUM STREQUAL EXPECTEDMD5HASH)
        set(${MD5MATCHES} 1 PARENT_SCOPE)
      else()
        set(${MD5MATCHES} 0 PARENT_SCOPE)
      endif()
  endfunction(CHECKMD5)

  # actual function for managing the download
  function(DOWNLOAD_IF_NOT_INSTALLED FILE_URL LOCALFILE TARGETPATH MD5HASH )
    find_file(FOUNDFILE ${LOCALFILE} ${TARGETPATH} )
    if(FOUNDFILE STREQUAL "FOUNDFILE-NOTFOUND")
        # set need download
        message(STATUS "need download of ${LOCALFILE} since not found")
        set( NEEDTODOWNLOAD 1 )
    else()
        # check md5
        message(STATUS "found existing file ${LOCALFILE}")
        CHECKMD5( ${FOUNDFILE} ${MD5HASH} MD5CORRECT )
        if( ${MD5CORRECT} STREQUAL "1" )
            # do not set download flag   
            set( NEEDTODOWNLOAD 0 )
        else( )
            # set need download
            message(STATUS "hash ${MD5HASH} not correct for file ${FOUNDFILE} ${MD5CORRECT}" )
            set( NEEDTODOWNLOAD 1 )
        endif( )
    endif()
  
    if( ${NEEDTODOWNLOAD} STREQUAL 1 )
        message(STATUS " downloading ... ")
        set(DOWNLOADLOCATION "${TARGETPATH}/${LOCALFILE}")
        FILE_DOWNLOAD( ${FILE_URL} ${DOWNLOADLOCATION} )
    else()
        message(STATUS " doing nothing ... ")
    endif()
    # in principle have to check now if download succeeded and has right MD5
    # TOBEDONE
    
    # this is annoying but we have to clear FOUNDFILE SINCE THIS IS TREATED LIKE A STATIC VARIABLE
    unset(FOUNDFILE CACHE)
  endfunction(DOWNLOAD_IF_NOT_INSTALLED)
  
  DOWNLOAD_IF_NOT_INSTALLED( "http://mnovak.web.cern.ch/mnovak/data/cms2015.root" "cms2015.root" 
                             "${CMAKE_SOURCE_DIR}/test/cmstestdata/" "fb3ede867532b372c5e6f7138d00c07e" )
  DOWNLOAD_IF_NOT_INSTALLED( "http://mnovak.web.cern.ch/mnovak/data/ExN03.root" "ExN03.root"
                             "${CMAKE_BINARY_DIR}"                   "b6b0cfdd5035117171bfe1b7f1f40c3f" )
  # add more data to download here
  #  
endif()


################################################################################

# Installation

export(PACKAGE VecGeom)

# Build
set(CONF_TYPE "build")
configure_file(VecGeomConfig.cmake.in "${CMAKE_BINARY_DIR}/VecGeomConfig.cmake" @ONLY)
# Installation
set(CONF_TYPE "install")
if (CUDA)
   set(CONF_LIBRARIES "${CONF_LIBRARY_DIR}/libvecgeomcuda.a" ${CONF_LIBRARIES}  "${CONF_LIBRARY_DIR}/libvecgeomcuda.a")
endif()
configure_file(VecGeomConfig.cmake.in "${CMAKE_BINARY_DIR}/lib/Cmake/VecGeomConfig.cmake" @ONLY)

# Install the VecGeomConfig.cmake
install(FILES "${CMAKE_BINARY_DIR}/lib/Cmake/VecGeomConfig.cmake" DESTINATION "lib/Cmake/" COMPONENT dev)

# Install headers and libraries
foreach(_folder base benchmarking management navigation volumes backend)
  install(DIRECTORY ${_folder} DESTINATION include)
endforeach()
install(TARGETS vecgeom DESTINATION lib)
if (CUDA)
  install(TARGETS vecgeomcuda DESTINATION lib)
endif()

