#!/bin/bash

# defaults

BACKEND=${BACKEND:-"-DVc=ON -DVECGEOM_BACKEND=Vc"}
BUILD_TYPE=${BUILD_TYPE:-"Release"}
SRCDIR=${SRCDIR:-${HOME}/src/vecgeom}
DESTDIR=${DESTDIR:-${PWD}}

BENCHMARK="ON"
VALIDATION="OFF"
BUILD_TESTING="ON"
GEANT4="ON"
ROOT="ON"
NO_SPECIALIZATION="ON"

# process options

for option in $@; do
case ${option} in
	# compilers
	icc|ICC|intel)
	export CC=icc CXX=icpc
	;;

	gcc|GCC|GNU)
	export CC=gcc CXX=g++
	;;

	clang|Clang)
	export CC=clang CXX=clang++
	;;

	# backends
	scalar|Scalar)
	BACKEND="-DVc=OFF -DVECGEOM_BACKEND=Scalar"
	;;

	vc|Vc|VC)
	BACKEND="-DVc=ON -DVECGEOM_BACKEND=Vc"
	;;

	# other options
	cuda|CUDA)
	CUDA="-DVECGEOM_ENABLE_CUDA=ON -DVECGEOM_NO_SPECIALIZATION=ON -DCUDA_VOLUME_SPECIALIZATION=OFF"
	;;

	test|ctest)     BUILD_TESTING="ON"  ;;
	notest|noctest) BUILD_TESTING="OFF" ;;

	bench|benchmark)     BENCHMARK="ON"  ;;
	nobench|nobenchmark) BENCHMARK="OFF" ;;

	validation)   VALIDATION="ON"  ;;
	novalidation) VALIDATION="OFF" ;;

	specialized)   NO_SPECIALIZATION="OFF" ;;
	unspecialized) NO_SPECIALIZATION="ON"  ;;

	geant4)    GEANT4="ON"  ;;
	nogeant4)  GEANT4="OFF" ;;

	root)      ROOT="ON"  ;;
	noroot)    ROOT="OFF" ;;
esac
done

echo -------------------------------------------------------------
echo
echo "Using CMake command:"
echo "cmake ${SRCDIR} -DCMAKE_INSTALL_PREFIX=${DESTDIR}          "
echo "    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${BACKEND} ${CUDA}    "
echo "    -DVECGEOM_ROOT=${ROOT} -DVECGEOM_GEANT4=${GEANT4}                      "
echo "    -DVECGEOM_TEST_BENCHMARK=${BENCHMARK} -DBUILD_TESTING=${BUILD_TESTING}"
echo "    -DVECGEOM_TEST_VALIDATION=${VALIDATION} -DVECGEOM_NO_SPECIALIZATION=${NO_SPECIALIZATION}"
echo
echo -------------------------------------------------------------

cmake ${SRCDIR} -DCMAKE_INSTALL_PREFIX=${DESTDIR}          \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${BACKEND} ${CUDA}    \
    -DVECGEOM_ROOT=${ROOT} -DVECGEOM_GEANT4=${GEANT4}                      \
    -DVECGEOM_TEST_BENCHMARK=${BENCHMARK} -DBUILD_TESTING=${BUILD_TESTING}              \
    -DVECGEOM_TEST_VALIDATION=${VALIDATION} -DVECGEOM_NO_SPECIALIZATION=${NO_SPECIALIZATION}
