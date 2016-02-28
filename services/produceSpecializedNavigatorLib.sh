# a shell script which steers the process to produces a shared lib of specialized navigators
# basic usage (in build directory) run script with a ROOT geometry file and a list of logical volume names 
# author: Sandro Wenzel; 26.2.2015

geomfile=$1;
shift;

assemblydir="navigatorAssemblyArea";
mkdir ${assemblydir};

# rest of arguments are interpreted as list of volumes
for i in $@;do

# for each volume make an tmp dir
mkdir $i

cd $i
echo "executing ../NavigationKernelBenchmarker ${geomfile} ${i}"
../NavigationKernelBenchmarker ../${geomfile} $i
cp simplenavoutpool.bin outpool.bin
../NavigationSpecializerTest ../${geomfile} $i
cp ${i}Navigator.h ../${assemblydir}
cd ..

done

cd ${assemblydir}
../LibraryGenerator $@
cmake -DVecGeom_DIR=/tmp/VecGeom/lib/CMake/VecGeom -DVc_DIR=~/local/Vc1.1/lib/cmake/Vc/
make