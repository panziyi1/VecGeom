#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedBox.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "volumes/UnplacedBox.h"
//#include "utilities/Visualizer.h"
#include <string>
#include <cmath>


using namespace vecgeom;


int main( int argc, char *argv[] ) {
/*
    if( argc < 9 )
    {
      std::cerr << "need to give geometry lib + logical volume name + local point + local dir\n";
      std::cerr << "example: " << argv[0] << " libCMSgeom.so CALO 10.0 0.8 -3.5 1 0 0\n";
      return 1;
    }
*/


  GeoManager &gGeoManager = GeoManager::Instance();
  gGeoManager.LoadGeometryFromSharedLib(argv[1]);
  //gGeoManager.LoadGeometryFromSharedLib("./libCMSgeom.so");



  VPlacedVolume * vecgeomplaced = gGeoManager.FindPlacedVolume(argv[2]);
  //std::string testvolume( argv[2] );

    double px = atof(argv[3]);
    double py = atof(argv[4]);
    double pz = atof(argv[5]);
    double dirx = atof(argv[6]);
    double diry = atof(argv[7]);
    double dirz = atof(argv[8]);
   
 

    // now get the VecGeom shape and benchmark it
    if( vecgeomplaced )
    {
/*
       Vector3D<Precision> minExtent,maxExtent;
       foundvolume->Extent(minExtent,maxExtent);
       Vector3D<Precision> delta = maxExtent -minExtent;
       foundvolume->Extent(minExtent,maxExtent);
       Vector3D<Precision> delta = maxExtent -minExtent;
       UnplacedBox worldUnplaced = UnplacedBox(delta*4);
       LogicalVolume newworld = LogicalVolume("newworld", &worldUnplaced);
       newworld.PlaceDaughter(foundvolume);
       VPlacedVolume *newworldPlaced = newworld.Place();
       gGeoManager.SetWorld(newworldPlaced);
       gGeoManager.CloseGeometry();
*/


  //      VPlacedVolume * vecgeomplaced = converted->Place();
        Vector3D<Precision> point(px,py,pz);
        Vector3D<Precision> dir(dirx,diry,dirz);

        std::cout<<"COMPARE DISTANCE point "<<point<<" distacne "<<dir<<std::endl;

	SOA3D<Precision> pointcontainer(4); pointcontainer.resize(4); 
	SOA3D<Precision> dircontainer(4); dircontainer.resize(4); 
	Precision * output = new Precision[4];
	Precision * steps = new Precision[4];

        if(argc > 9)
	{
          pointcontainer.set(0, point);
	  dircontainer.set(0, dir.x(), dir.y(), dir.z() );
          double px2 = atof(argv[9]);
          double py2 = atof(argv[10]);
          double pz2 = atof(argv[11]);
          double dirx2 = atof(argv[12]);
          double diry2 = atof(argv[13]);
          double dirz2 = atof(argv[14]);
          pointcontainer.set(1, px2,py2,pz2);
	  dircontainer.set(1, dirx2, diry2, dirz2 );
          pointcontainer.set(2, point);
	  dircontainer.set(2, dir.x(), dir.y(), dir.z() );
          pointcontainer.set(3, px2,py2,pz2);
	  dircontainer.set(3, dirx2, diry2, dirz2 );

          for( auto i=0;i<4;++i ){
	   steps[i] = vecgeom::kInfinity;
	  }

	}
        else{
	  for( auto i=0;i<4;++i ){
	  pointcontainer.set(i, point);
	  dircontainer.set(i, dir.x(), dir.y(), dir.z() );
	  steps[i] = vecgeom::kInfinity;
	  }
        }
	if( ! dir.IsNormalized() ){
	  std::cerr << "** Attention: Direction is not normalized **\n";
	  std::cerr << "** Direction differs from 1 by " << std::sqrt(dir.x()*dir.x() + dir.y()*dir.y()+ dir.z()*dir.z())-1.<< "\n";
	}
        double dist;
        bool contains = vecgeomplaced->Contains( point );
        std::cout << "VecGeom Capacity " << vecgeomplaced->Capacity( ) << "\n";
        std::cout << "VecGeom CONTAINS " << contains << "\n";
        dist = vecgeomplaced->DistanceToIn( point, dir );
        std::cout << "VecGeom DI " << dist << "\n";
        if(dist < vecgeom::kInfinity )
	  {
            double ddist = vecgeomplaced->DistanceToIn( point+dir*dist,-dir );
            EnumInside eInside = vecgeomplaced->Inside( point+dir*dist );
           std::cout << "VecGeom INSIDE(p=p+dist*dir) " << eInside << "\n";
           if(eInside == vecgeom::kOutside)
           std::cout << "VecGeom Distance seems to be to big  DI(p=p+dist*dir,-dir) " << ddist << "\n";
           if(eInside == vecgeom::kInside)
           std::cout << "VecGeom Distance seems to be to small DO(p=p+dist*dir,dir) " << ddist << "\n";
	  }
	vecgeomplaced->DistanceToIn( pointcontainer, dircontainer, steps, output );
        std::cout << "VecGeom DI-V " << output[0] << "\n";
        std::cout << "VecGeom DO " << vecgeomplaced->DistanceToOut( point, dir ) << "\n";
        vecgeomplaced->DistanceToOut( pointcontainer, dircontainer, steps, output );
	std::cout << "VecGeom DO-V " << output[0] << "\n";

        std::cout << "VecGeom SI " << vecgeomplaced->SafetyToIn( point ) << "\n";
        std::cout << "VecGeom SO " << vecgeomplaced->SafetyToOut( point ) << "\n";


#ifdef VECGEOM_USOLIDS
        VUSolid const * usolid = vecgeomplaced->ConvertToUSolids();
        if( usolid != NULL ){
          std::cout << "USolids Capacity " << const_cast<VUSolid*>(usolid)->Capacity(  ) << "\n";
          std::cout << "USolids CONTAINS " << usolid->Inside( point ) << "\n";
          std::cout << "USolids DI " << usolid->DistanceToIn( point, dir ) << "\n";

          Vector3D<Precision> norm; bool valid;
          std::cout << "USolids DO " << usolid->DistanceToOut( point, dir, norm, valid ) << "\n";
          std::cout << "USolids SI " << usolid->SafetyFromInside( point ) << "\n";
          std::cout << "USolids SO " << usolid->SafetyFromOutside( point ) << "\n";
        }
        else{
          std::cerr << "USOLID conversion failed\n";
        }
#endif

#ifdef VECGEOM_GEANT4
        G4ThreeVector g4p( point.x(), point.y(), point.z() );
        G4ThreeVector g4d( dir.x(), dir.y(), dir.z());

        G4VSolid const * g4solid = vecgeomplaced->ConvertToGeant4();
        if( g4solid != NULL ){
          std::cout << "G4 CONTAINS " << g4solid->Inside( g4p ) << "\n";
          std::cout << "G4 DI " << g4solid->DistanceToIn( g4p, g4d ) << "\n";
          std::cout << "G4 DO " << g4solid->DistanceToOut( g4p, g4d ) << "\n";
          std::cout << "G4 SI " << g4solid->DistanceToIn( g4p ) << "\n";
          std::cout << "G4 SO " << g4solid->DistanceToOut( g4p ) << "\n";
        }
        else{
          std::cerr << "G4 conversion failed\n";
        }
#endif
     /*
        double step=0;
        if( vecgeomplaced->Contains( &point[0] ) ){
            step = vecgeomplaced->DistFromInside( &point[0], &dir[0] );
        }
        else
        {
            step = vecgeomplaced->DistFromOutside( &point[0], &dir[0] );
        }
        Visualizer visualizer;
        visualizer.AddVolume( *vecgeomplaced );
        visualizer.AddPoint( point );
        visualizer.AddLine( point, point + step * dir );
        visualizer.Show();
     */
    }
    else
    {
        std::cerr << " NO SUCH VOLUME [" <<  argv[2] << "] FOUND ... EXITING \n";
        return 1;
    }
}

