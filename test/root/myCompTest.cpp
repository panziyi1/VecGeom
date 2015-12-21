#include <iostream>

#include "base/SOA3D.h"
#include "management/GeoManager.h"
#include "management/RootGeoManager.h"
#include "base/AOS3D.h"
#include "volumes/PlacedVolume.h"
#include "volumes/LogicalVolume.h"
#include "volumes/utilities/VolumeUtilities.h"
#include "navigation/NavigationState.h"
#include "navigation/SimpleNavigator.h"
#include "base/RNG.h"
#include "benchmarking/BenchmarkResult.h"
#include "navigation/ABBoxNavigator.h"

#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoBranchArray.h"

#include "utilities/Visualizer.h"
#include "base/Transformation3D.h"

#include <vector>
#include <set>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

using namespace VECGEOM_NAMESPACE;

// creates a four level box detector
// this modifies the global gGeoManager instance ( so no need for any return )
// void CreateRootGeom()
// {
//    double L = 10.;
//    double Lz = 10.;
//    const double Sqrt2 = sqrt(2.);
//    TGeoVolume * world =  ::gGeoManager->MakeBox("worldl",0, L, L, Lz );
//    TGeoVolume * boxlevel2 = ::gGeoManager->MakeBox("b2l",0, Sqrt2*L/2./2., Sqrt2*L/2./2., Lz );
//    TGeoVolume * boxlevel3 = ::gGeoManager->MakeBox("b3l",0, L/2./2., L/2./2., Lz);
//    TGeoVolume * boxlevel1 = ::gGeoManager->MakeBox("b1l",0, L/2., L/2., Lz );

//    boxlevel2->AddNode( boxlevel3, 0, new TGeoRotation("rot1",0,0,45));
//    boxlevel1->AddNode( boxlevel2, 0, new TGeoRotation("rot2",0,0,-45));
//    world->AddNode(boxlevel1, 0, new TGeoTranslation(-L/2.,0,0));
//    world->AddNode(boxlevel1, 1, new TGeoTranslation(+L/2.,0,0));
//    ::gGeoManager->SetTopVolume(world);
//    ::gGeoManager->CloseGeometry();
// }

void CreateSimpleRootGeom()
{
   double L = 10.;
   double Lz = 10.;
   const double Sqrt2 = sqrt(2.);
   TGeoVolume * world =  ::gGeoManager->MakeBox("worldl",0, L, L, Lz );
   // TGeoVolume * boxlevel2 = ::gGeoManager->MakeBox("b2l",0, Sqrt2*L/2./2., Sqrt2*L/2./2., Lz );
   // TGeoVolume * boxlevel3 = ::gGeoManager->MakeBox("b3l",0, L/2./2., L/2./2., Lz);
   TGeoVolume * boxlevel1 = ::gGeoManager->MakeBox("b1l",0, L/4., L, Lz );

   // boxlevel2->AddNode( boxlevel3, 0, new TGeoRotation("rot1",0,0,45));
   // boxlevel1->AddNode( boxlevel2, 0, new TGeoRotation("rot2",0,0,-45));
   world->AddNode(boxlevel1, 0, new TGeoTranslation(-3*L/4.,0.,0.));
   //world->AddNode(boxlevel1, 1, new TGeoTranslation(+L/2.,0,0));
   ::gGeoManager->SetTopVolume(world);
   ::gGeoManager->CloseGeometry();

   ::gGeoManager->Export("geom.root");
}

//void test_pointgenerationperlogicalvolume()
void test_pointgenerationperlogicalvolume(SOA3D<Precision> &globalpoints )
{
    int np = 1024;
    SOA3D<Precision> localpoints(np);
    //SOA3D<Precision> globalpoints(np);
    SOA3D<Precision> directions(np);

    // might need to resize this
    localpoints.resize(np);
    globalpoints.resize(np);
    directions.resize(np);

    volumeUtilities::FillGlobalPointsAndDirectionsForLogicalVolume(
        "worldl",
        localpoints,
        globalpoints,
        directions,
        0.5, np
    );

    assert( (int)localpoints.size() == np );
    assert( (int)globalpoints.size() == np );
    assert( (int)directions.size() == np );

    // test that points are really inside b1l; test also that they have to be in two different placed volumes
    std::set<VPlacedVolume const *> pvolumeset;
    SimpleNavigator nav;
    NavigationState * state = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth( ));
    for(int i=0;i<np;++i)
    {
        state->Clear();
        nav.LocatePoint( GeoManager::Instance().GetWorld(), globalpoints[i], *state, true );
    assert( std::strcmp( state->Top()->GetLogicalVolume()->GetLabel().c_str(), "worldl" ) == 0 );
        pvolumeset.insert( state->Top() );
    }
    // b1l should be placed two times
    assert( pvolumeset.size() == 1 );
    NavigationState::ReleaseInstance( state );
    std::cout << "test pointgenerationperlogicalvolume passed\n";
}

void Visualize(){
    RootGeoManager::Instance().LoadRootGeometry();
    VPlacedVolume const *myWorld = RootGeoManager::Instance().world();
	Visualizer visualizer;

	// for (Vector<Daughter>::const_iterator j = (*myWorld).GetDaughters().cbegin(),
 //             jEnd = (*myWorld).GetDaughters().cend(); j != jEnd; ++j) {
             
 //             visualizer.AddVolume(j.G);
 //             std::cout<< (*j)->GetLabel()<<std::endl;

	// }
int count=0;
for (auto i : (*myWorld).GetDaughters())
{   
	Transformation3D  const *trans = (*i).GetTransformation();
	std::cout<<*trans<<std::endl;
	//visualizer.AddVolume((*i).GetLogicalVolume()->Place(trans));
    //visualizer.AddVolume(*((*i).GetLogicalVolume())->Place(trans));
	visualizer.AddVolume(*i,*trans);
	//visualizer.AddVolume(*i);
	std::cout<<"Daugther no : "<< count <<std::endl;
	count++;
}
	// Vector<Daughter> daugtersVec = (*myWorld).GetDaughters();
	// for (auto i : daugtersVec)
 //        visualizer.AddVolume(i);
    visualizer.AddVolume(*myWorld);

    int np=1024;
    SOA3D<Precision> globalpoints(1024);
    test_pointgenerationperlogicalvolume(globalpoints);
    for(int i=0; i<np ;i++){
       
       visualizer.AddPoint(Vector3D<Precision>(globalpoints.x(i),globalpoints.y(i),globalpoints.z(i)));
    }
    visualizer.Show();
}
int main()
{
    CreateSimpleRootGeom();
    //ßRootGeoManager::Instance().LoadRootGeometry();
//    RootGeoManager::Instance().world()->PrintContent();
    RootGeoManager::Instance().PrintNodeTable();
    Visualize();

}