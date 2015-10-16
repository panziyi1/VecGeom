/// author: Andrei Gheata (andrei.gheata@cern.ch) Oct 14, 2015
#include <sstream>
#include <string>
#include "volumes/Tube.h"
#include "volumes/LogicalVolume.h"
#include "management/GeoManager.h"
#include "base/Stopwatch.h"
#include "SC15Navigators.h"
#include "navigation/SimpleNavigator.h"
#include "base/Vector3D.h"
#include <iostream>
using namespace vecgeom;

#define SETINNERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) InnerMostTubeNavigator<depth+1>::Instance() );

#define SETLAYERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) LayerNavigator<depth+1>::Instance() );

#define SETWORLDNAV(depth) \
  if(layer==0) vol->SetUserExtensionPtr( (void*) WorldNavigator<true>::Instance() );

#define WRITE_FILE_NAME "volumeImage.bmp" // output image name

static VNavigator const* GetNavigator(LogicalVolume const*lvol){
  return (VNavigator const*)lvol->GetUserExtensionPtr();
}


typedef struct tFILE_HEADER
{
  unsigned short bfType;
  unsigned long bfSize;
  unsigned short bfReserved1;
  unsigned short bfReserved2;
  unsigned long bfOffBits;
} FILE_HEADER;


typedef struct tINFO_HEADER
{
   unsigned long biSize;
   unsigned long biWidth;
   unsigned long biHeight;
   unsigned short biPlanes;
   unsigned short biBitCount;
   unsigned long biCompression;
   unsigned long biSizeImage;
   unsigned long biXPelsPerMeter;
   unsigned long biYPelsPerMeter;
   unsigned long biClrUsed;
   unsigned long biClrImportant;
} INFO_HEADER;

typedef struct tMY_BITMAP
{
  FILE_HEADER  bmpFileHeader;
  INFO_HEADER  bmpInfoHeader;
  unsigned char* bmpPalette;
  unsigned char* bmpRawData;
} MY_BITMAP;

typedef struct tBasket
{
  unsigned int pixel;        // pixel X/Y coordinates encoded in a single long
} Basket;

inline unsigned int EncodePixel(unsigned int ix, unsigned int iy) { return ( (iy << 16) | ix); }
inline void DecodePixel(unsigned int const pixel, unsigned int &ix, unsigned int &iy) { 
  ix = pixel & 0xFFFF; iy = pixel >> 16;
}  


// produce a bmp image out of pixel information given in volume_results
void make_bmp_header(MY_BITMAP * pBitmap, unsigned char * bmpBuf, int sizex, int sizey)
{
  int width_4= (sizex+ 3)&~3;
  unsigned int len= 0;

  // bitmap file header
  pBitmap->bmpFileHeader.bfType=0x4d42;
  pBitmap->bmpFileHeader.bfSize=sizey* width_4* 3+ 54;
  pBitmap->bmpFileHeader.bfReserved1= 0;
  pBitmap->bmpFileHeader.bfReserved2= 0;
  pBitmap->bmpFileHeader.bfOffBits= 54;

  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfType, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfSize, 4);
  len+= 4;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved1, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved2, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfOffBits, 4);
  len+= 4;

  // bitmap information header
  pBitmap->bmpInfoHeader.biSize= 40;
  pBitmap->bmpInfoHeader.biWidth= width_4;
  pBitmap->bmpInfoHeader.biHeight= sizey;
  pBitmap->bmpInfoHeader.biPlanes= 1;
  pBitmap->bmpInfoHeader.biBitCount= 24;
  pBitmap->bmpInfoHeader.biCompression= 0;
  pBitmap->bmpInfoHeader.biSizeImage= sizey* width_4* 3;
  pBitmap->bmpInfoHeader.biXPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biYPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biClrUsed= 0;
  pBitmap->bmpInfoHeader.biClrImportant=0;


  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSize, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biWidth, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biHeight, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biPlanes, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biBitCount, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biCompression, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSizeImage, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biXPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biYPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrUsed, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrImportant, 4);
  len+= 4;
}





void make_bmp(int const * volume_result, char const *name, int data_size_x, int data_size_y, bool linear = true)
{

  MY_BITMAP* pBitmap= new MY_BITMAP;
  FILE *pBitmapFile;
  int width_4= (data_size_x+ 3)&~3;
  unsigned char* bmpBuf;

  bmpBuf = (unsigned char*) new unsigned char[data_size_y* width_4* 3+ 54];
  printf("\n Write bitmap...\n");

  unsigned int len= 0;

  // bitmap file header
  pBitmap->bmpFileHeader.bfType=0x4d42;
  pBitmap->bmpFileHeader.bfSize=data_size_y* width_4* 3+ 54;
  pBitmap->bmpFileHeader.bfReserved1= 0;
  pBitmap->bmpFileHeader.bfReserved2= 0;
  pBitmap->bmpFileHeader.bfOffBits= 54;

  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfType, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfSize, 4);
  len+= 4;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved1, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved2, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfOffBits, 4);
  len+= 4;

  // bitmap information header
  pBitmap->bmpInfoHeader.biSize= 40;
  pBitmap->bmpInfoHeader.biWidth= width_4;
  pBitmap->bmpInfoHeader.biHeight= data_size_y;
  pBitmap->bmpInfoHeader.biPlanes= 1;
  pBitmap->bmpInfoHeader.biBitCount= 24;
  pBitmap->bmpInfoHeader.biCompression= 0;
  pBitmap->bmpInfoHeader.biSizeImage= data_size_y* width_4* 3;
  pBitmap->bmpInfoHeader.biXPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biYPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biClrUsed= 0;
  pBitmap->bmpInfoHeader.biClrImportant=0;

    memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSize, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biWidth, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biHeight, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biPlanes, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biBitCount, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biCompression, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSizeImage, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biXPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biYPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrUsed, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrImportant, 4);
  len+= 4;
  // find out maxcount before doing the picture
  int maxcount = 0;
  int x=0,y=0,origin_x=0;
  while( y< data_size_y )
  {
     while( origin_x< data_size_x )
     {
       int value = *(volume_result+y*data_size_x+origin_x);
       maxcount = ( value > maxcount )? value : maxcount;

       x++;
       origin_x++;
     }
     y++;
     x = 0;
     origin_x = 0;
  }
//  maxcount = std::log(maxcount + 1);

  x= 0;
  y= 0;
  origin_x= 0;

  int padding= width_4- data_size_x;
  int padding_idx= padding;
  unsigned char *imgdata= (unsigned char*) new unsigned char[data_size_y*width_4*3];

  int totalcount = 0;

  while( y< data_size_y )
  {
    while( origin_x< data_size_x )
    {
      int value = *(volume_result+y*data_size_x+origin_x);
      totalcount += value;
      if( linear ){
         *(imgdata+y*width_4*3+x*3+0)= (value/(1.*maxcount)) * 256;
         *(imgdata+y*width_4*3+x*3+1)= (value/(1.*maxcount)) * 256;
         *(imgdata+y*width_4*3+x*3+2)= (value/(1.*maxcount)) * 256;
      }
      else {
         *(imgdata+y*width_4*3+x*3+0)= (log(value+1))/(1.*(log(1+maxcount))) * 256;
         *(imgdata+y*width_4*3+x*3+1)= (log(value+1))/(1.*(log(1+maxcount))) * 256;
         *(imgdata+y*width_4*3+x*3+2)= (log(value+1))/(1.*(log(1+maxcount))) * 256;
      }
      x++;
      origin_x++;

      while( origin_x== data_size_x && padding_idx)
      {
      // padding 4-byte at bitmap image
        *(imgdata+y*width_4*3+x*3+0)= 0;
        *(imgdata+y*width_4*3+x*3+1)= 0;
        *(imgdata+y*width_4*3+x*3+2)= 0;
        x++;
        padding_idx--;
      }
      padding_idx= padding;
    }
    y++;
    x= 0;
    origin_x= 0;
  }

  memcpy(bmpBuf + 54, imgdata, width_4* data_size_y* 3);

  pBitmapFile = fopen(name, "wb");
  fwrite(bmpBuf, sizeof(char), width_4*data_size_y*3+54, pBitmapFile);


  fclose(pBitmapFile);
  delete[] imgdata;
  delete[] bmpBuf;
  delete pBitmap;

  std::cout << " wrote image file " << name <<  "\n";
  std::cout << " total count " << totalcount << "\n";
  std::cout << " max count " << maxcount << "\n";
}





void AssignNavigatorToVolume( LogicalVolume *vol, int layer, int maxlayers ){
  assert( maxlayers <= 10 ); // we are only listing 10 template specializations up to depth 10 here
  if( layer == maxlayers - 1 ) // it is the last layer
  {
    SETINNERNAV(0);
    SETINNERNAV(1);
    SETINNERNAV(2);
    SETINNERNAV(3);
    SETINNERNAV(4);
    SETINNERNAV(5);
    SETINNERNAV(6);
    SETINNERNAV(7);
    SETINNERNAV(8);
    SETINNERNAV(9);
    SETINNERNAV(10);
  }
  if( layer < maxlayers - 1 ){
    SETLAYERNAV(0);
    SETLAYERNAV(1);
    SETLAYERNAV(2);
    SETLAYERNAV(3);
    SETLAYERNAV(4);
    SETLAYERNAV(5);
    SETLAYERNAV(6);
    SETLAYERNAV(7);
    SETLAYERNAV(8);
    SETLAYERNAV(9);
    SETLAYERNAV(10);
  }
}


VPlacedVolume *CreateSimpleTracker(int nlayers) {
  std::cout << "Creating SimpleTracker geometry having " << nlayers << " layers"
            << std::endl;
  // World size
  const double world_size = 500.;

  // Top volume
  UnplacedBox *uTop = new UnplacedBox(world_size + 2, world_size + 2, world_size + 2);
  LogicalVolume *top = new LogicalVolume("world", uTop);
  top->SetUserExtensionPtr( (void *) WorldNavigator<true>::Instance() );
  
// Cylindrical layers
  double rmax, dz;
  std::string layerBase = "layer_";
  double deltaR = world_size/nlayers;

  LogicalVolume *mother = top;
  for (auto layer=0; layer<nlayers; ++layer) {
    rmax = world_size - layer * deltaR;
    dz = world_size - 0.1*layer*deltaR;
    std::ostringstream layerName;
    layerName << layerBase << layer;
    UnplacedTube *uLayer = new UnplacedTube(0, rmax, dz, 0, kTwoPi);
    LogicalVolume *layerVol = new LogicalVolume(layerName.str().c_str(), uLayer);
    
    AssignNavigatorToVolume(layerVol, layer, nlayers);

    // Place in mother
    mother->PlaceDaughter(layerName.str().c_str(), layerVol, &Transformation3D::kIdentity);
    // change mother to current lvol (to make a real hierarchy)
    mother = layerVol;
  }  

  VPlacedVolume *world = top->Place();
  GeoManager::Instance().SetWorld(world);
  GeoManager::Instance().CloseGeometry();
  return world;
}


int ScalarNavigation(Vector3D<Precision> p, Vector3D<Precision> dir) {
  // init navstates
  NavigationState * curnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  NavigationState * newnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  
  SimpleNavigator nav;
  nav.LocatePoint( GeoManager::Instance().GetWorld(), p, *curnavstate, true );

  int crossedvolumecount=0;
  while( ! curnavstate->IsOutside() ) {
    //
    //std::cout << "tracking in " << curnavstate->Top()->GetLogicalVolume()->GetName() << "\n";
      // get navigator object and move point
    VNavigator const *specialnav = GetNavigator(curnavstate->Top()->GetLogicalVolume());
    double step = specialnav->ComputeStepAndPropagatedState(p, dir, kInfinity, *curnavstate, *newnavstate);
    //std::cout << "step " << step << "\n";
    p = p + dir * (step + 1E-6);

    // pointer swap is enough
    auto *tmp = curnavstate;
    curnavstate = newnavstate;
    newnavstate = tmp;
    if (step>0.0) crossedvolumecount++;
  }

  // now test the vector progression (of coherent rays)
  // this is just testing the interface and makes sure that no trivial things go wrong
  NavigationState::ReleaseInstance(curnavstate);
  NavigationState::ReleaseInstance(newnavstate);
  return crossedvolumecount;
}

////////////////// VECTOR NAVIGATION

int VectorNavigation (SOA3D<Precision> points, SOA3D<Precision>dirs, int np) {

  NavStatePool *curnavstates = new NavStatePool(np, GeoManager::Instance().getMaxDepth());
  NavStatePool *newnavstates = new NavStatePool(np, GeoManager::Instance().getMaxDepth());
  for(auto i=0;i<np;++i){
    SimpleNavigator nav;
    nav.LocatePoint( GeoManager::Instance().GetWorld(), points[i], *(*curnavstates)[i], true );
  }

  double *steps    = (double*) _mm_malloc(np*sizeof(double),64);
  double *psteps    = (double*) _mm_malloc(np*sizeof(double),64);
  int crossedvolumecount=0;
  while(!(*curnavstates)[0]->IsOutside()) {
      //
     // std::cout << "tracking in " << (*curnavstates)[0]->Top()->GetLogicalVolume()->GetName() << "\n";
        // get navigator object and move point
      VNavigator const *specialnav = GetNavigator((*curnavstates)[0]->Top()->GetLogicalVolume());
      specialnav->ComputeStepsAndPropagatedStates(points, dirs, psteps, *curnavstates, *newnavstates, steps);

      // transport + crosschecks + pointerswap
      // We should vectorize this transport !!
      for(auto i=0; i<np; ++i){
        //assert(steps[i]==steps[0]);
        assert((*newnavstates)[i]->Top()==(*newnavstates)[0]->Top());

        points.set(i, points[i] + dirs[i]*(steps[i] + 1E-6));

      }
      // pointer swap to update navstates
      auto tmp=curnavstates;
      curnavstates = newnavstates;
      newnavstates = tmp;
      if (steps[0]>0.0) crossedvolumecount++;
  }

    delete curnavstates;
    delete newnavstates;
    return crossedvolumecount;
}

// a test function to verify correct functioning of the navigators
// for a simple test case
void TestScalarNavigation() {
  // setup point and direction in world
  Vector3D<Precision> p(-51.,0,0);
  Vector3D<Precision> dir(1.,0,0);
  std::cout<<"crossedvolumecount "<<ScalarNavigation(p,dir)<<std::endl;
}

void TestVectorNavigation() {
  // use vector navigation
  // choose 11 particles to also make sure we have tail treatment
  auto np=11;
  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);
  for(auto i=0;i<np;++i){
      points.set(i, -51,0,0);
      dirs.set(i, 1,0,0);
  }
  std::cout<<VectorNavigation(points, dirs,np)<<std::endl; 
}
// reproducing the pixel-by-pixel XRayBenchmark
// target: show speed gain from specialized navigators
// in comparision to current XRay-Benchmark
void XRayBenchmark(int axis, int pixel_width) {

  std::stringstream imagenamebase;
  imagenamebase << "simpleTrackerImage_";
  if(axis==1) imagenamebase << "x";
  if(axis==2) imagenamebase << "y";
  if(axis==3) imagenamebase << "z";

  Vector3D<Precision> minExtent,maxExtent;
  GeoManager::Instance().FindPlacedVolume(2)->Extent(minExtent,maxExtent);
  //GeoManager::Instance().GetWorld()->Extent(minExtent,maxExtent);
  
  double dx = (maxExtent - minExtent).x(); 
  double dy = (maxExtent - minExtent).y(); 
  double dz = (maxExtent - minExtent).z(); 

  Vector3D<Precision> orig = (maxExtent + minExtent)/2.;

  double axis1_start= 0.;
  double axis1_end= 0.;

  double axis2_start= 0.;
  double axis2_end= 0.;

  double pixel_axis= 1.;


  Vector3D<Precision> dir;
    if(axis== 1)
    {
      dir.Set(1., 0., 0.);
      axis1_start= orig.y()-dy;
      axis1_end= orig.y() + dy;
      axis2_start= orig.z()-dz;
      axis2_end= orig.z()+ dz;
      pixel_axis= (dy*2)/pixel_width;
    }
    else if(axis== 2)
    {
      dir.Set(0., 1., 0.);
      axis1_start= orig.x()-dx;
      axis1_end= orig.x()+ dx;
      axis2_start= orig.z()-dz;
      axis2_end= orig.z()+ dz;
      pixel_axis= (dx*2)/pixel_width;
    }
    else if(axis== 3)
    {
      dir.Set(0., 0., 1.);
      axis1_start= orig.x() -dx;
      axis1_end= orig.x()+ dx;
      axis2_start= orig.y()-dy;
      axis2_end= orig.y()+ dy;
      pixel_axis= (dx*2)/pixel_width;
    }

    // init data for image
    int data_size_x= (axis1_end-axis1_start)/pixel_axis;
    int data_size_y= (axis2_end-axis2_start)/pixel_axis;
    double pixel_width_1 = (axis1_end-axis1_start)/data_size_x;
    double pixel_width_2 = (axis2_end-axis2_start)/data_size_y;

    std::cout << "data_size_x = " << data_size_x << std::endl;
    std::cout << "data_size_y = " << data_size_y << std::endl;

    int *volume_result= (int*) new int[data_size_y * data_size_x*3];
   
    Stopwatch timer;
    timer.Start();

    for( int pixel_count_2 = 0; pixel_count_2 < data_size_y; ++pixel_count_2 ) {
       for( int pixel_count_1 = 0; pixel_count_1 < data_size_x; ++pixel_count_1 ) {
          double axis2_count = axis2_start + pixel_count_2 * pixel_width_2 + 1E-6;
          double axis1_count = axis1_start + pixel_count_1 * pixel_width_1 + 1E-6;

          Vector3D<Precision> p;
          if( axis== 1 )
              p.Set( orig.x(), axis1_count, axis2_count);
          else if( axis== 2)
              p.Set( axis1_count, orig.y(), axis2_count);
          else if( axis== 3)
              p.Set( axis1_count, axis2_count, orig.z());

          *(volume_result+pixel_count_2*data_size_x+pixel_count_1) = ScalarNavigation(p,dir);
      } // end inner loop
    } // end outer loop

   timer.Stop();
   std::cout << " XRay Elapsed time : "<< timer.Elapsed() << std::endl;

    std::stringstream VecGeomimage;
    VecGeomimage << imagenamebase.str();
    VecGeomimage << "_VecGeom.bmp";
    make_bmp(volume_result, VecGeomimage.str().c_str(), data_size_x, data_size_y);

}



void XRayBenchmarkVecNav(int axis, int pixel_width) {

  auto N=64;
  std::stringstream imagenamebase;
  imagenamebase << "simpleTrackerimage_";
  if(axis==1) imagenamebase << "x";
  if(axis==2) imagenamebase << "y";
  if(axis==3) imagenamebase << "z";

  Vector3D<Precision> minExtent,maxExtent;
  GeoManager::Instance().FindPlacedVolume(2)->Extent(minExtent,maxExtent);
  //GeoManager::Instance().GetWorld()->Extent(minExtent,maxExtent);
  
  double dx = (maxExtent - minExtent).x(); 
  double dy = (maxExtent - minExtent).y(); 
  double dz = (maxExtent - minExtent).z(); 

  Vector3D<Precision> orig = (maxExtent + minExtent)/2.;

  double axis1_start= 0.;
  double axis1_end= 0.;

  double axis2_start= 0.;
  double axis2_end= 0.;

  double pixel_axis= 1.;


  Vector3D<Precision> dir;
    if(axis== 1)
    {
      dir.Set(1., 0., 0.);
      axis1_start= orig.y()-dy;
      axis1_end= orig.y() + dy;
      axis2_start= orig.z()-dz;
      axis2_end= orig.z()+ dz;
      pixel_axis= (dy*2)/pixel_width;
    }
    else if(axis== 2)
    {
      dir.Set(0., 1., 0.);
      axis1_start= orig.x()-dx;
      axis1_end= orig.x()+ dx;
      axis2_start= orig.z()-dz;
      axis2_end= orig.z()+ dz;
      pixel_axis= (dx*2)/pixel_width;
    }
    else if(axis== 3)
    {
      dir.Set(0., 0., 1.);
      axis1_start= orig.x() -dx;
      axis1_end= orig.x()+ dx;
      axis2_start= orig.y()-dy;
      axis2_end= orig.y()+ dy;
      pixel_axis= (dx*2)/pixel_width;
    }

    // init data for image
    int data_size_x= (axis1_end-axis1_start)/pixel_axis;
    int data_size_y= (axis2_end-axis2_start)/pixel_axis;
    double pixel_width_1 = (axis1_end-axis1_start)/data_size_x;
    double pixel_width_2 = (axis2_end-axis2_start)/data_size_y;

    std::cout << "data_size_x = " << data_size_x << std::endl;
    std::cout << "data_size_y = " << data_size_y << std::endl;

    int *volume_result= (int*) new int[data_size_y * data_size_x*3];
   


  SOA3D<Precision> dirs(N);
  for(auto i=0;i<N;++i){
      dirs.set(i, dir.x(), dir.y(),dir.z());
  }

  Stopwatch timer;
  timer.Start();

    for( int pixel_count_2 = 0; pixel_count_2 < data_size_y; ++pixel_count_2 ) {
       for( int pixel_count_1 = 0; pixel_count_1 < data_size_x; ++pixel_count_1 ) {
          double axis2_count = axis2_start + pixel_count_2 * pixel_width_2 + 1E-6;
          double axis1_count = axis1_start + pixel_count_1 * pixel_width_1 + 1E-6;

          SOA3D<Precision> points(N);
          for(auto i=0;i<N;++i){
             if( axis== 1 )
                points.set( i, orig.x(), axis1_count, axis2_count );
             else if( axis== 2)
                points.set( i, axis1_count, orig.y(), axis2_count );
             else if( axis== 3)
                points.set( i, axis1_count, axis2_count, orig.z() );
          }

          *(volume_result+pixel_count_2*data_size_x+pixel_count_1) = VectorNavigation(points,dirs,N);
      } // end inner loop
   } // end outer loop 
   timer.Stop();
   std::cout << " XRayVecNav Elapsed time : "<< timer.Elapsed() << std::endl;

    std::stringstream VecGeomimage;
    VecGeomimage << imagenamebase.str();
    VecGeomimage << "_VecGeom.bmp";
    make_bmp(volume_result, VecGeomimage.str().c_str(), data_size_x, data_size_y);

}



void BasketBasedXRayBenchmark() {
  // to be filled in by Andrei

}


int main(int argc, char* argv[]) {
  int nlayers = 10;
  if (argc > 1) {
    nlayers = std::atoi(argv[1]);
  }
  CreateSimpleTracker(nlayers);

  int axis= 0;

  if( strcmp(argv[2], "x")==0 )
    axis= 1;
  else if( strcmp(argv[2], "y")==0 )
    axis= 2;
  else if( strcmp(argv[2], "z")==0 )
    axis= 3;
  else
  {
    std::cerr<< "Incorrect axis"<< std::endl<< std::endl;
    return 1;
  }

  double pixel_width= atof(argv[3]);
 
  XRayBenchmark(axis, pixel_width);
  XRayBenchmarkVecNav(axis, pixel_width);




/*
  // loop over all logical volumes and print navigator
  std::vector<LogicalVolume *> lvols;
  GeoManager::Instance().GetAllLogicalVolumes(lvols);

  for( auto v : lvols ){
    auto nav = (VNavigator*) v->GetUserExtensionPtr();
    std::cerr << v->GetName() << " has navigator " << nav->GetName() << "\n";
  }


  // test tracking
  TestScalarNavigation();
<<<<<<< HEAD
<<<<<<< HEAD
  std::cout << "start vector test\n";
  TestVectorNavigation();
=======
  std::cout << "start vector test\n";
  TestVectorNavigation();
*/
}
