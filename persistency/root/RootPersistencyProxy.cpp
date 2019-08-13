#include "persistency/root/include/RootPersistencyProxy.h"
#include <base/Array.h>

#include <TVirtualCollectionProxy.h>
#include <TGenCollectionStreamer.h>
#include <TCollectionProxyFactory.h>

template <class T>
void SetCollectionProxy(const char *altname)
{

  auto info = ROOT::Detail::TCollectionProxyInfo::Generate(ROOT::Detail::TCollectionProxyInfo::Pushback<T>());

  if (!info) {
    fprintf(stderr, "Could not create the collectionProxyInfo for %s\n", typeid(T).name());
  }
  auto cl = TClass::GetClass(typeid(T));

  TString keep = cl->GetName();

  R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

  // Trick the collection in getting the right type.
  cl->SetName(altname);
  if (cl) {
    cl->SetCollectionProxy(*info);
    auto proxy = cl->GetCollectionProxy();
    // Force initialization
    auto type = proxy->GetCollectionType();
    // fprintf(stderr, "Collection type is %d\n", type);

    auto streamer = dynamic_cast<TCollectionClassStreamer *>(cl->GetStreamer());
    if (streamer) {
      proxy = streamer->GetXYZ();
      type  = proxy->GetCollectionType();
      // fprintf(stderr, "Collection Streamer type is %d\n", type);
    }

  } else {
    fprintf(stderr, "Could not find the TClass for %s\n", typeid(T).name());
  }
  cl->SetName(keep);

  delete info;
}

void RootPersistencyProxy()
{
  // SetCollectionProxy<std::vector<const vecgeom::VPlacedVolume *>>();
  // SetCollectionProxy<std::vector<double>>();
  // SetCollectionProxy<std::vector<double *>>();
  // SetCollectionProxy<vecgeom::Array<double>>();

  SetCollectionProxy<vecgeom::Vector<const vecgeom::VPlacedVolume *>>("list<const vecgeom::VPlacedVolume *>");
  SetCollectionProxy<vecgeom::Vector<double>>("list<double>");
  SetCollectionProxy<vecgeom::Vector<float>>("list<float>");
  // This is is unlikely to be stored correctly as ROOT can't know the number of double pointed to by the double*.
  // (it is *unlikely* to be just 1 as it is assumed for objects)
  // SetCollectionProxy<vecgeom::Vector<double *>>("vector<double *>");
}
