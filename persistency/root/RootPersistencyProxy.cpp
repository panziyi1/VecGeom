#include "persistency/root/include/RootPersistencyProxy.h"

template <class T>
void SetCollectionProxy()
{

  auto info = ROOT::Detail::TCollectionProxyInfo::Generate(ROOT::Detail::TCollectionProxyInfo::Pushback<T>());

  if (!info) {
    fprintf(stderr, "Could not create the collectionProxyInfo for %s\n", typeid(T).name());
  }
  auto cl = TClass::GetClass(typeid(T));

  if (cl) {
    cl->SetCollectionProxy(*info);
  } else {
    fprintf(stderr, "Could not find the TClass for %s\n", typeid(T).name());
  }

  delete info;
}

void RootPersistencyProxy()
{
  SetCollectionProxy<vecgeom::Vector<const vecgeom::VPlacedVolume *>>();
  SetCollectionProxy<vecgeom::Vector<double>>();
  SetCollectionProxy<vecgeom::Vector<double *>>();
}