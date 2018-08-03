#include "test_fullgeo.hh"

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

void proxy()
{
  SetCollectionProxy<vecgeom::Vector<const vecgeom::VPlacedVolume *>>();
}