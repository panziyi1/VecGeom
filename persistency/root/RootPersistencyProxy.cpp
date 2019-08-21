#include "persistency/root/include/RootPersistencyProxy.h"
#include <base/Array.h>
#include <volumes/PolyhedronStruct.h>
#include <volumes/PolyconeSection.h>

#include <TVirtualCollectionProxy.h>
#include <TGenCollectionStreamer.h>
#include <TCollectionProxyFactory.h>
#include <TClassStreamer.h>

template <typename Value_t>
void Construct(void *m)
{
  ::new (m) Value_t();
}

template <>
void Construct<vecgeom::cxx::ZSegment>(void *m)
{
  ::new (m) vecgeom::cxx::ZSegment((TRootIOCtor *)nullptr);
}

template <class T>
struct ArrayType : public ROOT::Detail::TCollectionProxyInfo::Address<TYPENAME T::const_reference> {
  typedef T Cont_t;
  typedef typename T::iterator Iter_t;
  typedef typename T::value_type Value_t;
  typedef ROOT::Detail::TCollectionProxyInfo::Environ<Iter_t> Env_t;
  typedef Env_t *PEnv_t;
  typedef Cont_t *PCont_t;
  typedef Value_t *PValue_t;

  virtual ~ArrayType() {}

  static inline PCont_t object(void *ptr) { return PCont_t(PEnv_t(ptr)->fObject); }
  static void *size(void *env)
  {
    PEnv_t e = PEnv_t(env);
    e->fSize = ROOT::Detail::TCollectionProxyInfo::SfinaeHelper::GetContainerSize(*PCont_t(e->fObject));
    return &e->fSize;
  }
  static void *clear(void *env)
  {
    // object(env)->clear();
    return 0;
  }
  static void *first(void *env)
  {
    PEnv_t e  = PEnv_t(env);
    PCont_t c = PCont_t(e->fObject);
#if 0
         // Assume iterators do not need destruction
         ::new(e->buff) Iter_t(c->begin());
#endif
    e->fIterator = c->begin();
    e->fSize     = ROOT::Detail::TCollectionProxyInfo::SfinaeHelper::GetContainerSize(*c);
    if (0 == e->fSize) return e->fStart = 0;
    TYPENAME T::const_reference ref = *(e->iter());
    return e->fStart                = ArrayType<T>::address(ref);
  }
  static void *next(void *env)
  {
    PEnv_t e  = PEnv_t(env);
    PCont_t c = PCont_t(e->fObject);
    for (; e->fIdx > 0 && e->iter() != c->end(); ++(e->iter()), --e->fIdx) {
    }
    // TODO: Need to find something for going backwards....
    if (e->iter() == c->end()) return 0;
    TYPENAME T::const_reference ref = *(e->iter());
    return ArrayType<T>::address(ref);
  }
  static void *construct(void *what, size_t size)
  {
    PValue_t m = PValue_t(what);
    for (size_t i = 0; i < size; ++i, ++m) {
      // SOMETIMES pass nothings, sometimes pass ((TRootIOCtor *)nullptr)
      //::new(m) Value_t();
      Construct<Value_t>(m);
    }
    return 0;
  }
  static void *collect(void *coll, void *array)
  {
    PCont_t c  = PCont_t(coll);
    PValue_t m = PValue_t(array);
    for (Iter_t i = c->begin(); i != c->end(); ++i, ++m)
      ::new (m) Value_t(*i);
    return 0;
  }
  static void destruct(void *what, size_t size)
  {
    PValue_t m = PValue_t(what);
    for (size_t i = 0; i < size; ++i, ++m)
      m->~Value_t();
  }

  static const bool fgLargeIterator = sizeof(typename Cont_t::iterator) > TGenCollectionProxy::fgIteratorArenaSize;
  typedef ROOT::Detail::TCollectionProxyInfo::Iterators<Cont_t, fgLargeIterator> Iterators_t;
};

template <class T>
struct ForArray : public ArrayType<T> {
  typedef T Cont_t;
  typedef typename T::iterator Iter_t;
  typedef typename T::value_type Value_t;
  typedef ROOT::Detail::TCollectionProxyInfo::Environ<Iter_t> Env_t;
  typedef Env_t *PEnv_t;
  typedef Cont_t *PCont_t;
  typedef Value_t *PValue_t;

  static void resize(void *obj, size_t n)
  {
    PCont_t c = PCont_t(obj);
    c->Allocate(n);
    ArrayType<T>::construct(c->begin(), n);
  }

  static void *feed(void *from, void *to, size_t size)
  {
    PCont_t c  = PCont_t(to);
    PValue_t m = PValue_t(from);
    if ((size_t)c->size() != size) c->Allocate(size);
    for (size_t i = 0; i < size; ++i, ++m)
      (*c)[i] = *m;
    return 0;
  }
  static int value_offset() { return 0; }
};

template <class T, template <typename> class Q>
void SetCollectionProxyImpl(const char *altname)
{

  auto info = ROOT::Detail::TCollectionProxyInfo::Generate(Q<T>());

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
    (void)type; // silence unused variable warning.
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

template <class T>
void SetCollectionProxyPushBack(const char *altname)
{
  SetCollectionProxyImpl<T, ROOT::Detail::TCollectionProxyInfo::Pushback>(altname);
}

template <class T>
void SetCollectionProxyArray(const char *altname)
{
  SetCollectionProxyImpl<T, ForArray>(altname);
}

class EmptyClassStreamer : public TClassStreamer {
public:
  TClassStreamer *Generate() const override
  {
    // Virtual copy constructor.
    return new EmptyClassStreamer(*this);
  }

  ~EmptyClassStreamer() override{};

  virtual void operator()(TBuffer &, void *) override {}
  virtual void Stream(TBuffer &, void *, const TClass *) override {}
};

void RootPersistencyProxy()
{
  // SetCollectionProxy<std::vector<const vecgeom::VPlacedVolume *>>();
  // SetCollectionProxy<std::vector<double>>();
  // SetCollectionProxy<std::vector<double *>>();
  // SetCollectionProxy<vecgeom::Array<double>>();

  SetCollectionProxyPushBack<vecgeom::Vector<const vecgeom::VPlacedVolume *>>("list<const vecgeom::VPlacedVolume *>");
  SetCollectionProxyPushBack<vecgeom::Vector<double>>("list<double>");
  SetCollectionProxyPushBack<vecgeom::Vector<float>>("list<float>");
  SetCollectionProxyPushBack<vecgeom::Vector<vecgeom::PolyconeSection>>("list<vecgeom::PolyconeSection>");
  // This is is unlikely to be stored correctly as ROOT can't know the number of double pointed to by the double*.
  // (it is *unlikely* to be just 1 as it is assumed for objects)
  // SetCollectionProxy<vecgeom::Vector<double *>>("vector<double *>");

  // auto info = ROOT::Detail::TCollectionProxyInfo::Generate(ForArray<vecgeom::Array<double>>());
  SetCollectionProxyArray<vecgeom::Array<float>>("list<float>");
  SetCollectionProxyArray<vecgeom::Array<double>>("list<double>");
  SetCollectionProxyArray<vecgeom::Array<vecgeom::ZSegment>>("list<vecgeom::ZSegment>");
  SetCollectionProxyArray<vecgeom::Array<bool>>("list<bool>");

  TClass *cl = TClass::GetClass("vecgeom::AlignedBase");
  if (cl) cl->AdoptStreamer(new EmptyClassStreamer());
}
