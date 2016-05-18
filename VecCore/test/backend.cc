#include <VecCore/VecCore>

#include <type_traits>
#include <gtest/gtest.h>

using namespace testing;

#if defined(GTEST_HAS_TYPED_TEST) && defined(GTEST_HAS_TYPED_TEST_P)

template <class Backend>
using FloatTypes = Types
<
  typename Backend::Float_v,
  typename Backend::Double_v
>;

template <class Backend>
using IntTypes = Types
<
  typename Backend::Int16_v,
  typename Backend::UInt16_v,
  typename Backend::Int32_v,
  typename Backend::UInt32_v
>;

template <class Backend>
using VectorTypes = Types
<
  typename Backend::Float_v,
  typename Backend::Double_v,
  typename Backend::Int16_v,
  typename Backend::UInt16_v,
  typename Backend::Int32_v,
  typename Backend::UInt32_v
>;

///////////////////////////////////////////////////////////////////////////////

template <class T> class VectorTypeTest : public Test {
public:
  using Scalar_t = typename vecCore::ScalarType<T>::Type;
  using Vector_t = T;
};

///////////////////////////////////////////////////////////////////////////////

template <class T> class ConstructorTest : public VectorTypeTest<T> {};

TYPED_TEST_CASE_P(ConstructorTest);

TYPED_TEST_P(ConstructorTest, Default)
{
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE((std::is_constructible<Vector_t>::value));
}

TYPED_TEST_P(ConstructorTest, Copy)
{
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE((std::is_copy_constructible<Vector_t>::value));
}

TYPED_TEST_P(ConstructorTest, Move)
{
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE((std::is_move_constructible<Vector_t>::value));
}

TYPED_TEST_P(ConstructorTest, FromScalar)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE((std::is_constructible<Vector_t, Scalar_t>::value));
}

TYPED_TEST_P(ConstructorTest, FromRefToScalar)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE((std::is_constructible<Vector_t, const Scalar_t &>::value));
}

TYPED_TEST_P(ConstructorTest, FromPtrToScalar)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  Scalar_t tmp[vecCore::VectorSize<Vector_t>()];
  Scalar_t *addr = &tmp[0];

  for (vecCore::UInt_s i = 0; i < vecCore::VectorSize<Vector_t>(); i++)
    tmp[i] = static_cast<Scalar_t>(i);

  Vector_t x = vecCore::FromPtr<Vector_t>(addr);

  for (vecCore::UInt_s i = 0; i < vecCore::VectorSize<Vector_t>(); i++)
    EXPECT_TRUE(!vecCore::MaskEmpty(x == Vector_t(Scalar_t(i))));
}

REGISTER_TYPED_TEST_CASE_P(ConstructorTest, Default, Copy, Move, FromScalar, FromRefToScalar, FromPtrToScalar);

///////////////////////////////////////////////////////////////////////////////

template <class T> class ArithmeticsTest : public VectorTypeTest<T> {};

TYPED_TEST_CASE_P(ArithmeticsTest);

TYPED_TEST_P(ArithmeticsTest, Addition)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t lhs(static_cast<Scalar_t>(1));
  Vector_t rhs(static_cast<Scalar_t>(1));
  Vector_t res(static_cast<Scalar_t>(2));

  EXPECT_TRUE(vecCore::MaskFull(res == (lhs + rhs)));
}

TYPED_TEST_P(ArithmeticsTest, Subtraction)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t lhs(static_cast<Scalar_t>(3));
  Vector_t rhs(static_cast<Scalar_t>(1));
  Vector_t res(static_cast<Scalar_t>(2));

  EXPECT_TRUE(vecCore::MaskFull(res == (lhs - rhs)));
}

TYPED_TEST_P(ArithmeticsTest, Multiplication)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t lhs(static_cast<Scalar_t>(2));
  Vector_t rhs(static_cast<Scalar_t>(2));
  Vector_t res(static_cast<Scalar_t>(4));

  EXPECT_TRUE(vecCore::MaskFull(res == (lhs * rhs)));
}

TYPED_TEST_P(ArithmeticsTest, Division)
{
  using Scalar_t = typename TestFixture::Scalar_t;
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t lhs(static_cast<Scalar_t>(10));
  Vector_t rhs(static_cast<Scalar_t>(5));
  Vector_t res(static_cast<Scalar_t>(2));

  EXPECT_TRUE(vecCore::MaskFull(res == (lhs / rhs)));
}

REGISTER_TYPED_TEST_CASE_P(ArithmeticsTest, Addition, Subtraction, Multiplication, Division);

///////////////////////////////////////////////////////////////////////////////

template <class T> class VectorInterfaceTest : public VectorTypeTest<T> {};

TYPED_TEST_CASE_P(VectorInterfaceTest);

TYPED_TEST_P(VectorInterfaceTest, VectorSize)
{
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE(vecCore::VectorSize<Vector_t>() > 0);
}

TYPED_TEST_P(VectorInterfaceTest, VectorSizeVariable)
{
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t x;

  EXPECT_TRUE(vecCore::VectorSize(x) > 0);
}

TYPED_TEST_P(VectorInterfaceTest, StoreToPtr) {
  using Vector_t = typename TestFixture::Vector_t;
  using Scalar_t = typename TestFixture::Scalar_t;

  auto kVS = vecCore::VectorSize<Vector_t>();
  auto N = 2 * kVS;
  Scalar_t input[N];
  Scalar_t output[N];

  // init input; output
  for (vecCore::UInt_s i = 0; i < N; ++i) {
    input[i] = 2 * i;
    output[i] = 0;
  }

  // transfer to output via FromPtr; Store sequence
  for (vecCore::UInt_s i = 0; i < 2; ++i) {
    Vector_t tmp(vecCore::FromPtr<Vector_t>(&input[i * kVS]));
    vecCore::Store<Vector_t>(tmp, &output[i * kVS]);
  }

  // assert input == output
  for (vecCore::UInt_s i = 0; i < N; ++i)
    EXPECT_EQ(input[i], output[i]);
}

TYPED_TEST_P(VectorInterfaceTest, VectorLaneRead) {
  using Vector_t = typename TestFixture::Vector_t;
  using Scalar_t = typename TestFixture::Scalar_t;

  auto kVS = vecCore::VectorSize<Vector_t>();
  Scalar_t input[kVS];

  for (vecCore::UInt_s i = 0; i < kVS; ++i) {
    input[i] = i;
  }

  Vector_t tmp(vecCore::FromPtr<Vector_t>(&input[0]));

  for (vecCore::UInt_s i = 0; i < kVS; ++i)
    EXPECT_EQ(input[i], vecCore::LaneAt<Vector_t>(tmp, i));
}

TYPED_TEST_P(VectorInterfaceTest, MaskLaneRead) {
  using Vector_t = typename TestFixture::Vector_t;
  using Scalar_t = typename TestFixture::Scalar_t;

  auto kVS = vecCore::VectorSize<Vector_t>();
  Scalar_t input[kVS];

  for (vecCore::UInt_s i = 0; i < kVS; ++i) {
    input[i] = (i % 2 == 0) ? i : -i;
  }

  Vector_t tmp(vecCore::FromPtr<Vector_t>(&input[0]));

  auto mask = tmp > Vector_t(0);

  for (vecCore::UInt_s i = 0; i < kVS; ++i)
    EXPECT_EQ(input[i] > Scalar_t(0), vecCore::MaskLaneAt(mask, i));
}

TYPED_TEST_P(VectorInterfaceTest, Gather) {
  using Vector_t = typename TestFixture::Vector_t;
  using Scalar_t = typename TestFixture::Scalar_t;
  using Index_v = typename vecCore::Index_v<Vector_t>;
  using Index_t = typename vecCore::ScalarType<Index_v>::Type;

  size_t N = vecCore::VectorSize<Vector_t>();

  Scalar_t input[16 * N];
  Index_t  index[16 * N];

  for (vecCore::UInt_s i = 0; i < 16*N; ++i) {
    input[i] = (i % 2 == 0) ? 0 : 1, index[i] = i;
  }

  Vector_t x, y;
  Index_v idx1(vecCore::FromPtr<Index_v>(&index[0]));
  Index_v idx2(vecCore::FromPtr<Index_v>(&index[N]));

  x = vecCore::FromPtr<Vector_t>(&input[0]);
  y = vecCore::Gather<Vector_t>(&input[0], idx1);

  EXPECT_TRUE(vecCore::MaskFull(x == y));

  x = vecCore::FromPtr<Vector_t>(&input[N]);
  y = vecCore::Gather<Vector_t>(&input[0], idx2);

  EXPECT_TRUE(vecCore::MaskFull(x == y));
  x = Vector_t(0);
  y = vecCore::Gather<Vector_t>(&input[0], 2*idx1);

  EXPECT_TRUE(vecCore::MaskFull(x == y));

  idx2 = 2*idx1+1;
  y = vecCore::Gather<Vector_t>(&input[0], idx2);

  EXPECT_TRUE(vecCore::MaskFull(x < y));

  idx2 = 8*idx1+1;
  y = vecCore::Gather<Vector_t>(&input[0], idx2);

  EXPECT_TRUE(vecCore::MaskFull(x < y));
}

REGISTER_TYPED_TEST_CASE_P(VectorInterfaceTest, VectorSize, VectorSizeVariable, StoreToPtr, VectorLaneRead,
                           MaskLaneRead, Gather);

///////////////////////////////////////////////////////////////////////////////

template <class T> class VectorMaskTest : public VectorTypeTest<T> {};

TYPED_TEST_CASE_P(VectorMaskTest);

TYPED_TEST_P(VectorMaskTest, Constructor)
{
  using Vector_t = typename TestFixture::Vector_t;

  EXPECT_TRUE(std::is_constructible<vecCore::Mask_v<Vector_t>>::value);
  EXPECT_TRUE(std::is_copy_constructible<vecCore::Mask_v<Vector_t>>::value);
  EXPECT_TRUE(std::is_move_constructible<vecCore::Mask_v<Vector_t>>::value);
}

TYPED_TEST_P(VectorMaskTest, MaskFull)
{
  using Vector_t = typename TestFixture::Vector_t;

  vecCore::Mask_v<Vector_t> tmask(true), fmask(false);

  EXPECT_TRUE(vecCore::MaskFull(tmask));
  EXPECT_TRUE(vecCore::MaskFull(!fmask));
}

TYPED_TEST_P(VectorMaskTest, MaskEmpty)
{
  using Vector_t = typename TestFixture::Vector_t;

  vecCore::Mask_v<Vector_t> tmask(true), fmask(false);

  EXPECT_TRUE(vecCore::MaskEmpty(fmask));
  EXPECT_TRUE(vecCore::MaskEmpty(!tmask));
}

TYPED_TEST_P(VectorMaskTest, MaskAssign) {
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t a(0), b(1);

  vecCore::MaskedAssign(a, a > b, b);
  EXPECT_TRUE(vecCore::MaskFull(a == Vector_t(0)));

  vecCore::MaskedAssign(a, b > a, b);
  EXPECT_TRUE(vecCore::MaskFull(a == b));
}

TYPED_TEST_P(VectorMaskTest, MaskAssign2) {
  using Vector_t = typename TestFixture::Vector_t;
  using Scalar_t = typename TestFixture::Scalar_t;

  auto kVS = vecCore::VectorSize<Vector_t>();
  Scalar_t input[kVS];
  Scalar_t output[kVS];

  for (vecCore::UInt_s i = 0; i < kVS; ++i) {
    input[i] = (i % 2 == 0) ? i : -i;
    output[i] = (input[i] > 0) ? input[i] : 0;
  }

  Vector_t c(vecCore::FromPtr<Vector_t>(&input[0])), dest(0);
  vecCore::MaskedAssign(dest, c > Vector_t(0), c);

  for (vecCore::UInt_s i = 0; i < kVS; ++i)
    EXPECT_EQ(vecCore::LaneAt<Vector_t>(dest, i), output[i]);
}

TYPED_TEST_P(VectorMaskTest, Blend)
{
  using Vector_t = typename TestFixture::Vector_t;

  Vector_t a(0), b(1);

  a = vecCore::Blend(a > b, a, b);
}

REGISTER_TYPED_TEST_CASE_P(VectorMaskTest, Constructor, MaskFull, MaskEmpty, MaskAssign, MaskAssign2, Blend);

#define TEST_BACKEND_P(name, x)                                                            \
  INSTANTIATE_TYPED_TEST_CASE_P(name, ConstructorTest, VectorTypes<vecCore::backend::x>);  \
  INSTANTIATE_TYPED_TEST_CASE_P(name, ArithmeticsTest, VectorTypes<vecCore::backend::x>);  \
  INSTANTIATE_TYPED_TEST_CASE_P(name, VectorMaskTest, VectorTypes<vecCore::backend::x>);   \
  INSTANTIATE_TYPED_TEST_CASE_P(name, VectorInterfaceTest, VectorTypes<vecCore::backend::x>)

#define TEST_BACKEND(x) TEST_BACKEND_P(x, x)

///////////////////////////////////////////////////////////////////////////////

TEST_BACKEND(Scalar);
TEST_BACKEND(ScalarWrapper);
TEST_BACKEND_P(SIMD, SIMD<16>);

#ifdef VECCORE_ENABLE_VC
TEST_BACKEND(VcScalar);
TEST_BACKEND(VcVector);
TEST_BACKEND_P(VcSimdArray, VcSimdArray<16>);
#endif

#else // if !GTEST_HAS_TYPED_TEST
TEST(DummyTest, TypedTestsAreNotSupportedOnThisPlatform) {}
#endif

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
