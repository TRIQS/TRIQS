#include <triqs/test_tools/gfs.hpp>

TEST(Gf, Base) {

  triqs::clef::placeholder<0> t_;
  double beta = 1;

  auto Gt = gf<imtime>{{beta, Fermion, 100}, {2, 2}};

  Gt(t_) << 8 * t_;
  auto gtr           = real(Gt);
  array<double, 3> D = real(Gt.data());
  EXPECT_ARRAY_NEAR(D, gtr.data());
  EXPECT_TRUE(is_gf_real(Gt));

  Gt(t_) << 8i * t_;
  EXPECT_FALSE(is_gf_real(Gt));
}
MAKE_MAIN;
