#include <triqs/test_tools/gfs.hpp>
#include <triqs/gfs.hpp>

using namespace triqs::clef;
using namespace triqs::lattice;

TEST(Gf, ProdGfTail) {

  double beta = 1;
  placeholder<0> iw_; 

  auto G = gf<imfreq, matrix_valued>{{beta, Fermion, 50}, {2,2}};
  G(iw_) << 1.0 / iw_ + 0.25 / ( iw_ * iw_ * iw_ ); 

  auto GG = G; 

  std::cout << " Tail of G \n " << GG.singularity() << "\n"; 

  GG *= G; 

  std::cout << " Tail of product G*G \n " << GG.singularity() << "\n"; 

  auto GG_res = G; 
  GG_res(iw_) << ( 1.0 / iw_ + 0.25 / ( iw_ * iw_ * iw_ ) ) * ( 1.0 / iw_ + 0.25 / ( iw_ * iw_ * iw_ ) );

  std::cout << " Expected Tail \n " << GG_res.singularity() << "\n"; 

  EXPECT_GF_NEAR( GG, GG_res ); 
}

TEST(Gf, AutoAssignMatrixGf2_product) {

 double beta = 2.3;
 auto g2 = gf<cartesian_product<imfreq, imfreq>, matrix_valued>{{{beta, Fermion, 10}, {beta, Fermion, 10}}, {2, 2}};

 placeholder<0> i_;
 placeholder<1> j_;
 placeholder<3> om_;
 placeholder<4> nu_;

 g2(om_, nu_)(i_, j_) << i_ + j_ + om_ * nu_;

 // CHECK
 for (int i = 0; i < 2; ++i)
  for (int j = 0; j < 2; ++j)
   for (int om = 0; om < 10; ++om)
    for (int nu = 0; nu < 10; ++nu) {
     auto xom = ((2 * om + 1) * M_PI * 1_j / beta);
     auto xnu = ((2 * nu + 1) * M_PI * 1_j / beta);
     EXPECT_CLOSE(g2.data()(10+om, 10 +nu, i, j), i + j + xom * xnu);
    }
}

TEST(Gf, AutoAssignMatrixGf2_sum_product) {

 double beta = 2.3;
 auto g2 = gf<cartesian_product<imfreq, imfreq>, matrix_valued>{{{beta, Fermion, 10}, {beta, Fermion, 10}}, {2, 2}};

 placeholder<0> i_;
 placeholder<1> j_;
 placeholder<3> om_;
 placeholder<4> nu_;

 g2(om_, nu_)(i_, j_) << i_ + j_ + om_ * (nu_ + om_);

 // CHECK
 for (int i = 0; i < 2; ++i)
  for (int j = 0; j < 2; ++j)
   for (int om = 0; om < 10; ++om)
    for (int nu = 0; nu < 10; ++nu) {
     auto xom = ((2 * om + 1) * M_PI * 1_j / beta);
     auto xnu = ((2 * nu + 1) * M_PI * 1_j / beta);
     EXPECT_CLOSE(g2.data()(10+om, 10 +nu, i, j), i + j + xom * (xnu + xom));
    }
}


MAKE_MAIN; 
