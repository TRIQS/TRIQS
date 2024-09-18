// Copyright (c) 2017 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2017 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2019-2020 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Olivier Parcollet, Nils Wentzell

#include <triqs/test_tools/gfs.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

using namespace nda::clef;
using namespace triqs::lattice;

TEST(Gf, ProdGfTail) {

  double beta = 1;
  placeholder<0> iw_;

  auto G = gf<imfreq, matrix_valued>{{beta, Fermion, 50}, {2, 2}};
  G(iw_) << 1.0 / iw_ + 0.25 / (iw_ * iw_ * iw_);

  auto GG = G;

  GG *= G;

  auto GG_res = G;
  GG_res(iw_) << (1.0 / iw_ + 0.25 / (iw_ * iw_ * iw_)) * (1.0 / iw_ + 0.25 / (iw_ * iw_ * iw_));

  EXPECT_GF_NEAR(GG, GG_res);
}

TEST(Gf, AutoAssignMatrixGf2_product) {

  double beta = 2.3;
  auto g2     = gf<prod<imfreq, imfreq>, matrix_valued>{{{beta, Fermion, 10}, {beta, Fermion, 10}}, {2, 2}};

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
          auto xom = ((2 * om + 1) * M_PI * 1i / beta);
          auto xnu = ((2 * nu + 1) * M_PI * 1i / beta);
          EXPECT_CLOSE(g2.data()(10 + om, 10 + nu, i, j), i + j + xom * xnu);
        }
}

TEST(Gf, AutoAssignMatrixGf2_sum_product) {

  double beta = 2.3;
  auto g2     = gf<prod<imfreq, imfreq>, matrix_valued>{{{beta, Fermion, 10}, {beta, Fermion, 10}}, {2, 2}};

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
          auto xom = ((2 * om + 1) * M_PI * 1i / beta);
          auto xnu = ((2 * nu + 1) * M_PI * 1i / beta);
          EXPECT_CLOSE(g2.data()(10 + om, 10 + nu, i, j), i + j + xom * (xnu + xom));
        }
}

MAKE_MAIN;
