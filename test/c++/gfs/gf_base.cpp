// Copyright (c) 2015-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2015-2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2023 Simons Foundation
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
// Authors: Michel Ferrero, Olivier Parcollet, Nils Wentzell

#undef NDEBUG
#include <triqs/test_tools/gfs.hpp>

TEST(Gf, Base) {

  nda::clef::placeholder<0> om_;
  double beta = 1;

  auto G  = gf<imfreq>{{beta, Fermion}, {2, 2}};
  auto Gc = G;
  auto G3 = G;

  auto Gt = gf<imtime>{{beta, Fermion, 100}, {2, 2}};

  auto Gv = G();

  EXPECT_ARRAY_ZERO(G(0));

  Gv[0] = 20;

  EXPECT_ARRAY_NEAR(Gv(0), G(0));
  EXPECT_ARRAY_NEAR(G(0), matrix<double>{{20, 0}, {0, 20}});

  Gv[0] = 0;

  auto Gv2 = slice_target(G(), range(0, 1), range(0, 1));
  Gv2[0]   = 10;
  EXPECT_ARRAY_NEAR(G(0), matrix<double>{{10, 0}, {0, 0}});

  EXPECT_ARRAY_NEAR(eval(G(om_), om_ = 0), G(0));
  EXPECT_ARRAY_NEAR(eval(Gv(om_), om_ = 0), G(0));

  // Lazy assig
  G(om_) << (0.2 + om_ + 2.1);
  EXPECT_ARRAY_NEAR(G(0), matrix<dcomplex>{{2.3 + 3.1415926535897931i, 0i}, {0i, 2.3 + 3.1415926535897931i}});

  Gv(om_) << 1 / (om_ + 2.3);

  dcomplex z = 1 / (M_PI * 1i + 2.3);
  EXPECT_ARRAY_NEAR(G(0), matrix<dcomplex>{{z, 0i}, {0i, z}});

  // copy
  Gc = G;
  EXPECT_GF_NEAR(G, Gc);

  // operations on gf
  G3 = G + Gc;

  // test for density
  EXPECT_ARRAY_NEAR(density(G3), matrix<double>{{1.8177540779781256, 0.0}, {0.0, 1.8177540779781256}});

  rw_h5(G, "G");

  {
    auto G0w     = gf<imfreq, scalar_valued>{{beta, Fermion, 100}};
    auto D0w     = gf<imfreq, scalar_valued>{{beta, Boson, 100}};
    auto Sigma_w = gf<imfreq, scalar_valued>{{beta, Fermion, 100}};
    G0w(om_) << 1 / (om_ - 3);

    for (auto nu : D0w.mesh()) Sigma_w(om_) << 2 * G0w(om_ - nu);
  }

  {
    auto G0w     = gf<imfreq, matrix_valued>{{beta, Fermion, 100}, {1, 1}};
    auto D0w     = gf<imfreq, matrix_valued>{{beta, Boson, 100}, {1, 1}};
    auto Sigma_w = gf<imfreq, matrix_valued>{{beta, Fermion, 100}, {1, 1}};
    G0w(om_) << 1 / (om_ - 3);

    for (auto nu : D0w.mesh()) Sigma_w(om_) << 2 * G0w(om_ - nu);
  }
}

TEST(Gf, RealImag) {

  double beta = 1;
  auto g      = gf<imfreq, matrix_valued>{{beta, Fermion, 100}, {2, 2}};

  nda::clef::placeholder<0> iw_;
  g(iw_) << 1 / (iw_ - 2.0);

  auto g_re_im = gf<imfreq, matrix_valued>{real(g) + 1i * imag(g)};
  EXPECT_GF_NEAR(g, g_re_im);
}

// Test the technique to avoid the infinity
TEST(Gf, EvaluatorMatrix) {

  double beta = 1;
  auto g      = gf<imfreq, matrix_valued>{{beta, Fermion, 100}, {2, 2}};

  nda::clef::placeholder<0> iw_;
  g(iw_) << 1 / (iw_ + 2.3);
  std::cout << g[0] << std::endl;

  auto f     = matsubara_freq{120, beta, Fermion};
  auto exact = matrix<dcomplex>{{1 / (f + 2.3), 0i}, {0i, 1 / (f + 2.3)}};

  EXPECT_ARRAY_NEAR(g(f), exact, 1e-14);
}

// Test the technique to avoid the infinity
TEST(Gf, EvaluatorScalar) {

  double beta = 1;
  auto g      = gf<imfreq, scalar_valued>{{beta, Fermion, 100}};
  static_assert(std::is_reference_v<decltype(g[*g.mesh().begin()])>);

  nda::clef::placeholder<0> om_;
  g(om_) << 1 / (om_ + 2.3);

  auto f     = matsubara_freq{120, beta, Fermion};
  auto exact = 1 / (f + 2.3);

  EXPECT_COMPLEX_NEAR(g(f), exact, 1e-14);
}

// Test the technique to avoid the infinity
TEST(Gf, PhNoInfinity) {

  double beta = 1;
  auto g      = gf<imfreq, matrix_valued>{{beta, Fermion}, {2, 2}};

  nda::clef::placeholder<0> om_;
  g(om_) << 1 / (om_ + 2.3);
}

// Test the technique to avoid the infinity
TEST(Gf, PhNoInfinity_tau) {

  double beta = 1, a = 1;
  auto g = gf<imtime, matrix_valued>{{beta, Fermion, 10000}, {2, 2}};

  using std::exp;
  nda::clef::placeholder<0> tau_;
  g(tau_) << exp(-a * tau_) / (1 + exp(-beta * a));
}

TEST(Gf, ZeroM) {

  double beta = 1;
  auto g      = gf<imfreq>{{beta, Fermion}, {2, 2}};
  EXPECT_ARRAY_NEAR(zeros<dcomplex>(g.target_shape()), matrix<double>{{0, 0}, {0, 0}});
}
TEST(Gf, ZeroS) {

  double beta = 1;
  auto g      = gf<imfreq, scalar_valued>{{beta, Fermion}};
  EXPECT_COMPLEX_NEAR(zeros<dcomplex>(g.target_shape()), 0);
}

TEST(Gf, SliceTargetScalar) {
  double beta = 1;
  auto g      = gf<imfreq>{{beta, Fermion}, {2, 2}};

  auto gs = slice_target_to_scalar(g, 0, 0);
}

TEST(Gf, TargetSpaceLoop) {
  double beta = 1;
  auto g1     = gf<imfreq>{{beta, Fermion}, {2, 2}};
  auto g2     = gf{g1};

  nda::clef::placeholder<0> iw_;
  nda::clef::placeholder<1> i_;
  nda::clef::placeholder<2> j_;

  static_assert(!mesh::is_product<decltype(auto{g1.mesh()})>);
  g1[iw_](i_, j_) << 1. / (iw_ + i_ + j_ * 0.1);

  for (auto iw : g2.mesh())
    for (auto [i, j] : g2.target_indices()) g2[iw](i, j) = 1. / (iw + i + j * 0.1);

  EXPECT_GF_NEAR(g1, g2);
}

TEST(Gf, MeshCheck) {
  double beta = 1;

  auto iw_mesh     = mesh::imfreq{beta, Fermion, 5};
  auto iw_mesh_big = mesh::imfreq{beta, Fermion, 10};

  auto g1 = gf<imfreq>{iw_mesh, {1, 1}};
  auto g2 = gf<imfreq>{iw_mesh_big, {1, 1}};

  EXPECT_DEATH(g1[*g2.mesh().begin()], "mesh_hash.* violated");
  EXPECT_DEATH(g2[*g1.mesh().begin()], "mesh_hash.* violated");
}

TEST(Gf, EvalSlice) {
  auto t_mesh = mesh::refreq({-10., 10., 100});
  gf<refreq, matrix_valued> g(t_mesh, {2, 2});

  auto g5 = g(5.0);
  g5(range::all, 0);
  auto g_eval_slice = array<dcomplex, 1>{g(5.0)(range::all, 0)};
  EXPECT_EQ(g_eval_slice, (array<dcomplex, 1>{0.0, 0.0}));
}

TEST(Gf, TransposeComparison) {

  auto m    = mesh::retime{0, 10, 99};
  auto Ginv = gf<retime>{m, {2, 2}};

  placeholder<0> w;
  placeholder<1> i;
  placeholder<2> j;

  Ginv[w](i, j) << w - (i + j) + 0.001i;

  EXPECT_EQ(Ginv[0], transpose(Ginv[0]));
  EXPECT_EQ(Ginv[10], transpose(Ginv[10]));
}

MAKE_MAIN;
