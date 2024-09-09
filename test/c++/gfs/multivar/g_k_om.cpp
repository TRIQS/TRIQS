// Copyright (c) 2014-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2014-2018 Centre national de la recherche scientifique (CNRS)
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

#ifndef NDA_DEBUG
#define NDA_DEBUG
#endif

#include <triqs/test_tools/gfs.hpp>
using namespace triqs::clef;
using namespace triqs::lattice;
using triqs::clef::placeholder;

template <typename Function, typename Mesh> auto sum_gf(Function const &f, Mesh const &m) {
  auto res = make_regular(0 * f(*(m.begin())));
  for (auto const &x : m) //
    res = res + f(x);
  return res;
}

namespace nda::clef {
  CLEF_MAKE_FNT_LAZY(sum_gf);
}

placeholder<0> k_;
placeholder<1> w_;
placeholder<2> wp_;
using gf_bz_imfreq_mat = gf<prod<brzone, imfreq>, matrix_valued>;
using gf_bz_imtime_mat = gf<prod<brzone, imtime>, matrix_valued>;

auto _ = all_t{};

// ------------------------------------------------------------

TEST(Gf, GkOm) {

  double beta = 1;
  int n_bz    = 20;

  auto bz    = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
  auto g_eps = gf<brzone>{{bz, n_bz}, {1, 1}};

  auto G     = gf_bz_imfreq_mat{{{bz, n_bz}, {beta, Fermion, 100}}, {1, 1}};
  auto Sigma = gf<imfreq, matrix_valued>{{beta, Fermion, 100}, {1, 1}};
  Sigma()    = 0;

  auto eps_k = -2 * (cos(k_[0]) + cos(k_[1]));

  G(k_, w_) << 1 / (w_ - eps_k - 1 / (w_ + 2));

  G[k_, w_] << G[k_, w_] + 1 / (w_ - eps_k - 1 / (w_ + 2));

  G[k_, wp_] << 2 * G[k_, wp_];

  G[k_, w_] << 2 * Sigma(w_);

  G[k_, w_] << 2 * Sigma[w_];

  G[k_, w_] << 1 / (w_ - eps_k - 1 / (w_ + 2));

  rw_h5(G, "ess_g_k_om.h5", "g");
}

// -------------------------------------------------------------
using K_t = std::array<double, 3>;

TEST(Gkom, Eval) {
  double beta = 5;

  int n_k = 40;
  int n_w = 5;

  auto bz = brillouin_zone{bravais_lattice{nda::eye<double>(2)}};
  auto g  = gf<prod<brzone, imfreq>, scalar_valued>{{{bz, n_k}, {beta, Fermion, n_w}}};

  static_assert(std::tuple_size_v<prod<brzone, imfreq>> == 2);
  static_assert(std::tuple_size_v<prod<brzone, imfreq>::mesh_point_t> == 2);
  //Gk(k_) << -2 * (cos(k_[0]) + cos(k_[1]));
  for (auto [k, w] : g.mesh()) g[k, w] = 1 / (w + 3 + 2 * (cos(k[0]) + cos(k[1])));

  // reevaluate on the mesh itself.
  for (auto [k, w] : g.mesh()) {
    EXPECT_NEAR(k[0], k.index()[0] / double(n_k) * 2 * M_PI, 1.e-14);
    EXPECT_NEAR(k[1], k.index()[1] / double(n_k) * 2 * M_PI, 1.e-14);
    dcomplex res = 1 / (w + 3 + 2 * (cos(k[0]) + cos(k[1])));
    EXPECT_COMPLEX_NEAR(g(K_t{k[0], k[1], 0}, w), res, 1.e-14);
  }

  //// evaluate on a larger grid
  int n_k2 = 3 * n_k;
  for (int nkx = 0; nkx < n_k2; ++nkx)
    for (int nky = 0; nky < n_k2; ++nky)
      for (auto w : std::get<1>(g.mesh())) {
        double kx    = nkx / double(n_k2) * 2 * M_PI;
        double ky    = nky / double(n_k2) * 2 * M_PI;
        dcomplex res = 1 / (w + 3 + 2 * (cos(kx) + cos(ky)));
        EXPECT_COMPLEX_NEAR(g(K_t{kx, ky, 0}, w), res, 0.05 * std::abs(res));
      }
}

// -------------------------------------------------------------

// test g(k, all_t)
TEST(Gkom, EvalSlice) {
  double beta = 5;

  int n_k = 40;
  int n_w = 5;

  auto bz = brillouin_zone{bravais_lattice{nda::eye<double>(2)}};
  auto g  = gf<prod<brzone, imfreq>, scalar_valued>{{{bz, n_k}, {beta, Fermion, n_w}}};

  //Gk(k_) << -2 * (cos(k_[0]) + cos(k_[1]));
  for (auto [k, w] : g.mesh()) g[k, w] = 1 / (w + 3 + 2 * (cos(k[0]) + cos(k[1])));

  auto pi_pi = K_t{M_PI, M_PI, 0};

  auto gw = g(pi_pi, all_t{});

  EXPECT_TRUE(std::get<1>(g.mesh()) == gw.mesh());
  for (auto w : gw.mesh()) EXPECT_COMPLEX_NEAR(gw[w], g(pi_pi, w), 1.e-12);
}

// ----------------------------------------------------------------

TEST(Gkom, DLR) {
  auto beta    = 5.0;
  auto mu      = 0.5;
  auto wmax    = 10 * beta;
  auto eps_DLR = 1e-10;
  auto const iw_mesh = dlr_imfreq{beta, Fermion, wmax, eps_DLR};

  auto bz = brillouin_zone{bravais_lattice{nda::eye<double>(2)}};
  auto k_mesh = mesh::brzone{bz, /*n_k =*/ 40};

  auto g = gf<prod<brzone, dlr_imfreq>, scalar_valued>{{k_mesh, iw_mesh}};
  auto eps = [](auto k) { return -2 * (cos(k[0]) + cos(k[1])); };
  for (auto [k, iw] : g.mesh()) g[k, iw] = 1 / (iw + mu - eps(k));

  auto gktau = make_gf_from_fourier<1>(g);

  auto onefermion = [&](auto tau, auto w) { return -std::exp(-w * tau) / (1.0 + std::exp(-beta * w)); };
  for (auto [k, tau] : gktau.mesh()) { EXPECT_COMPLEX_NEAR(gktau[k, tau], onefermion(tau, -mu + eps(k)), 100*eps_DLR); }
}

MAKE_MAIN;
