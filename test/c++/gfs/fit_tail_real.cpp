// Copyright (c) 2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2018 Centre national de la recherche scientifique (CNRS)
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
// Authors: Olivier Parcollet, Nils Wentzell

#include <triqs/test_tools/gfs.hpp>

using namespace nda;

TEST(FitTailReal, Basic) { // NOLINT

  nda::clef::placeholder<0> w_;
  double omega_max = 10;
  int L            = 201;

  auto w_mesh = mesh::refreq{-omega_max, omega_max, L};

  // Set the fraction of mesh points to use for the tail fit
  double tail_fraction = 0.3;
  w_mesh.get_tail_fitter(tail_fraction);

  auto gw   = gf<refreq>{w_mesh, {1, 1}};
  auto gw_s = gf<refreq, scalar_valued>{w_mesh};

  // Initialize the Green functions
  array<dcomplex, 1> c{0, 1, 3, 5};
  gw(w_) << c(0) + c(1) / (w_ + 1e-14) + c(2) / (w_ * w_ + 1e-14) + c(3) / (w_ * w_ * w_ + 1e-14);
  gw_s(w_) << c(0) + c(1) / (w_ + 1e-14) + c(2) / (w_ * w_ + 1e-14) + c(3) / (w_ * w_ * w_ + 1e-14);

  // ==== Fix only the 0th moment to 0
  {
    auto known_moments   = array<dcomplex, 3>{{{0.0}}};
    auto known_moments_s = array<dcomplex, 1>{0.0};

    // Get the tail using least square fitting
    auto [tail, err]     = fit_tail(gw, known_moments);
    auto [tail_s, err_s] = fit_tail(gw_s, known_moments_s);

    EXPECT_ARRAY_NEAR(c, tail(range(0, 4), 0, 0), 1e-9);
    EXPECT_ARRAY_NEAR(c, tail_s(range(0, 4)), 1e-9);
  }

  // ==== Now fix both the 0th and 1st moment
  {
    auto known_moments   = array<dcomplex, 3>{{{0.0}}, {{1.0}}};
    auto known_moments_s = array<dcomplex, 1>{0.0, 1.0};

    // Get the tail using least square fitting
    auto [tail, err]     = fit_tail(gw, known_moments);
    auto [tail_s, err_s] = fit_tail(gw_s, known_moments_s);

    EXPECT_ARRAY_NEAR(c, tail(range(4), 0, 0), 1e-9);
    EXPECT_ARRAY_NEAR(c, tail_s(range(4)), 1e-9);
  }
}

// ------------------------------------------------------------------------------

TEST(FitTailReal, Complex) { // NOLINT

  nda::clef::placeholder<0> w_;
  double omega_max = 100;
  int L            = 201;

  auto gw = gf<refreq>{{-omega_max, omega_max, L}, {1, 1}};

  // Initialize the Green functions
  auto a = dcomplex(1.0, 0.4);
  gw(w_) << 1 / (w_ - a);

  // Fix both the 0th and 1st moment
  auto known_moments = array<dcomplex, 3>{{{0.0}}, {{1.0}}};

  auto [tail, err] = fit_tail(gw, known_moments);

  auto tail_exact = array<dcomplex, 1>{dcomplex(0, 0), dcomplex(1, 0), a, std::pow(a, 2), std::pow(a, 3)};

  EXPECT_ARRAY_NEAR(tail_exact, tail(range(5), 0, 0), 1e-7);
}

// ------------------------------------------------------------------------------

TEST(FitTailReal, Multivar) { // NOLINT

  nda::clef::placeholder<0> k_;
  nda::clef::placeholder<1> w_;

  int N_k          = 4;
  double omega_max = 100;
  int L            = 201;

  auto BL        = bravais_lattice{matrix<double>{{1, 0}, {0, 1}}};
  auto BZ        = brillouin_zone{BL};
  auto k_mesh    = mesh::brzone(BZ, N_k);
  auto w_mesh    = mesh::refreq{-omega_max, omega_max, L};
  auto prod_mesh = k_mesh * w_mesh;

  auto g = gf{prod_mesh, {1, 1}};

  // Initialize the Multivariable Green functions
  auto ieta = dcomplex(0.0, 0.001);
  g(k_, w_) << 1 / (w_ - cos(k_[0]) * cos(k_[1]) - ieta);

  // Fix both the 0th and 1st moment
  auto known_moments                 = array<dcomplex, 4>(2, N_k * N_k, 1, 1);
  known_moments(0, range::all, 0, 0) = 0.0;
  known_moments(1, range::all, 0, 0) = 1.0;

  // Fit for all k-points. Resulting shape is (N_orders, N_k * N_k, 1, 1)
  auto [tail, err] = fit_tail<1>(g, known_moments);

  // Calculate the exact tail
  auto tail_exact = array<dcomplex, 2>(5, N_k * N_k);
  for (auto k : k_mesh) {
    dcomplex pole                          = cos(k[0]) * cos(k[1]) + ieta;
    tail_exact(range::all, k.data_index()) = array<dcomplex, 1>{dcomplex(0.0, 0.0), dcomplex(1.0, 0.0), pole, std::pow(pole, 2), std::pow(pole, 3)};
  }

  EXPECT_ARRAY_NEAR(tail_exact, tail(range(5), range::all, 0, 0), 1e-7);
}

MAKE_MAIN;
