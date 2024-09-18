// Copyright (c) 2013-2014 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2013-2014 Centre national de la recherche scientifique (CNRS)
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

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs;
using namespace triqs::gfs;
using namespace std::complex_literals;

template <typename T> struct check {
  //  static_assert(std::is_same<typename T::view_type, view_or_type_t<T>>::value, "err");
  static_assert(std::is_same<typename T::regular_type, regular_t<T>>::value, "err");
};

TEST(ViewTools, TypeTraits) {
  check<array<int, 1>>();
  check<array_view<int, 1>>();
  check<matrix<int>>();
  check<matrix_view<int>>();
}

TEST(ViewTools, PosFreqView) {

  double beta  = 1;
  int n_iw     = 100;
  auto iw_mesh = mesh::imfreq{beta, Fermion, n_iw};

  // -- Build and Init Green function
  auto g_iw = gf<imfreq, matrix_valued>{iw_mesh, {1, 1}};
  nda::clef::placeholder<0> iw_;
  g_iw(iw_) << 1.0 / (iw_ + 1.0 + 1i);

  // -- Test positive_freq_view()
  auto g_iw_cv = gf_const_view(g_iw);

  auto g_iw_pfv   = positive_freq_view(g_iw);
  auto g_iw_pfv_v = positive_freq_view(g_iw()); // From rvalue view
  //auto g_iw_pfv_v = positive_freq_view(gf{g_iw}); // From temporary gf: This should fail compilation
  auto g_iw_pfv_cv = positive_freq_view(g_iw_cv); // From lvalue view

  auto pos_freq_data_view = g_iw.data()(range(n_iw, 2 * n_iw), ellipsis());
  EXPECT_ARRAY_NEAR(pos_freq_data_view, g_iw_pfv.data());
  EXPECT_ARRAY_NEAR(pos_freq_data_view, g_iw_pfv_v.data());
  EXPECT_ARRAY_NEAR(pos_freq_data_view, g_iw_pfv_cv.data());
}

MAKE_MAIN;
