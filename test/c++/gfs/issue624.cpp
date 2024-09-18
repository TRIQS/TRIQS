// Copyright (c) 2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018 Simons Foundation
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
// Authors: Olivier Parcollet

#include <triqs/test_tools/arrays.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

using namespace triqs::gfs;
using namespace triqs;

TEST(gf, t2) {

  using g_t = gf<imtime, matrix_valued>;

  auto zero = g_t{{1, Fermion, 10}, {1, 1}};
  zero()    = 0;

  auto gg = zero;

  gg = (zero + zero) / 2.0;
}
MAKE_MAIN;
