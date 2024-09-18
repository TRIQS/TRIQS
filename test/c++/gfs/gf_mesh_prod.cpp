// Copyright (c) 2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2019 Simons Foundation
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

#include <triqs/test_tools/gfs.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

using namespace nda::clef;
using namespace triqs::lattice;

TEST(GfMeshProd, Prod) {

  mesh::imfreq m{1, Fermion, 10};
  mesh::imtime p{1, Fermion, 10};

  mesh::prod<imfreq, imfreq> m2{m, m};
  mesh::prod<imfreq, imfreq, imfreq> m3{m, m, m};
  mesh::prod<imfreq, imfreq, imfreq, imfreq> m4{m, m, m, m};
  mesh::prod<imfreq, imfreq, imfreq, imfreq, imfreq> m5{m, m, m, m, m};

  mesh::prod<imfreq, imtime> mp{m, p};
  mesh::prod<imtime, imfreq> pm{p, m};
  mesh::prod<imfreq, imtime, imfreq> mpm{m, p, m};
  mesh::prod<imfreq, imfreq, imtime, imfreq> mmpm{m, m, p, m};
  mesh::prod<imfreq, imfreq, imtime, imfreq, imtime> mmpmp{m, m, p, m, p};

  EXPECT_EQ(m2, m * m);
  EXPECT_EQ(m4, m2 * m2);
  EXPECT_EQ(m3, m * m2);
  EXPECT_EQ(m3, m2 * m);
  EXPECT_EQ(m5, m2 * m * m * m);
  EXPECT_EQ(m5, m * m2 * m2);

  EXPECT_EQ(mp, m * p);
  EXPECT_EQ(mmpm, m2 * pm);
  EXPECT_EQ(mpm, mp * m);
  EXPECT_EQ(mpm, m * pm);
  EXPECT_EQ(mmpmp, m * mp * mp);
}
MAKE_MAIN;
