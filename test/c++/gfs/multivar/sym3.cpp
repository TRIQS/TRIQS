// Copyright (c) 2014-2016 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2014-2016 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2019 Simons Foundation
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
// Authors: Michel Ferrero, Olivier Parcollet

#include <triqs/test_tools/gfs.hpp>
using namespace triqs::lattice;
/*
// In this test, we define a new evaluator for a multivar gf
// to implement the symmetry g(nu, nup, -omega) = conj(g(-nu,-nup,omega))
// The principle is simple : just write an evaluator with a view of the gf with the default evaluator.

struct my_evaluator {
 static constexpr int arity = 3;

 // keep a view of the gf with the default evaluator
 gf_const_view<prod<imfreq, imfreq, imfreq>, scalar_valued, nothing> gv;

 template <typename G> my_evaluator(G *g) : gv(*g) {}

 // all evaluation are rewritten with in terms of the default evaluator.
 template <typename G, typename Nu, typename Nup, typename Om>
 auto operator()(G const &g, Nu const &nu, Nup const &nup, Om const &omega) const {
  if (omega < 0)
   return conj(gv(-nu, -nup, -omega));
  else
   return gv(nu, nup, omega);
 }
};

using gf3_s = gf<prod<imfreq, imfreq, imfreq>, scalar_valued, nothing>;
using gf3_s_modif = gf<prod<imfreq, imfreq, imfreq>, scalar_valued, nothing, my_evaluator>;

TEST(Gf, EvaluatorWithSymmetry) {
 int nw = 5;
 double beta = 10;
 clef::placeholder<0> k_;
 clef::placeholder<1> q_;
 clef::placeholder<2> r_;
 clef::placeholder<3> iw_;
 clef::placeholder<4> inu_;
 clef::placeholder<5> inup_;

 // This would throw
 auto g_classic = gf3_s{{{beta, Fermion, 2 * nw}, {beta, Fermion, 2 * nw}, {beta, Boson, nw, imfreq::option::positive_frequencies_only}}};

 auto g = gf3_s_modif{{{beta, Fermion, 2 * nw}, {beta, Fermion, 2 * nw}, {beta, Boson, nw, imfreq::option::positive_frequencies_only}}};

 g(inu_, inup_, iw_) << inu_ + 10 * inup_ + 100 * iw_;

 auto nu = matsubara_freq{2, beta, Fermion};
 auto nup = matsubara_freq{2, beta, Fermion};
 auto omega = matsubara_freq{3, beta, Boson};

 EXPECT_CLOSE(g(nu, nup, omega), nu + 10 * nup + 100 * omega);
 EXPECT_CLOSE(g(-nu, -nup, omega), - nu - 10 * nup + 100 * omega);
 EXPECT_CLOSE(g(nu, nup, -omega), nu + 10 * nup - 100 * omega);

 EXPECT_THROW(g_classic(nu, nup, -omega), std::exception);
}
*/
MAKE_MAIN;
