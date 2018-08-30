/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 * Copyright (C) 2018, The Simons Foundation
 * Author: H. U.R. Strand
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "brillouin_zone.hpp"
#include <triqs/utility/complex_ops.hpp>

namespace triqs {
  namespace lattice {

    /**
   For tightbinding Hamiltonian with fully localised orbitals
   Overlap between orbital is taken as unit matrix.
  */
    class tight_binding {

      using displs_t = std::vector<std::vector<long>>;
      using matrices_t = std::vector<matrix<dcomplex>>;
      
      bravais_lattice bl_;
      displs_t all_disp;
      matrices_t all_matrices;

      public:
      ///
      tight_binding(bravais_lattice const &bl, std::vector<std::vector<long>> all_disp, std::vector<matrix<dcomplex>> all_matrices);

      /// Underlying lattice
      bravais_lattice const &lattice() const { return bl_; }

      /// Number of bands, i.e. size of the matrix t(k)
      int n_bands() const { return bl_.n_orbitals(); }

      // calls F(R, t(R)) for all R
      template <typename F> friend void foreach (tight_binding const &tb, F f) {
        int n = tb.all_disp.size();
        for (int i = 0; i < n; ++i) f(tb.all_disp[i], tb.all_matrices[i]);
      }

      const displs_t & displacements() const { return all_disp; }
      const matrices_t & matrices() const { return all_matrices; }

      // a simple function on the domain brillouin_zone
      struct fourier_impl {
        tight_binding const &tb;
        const int nb;

        using domain_t = brillouin_zone;
        brillouin_zone domain() const { return {tb.lattice()}; }

        TRIQS_CLEF_IMPLEMENT_LAZY_CALL();

        template <typename K> std14::enable_if_t<!clef::is_clef_expression<K>::value, matrix<dcomplex>> operator()(K const &k) const {
          matrix<dcomplex> res(nb, nb);
          res() = 0;

	  /*

	  // Previous Implementation
          foreach (tb, [&](std::vector<long> const &displ, matrix<dcomplex> const &m) {
            double dot_prod = 0;
            int imax        = displ.size();
            for (int i = 0; i < imax; ++i) dot_prod += k(i) * displ[i];
            //double dot_prod = k(0) * displ[0] + k(1) * displ[1];
            res += m * exp(2_j * M_PI * dot_prod);
          })
            ;
	  */

#pragma omp declare reduction (+ : triqs::arrays::matrix<dcomplex> : omp_out += omp_in ) \
  initializer(omp_priv = triqs::arrays::matrix<dcomplex>(omp_orig.shape()))

#pragma omp parallel for reduction(+:res)
	  for(int idx = 0; idx < tb.all_disp.size(); idx++) {

	    auto const & displ = tb.displacements()[idx];
	    auto const & m = tb.matrices()[idx];
	    
            double dot_prod = 0;
            int imax        = displ.size();
            for (int i = 0; i < imax; ++i) dot_prod += k(i) * displ[i];
	    matrix<dcomplex> tmp = m * exp(2_j * M_PI * dot_prod);

            res += tmp;
	  } 

          return res;
        }
      };

      fourier_impl friend fourier(tight_binding const &tb) { return {tb, tb.n_bands()}; }

    }; // tight_binding

    /**
   Factorized version of hopping (for speed)
   k_in[:,n] is the nth vector
   In the result, R[:,:,n] is the corresponding hopping t(k)
   */
    array<dcomplex, 3> hopping_stack(tight_binding const &TB, arrays::array_const_view<double, 2> k_stack);
    // not optimal ordering here

    std::pair<array<double, 1>, array<double, 2>> dos(tight_binding const &TB, int nkpts, int neps);
    std::pair<array<double, 1>, array<double, 1>> dos_patch(tight_binding const &TB, const array<double, 2> &triangles, int neps, int ndiv);
    array<double, 2> energies_on_bz_path(tight_binding const &TB, k_t const &K1, k_t const &K2, int n_pts);
    array<dcomplex, 3> energy_matrix_on_bz_path(tight_binding const &TB, k_t const &K1, k_t const &K2, int n_pts);
    array<double, 2> energies_on_bz_grid(tight_binding const &TB, int n_pts);
  } // namespace lattice
} // namespace triqs
