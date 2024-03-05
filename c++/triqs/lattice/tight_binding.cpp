// Copyright (c) 2013-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2013-2018 Centre national de la recherche scientifique (CNRS)
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
// Authors: Thomas Ayral, Michel Ferrero, Alexander Hampel, Olivier Parcollet, Nils Wentzell

#include "tight_binding.hpp"
#include <nda/algorithms.hpp>
#include <nda/linalg/eigenelements.hpp>
#include "grid_generator.hpp"
namespace triqs {
  namespace lattice {

    using namespace arrays;

    tight_binding::tight_binding(bravais_lattice bl, std::vector<nda::vector<long>> displ_vec, std::vector<nda::matrix<dcomplex>> overlap_mat_vec)
       : bl_(std::move(bl)), displ_vec_(std::move(displ_vec)), overlap_mat_vec_(std::move(overlap_mat_vec)) {

      // checking inputs
      if (displ_vec_.size() != overlap_mat_vec_.size()) TRIQS_RUNTIME_ERROR << " Number of displacements != Number of matrices";
      for (int i = 0; i < displ_vec_.size(); ++i) {
        if (displ_vec_[i].size() != bl_.ndim())
          TRIQS_RUNTIME_ERROR << "displacement of incorrect size : got " << displ_vec_[i].size() << "instead of " << bl_.ndim();
        if (first_dim(overlap_mat_vec_[i]) != n_orbitals())
          TRIQS_RUNTIME_ERROR << "the first dim matrix is of size " << first_dim(overlap_mat_vec_[i]) << " instead of " << n_orbitals();
        if (second_dim(overlap_mat_vec_[i]) != n_orbitals())
          TRIQS_RUNTIME_ERROR << "the second dim matrix is of size " << second_dim(overlap_mat_vec_[i]) << " instead of " << n_orbitals();

        // check hermiticity of hoppings: Hij(+R)= Hji(-R)* by looping of all displacements again
        bool found = false;
        for (int j = 0; j < displ_vec_.size(); ++j) {
          if (displ_vec_[i] == -displ_vec_[j]) {
            found = true;
            if (max_element(abs(overlap_mat_vec_[i] - dagger(overlap_mat_vec_[j]))) > 1.e-12)
              TRIQS_RUNTIME_ERROR << "hopping matrix of displacement " << displ_vec_[i] << overlap_mat_vec_[i]
                                  << "\nis not hermitian conjugate of matrix for displacement " << displ_vec_[j] << overlap_mat_vec_[j] << "\n";
            break;
          }
        }
        if (not found) TRIQS_RUNTIME_ERROR << "opposite hopping vector of " << displ_vec_[i] << " cannot be found";
      }
    }

    //------------------------------------------------------

    std::pair<array<double, 1>, array<double, 2>> dos(tight_binding const &TB, int nkpts, int neps) {

      // loop on the BZ
      int ndim = TB.lattice().ndim();
      int norb = TB.lattice().n_orbitals();
      grid_generator grid(ndim, nkpts);
      array<double, 1> tempeval(norb);
      array<dcomplex, 3> evec(norb, norb, grid.size());
      array<double, 2> eval(norb, grid.size());
      if (norb == 1)
        for (; grid; ++grid) {
          double ee                = real(TB.fourier((*grid)(range(ndim)))(0, 0));
          eval(0, grid.index())    = ee;
          evec(0, 0, grid.index()) = 1;
        }
      else
        for (; grid; ++grid) {
          // cerr<<" index = "<<grid.index()<<endl;
          array_view<double, 1> eval_sl   = eval(range::all, grid.index());
          array_view<dcomplex, 2> evec_sl = evec(range::all, range::all, grid.index());
          std::tie(eval_sl, evec_sl)      = linalg::eigenelements(TB.fourier((*grid)(range(ndim))));
          // cerr<< " point "<< *grid <<  " value "<< eval_sl<< endl; //" "<< (*grid) (range(ndim)) << endl;
        }

      // define the epsilon mesh, etc.
      array<double, 1> epsilon(neps);
      double epsmax = max_element(eval);
      double epsmin = min_element(eval);
      double deps   = (epsmax - epsmin) / neps;
      // for (int i =0; i< neps; ++i) epsilon(i)= epsmin+i/(neps-1.0)*(epsmax-epsmin);
      for (int i = 0; i < neps; ++i) epsilon(i) = epsmin + (i + 0.5) * deps;

      // bin the eigenvalues according to their energy
      // NOTE: a is defined as an integer. it is the index for the DOS.
      // REPORT <<"Starting Binning ...."<<endl;
      array<double, 2> rho(neps, norb);
      rho() = 0;
      for (int l = 0; l < norb; l++) {
        for (int j = 0; j < grid.size(); j++) {
          int a = int((eval(l, j) - epsmin) / deps);
          if (a == int(neps)) a = a - 1;
          for (int k = 0; k < norb; k++) { rho(a, k) += real(conj(evec(k, l, j)) * evec(k, l, j)); }
        }
      }
      rho /= grid.size() * deps;
      return std::make_pair(epsilon, rho);
    }

    //----------------------------------------------------------------------------------

    std::pair<array<double, 1>, array<double, 1>> dos_patch(tight_binding const &TB, const array<double, 2> &triangles, int neps, int ndiv) {
      // WARNING: This version only works for a single band Hamiltonian in 2 dimensions!!!!
      // triangles is an array of points defining the triangles of the patch
      // neps in the number of bins in energy
      // ndiv in the number of divisions used to divide the triangles

      // int ndim=TB.lattice().dim();
      // int norb=TB.lattice().n_orbitals();
      int ntri = triangles.shape(0) / 3;
      array<double, 1> dos(neps);

      // Check consistency
      int ndim = TB.lattice().ndim();
      // int norb=TB.lattice().n_orbitals();
      if (ndim != 2) TRIQS_RUNTIME_ERROR << "dos_patch : dimension 2 only !";
      if (triangles.shape(1) != ndim) TRIQS_RUNTIME_ERROR << "dos_patch : the second dimension of the 'triangle' array in not " << ndim;

      // Every triangle has ndiv*ndiv k points
      int nk        = ntri * ndiv * ndiv;
      int k_index   = 0;
      double epsmax = -100000, epsmin = 100000;
      array<dcomplex, 2> thop(1, 1);
      array<double, 1> energ(nk), weight(nk);

      // a, b, c are the corners of the triangle
      // g the center of gravity taken from a
      array<double, 1> a(ndim), b(ndim), c(ndim), g(ndim), rv(ndim);
      int pt = 0;
      double s, t;

      // loop over the triangles
      for (int tri = 0; tri < ntri; tri++) {
        a = triangles(pt, range::all);
        pt++;
        b = triangles(pt, range::all);
        pt++;
        c = triangles(pt, range::all);
        pt++;
        g = ((a + b + c) / 3.0 - a) / double(ndiv);

        // the area around a k point might be different from one triangle to the other
        // so I use it to weight the sum in the dos
        double area = abs(0.5 * ((b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0))) / (ndiv * ndiv));

        for (int i = 0; i < ndiv; i++) {
          s = i / double(ndiv);
          for (int j = 0; j < ndiv - i; j++) {
            t = j / double(ndiv);
            for (int k = 0; k < 2; k++) {

              rv = a + s * (b - a) + t * (c - a) + (k + 1.0) * g;

              if (k == 0 || j < ndiv - i - 1) {

                energ(k_index) = real(TB.fourier(rv)(0, 0));
                // compute(rv);
                // energ(k_index) = real(tk_for_eval(1,1)); //tk_for_eval is Fortran array
                weight(k_index) = area;
                if (energ(k_index) > epsmax) epsmax = energ(k_index);
                if (energ(k_index) < epsmin) epsmin = energ(k_index);
                k_index++;
              }
            }
          }
        }
      }
      // check consistency
      assert(k_index == nk);

      // define the epsilon mesh, etc.
      array<double, 1> epsilon(neps);
      double deps = (epsmax - epsmin) / neps;
      for (int i = 0; i < neps; ++i) epsilon(i) = epsmin + i / (neps - 1.0) * (epsmax - epsmin);

      // bin the eigenvalues according to their energy
      int ind;
      dos() = 0.0;
      for (int j = 0; j < nk; j++) {
        ind = int((energ(j) - epsmin) / deps);
        if (ind == int(neps)) ind--;
        dos(ind) += weight(j);
      }
      dos /= deps; // Normalize the DOS
      return {std::move(epsilon), std::move(dos)};
    }

    //------------------------------------------------------
    array<dcomplex, 3> hopping_stack(tight_binding const &TB, nda::array_const_view<double, 2> k_stack) {
      array<dcomplex, 3> res(TB.n_orbitals(), TB.n_orbitals(), k_stack.shape(1));
      for (int i = 0; i < k_stack.shape(1); ++i) res(range::all, range::all, i) = TB.fourier(k_stack(range::all, i));
      return res;
    }

    //------------------------------------------------------
    array<double, 2> energies_on_bz_path(tight_binding const &TB, k_t const &K1, k_t const &K2, int n_pts) {
      int norb = TB.lattice().n_orbitals();
      int ndim = TB.lattice().ndim();
      array<double, 2> eval(norb, n_pts);
      k_t dk = (K2 - K1) / double(n_pts), k = K1;
      for (int i = 0; i < n_pts; ++i, k += dk) { eval(range::all, i) = linalg::eigenvalues(TB.fourier(k(range(ndim)))()); }
      return eval;
    }

    //------------------------------------------------------
    array<dcomplex, 3> energy_matrix_on_bz_path(tight_binding const &TB, k_t const &K1, k_t const &K2, int n_pts) {
      int norb = TB.lattice().n_orbitals();
      int ndim = TB.lattice().ndim();
      array<dcomplex, 3> eval(norb, norb, n_pts);
      k_t dk = (K2 - K1) / double(n_pts), k = K1;
      for (int i = 0; i < n_pts; ++i, k += dk) { eval(range::all, range::all, i) = TB.fourier(k(range(ndim)))(); }
      return eval;
    }

    //------------------------------------------------------
    array<double, 2> energies_on_bz_grid(tight_binding const &TB, int n_pts) {

      int norb = TB.lattice().n_orbitals();
      int ndim = TB.lattice().ndim();
      grid_generator grid(ndim, n_pts);
      array<double, 2> eval(norb, grid.size());
      for (; grid; ++grid) { eval(range::all, grid.index()) = linalg::eigenvalues(TB.fourier((*grid)(range(ndim)))()); }
      return eval;
    }

  } // namespace lattice
} // namespace triqs
