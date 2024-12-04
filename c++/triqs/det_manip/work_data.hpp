// Copyright (c) 2013-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2013-2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2023 Simons Foundation
// Copyright (c) 2016 Igor Krivenko
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
// Authors: Michel Ferrero, JaksaVucicevic, Igor Krivenko, Henri Menke, Laura Messio, Olivier Parcollet, Priyanka Seth, Hugo U. R. Strand, Nils Wentzell

#pragma once

#include <nda/nda.hpp>

#include <vector>

namespace triqs::det_manip::detail {

  // ================ Work Data Types =====================

  // For single-row/column operations
  template <typename x_type, typename y_type, typename value_type> struct work_data_type1 {
    x_type x;
    y_type y;
    long i, j, ireal, jreal;
    // MB = A^(-1)*B,
    // MC = C*A^(-1)
    nda::vector<value_type> MB, MC, B, C;
    // ksi = newdet/det
    value_type ksi;
    void resize(long N) {
      MB.resize(N);
      MC.resize(N);
      B.resize(N);
      C.resize(N);
    }
  };

  // For multiple-row/column operations
  template <typename x_type, typename y_type, typename value_type> struct work_data_typek {
    std::vector<x_type> x;
    std::vector<y_type> y;
    std::vector<long> i, j, ireal, jreal;
    // MB = A^(-1)*B,
    // MC = C*A^(-1)
    nda::matrix<value_type> MB, MC, B, C, ksi;
    void resize(long N, long k) {
      if (k < 2) return;
      x.resize(k);
      y.resize(k);
      i.resize(k);
      j.resize(k);
      ireal.resize(k);
      jreal.resize(k);
      MB.resize(N, k);
      MC.resize(k, N);
      B.resize(N, k);
      C.resize(k, N);
      ksi.resize(k, k);
    }
    value_type det_ksi(long k) const {
      if (k == 2) {
        return ksi(0, 0) * ksi(1, 1) - ksi(1, 0) * ksi(0, 1);
      } else if (k == 3) {
        return                                 // Rule of Sarrus
           ksi(0, 0) * ksi(1, 1) * ksi(2, 2) + //
           ksi(0, 1) * ksi(1, 2) * ksi(2, 0) + //
           ksi(0, 2) * ksi(1, 0) * ksi(2, 1) - //
           ksi(2, 0) * ksi(1, 1) * ksi(0, 2) - //
           ksi(2, 1) * ksi(1, 2) * ksi(0, 0) - //
           ksi(2, 2) * ksi(1, 0) * ksi(0, 1);  //
      } else {
        auto Rk = nda::range(k);
        return nda::determinant(ksi(Rk, Rk));
      };
    }
  };

  // For refill operations
  template <typename x_type, typename y_type, typename value_type> struct work_data_type_refill {
    std::vector<x_type> x_values;
    std::vector<y_type> y_values;
    nda::matrix<value_type> M;
    void reserve(long N) {
      x_values.reserve(N);
      y_values.reserve(N);
      M.resize(N, N);
    }
  };

  // Data storage for temporary data used in the det_manip class when inserting a new row and column.
  //
  // - x and y: MatrixBuilder arguments for the new row and column.
  // - i and j: Positions of the new row and column in the original matrix F.
  // - S_inv: Inverse of the new diagonal element of the inverse matrix M, i.e. S^{-1} = D - C M^{(n)} B
  // - B and C: New column and row of the matrix G excluding the diagonal element D.
  // - MB and CM: Products of the current inverse matrix M and the new column B and row C, respectively.
  template <typename X, typename Y, typename T> struct work_data_insert {
    X x;
    Y y;
    long i;
    long j;
    T S_inv;
    nda::vector<T> B;
    nda::vector<T> C;
    nda::vector<T> MB;
    nda::vector<T> CM;

    // Get current capacity of the data storages.
    auto capacity() const { return B.size(); }

    // Reserve memory and resize the data storages if needed.
    void reserve(long cap) {
      if (cap > capacity()) {
        B.resize(cap);
        C.resize(cap);
        MB.resize(cap);
        CM.resize(cap);
      }
    }
  };

  // Data storage for temporary data used in the det_manip class when inserting k new rows and columns.
  //
  // - x and y: MatrixBuilder arguments for the new rows and columns.
  // - i and j: Positions of the new rows and columns in the original matrix F.
  // - S_inv: Inverse of the block matrix S, i.e. S^{-1} = D - C M^{(n)} B (= D for M^{(0)}).
  // - B and C: New columns and rows of the matrix G excluding the block matrix D.
  // - MB and CM: Products of the current inverse matrix M and B and C, respectively.
  // - det: Determinant of the matrix G^{(n+k)} including the new rows and columns.
  // - sign: Sign associated with the permutation matrices P^{(n+k)}_r and P^{(n+k)}_c including the new rows and
  //   columns.
  template <typename X, typename Y, typename T> struct work_data_insert_k {
    std::vector<X> x;
    std::vector<Y> y;
    std::vector<long> i;
    std::vector<long> j;
    nda::matrix<T> S_inv;
    nda::matrix<T> B;
    nda::matrix<T> C;
    nda::matrix<T> MB;
    nda::matrix<T> CM;

    // Get current capacity of the data storages.
    auto capacity() const { return std::make_pair(B.shape()[0], B.shape()[1]); }

    // Reserve memory and resize the data storages if needed.
    void reserve(long cap, long k) {
      auto const [n_cap, k_cap] = capacity();
      if (cap > n_cap || k > k_cap) {
        S_inv.resize(k, k);
        B.resize(cap, k);
        C.resize(k, cap);
        MB.resize(cap, k);
        CM.resize(k, cap);
      }
    }
  };

  // Data storage for temporary data used in the det_manip class when removing 1 row and column.
  //
  // - i and j: Positions of the row and column in the original matrix F.
  // - ip and jp: Positions of the row and column in the matrix G.
  // - S: Diagonal element of \widetilde{M}^{(n)}.
  template <typename T> struct work_data_remove {
    long i;
    long j;
    long ip;
    long jp;
    T S;
  };

  // Data storage for temporary data used in the det_manip class when removing k rows and columns.
  //
  // - i and j: Positions of the rows and columns in the original matrix F.
  // - ip and jp: Positions of the rows and columns in the matrix G.
  // - S: Block matrix of \widetilde{M}^{(n)}.
  template <typename T> struct work_data_remove_k {
    std::vector<long> i;
    std::vector<long> j;
    std::vector<long> ip;
    std::vector<long> jp;
    nda::matrix<T> S;

    // Get current capacity of the data storages.
    auto capacity() const { return S.shape()[0]; }

    // Reserve memory and resize the data storages if needed.
    void reserve(long cap) {
      if (cap > capacity()) {
        ip.resize(cap);
        jp.resize(cap);
        S.resize(cap, cap);
      }
    }
  };

  // Calculate the determinant of the matrix M(nda::range(k), nda::range(k)).
  template <typename T> T determinant(nda::matrix<T> const &M, long k) {
    switch (k) {
      case 0: return 1;
      case 1: return M(0, 0);
      case 2: return M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1);
      case 3:
        return M(0, 0) * M(1, 1) * M(2, 2) + M(0, 1) * M(1, 2) * M(2, 0) + M(0, 2) * M(1, 0) * M(2, 1) - M(2, 0) * M(1, 1) * M(0, 2)
           - M(2, 1) * M(1, 2) * M(0, 0) - M(2, 2) * M(1, 0) * M(0, 1);
      default: return nda::determinant(M(nda::range(k), nda::range(k)));
    }
  }

} // namespace triqs::det_manip::detail
