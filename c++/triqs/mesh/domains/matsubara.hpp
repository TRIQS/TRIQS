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
// Authors: Thomas Ayral, Philipp Dumitrescu, Laura Messio, Olivier Parcollet, Nils Wentzell

#pragma once

#include <numbers>
#include <fmt/core.h>
#include <h5/h5.hpp>

#include "../utils.hpp"
#include "../../utility/exceptions.hpp"
#include "../../utility/kronecker.hpp"

namespace triqs::mesh {

  /**
  * A matsubara frequency, i.e.
  *   * n : int, the index
  *   * beta : double, the temperature inverse
  *   * statistic : Fermion or Boson
  *
  * * Can be cast into a complex.
  *
  * * Every operation is done by casting to complex, except addition and substraction of matsubara_freq, which return matsubara_freq
  *   and work on the index
  **/
  struct matsubara_freq {
    long n                   = 0; /// The matsubara index
    double beta              = 0.0;
    statistic_enum statistic = Fermion;

    matsubara_freq() = default;
    matsubara_freq(long n_, double beta_, statistic_enum stat_) : n(n_), beta(beta_), statistic(stat_) {}

    /// Matsubara index accessor (consistency with imfreq::mesh_point_t interface)
    [[nodiscard]] long index() const { return n; }

    using cast_t = std::complex<double>;
    operator cast_t() const { return std::complex<double>{0, std::numbers::pi * (2 * n + statistic) / beta}; }
  };

  inline std::ostream &operator<<(std::ostream &out, matsubara_freq const &y) { return out << std::complex<double>(y); }

  inline matsubara_freq operator+(matsubara_freq const &x, matsubara_freq const &y) {
    return {x.n + y.n + (x.statistic & y.statistic), x.beta, ((x.statistic ^ y.statistic) == 1 ? Fermion : Boson)};
  }

  inline matsubara_freq operator-(matsubara_freq const &x, matsubara_freq const &y) {
    return {x.n - y.n - (~x.statistic & y.statistic), x.beta, ((x.statistic ^ y.statistic) == 1 ? Fermion : Boson)};
  }

  inline matsubara_freq operator-(matsubara_freq const &mp) { return {-(mp.n + (mp.statistic == Fermion ? 1 : 0)), mp.beta, mp.statistic}; }

  inline std::complex<double> operator/(matsubara_freq const &x, matsubara_freq const &y) {
    return std::complex<double>{x} * std::complex<double>{y};
  }

  inline std::complex<double> operator*(matsubara_freq const &x, matsubara_freq const &y) {
    return std::complex<double>{x} * std::complex<double>{y};
  }

#define IMPL_OP(OP)                                                                                                                                  \
  template <typename T>                                                                                                                              \
    requires(not std::is_base_of_v<matsubara_freq, std::decay_t<T>>)                                                                                 \
  inline auto operator OP(matsubara_freq const &iw, T &&y) {                                                                                         \
    return dcomplex(iw) OP std::forward<T>(y);                                                                                                       \
  }                                                                                                                                                  \
  template <typename T>                                                                                                                              \
    requires(not std::is_base_of_v<matsubara_freq, std::decay_t<T>>)                                                                                 \
  inline auto operator OP(T &&x, matsubara_freq const &iw) {                                                                                         \
    return std::forward<T>(x) OP dcomplex(iw);                                                                                                       \
  }
  IMPL_OP(+);
  IMPL_OP(-);
  IMPL_OP(*);
  IMPL_OP(/);
#undef IMPL_OP

  inline bool kronecker(matsubara_freq const &freq) { return freq.n == 0; }
  inline bool kronecker(matsubara_freq const &f1, matsubara_freq const &f2) { return f1.n == f2.n; }

  //---------------------------------------------------------------------------------------------------------
  // Domains

  struct matsubara_time_domain; // Forward Declaration

  struct [[deprecated("matsubara_freq_domain is deprecated")]] matsubara_freq_domain {
    using value_t = matsubara_freq;

    double beta              = 0.0;
    statistic_enum statistic = Fermion;

    [[nodiscard]] bool contains(value_t const &pt) const { return (pt.beta == beta) && (pt.statistic == statistic); }

    matsubara_freq_domain() = default;
    matsubara_freq_domain(double beta_, statistic_enum statistic_) : beta{beta_}, statistic(statistic_) {
      if (beta < 0) TRIQS_RUNTIME_ERROR << "Matsubara domain construction : beta < 0 : beta =" << beta << "\n";
    }
    explicit matsubara_freq_domain(matsubara_time_domain const &x);

    bool operator==(matsubara_freq_domain const &) const = default;
    bool operator!=(matsubara_freq_domain const &) const = default;

    [[nodiscard]] static std::string hdf5_format() { return "MatsubaraFreqDomain"; }

    /// Write into HDF5
    friend void h5_write(h5::group fg, std::string const &subgroup_name, matsubara_freq_domain const &d) {
      h5::group gr = fg.create_group(subgroup_name);
      h5_write(gr, "beta", d.beta);
      h5_write(gr, "statistic", (d.statistic == Fermion) ? "F" : "B");
    }

    /// Read from HDF5
    friend void h5_read(h5::group fg, std::string const &subgroup_name, matsubara_freq_domain &d) {
      h5::group gr      = fg.open_group(subgroup_name);
      std::string stats = " ";
      h5_read(gr, "beta", d.beta);
      h5_read(gr, "statistic", stats);
      d.statistic = stats == "F" ? Fermion : Boson;
    }

    friend std::ostream &operator<<(std::ostream &sout, matsubara_freq_domain const &d) {
      auto stat_cstr = (d.statistic == Boson ? "Boson" : "Fermion");
      return sout << fmt::format("Matsubara frequency domain with beta = {}, statistic = {}", d.beta, stat_cstr);
    }
  };

  struct [[deprecated("matsubara_time_domain is deprecated")]] matsubara_time_domain {
    using value_t = double;

    double beta              = 0.0;
    statistic_enum statistic = Fermion;

    [[nodiscard]] bool contains(value_t const &pt) const { return (0 <= pt) && (pt <= beta); }

    [[nodiscard]] constexpr value_t min() const { return 0.0; }
    [[nodiscard]] value_t max() const { return beta; }

    matsubara_time_domain() = default;
    matsubara_time_domain(double beta_, statistic_enum statistic_) : beta{beta_}, statistic(statistic_) {
      if (beta < 0) TRIQS_RUNTIME_ERROR << "Matsubara domain construction : beta < 0 : beta =" << beta << "\n";
    }
    explicit matsubara_time_domain(matsubara_freq_domain const &x);

    bool operator==(matsubara_time_domain const &) const = default;
    bool operator!=(matsubara_time_domain const &) const = default;

    static std::string hdf5_format() { return "MatsubaraTimeDomain"; }

    /// Write into HDF5
    friend void h5_write(h5::group fg, std::string const &subgroup_name, matsubara_time_domain const &d) {
      h5::group gr = fg.create_group(subgroup_name);
      h5_write(gr, "beta", d.beta);
      h5_write(gr, "statistic", (d.statistic == Fermion ? "F" : "B"));
    }

    /// Read from HDF5
    friend void h5_read(h5::group fg, std::string const &subgroup_name, matsubara_time_domain &d) {
      h5::group gr = fg.open_group(subgroup_name);
      double b;
      std::string stats{};
      h5_read(gr, "beta", b);
      h5_read(gr, "statistic", stats);
      d = matsubara_time_domain(b, (stats == "F" ? Fermion : Boson));
    }

    friend std::ostream &operator<<(std::ostream &sout, matsubara_time_domain const &d) {
      auto stat_cstr = (d.statistic == Boson ? "Boson" : "Fermion");
      return sout << fmt::format("Matsubara time domain with beta = {}, statistic = {}", d.beta, stat_cstr);
    }
  };

  inline matsubara_freq_domain::matsubara_freq_domain(matsubara_time_domain const &x) : matsubara_freq_domain{x.beta, x.statistic} {};
  inline matsubara_time_domain::matsubara_time_domain(matsubara_freq_domain const &x) : matsubara_time_domain{x.beta, x.statistic} {};

} // namespace triqs::mesh
