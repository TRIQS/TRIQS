/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012-2015 by O. Parcollet
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

namespace triqs {
namespace gfs {

 /*------------------------------------------------------------------------------------------------------
  *                              HDF5
  *-----------------------------------------------------------------------------------------------------*/

 template <typename G, typename Target> constexpr bool gf_has_target() {
  return std::is_same<typename G::target_t, Target>::value;
 }
/*
 /// ---------------------------  

 // Some work that may be necessary before writing (some compression, see imfreq)
 // Default : do nothing
 template <typename M, typename T> struct gf_h5_before_write {
  template <typename G> static G const &invoke(h5::group gr, G const &g) { return g; }
 };

 // Before writing to h5, check if I can save the positive freq only
 template <typename T> struct gf_h5_before_write<imfreq, T> {
  template <typename G> static gf_const_view<imfreq, T> invoke(h5::group gr, G const &g) {
   if (is_gf_real_in_tau(g, 1.e-13)) return positive_freq_view(g);
   return g;
  }
 };
*/
 /// ---------------------------  

 // FIXME : C17 : REMOVE THIS dispatch with a constexpr if
 // Some work that may be necessary after the read (for backward compatibility e.g.)
 // Default : do nothing
 template <typename M, typename T> struct gf_h5_after_read {
  template <typename G> static void invoke(h5::group gr, G &g) {}
 };

 // After reading from h5, is the function is for freq >0, unfold it to the full mesh
 template <typename T> struct gf_h5_after_read<imfreq, T> {
  template <typename G> static void invoke(h5::group gr, G &g) {
   if (g.mesh().positive_only()) g = make_gf_from_real_gf(make_const_view(g));
  }
 };
 // same, for python interface
 template <typename T> gf<imfreq, T>  _gf_h5_after_read(gf_view<imfreq, T> g) { 
   if (g.mesh().positive_only()) return make_gf_from_real_gf(make_const_view(g));
   else return g;
  }

 /// ---------------------------

 // the h5 write and read of gf members, so that we can specialize it e.g. for block gf
 template <typename V, typename T> struct gf_h5_rw {

  //template <typename G> static void write(h5::group gr, G const &g) { _write(gr, gf_h5_before_write<V, T>::invoke(gr, g)); }

  template <typename G> static void write(h5::group gr, G const &g) {
   // write the data
   //constexpr bool _can_compress = (gf_has_target<G, imtime>() or gf_has_target<G, legendre>());
   //if (_can_compress and is_gf_real(g))
   // h5_write(gr, "data", array<double, G::data_t::rank>(real(g.data())));
   //else
   h5_write(gr, "data", g.data());
   h5_write(gr, "singularity", g._singularity);
   h5_write(gr, "mesh", g._mesh);
   h5_write(gr, "indices", g._indices);
  }

  template <typename G> static void read(h5::group gr, G &g) {
   h5_read(gr, "data", g._data);
   h5_read(gr, "singularity", g._singularity);
   h5_read(gr, "mesh", g._mesh);
   h5_read(gr, "indices", g._indices);
   gf_h5_after_read<V, T>::invoke(gr, g);
  }
 };

 /// ---------------------------

 // FIXME The mesh info is duplicated with the function it reads ... Need to pass a ref to the parent function ?
 template <typename Var> struct gf_h5_rw<Var, tail_valued<matrix_valued>> {
  template <typename G> static void write(h5::group gr, G const &g) {
   h5_write(gr, "data", g._data);
   h5_write(gr, "mesh", g._mesh);
  }
  template <typename G> static void read(h5::group gr, G &g) {
   h5_read(gr, "data", g._data);
   h5_read(gr, "mesh", g._mesh);
  }
 };

} // triqs::gfs
}

