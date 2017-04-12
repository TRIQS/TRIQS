/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012-2013 by O. Parcollet
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
#include "./segment.hpp"

namespace triqs {
namespace gfs {

 struct refreq {};

 template <> struct gf_mesh<refreq> : segment_mesh {
  using var_t = refreq;
  template <typename... T> gf_mesh(T &&... x) : segment_mesh(std::forward<T>(x)...) {}
  // using segment_mesh::segment_mesh;
  
  friend std::string get_triqs_hdf5_data_scheme(gf_mesh const &) { return "MeshReFreq";}
 
  friend void h5_write(h5::group fg, std::string const &subgroup_name, gf_mesh const &m) {
   h5_write_impl(fg, subgroup_name, m,"MeshReFreq");
  }
  
  friend void h5_read(h5::group fg, std::string const & subgroup_name, gf_mesh &m) {
   h5_read_impl(fg, subgroup_name, m,"MeshReFreq");
  }   
 };

}
}

