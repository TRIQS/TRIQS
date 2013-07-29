/*******************************************************************************
 * 
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_RE_IM_FREQ_H
#define TRIQS_GF_RE_IM_FREQ_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./refreq.hpp"
#include "./imfreq.hpp"
#include "./meshes/product.hpp"

namespace triqs { namespace gfs { 
  
  struct re_im_freq {};
  
  namespace gfs_implementation {
    
    // the mesh
    template<typename Opt> struct mesh<re_im_freq,Opt>  { 
      typedef typename mesh<refreq,Opt>::type m1_t; 
      typedef typename mesh<imfreq,Opt>::type m2_t; 
      typedef mesh_product<m1_t,m2_t> type;
      static type make (double wmin, double wmax, size_t n_freq_re, double beta, statistic_enum S, size_t n_freq_im) { 
        return {gfs::make_gf_mesh<refreq,Opt>(wmin,wmax,n_freq_re,full_bins), make_gf_mesh<imfreq,Opt>(beta, S, n_freq_im)};
      }
    };
    
    // singularity
    template<typename Opt> struct singularity<re_im_freq,scalar_valued,Opt>  { typedef gf<refreq,scalar_valued> type;};
    
    // h5 name
    template<typename Opt> struct h5_name<re_im_freq,scalar_valued,Opt> { static std::string invoke(){ return  "GfReImFreq";}};
    
    /// ---------------------------  data access  ---------------------------------
    
    template<typename Opt> struct data_proxy<re_im_freq,scalar_valued,Opt> : data_proxy_array<std::complex<double>,1> {};
    
    /// ---------------------------  evaluator ---------------------------------  
    
    template<typename Opt>
    struct evaluator<re_im_freq,scalar_valued,Opt> {
      static constexpr int arity = 2;
      template<typename G>
      std::complex<double> operator() (G const * g, double w, long n)  const {
        auto & data = g->data();
        auto & mesh = g->mesh();
        size_t nr; double wr; bool in;
        std::tie(in, nr, wr) = windowing( std::get<0>(g->mesh().components()), w);
        if (!in) TRIQS_RUNTIME_ERROR <<" Evaluation out of bounds";
        auto gg = [g,data,mesh]( size_t nr, size_t n) {return data(mesh.index_to_linear(std::tuple<size_t,size_t>{nr,n}));};
        return wr * gg(nr,n) + (1-wr) * gg(nr+1,n) ;
      } 
    };
    
    // -------------------------------   Factories  --------------------------------------------------
    
    template<typename Opt> struct factories<re_im_freq, scalar_valued,Opt> {
      typedef gf<re_im_freq, scalar_valued,Opt> gf_t;
//       typedef typename mesh<re_im_freq, Opt>::type mesh_t;
      
      static gf_t make_gf(double wmin, double wmax, size_t nw, double beta, statistic_enum S, size_t nwn) { 
        auto m =  make_gf_mesh<re_im_freq,Opt>(wmin, wmax, nw, beta, S, nwn);
        typename gf_t::data_non_view_t A(m.size()); 
        A() =0;
        return gf_t (m, std::move(A), gfs::make_gf<refreq,scalar_valued>(wmin, wmax, nw), nothing() ) ;
      }
    };
    
} // gfs_implementation

}}
#endif

