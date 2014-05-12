/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by O. Parcollet
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
#include "./tools.hpp"
#include "./gf.hpp"
#include "./meshes/product.hpp"
#include "./evaluators.hpp"

namespace triqs {
namespace gfs {

 template <typename... Ms> struct cartesian_product {
  using type = std::tuple<Ms...>;
  static constexpr size_t size = sizeof...(Ms);
 };

 // use alias
 template <typename... Ms> struct cartesian_product<std::tuple<Ms...>> : cartesian_product<Ms...> {};

 // the mesh is simply a cartesian product
 template <typename Opt, typename... Ms> struct gf_mesh<cartesian_product<Ms...>, Opt> : mesh_product<gf_mesh<Ms, Opt>...> {
  // using mesh_product< gf_mesh<Ms,Opt> ... >::mesh_product< gf_mesh<Ms,Opt> ... > ;
  using B = mesh_product<gf_mesh<Ms, Opt>...>;
  gf_mesh() = default;
  gf_mesh(gf_mesh<Ms, Opt>... ms) : B{std::move(ms)...} {}
 };

 namespace gfs_implementation {

  /// ---------------------------  hdf5 ---------------------------------

  // h5 name : name1_x_name2_.....
  template <typename Opt, typename... Ms> struct h5_name<cartesian_product<Ms...>, matrix_valued, Opt> {
   static std::string invoke() {
    return triqs::tuple::fold([](std::string a, std::string b) { return a + std::string(b.empty() ? "" : "_x_") + b; },
                              std::make_tuple(h5_name<Ms, matrix_valued, Opt>::invoke()...), std::string());
   }
  };
  template <typename Opt, int R, typename... Ms>
  struct h5_name<cartesian_product<Ms...>, tensor_valued<R>, Opt> : h5_name<cartesian_product<Ms...>, matrix_valued, Opt> {};

  // a slight difference with the generic case : reinterpret the data array to avoid flattening the variables
  template <typename Opt, int R, typename... Ms> struct h5_rw<cartesian_product<Ms...>, tensor_valued<R>, Opt> {
   using g_t = gf<cartesian_product<Ms...>, tensor_valued<R>, Opt>;

   static void write(h5::group gr, typename g_t::const_view_type g) {
    h5_write(gr, "data", reinterpret_linear_array(g.mesh(), g().data()));
    h5_write(gr, "singularity", g._singularity);
    h5_write(gr, "mesh", g._mesh);
    h5_write(gr, "symmetry", g._symmetry);
   }

   template <bool IsView>
   static void read(h5::group gr, gf_impl<cartesian_product<Ms...>, tensor_valued<R>, Opt, IsView, false> &g) {
    using G_t = gf_impl<cartesian_product<Ms...>, tensor_valued<R>, Opt, IsView, false>;
    h5_read(gr, "mesh", g._mesh);
    auto arr = arrays::array<typename G_t::data_t::value_type, sizeof...(Ms) + R>{};
    h5_read(gr, "data", arr);
    auto sh = arr.shape();
    arrays::mini_vector<size_t, R + 1> sh2;
    sh2[0] = g._mesh.size();
    for (int u = 1; u < R + 1; ++u) sh2[u] = sh[sizeof...(Ms) - 1 + u];
    g._data = arrays::array<typename G_t::data_t::value_type, R + 1>{sh2, std::move(arr.storage())};
    h5_read(gr, "singularity", g._singularity);
    h5_read(gr, "symmetry", g._symmetry);
   }
  };

  /// ---------------------------  data access  ---------------------------------

  template <typename Opt, typename... Ms>
  struct data_proxy<cartesian_product<Ms...>, scalar_valued, Opt> : data_proxy_array<std::complex<double>, 1> {};
  template <typename Opt, typename... Ms>
  struct data_proxy<cartesian_product<Ms...>, matrix_valued, Opt> : data_proxy_array<std::complex<double>, 3> {};
  template <int R, typename Opt, typename... Ms>
  struct data_proxy<cartesian_product<Ms...>, tensor_valued<R>, Opt> : data_proxy_array<std::complex<double>, R + 1> {};
 
  template <typename Opt, typename M0>
  struct data_proxy<cartesian_product<M0,imtime>, matrix_valued, Opt> : data_proxy_array<double, 3> {};
  
  /// ---------------------------  evaluator ---------------------------------

  /**
   * This the multi-dimensional evaluator.
   * It combine the evaluator of each components, as long as they are  a linear form
   * eval(g, x) = \sum_i w_i g( n_i(x)) , with w some weight and n_i some points on the grid.
   * Mathematically, it is written as (example of evaluating g(x1,x2,x3,x4)).
   * Notation : eval(X)  : g -> g(X)
   * eval(x1,x2,x3,x4) (g) = eval (x1) ( binder ( g, (),  (x2,x3,x4)) )
   * binder( g, (), (x2,x3,x4)) (p1) = eval(x2)(binder (g,(p1),(x3,x4)))
   * binder( g, (p1), (x3,x4))  (p2) = eval(x3)(binder (g,(p1,p2),(x4)))
   * binder( g, (p1,p2), (x4))  (p3) = eval(x4)(binder (g,(p1,p2,p3),()))
   * binder( g, (p1,p2,p3),())  (p4) = g[p1,p2,p3,p4]
   *
   * p_i are points on the grids, x_i points in the domain.
   *
   * Unrolling the formula gives (for 2 variables, with 2 points interpolation)
   * eval(xa,xb) (g) = eval (xa) ( binder ( g, (),  (xb)) ) =
   *    w_1(xa)  binder ( g, (),  (xb))( n_1(xa)) + w_2(xa)  binder ( g, (), (xb))( n_2(xa))
   *  =  w_1(xa) ( eval(xb)(  binder ( g, (n_1(xa) ),  ())))  + 1 <-> 2
   *  =  w_1(xa) ( W_1(xb) *  binder ( g, (n_1(xa) ),  ())(N_1(xb)) + 1<->2 )  + 1 <-> 2
   *  =  w_1(xa) ( W_1(xb) *  g[n_1(xa), N_1(xb)] + 1<->2 )  + 1 <-> 2
   *  =  w_1(xa) ( W_1(xb) *  g[n_1(xa), N_1(xb)] + W_2(xb) *  g[n_1(xa), N_2(xb)] )  + 1 <-> 2
   *  which is the expected formula
   */
  // implementation : G = gf, Tn : tuple of n points, Ev : tuple of evaluators (the evals functions),
  // pos = counter from #args-1 =>0
  // NB : the tuple is build in reverse with respect to the previous comment.
  template <typename G, typename Tn, typename Ev, int pos> struct binder;

  template <int pos, typename G, typename Tn, typename Ev> binder<G, Tn, Ev, pos> make_binder(G const *g, Tn tn, Ev const &ev) {
   return binder<G, Tn, Ev, pos>{g, std::move(tn), ev};
  }

  template <typename G, typename Tn, typename Ev, int pos> struct binder {
   G const *g;
   Tn tn;
   Ev const &evals;
   auto operator()(size_t p) const
       DECL_AND_RETURN(std::get<pos>(evals)(make_binder<pos - 1>(g, triqs::tuple::push_front(tn, p), evals)));
  };

  template <typename G, typename Tn, typename Ev> struct binder<G, Tn, Ev, -1> {
   G const *g;
   Tn tn;
   Ev const &evals;
   auto operator()(size_t p) const DECL_AND_RETURN(triqs::tuple::apply(on_mesh(*g), triqs::tuple::push_front(tn, p)));
  };

  // now the multi d evaluator itself.
  template <typename Target, typename Opt, typename... Ms> struct evaluator<cartesian_product<Ms...>, Target, Opt> {
   static constexpr int arity = sizeof...(Ms);
   mutable std::tuple<evaluator_fnt_on_mesh<Ms>...> evals;

   struct _poly_lambda { // replace by a polymorphic lambda in C++14
    template <typename A, typename B, typename C> void operator()(A &a, B const &b, C const &c) const {
     a = A{b, c};
    }
   };

   template <typename G, typename... Args>
   // std::complex<double> operator() (G const * g, Args && ... args) const {
   auto operator()(G const *g, Args &&... args)
       const -> decltype(std::get<sizeof...(Args) - 1>(evals)(make_binder<sizeof...(Args) - 2>(g, std::make_tuple(), evals)))
       // when do we get C++14 decltype(auto) ...!?
   {
    static constexpr int R = sizeof...(Args);
    // build the evaluators, as a tuple of ( evaluator<Ms> ( mesh_component, args))
    triqs::tuple::call_on_zip(_poly_lambda(), evals, g->mesh().components(), std::make_tuple(args...));
    return std::get<R - 1>(evals)(make_binder<R - 2>(g, std::make_tuple(), evals));
   }
  };

 } // gf_implementation
}
}

