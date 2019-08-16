
namespace triqs::mesh {

  /*----------------------------------------------------------
  *  closest_point mechanism
  *  This trait will contain the specialisation to get
  *  the closest point ...
  *--------------------------------------------------------*/
  template <typename Mesh, typename Target> struct closest_point;

  // implementation
  template <typename... T> struct closest_pt_wrap;

  template <typename T> struct closest_pt_wrap<T> : tag::mesh_point {
    T value;
    template <typename U> explicit closest_pt_wrap(U &&x) : value(std::forward<U>(x)) {}
  };

  template <typename T1, typename T2, typename... Ts> struct closest_pt_wrap<T1, T2, Ts...> : tag::mesh_point {
    std::tuple<T1, T2, Ts...> value_tuple;
    template <typename... U> explicit closest_pt_wrap(U &&... x) : value_tuple(std::forward<U>(x)...) {}
  };

  /// DOC ?
  template <typename... T> closest_pt_wrap<T...> closest_mesh_pt(T &&... x) { return closest_pt_wrap<T...>{std::forward<T>(x)...}; }

  //-------------------------------------------------------
  // closest mesh point on the grid
  // ------------------------------------------------------

  template <typename Target> struct closest_point<imtime, Target> {
    // index_t is int
    template <typename M, typename T> static int invoke(M const &mesh, closest_pt_wrap<T> const &p) {
      double x = double(p.value) + 0.5 * mesh.delta();
      int n    = std::floor(x / mesh.delta());
      return n;
    }
  };

  //-------------------------------------------------------
  // closest mesh point on the grid
  // ------------------------------------------------------

  struct closest_point_linear_mesh {
    // index_t is int
    template <typename M, typename T> static int invoke(M const &mesh, closest_pt_wrap<T> const &p) {
      double x = double(p.value) - mesh.x_min() + 0.5 * mesh.delta();
      int n    = std::floor(x / mesh.delta());
      return n;
    }
  };

  //-------------------------------------------------------
  // closest mesh point on the grid
  // ------------------------------------------------------

  template <typename Target> struct closest_point<b_zone, Target> {

    template <typename T> static b_zone::index_t invoke(b_zone const &m, closest_pt_wrap<T> const &p) {

      b_zone::index_t result; // array<Int, 3>
      // FIXME : works only on square lattice
      auto id = m.get_interpolation_data(p.value);
      for (int u = 0; u < 3; ++u) {
        int i     = (id.wa[0][u] > id.wa[1][u] ? 0 : 1);
        result[u] = id.ia[i][u];
      }
      return result;
    }
  };
  //-------------------------------------------------------
  // For all mesh represented by a linear grid, the code is the same

  //template <typename Target> struct closest_point<imtime, Target> : closest_point_linear_mesh{};
  template <typename Target> struct closest_point<retime, Target> : closest_point_linear_mesh {};
  template <typename Target> struct closest_point<refreq, Target> : closest_point_linear_mesh {};

  //-------------------------------------------------------
  // closest mesh point on the grid
  // ------------------------------------------------------

  template <typename... Ms, typename Target> struct closest_point<cartesian_product<Ms...>, Target> {
    using index_t = typename cartesian_product<Ms...>::index_t;

    template <typename M, typename... T, size_t... Is> static index_t _impl(M const &m, closest_pt_wrap<T...> const &p, std::index_sequence<Is...>) {
      return index_t(closest_point<Ms, Target>::invoke(std::get<Is>(m), closest_pt_wrap<T>{std::get<Is>(p.value_tuple)})...);
    }

    template <typename M, typename... T> static index_t invoke(M const &m, closest_pt_wrap<T...> const &p) {
      return _impl(m, p, std::index_sequence_for<T...>{});
    }
  };

} // namespace triqs::mesh
