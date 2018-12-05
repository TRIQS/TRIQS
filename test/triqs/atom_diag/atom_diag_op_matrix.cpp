// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

#include <triqs/test_tools/gfs.hpp>

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp> // gf_struct_t
using gf_struct_t = triqs::hilbert_space::gf_struct_t;

using namespace triqs::arrays;
using namespace triqs::hilbert_space;
using namespace triqs::atom_diag;
using namespace triqs::operators;

// -----------------------------------------------------------------------------

using h_scalar_t = double;

using linindex_t = std::map<std::pair<int, int>, int>;

using block_index_map_t = std::map<std::variant<int, std::string>, int>;

// -----------------------------------------------------------------------------
linindex_t make_linear_index(const gf_struct_t &gf_struct, const fundamental_operator_set &fops) {
  linindex_t linindex;
  int block_index = 0;
  for (auto const &bl : gf_struct) {
    int inner_index = 0;
    for (auto const &a : bl.second) {
      linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first, a}];
      inner_index++;
    }
    block_index++;
  }
  return linindex;
}

// -----------------------------------------------------------------------------
TEST(atom_diag, op_matrix) {

  gf_struct_t gf_struct{{"up", {0}}, {"dn", {0}}};
  fundamental_operator_set fops(gf_struct);
  auto linindex = make_linear_index(gf_struct, fops);

  block_index_map_t bidx_map;
  int bidx = 0;
  for (auto const &bl : gf_struct) {
    bidx_map[bl.first] = bidx;
    bidx++;
  }

  // -----------------------------------------------------------------------------
  // atom_diag

  double U  = 1.0;
  double mu = 0.5 * U;

  many_body_operator_real H;
  H += -mu * (n("up", 0) + n("dn", 0)) + U * n("up", 0) * n("dn", 0);

  auto ad = triqs::atom_diag::atom_diag<false>(H, fops);
  std::cout << "Found " << ad.n_subspaces() << " subspaces." << std::endl;

  EXPECT_EQ(ad.n_subspaces(), 4);

  // -----------------------------------------------------------------------------

  {
    auto op = ad.get_op_mat(n("up", 0));

    std::cout << "op =\n" << op << "\n";

    EXPECT_EQ(op.n_blocks(), 4);
    for (auto b : range(op.n_blocks())) {
      if (op.connection(b) != -1) EXPECT_EQ(op.block_mat[b], matrix<double>({{1.}}));
    }
  }

  // -----------------------------------------------------------------------------

  {
    auto op = ad.get_op_mat(n("dn", 0));

    std::cout << "op =\n" << op << "\n";

    EXPECT_EQ(op.n_blocks(), 4);
    for (auto b : range(op.n_blocks())) {
      if (op.connection(b) != -1) EXPECT_EQ(op.block_mat[b], matrix<double>({{1.}}));
    }
  }
}

MAKE_MAIN;
