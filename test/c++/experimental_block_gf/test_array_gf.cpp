#include <triqs/gfs.hpp>
#include <nda/nda.hpp>

int main() {
    using namespace triqs::gfs;
    using namespace triqs::mesh;
    using namespace nda::clef::literals;

    // imfreq mesh
    auto mesh = imfreq{10, Fermion, 5};

    // gf
    auto g = gf<imfreq, scalar_real_valued>{mesh};
    g(iw_) << 1.0;

    // array of gfs
    auto g_arr = nda::array<gf<imfreq, scalar_real_valued>, 2>{{2, 2}};
    g_arr = g;

    // sum of gf arrays
    auto g_sum = nda::array<gf<imfreq, scalar_real_valued>, 2>{g_arr + g_arr};
    std::cout << g_sum(0, 1).data() << std::endl;
}