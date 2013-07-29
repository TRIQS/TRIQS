//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gfs/imfreq.hpp> 
#include <triqs/gfs/imtime.hpp> 
#include <triqs/gfs/block.hpp> 

#include <triqs/gfs/local/functions.hpp> 
#include <boost/mpi.hpp>

namespace tql= triqs::clef;
namespace tqa= triqs::arrays;
using tqa::range;
using triqs::arrays::make_shape;
using triqs::gfs::Fermion;
using triqs::gfs::imfreq;
using triqs::gfs::imtime;
using triqs::gfs::make_gf;
using triqs::gfs::gf;
using triqs::gfs::block_index;



#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main(int argc, char* argv[]) {
 boost::mpi::environment env(argc,argv);
 boost::mpi::communicator c;

 double beta =1;
 auto G = make_gf<imfreq> (beta, Fermion, make_shape(2,2));
 triqs::clef::placeholder<0> om_;
 //G(om_) << (om_ - 2.1);
 auto G2 = G;
 std::cout << c.rank() << "\t"<<  G.singularity()<< std::endl;
 boost::mpi::reduce(c, G, G2, std::plus<gf<imfreq>>(), 0);
 std::cout << c.rank() << "\t" << G2.singularity()<< std::endl;

 G(om_) << (om_ - 2.1);

 //std::cout << c.rank() << "\t"<<  G.singularity()<< std::endl;
 //boost::mpi::reduce(c, G, G2, std::plus<gf<imfreq>>(), 0);
// std::cout << c.rank() << "\t" << G2.singularity()<< std::endl;

 
 auto g3 = G2 + G;

 std::cout << c.rank() << "\t" << g3.singularity()<< std::endl;


 //auto Gi = make_gf<imfreq> (beta, Fermion, make_shape(2,2));
 //G(om_) << (om_ - 2.1);
 auto g4 = g3 + G;
 //std::cout << c.rank() << "\t" << Gi.singularity()<< std::endl;
 std::cout << c.rank() << "\t" << g4.singularity()<< std::endl;
 std::cout << c.rank() << "\t" << g3.singularity() + g4.singularity()<< std::endl;
 
}
