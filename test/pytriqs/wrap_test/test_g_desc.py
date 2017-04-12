from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.wrap_test.test_g", doc = " Doc of my_module ")
module.add_include("<triqs/../test/pytriqs/wrap_test/g.hpp>") 
module.add_include("<triqs/python_tools/converters/block_gf.hpp>")
#module.add_include("<triqs/python_tools/converters/gf.hpp>")

module.add_function (name = "make_bgf", signature = "block_gf_view<imfreq> (double a)", doc = "DOC of print_a")
module.add_function (name = "pass_bgf", signature = "void (block_gf_view<imfreq> g) ", doc = "DOC of print_a")

module.add_function (name = "make_sgf", signature = "gf_view<imfreq,scalar_valued> (double a)", doc = "DOC of print_a")
module.add_function (name = "pass_sgf", signature = "void (gf_view<imfreq,scalar_valued> g)", doc = "DOC of print_a")


if __name__ == '__main__' : 
   module.generate_code()

