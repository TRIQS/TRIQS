
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
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

#include "MersenneRNG.hpp"

namespace triqs {
  namespace mc_tools {
    namespace RandomGenerators {

      double RandMT::eval() {
#ifdef DEBUG
        double r = (double)(randomMT()) / 0xFFFFFFFFU;
        std::cout << "RANDOM " << r << std::endl;
        return r;
#else
        return ((double)(randomMT()) / 0xFFFFFFFFU);
#endif
      }
    } // namespace RandomGenerators
  }   // namespace mc_tools
} // namespace triqs
