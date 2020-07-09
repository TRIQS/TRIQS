/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation
 *   author : N. Wentzell
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

#warning The use of the triqs/mpi/base.hpp header deprecated. \
  The mpi functionality has been moved into the mpi/mpi.hpp header \
  in order to make it usable also outside of the triqs project. \
  The namespace was renamed from triqs::mpi to just mpi. \
  Please adjust your application accordingly.

#include <mpi/mpi.hpp>
#include <mpi/vector.hpp>
#include <mpi/string.hpp>

namespace triqs::mpi {
  using namespace ::mpi;
}
