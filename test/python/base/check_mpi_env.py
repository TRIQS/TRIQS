# Copyright (c) 2013-2014 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
# Copyright (c) 2013-2014 Centre national de la recherche scientifique (CNRS)
# Copyright (c) 2019-2020 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https:#www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: Olivier Parcollet, Nils Wentzell

import os
import sys


from triqs.utility.mpi import check_for_mpi


is_mpi = check_for_mpi()

assert is_mpi  == True, 'MPI environment could not be properly deterimend, check your MPI env'





