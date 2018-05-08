################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2017 by N. Wentzell, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import unittest
import numpy as np

from pytriqs.gf import *

class TestInitGfMultivar(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        # Call base class constructor
        super(TestInitGfMultivar, self).__init__(*args, **kwargs)

        # Set up shared variables
        self.beta = 10.0
        self.n_iw = 10
        iw_mesh = MeshImFreq(self.beta, 'Fermion', self.n_iw)
        self.prod_mesh = MeshProduct(iw_mesh, iw_mesh)

    def test_init_data(self):
        data=np.zeros( (2*self.n_iw, 2*self.n_iw, 1,2,3,4) )
        G = Gf(mesh=self.prod_mesh, data=data)

        self.assertEqual([i for i in G.mesh.size_of_components()], [20,20])
        self.assertEqual(G.target_shape, (1,2,3,4))

        # Check that construction with incompatible indices throws
        with self.assertRaises(AssertionError):
            Gf(mesh=self.prod_mesh, data=data, indices=['a','b'])

    def test_init_indices(self):
        indices = ['a','b','c','d'] # Assumes matrix for target_shape
        G = Gf(mesh=self.prod_mesh, indices=indices)
        self.assertEqual(G.target_shape, (4,4))

        indices = [['a','b','c'],['d']] # Defines target_shape
        G = Gf(mesh=self.prod_mesh, indices=indices)
        self.assertEqual(G.target_shape, (3,1))

        # Check that construction with incompatible indices throws
        with self.assertRaises(AssertionError):
            Gf(mesh=self.prod_mesh, indices=['a'], target_shape=(1,1,1,1))

if __name__ == '__main__':
    unittest.main()
