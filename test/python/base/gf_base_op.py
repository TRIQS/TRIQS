# Copyright (c) 2013-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
# Copyright (c) 2013-2018 Centre national de la recherche scientifique (CNRS)
# Copyright (c) 2018-2023 Simons Foundation
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
# Authors: Michel Ferrero, Alexander Hampel, Jonathan Karp, Olivier Parcollet, Nils Wentzell

from h5 import *
from triqs.gf import *
from triqs.utility.comparison_tests import *

import numpy as np, copy
from numpy import linalg
from math import pi
import unittest

def max_abs(a) :
    return np.amax(np.abs(a))

def onefermion(tau, eps, beta):
    from math import exp
    return -exp(-eps * tau) / (1 + exp(-beta * eps));

class test_Gf_Base_Op(unittest.TestCase):

    def setUp(self):

        self.precision = 1.e-6
        self.beta = 50.0

    def test_Base_Op(self):

        iw_mesh = MeshImFreq(beta=self.beta, S = "Fermion", n_iw = 1000)

        ga = Gf(mesh=iw_mesh, target_shape=(2,2), name = "a1Block")
        gb = Gf(mesh=iw_mesh, target_shape=(2,2), name = "b1Block")

        G = BlockGf(name_list = ('a','b'), block_list = (ga,gb), make_copies = False)
        G << iOmega_n + 2.0

        def matsu(n) :
            return (2*n+1)*pi/self.beta * 1j

        for ii, g in G :
            N = g.data.shape[0]
            for n in range(N//2) :
                assert_array_close_to_scalar( g.data[n+N//2],  matsu(n) + 2.0)

        # Arithmetic operations
        G2 = G.copy()
        G2 << G * G + 1.5 * G

        G2 += G
        G2['b'] += np.eye(G['a'].data.shape[0])
        G2['a'] -= 1.3-1.0j

        # inverse:
        G << inverse(G)

        #  Density:
        dens = G.total_density()
        assert abs(dens - 4.000001283004012) < self.precision, "oops dens =  %s"%dens

        # FT:
        f = lambda g,L : Gf(mesh = MeshImTime(beta=self.beta, S="Fermion", n_tau = L), target_shape = g.target_shape)
        gt = BlockGf(name_block_generator = [ (n,f(g,2001) ) for n,g in G], make_copies=False, name='gt')
        for (i,gtt) in gt : gtt.set_from_fourier(G[i])

        res = np.array([[[  3.14815470e-04,   0.00000000e+00],
                [  0.00000000e+00,   3.14815470e-04]],

               [[ -1.48721028e-04,   0.00000000e+00],
                [  0.00000000e+00,  -1.48721028e-04]],

               [[  7.46732524e-05,   0.00000000e+00],
                [  0.00000000e+00,   7.46732524e-05]]])

        assert_arrays_are_close(gt['a'].data[:3], res, 1.e-3)

        # Matrix operations:
        ga2 = Gf(indices = [1,2,3], mesh=MeshImFreq(beta=self.beta, S="Fermion", n_iw=1000), name = "a1Block")
        mat = np.array([[1.0,0.0,1.0],[-1.0,1.0,0.0]], complex)

        ga2.from_L_G_R(mat.transpose(),ga,mat)

        res = np.array([[[ 0.99901401-0.03138495j, -0.49950701+0.01569248j,
                  0.49950701-0.01569248j],
                [-0.49950701+0.01569248j,  0.49950701-0.01569248j,  0.00000000+0.j        ],
                [ 0.49950701-0.01569248j,  0.00000000+0.j        ,
                  0.49950701-0.01569248j]],

               [[ 0.99119556-0.09341798j, -0.49559778+0.04670899j,
                  0.49559778-0.04670899j],
                [-0.49559778+0.04670899j,  0.49559778-0.04670899j,  0.00000000+0.j        ],
                [ 0.49559778-0.04670899j,  0.00000000+0.j        ,
                  0.49559778-0.04670899j]],

               [[ 0.97592014-0.15329718j, -0.48796007+0.07664859j,
                  0.48796007-0.07664859j],
                [-0.48796007+0.07664859j,  0.48796007-0.07664859j,  0.00000000+0.j        ],
                [ 0.48796007-0.07664859j,  0.00000000+0.j        ,
                  0.48796007-0.07664859j]]])
        shift = ga2.data.shape[0]//2
        assert_arrays_are_close(ga2.data[shift:shift+3], res)

        # conjugate:
        Gc = G.conjugate()

        res = np.array([[[ 0.49950701+0.01569248j,  0.00000000-0.j        ],
                [-0.00000000-0.j        ,  0.49950701+0.01569248j]],

               [[ 0.49559778+0.04670899j,  0.00000000-0.j        ],
                [-0.00000000-0.j        ,  0.49559778+0.04670899j]],

               [[ 0.48796007+0.07664859j,  0.00000000-0.j        ],
                [-0.00000000-0.j        ,  0.48796007+0.07664859j]]])
        shift = Gc['a'].data.shape[0]//2
        assert_arrays_are_close(Gc['a'].data[shift:shift+3], res)

        # to be continued
        # tranpose

        with HDFArchive('gf_base_op_test.h5','w') as h :
            h['g'] = G
            h['gt'] = gt

        # Check reading out of archive
        with HDFArchive('gf_base_op_test.h5','r') as h :
            assert_block_gfs_are_close(h['g'], G)
            assert_block_gfs_are_close(h['gt'], gt)

        # Pickle (also use in mpi, etc...)
        import pickle

        def check_pickle(g):
            s = pickle.dumps(g)
            g2 = pickle.loads(s)
            assert_block_gfs_are_close(g,g2)

        check_pickle(G)
        check_pickle(gt)

        # some basic checks for MeshImTime
        tau_mesh = MeshImTime(beta=self.beta, S = "Fermion", n_tau = 1000)

        ga_tau = Gf(mesh=tau_mesh, target_shape=(2,2), name = "a1Block")
        gb_tau = Gf(mesh=tau_mesh, target_shape=(2,2), name = "b1Block")

        G_tau = BlockGf(name_list = ('a','b'), block_list = (ga_tau,gb_tau), make_copies = False)
        for block, g in G_tau:
            for tau in g.mesh:
                g[tau] = onefermion(tau, 1.413, 1e-10)

        # Arithmetic operations
        G2_tau = G_tau.copy()
        G2_tau << 1.5 * G_tau

        G2_tau += G_tau
        G2_tau['b'] += np.eye(G_tau['b'].target_shape[0])
        G2_tau['a'] -= 1.3-1.0j

    def test_Mat_Prod(self):

        iw_mesh = MeshImFreq(beta=self.beta, S = "Fermion", n_iw = 50)

        Mat = np.array([[1, 2], [3, 4]])

        G = Gf(mesh=iw_mesh, target_shape=(2,2), name = "G_iw")
        G << iOmega_n * Mat

        G_exact = G.copy()
        G_exact[0, 0] << Mat[0, 0] * iOmega_n
        G_exact[1, 0] << Mat[1, 0] * iOmega_n
        G_exact[0, 1] << Mat[0, 1] * iOmega_n
        G_exact[1, 1] << Mat[1, 1] * iOmega_n

        assert_gfs_are_close(G, G_exact)
        assert_gfs_are_close(Mat * G * linalg.inv(Mat), G_exact)

        # ======

        G = Gf(mesh=MeshProduct(iw_mesh, iw_mesh), target_shape=(2,2), name = "G_iw_iw")
        for iw1, iw2 in G.mesh:
            G[iw1, iw2] = np.identity(2) / (iw1 + 2.0 * iw2 + 4.0)

        G_exact = G.copy()
        for iw1, iw2 in G.mesh:
            G_exact[iw1, iw2][0, 0] = Mat[0, 0] / (iw1 + 2.0 * iw2 + 4.0)
            G_exact[iw1, iw2][1, 0] = Mat[1, 0] / (iw1 + 2.0 * iw2 + 4.0)
            G_exact[iw1, iw2][0, 1] = Mat[0, 1] / (iw1 + 2.0 * iw2 + 4.0)
            G_exact[iw1, iw2][1, 1] = Mat[1, 1] / (iw1 + 2.0 * iw2 + 4.0)

        assert_gfs_are_close(G * Mat, G_exact)
        assert_gfs_are_close(Mat * G, G_exact)
        G *= Mat
        assert_gfs_are_close(G, G_exact)
        assert_gfs_are_close(Mat * G * linalg.inv(Mat), G_exact)

    def test_different_rank_prod(self):
        mesh = MeshImFreq(beta=self.beta, S="Fermion", n_iw=10)
        G1 = Gf(mesh=mesh, target_shape=[2,2])
        G1 << inverse(iOmega_n + 2)
        G2 = Gf(mesh=mesh, target_shape=[])
        G2 << inverse(iOmega_n - 2)

        G_loop = G1.copy()
        for i in range(len(G1.mesh)):
            G_loop.data[i] = G1.data[i] * G2.data[i]

        G3 = G1.copy()
        G3 *= G2

        assert_gfs_are_close(G1*G2, G_loop)
        assert_gfs_are_close(G2*G1, G_loop)
        assert_gfs_are_close(G3, G_loop)

if __name__ == '__main__':
    unittest.main()
