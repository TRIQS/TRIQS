from pytriqs.gf import *

g = GfImFreq(indices = [range(3), range(5)], beta = 40, n_points = 1000)
gt = g.transpose()

assert g.data.shape == (2000, 3, 5)
assert g.tail.data.shape[1:] == (3, 5)
assert gt.data.shape == (2000, 5, 3)
assert gt.tail.data.shape[1:] == (5, 3)
