import numpy as np
import os
import sys


def LoadDist(type, size, pt_flag, block=None):
    if pt_flag:
        pt = 'PT'
    else:
        pt = f'block{block}'

    filename = f'N{size}/dis_ave/dis_ave_{type}_dist_{pt}_size{size}.dat'
    p = np.loadtxt(filename, dtype=np.float64, delimiter="\t", usecols=[0, 1])
    return p[:, 0], p[:, 1]


class Minimizer:
    def __init__(self, pq, pc, phi, c, bins_c, bins_q, bin_size_c, bin_size_q):
        self.pq = pq
        self.pc = pc
        self.c = c
        self.bins_q = bins_q
        self.bins_c = bins_c
        self.bin_size_q = bin_size_q
        self.bin_size_c = bin_size_c


size = int(sys.argv[1])
pt = int(sys.argv[2])
npt = int(sys.argv[3])
os.system('clear && pfetch')
phi, pq = LoadDist('parisi', size, pt)
c, pc = LoadDist('theo_ifo', size, pt)
phi = np.reshape(phi, newshape=(npt, -1))
c = np.reshape(c, newshape=(npt, -1))
pq = np.reshape(pq, newshape=(npt, -1))
pc = np.reshape(pc, newshape=(npt, -1))

c = c[0]
phi = phi[0]

bins_q = len(phi)
bins_c = len(c)

bin_size_q = 2./bins_q
bin_size_c = 2./bins_c

for k in range(npt):
    norm_q = np.sum(pq[k]*bin_size_q)
    norm_c = np.sum(pc[k]*bin_size_c)
    pq[k] /= norm_q
    pc[k] /= norm_c

print(bins_q)
print(bins_c)
print(bin_size_q)
print(bin_size_c)
print(np.sum(pq[0]*bin_size_q))
print(np.sum(pc[0]*bin_size_c))


