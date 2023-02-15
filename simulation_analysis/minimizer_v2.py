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
    def __init__(self, pq, pc, phi, c, bins_c, bins_q,
                 bin_size_c, bin_size_q, temperature_index):
        self.pq = pq[temperature_index]
        self.pc = pc[temperature_index]
        self.c = c
        self.phi = phi
        self.bins_q = bins_q
        self.bins_c = bins_c
        self.bin_size_q = bin_size_q
        self.bin_size_c = bin_size_c

    def Entropy(self, distribution=None):
        e = 0.
        if distribution:
            for k in range(len(distribution)):
                if distribution[k] > 0:
                    e += distribution[k] * np.log(distribution[k])
            return -e
        else:
            for k in range(len(self.pc)):
                if self.pc[k] != 0:
                    e += self.pc[k] * np.log(self.pc[k])
            return -e

    def Mean(self, distribution, x):
        return np.sum(distribution * x) / np.sum(distribution)

    def g(self, q, a, b, c):
        return -1./(2 * c**2) * (a*self.phi**4 - 2*a*phi**2*(q-b))

    def Zeta(self, q, a, b, c):
        z = np.sum(self.pq * self.bin_size_q * np.exp(-1./(2*c**2) *
                   (a**2*self.phi**4-2*a*self.phi**2*(q-b))))
        return z

    def paux(self, a, b, c):
        p = []

        for i in self.c:
            p.append(1./np.sqrt(2*np.pi*c**2)*np.exp(-1./(2*c**2)*(i-b)**2) *
                     self.Zeta(i, a, b, c))
        return np.array(p, dtype=np.float64)

    def p(self, q, a, b, c):
        p = self.bin_size_q*self.pq*np.exp(-1./(2*c**2)*(a*self.phi**4-2*a *
                                                         self.phi**2*(q-b)))
        return p

    def KL(self, a, b, c):
        paux = self.paux(a, b, c)

        k = 0.
        for i in range(len(self.c)):
            if paux[i] > 0 and self.pc[i] > 0:
                k += self.pc[i] * np.log(self.pc[i]/paux[i]) * bin_size_c
        return k

    def Gradient(self, a, b, c):
        ma = []
        mb = []
        mc = []
        for i in self.c:
            ma.append(self.Mean(self.p(i, a, b, c), (a*self.phi**4 - self.phi**2*(i-b))))
            mb.append(self.Mean(self.p(i, a, b, c), self.phi**2))
            mc.append(self.Mean(self.p(i, a, b, c), (a**2*self.phi**4-2*a*self.phi**2*(i-b))))
        ma = np.array(ma, dtype=np.float64)
        mb = np.array(mb, dtype=np.float64)
        mc = np.array(mc, dtype=np.float64)

        da = 1./c**2 * self.Mean(self.pc, ma)
        db = 1./c**2 * (b - self.Mean(self.pc, self.c)) + a/c**2 * self.Mean(self.pc, mb)
        dc = 1./c - 1./c**3 * (self.Mean(self.pc, self.c**2)+b**2-2*b*self.Mean(self.pc, self.c)) + 1./c**3 * (self.Mean(self.pc, mc))

        return da, db, dc


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

npt = 1
for k in range(npt):
    minimizer = Minimizer(pq=pq, pc=pc,
                          phi=phi, c=c,
                          bins_c=bins_c,
                          bins_q=bins_q,
                          bin_size_c=bin_size_c,
                          bin_size_q=bin_size_q,
                          temperature_index=k)

    print(minimizer.KL(0.5, 0.5, 0.5), minimizer.Gradient(0.5, 0.5, 0.5))
    # print(minimizer.Zeta(0.5, 0.5, 0.5))
