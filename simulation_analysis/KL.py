import numpy as np
import sys
import math


class Minimizer:
    def __init__(self, pq, pc, phi, c, bins_c, bins_q,
                 bin_size_c, bin_size_q, temperature_index, temperature):
        self.pq = pq[temperature_index]
        self.pc = pc[temperature_index]
        self.c = c
        self.phi = phi
        self.bins_q = bins_q
        self.bins_c = bins_c
        self.bin_size_q = bin_size_q
        self.bin_size_c = bin_size_c
        self.temperature_index = temperature_index
        self.temperature = temperature

    def Mean(self, distribution, x):
        return np.sum(distribution * x) / np.sum(distribution)

    def g(self, q, a, b, c):
        return 1./(2 * c**2) * (a*self.phi**4 - 2*a*self.phi**2*(q-b))

    def paux(self, a, b, c):
        p = []

        for q in self.c:
            p.append(1./np.sqrt(2*np.pi*c**2) * np.exp(-1./(2*c**2) * (q - b)**2) * self.Mean(self.pq, np.exp(-self.g(q, a, b, c))))
        return np.array(p, dtype=np.float64)

    def KL(self, a, b, c, model):
        if model == 0:
            paux = self.paux(a, b, c)
        elif model == 1:
            paux = self.paux2(a, b, c)
        else:
            sys.exit(-1)
        # paux /= np.sum(paux*self.bin_size_c)
        k = 0.
        for i in range(len(self.c)):
            if paux[i] > 0 and self.pc[i] > 0:
                k += self.pc[i] * np.log(self.pc[i]/paux[i]) * self.bin_size_c
        return k

    def GradientG(self, q, a, b, c):
        da = 1./c**2 * (a * self.phi**4 - self.phi**2 * (q-b))
        db = 1./c**2 * a * self.phi**2
        dc = 1./c**3 * (a**2 * self.phi**4 - 2*a*self.phi**2*(q - b))

        return da, db, dc

    def g2(self, q, a, b, c):
        return (q - a * self.phi**2 - b)**2/(2*c**2*(1-self.phi**2))

    def GradientG2(self, q, a, b, c):
        da = -self.phi**2*(q-a*self.phi**2-b)/(c**2*(1-self.phi**2))
        db = -(q-a*self.phi**2-b)/(c**2*(1-self.phi**2))
        dc = -(q-a*self.phi**2-b)**2/(c**3*(1-self.phi**2))

        return da, db, dc

    def f2(self, a, b, c):
        f = c*np.sqrt(2*np.pi*(1-self.phi**2))
        # print(1./f)
        return f

    def paux2(self, a, b, c):
        p = []

        for q in self.c:
            p.append(self.Mean(self.pq, 1./self.f2(a, b, c) * np.exp(-self.g2(q, a, b, c))))
        p = np.array(p, dtype=np.float64)
        # print(p)
        return p

    def GradientKL2(self, a, b, c):
        ma = []
        mb = []
        mc = []

        for q in self.c:
            grad = self.GradientG2(q, a, b, c)
            #print(grad)
            norm = self.Mean(self.pq, 1./self.f2(a, b, c) * np.exp(-self.g2(q, a, b, c)))
            # print(norm)
            ma.append(self.Mean(self.pq, 1./self.f2(a, b, c) * np.exp(-self.g2(q, a, b, c)) * grad[0])/norm)
            mb.append(self.Mean(self.pq, 1./self.f2(a, b, c) * np.exp(-self.g2(q, a, b, c)) * grad[1])/norm)
            mc.append(self.Mean(self.pq, 1./self.f2(a, b, c) * np.exp(-self.g2(q, a, b, c)) * grad[2])/norm)
        ma = np.array(ma, dtype=np.float64)
        mb = np.array(mb, dtype=np.float64)
        mc = np.array(mc, dtype=np.float64)

        da = self.Mean(self.pc, ma)
        db = self.Mean(self.pc, mb)
        dc = 1./c + self.Mean(self.pc, mc)

        return da, db, dc

    def GradientKL(self, a, b, c):
        ma = []
        mb = []
        mc = []

        for q in self.c:
            norm = self.Mean(self.pq, np.exp(-self.g(q, a, b, c)))
            grad = self.GradientG(q, a, b, c)
            na = self.Mean(self.pq, np.exp(-self.g(q, a, b, c)) * grad[0])
            nb = self.Mean(self.pq, np.exp(-self.g(q, a, b, c)) * grad[1])
            nc = self.Mean(self.pq, np.exp(-self.g(q, a, b, c)) * grad[2])

            ma.append(na/norm)
            mb.append(nb/norm)
            mc.append(nc/norm)
        ma = np.array(ma, dtype=np.float64)
        mb = np.array(mb, dtype=np.float64)
        mc = np.array(mc, dtype=np.float64)

        ma = self.Mean(self.pc, ma)
        mb = self.Mean(self.pc, mb)
        mc = self.Mean(self.pc, mc)

        da = ma
        db = mb - 1./c**2 * self.Mean(self.pc, (self.c - b))
        dc = mc - 1./c**3 * self.Mean(self.pc, (self.c - b)**2) + 1./c

        return da, db, dc

    def GD(self, loss_file, model, learning_rate=0.01, epochs=1000,
           tolerance=1e-3, beta=0.9,
           print_rate=1, parameters_range=2,
           momentum=True):
        a = np.random.uniform(1, parameters_range)
        b = np.random.uniform(1, parameters_range)
        c = np.random.uniform(1, parameters_range)
        if momentum:
            M = np.zeros(3)

        for i in range(epochs):
            k = self.KL(a, b, c, model)

            if model == 0:
                grad = np.array(self.GradientKL(a, b, c), dtype=np.float64)
            elif model == 1:
                grad = np.array(self.GradientKL2(a, b, c), dtype=np.float64)
            else:
                sys.exit(-1)
            h = grad
            loss_file.write(f'{i}\t{a:.4e}\t{b:.4e}\t{c:.4e}\t{k:.4e}\t{grad[0]:.4e}\t{grad[1]:.4e}\t{grad[2]:.4e}\t{self.temperature}\n')
            if i % print_rate == 0:
                if momentum:
                    print(f'[{self.temperature_index+1}/{i}]{k:.3e} [{grad[0]:.2e} {grad[1]:.2e} {grad[2]:.2e}] [{a:.2e} {b:.2e} {c:.2e}][{M[0]:.2e} {M[1]:.2e} {M[2]:.2e}]')
                else:
                    print(f'[{self.temperature_index+1}/{i}]{k:.3e} [{grad[0]:.2e} {grad[1]:.2e} {grad[2]:.2e}] [{a:.2e} {b:.2e} {c:.2e}]')
            if np.max(np.abs(grad)) < tolerance:
                break
            else:
                if momentum:
                    h = beta*M + (1-beta)*h
                    M = h
                a -= learning_rate * h[0]
                b -= learning_rate * h[1]
                c -= learning_rate * h[2]
                if math.isnan(a):
                    sys.exit(-1)
                    break

        if model == 0:
            paux = self.paux(a, b, c)
        elif model == 1:
            paux = self.paux2(a, b, c)
        paux /= np.sum(paux*self.bin_size_c)

        return a, b, c, paux

    def Mapper(self, range_a, range_b, range_c, model, points=1000, file=None):
        a = np.linspace(0.0001, range_a, num=points, dtype=np.float64)
        b = np.linspace(-range_b, range_b, num=points, dtype=np.float64)
        c = np.linspace(-range_c, range_c, num=points, dtype=np.float64)
        tot = len(a)*len(b)*len(c)
        count = 1.
        for ai in a:
            for bi in b:
                for ci in c:
                    if file:
                        print(f'[{count/tot:.2%}]{ai:.4e}\t{bi:.4e}\t{ci:.4e}\t{self.KL(ai, bi, ci, model):.4e}')
                        file.write(f'{ai:.4e}\t{bi:.4e}\t{ci:.4e}\t{self.KL(ai, bi, ci, model):.4e}\n')
                    else:
                        print(f'[{count/tot:.2%}]{ai:.4e}\t{bi:.4e}\t{ci:.4e}\t{self.KL(ai, bi, ci,model):.4e}')

                    count += 1.

        file.write("\n\n")
