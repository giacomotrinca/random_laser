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

    def f(self):
        return np.sqrt(1 - self.phi**2)
        #return 1.

    def g3(self, q, a, b, c, vec_flag=False):
        
        if vec_flag:
            g = []
            for qi in self.c:
                g.append(1./(np.sqrt(2.)*c) * (qi - a * self.phi**2 - b))
            return np.array(g, dtype=np.float64)
        else:
            return 1./(np.sqrt(2.)*c) * (self.c[q] - a * self.phi**2 - b)

    def pauxGeneral(self, a, b, c):
        p = []
        for q in range(len(self.pc)):
            p.append(1./(c*np.sqrt(2.))*self.Mean(self.pq, np.exp(-self.g3(q, a, b, c)**2/self.f()**2)/self.f()))
        return np.array(p, dtype=np.float64)

    def NormInf(self, a, b, c):
        paux = self.pauxGeneral(a, b, c)
        L = np.abs(paux - self.pc)
        q = np.argmax(L)
        return q, np.amax(L)

    def GradientG3(self, q, a, b, c):
        da = -1./(c*np.sqrt(2.)) * self.phi**2
        db = -1./(c*np.sqrt(2.))
        dc = -1./c * self.g3(q, a, b, c)[q]

        return da, db, dc

    def GradientPaux(self, q, a, b, c):
        
        grad = self.GradientG3(q, a, b, c)
        paux = self.pauxGeneral(a, b, c)
        da = -np.sqrt(2./np.pi) * 1./c * self.Mean(self.pq, self.g3(q, a, b, c)/self.f() * np.exp(-self.g3(q, a, b, c)**2/self.f()**2)* grad[0])
        db = -np.sqrt(2./np.pi) * 1./c * self.Mean(self.pq, self.g3(q, a, b, c)/self.f() * np.exp(-self.g3(q, a, b, c)**2/self.f()**2)* grad[1])
        dc = -np.sqrt(2./np.pi) * 1./c * self.Mean(self.pq, self.g3(q, a, b, c)/self.f() * np.exp(-self.g3(q, a, b, c)**2/self.f()**2)* grad[2]) - 1./c * paux[q]
     
        return da, db, dc

    def GradientNormInf(self, q, a, b, c):
        # print(q)
        paux = self.pauxGeneral(a, b, c)
        grad = self.GradientPaux(q, a, b, c)

        da = np.sign(paux[q] - self.pc[q]) * np.abs(grad[0])
        db = np.sign(paux[q] - self.pc[q]) * np.abs(grad[1])
        dc = np.sign(paux[q] - self.pc[q]) * np.abs(grad[2])
        return da, db, dc

    def Norm2(self, a, b, c):
        paux = self.pauxGeneral(a, b, c)
        return np.sum((paux - self.pc)**2 * self.bin_size_c), paux
    
    def NumericalGradient(self, a, b, c, model):

        h = 1e-6

        if model==0 or model == 1:
            loss = self.KL(a, b, c, model)
            loss_a = self.KL(a+h, b, c, model)
            loss_b = self.KL(a, b+h, c, model)
            loss_c = self.KL(a, b, c+h, model)
        elif model == 2:
            max_index, loss = self.NormInf(a, b, c)
            max_index_a, loss_a = self.NormInf(a+h, b, c)
            max_index_b, loss_b = self.NormInf(a, b+h, c)
            max_index_c, loss_c = self.NormInf(a, b, c+h)
            del max_index
            del max_index_a
            del max_index_b
            del max_index_c

        elif model == 3:
            loss, paux = self.Norm2(a, b, c)
            loss_a, paux_a = self.Norm2(a+h, b, c)
            loss_b, paux_b = self.Norm2(a, b+h, c)
            loss_c, paux_c = self.Norm2(a, b, c+h)
            del paux
            del paux_a
            del paux_b
            del paux_c
        else:
            sys.exit(-1)
        
        da = (loss_a - loss)/h
        db = (loss_b - loss)/h
        dc = (loss_c - loss)/h

        return da, db, dc

           
    def GradientNorm2(self, paux, a, b, c):
        f = self.f()
        #print(np.shape(f))
        g = self.g3(None, a, b, c, True)
        g2=g**2
        g22 = g2
        da = np.zeros(shape=(len(self.pc), len(self.pq)))
        db = np.zeros(shape=(len(self.pc), len(self.pq)))
        dc = np.zeros(shape=(len(self.pc), len(self.pq)))
        for q in range(len(self.pc)):
            g[q, :] = g[q, :]/f
            g2[q, :]=g2[q, :]/f**2
            g22[q, :] = g22[q,:] / f
            da[q] = g[q, :] * self.phi**2 * np.exp(-g2[q, :]) * (paux[q] - self.pc[q])
            db[q] = g[q, :] * np.exp(-g2[q, :]) * (paux[q] - self.pc[q])
            dc[q] = g22[q, :] * np.exp(-g2[q, :]) * (paux[q] - self.pc[q])


        
        da = 2./(np.sqrt(np.pi) * c**2) * np.sum(da, axis = 0)
        db = 2./(np.sqrt(np.pi) * c**2) * np.sum(db, axis = 0)
        dc = np.sqrt(2./np.pi)/c**2 * np.sum(dc, axis=0)
        da = self.Mean(self.pq, da)
        db = self.Mean(self.pq, db)
        dc = self.Mean(self.pq, dc)  - 1./c * np.sum(paux * (paux - self.pc))
        
        
        return da, db, dc

    def GD(self, loss_file, model, learning_rate=0.01, epochs=1000,
           tolerance=1e-3, beta=0.9,
           print_rate=1, parameters_range=2,
           momentum=True, numerical_gradient=True):
        a = np.random.uniform(1, parameters_range)
        b = np.random.uniform(1, parameters_range)
        c = np.random.uniform(1, parameters_range)
        if momentum:
            M = np.zeros(3)
            #print(M)
        k0 = 0.
        for i in range(epochs):
            if model == 0 or model == 1:
                k = self.KL(a, b, c, model)
                if model == 0:
                    if numerical_gradient:
                        grad = np.array(self.NumericalGradient(a, b, c, model), dtype=np.float64)
                    else:
                        grad = np.array(self.GradientKL(a, b, c), dtype=np.float64)
                elif model == 1:
                    if numerical_gradient:
                        grad = np.array(self.NumericalGradient(a, b, c, model), dtype=np.float64)
                    else:
                        grad = np.array(self.GradientKL2(a, b, c), dtype=np.float64)
            elif model == 2:
                max_index, k = self.NormInf(a, b, c)
                # print(max_index)
                if numerical_gradient:
                    grad = np.array(self.NumericalGradient(a, b, c, model), dtype=np.float64)
                else:
                    grad = np.array(self.GradientNormInf(max_index, a, b, c), dtype=np.float64)
            elif model == 3:
                k, paux = self.Norm2(a, b, c)
                if numerical_gradient:
                    grad = np.array(self.NumericalGradient(a, b, c, model), dtype=np.float64)
                else:
                    grad = np.array(self.GradientNorm2(paux, a, b, c), dtype =np.float64)
            else:
                sys.exit(-1)
            h = grad
            #print(np.shape(h))
            loss_file.write(f'{i}\t{k:.4e}\t{a:.4e}\t{b:.4e}\t{c:.4e}\t{grad[0]:.4e}\t{grad[1]:.4e}\t{grad[2]:.4e}\t{self.temperature}\n')
            if i % print_rate == 0:
                if momentum:
                    print(f'[{self.temperature_index+1}/{i}]{k:.3e} [{grad[0]:.2e} {grad[1]:.2e} {grad[2]:.2e}] [{a:.2e} {b:.2e} {c:.2e}][{M[0]:.2e} {M[1]:.2e} {M[2]:.2e}]')
                else:
                    print(f'[{self.temperature_index+1}/{i}]{k:.3e} [{np.abs(grad[0]):.2e} {np.abs(grad[1]):.2e} {np.abs(grad[2]):.2e}] [{a:.2e} {b:.2e} {c:.2e}]')
            if np.max(np.abs(grad)) < tolerance or np.abs(k - k0) < tolerance:
                break
            else:
                k0 = k
                if momentum:
                    #print(h)
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
        elif model == 2:
            paux = self.pauxGeneral(a, b, c)
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
