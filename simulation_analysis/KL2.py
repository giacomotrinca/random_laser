import numpy as np
from scipy.optimize import minimize
from scipy.special import kl_div
import sys
import os
import math

class Minimizer:
    def __init__(self, pq, pc, phi, c, bins_c, bins_q,
                 bin_size_c, bin_size_q, temperature_index, 
                 temperature, model, var_flag, edge,attempt,
                 c_threshold=1e-5):
        self.pq = pq
        self.pc = pc
        self.c = c
        self.phi = phi
        self.bins_q = bins_q
        self.bins_c = bins_c
        self.bin_size_q = bin_size_q
        self.bin_size_c = bin_size_c
        self.temperature_index = temperature_index
        self.temperature = temperature
        self.model = model
        self.SpinModel = var_flag
        self.attempt = attempt
        self.c_threshold=c_threshold
        self.edge=edge

        self.pc = self.Normalize(self.pc, self.bin_size_c)
        self.pq = self.Normalize(self.pq, self.bin_size_q)
        #print(self.pc)
    
    '''This is an auxiliary function that depends on the specific model we are considering
    for the sherrington-kirkpatrick model we have sqrt(1-q**2)
    we don't know the f for or random laser model. Initially we can set it to be identically 1'''
    def f(self):
        if self.SpinModel:
            return np.sqrt(self.edge**2 - self.phi**2)
        else:
            return 1.

    '''This is an auxiliary function used for the seek of notation'''
    def g(self, q, a, b, c, vec_flag=False):
        
        if vec_flag:
            g = []
            for qi in self.c:
                g.append(1./(np.sqrt(2.)*c) * (qi - a * self.phi**2 - b))
            return np.array(g, dtype=np.float64)
        else:
            return 1./(np.sqrt(2.)*c) * (self.c[q] - a * self.phi**2 - b)
    
    '''This function computes the mean over the ensable given by distribution array'''
    def Mean(self, distribution, x):
        return np.sum(distribution * x) / np.sum(distribution)
    
    def Normalize(self, distribution, delta):
        sum = 0.
        for i in range(len(distribution)):
            sum += distribution[i] * delta
        
        return distribution/sum
            
    '''This function computes the auxiliary distribution'''
    def paux(self, params):
        a, b, c = params
        p = []
        if c <= self.c_threshold:
            
            for q in self.c:
                #print(q, b)
                k = np.sqrt((np.abs(q - b))/a)
                #print(q, b)
                #print(k)
                if math.isnan(k):
                    print('complex sqrt!')
                    sys.exit(-1)
                index = np.abs(self.phi-k).argmin()
                p.append(self.pq[index]/(a*k))
        else:
            for q in range(len(self.pc)):
                p.append(1./(c*np.sqrt(2.))*self.Mean(self.pq, np.exp(-self.g(q, a, b, c)**2/self.f()**2)/self.f()))
        paux = np.array(p, dtype=np.float64)
        
        return self.Normalize(paux, self.bin_size_c)
    
    def loss(self, params):
        
        paux = self.paux(params)
        # Kullback Leibler KL(IFO||AUX)
        if self.model == 0:
            
            k = 0.
            for i in range(len(self.c)):
                if paux[i] > 0. and self.pc[i] > 0.:
                    k += self.pc[i] * np.log(self.pc[i]/paux[i]) * self.bin_size_c
                
                elif paux[i] > 0. and self.pc[0] == 0.:
                    k += 1e-10*np.log(1e-10/paux[i]) * self.bin_size_c
                elif paux[i] == 0. and self.pc[i] > 0.:
                    k += self.pc[i] * np.log(self.pc[i]/1e-10) * self.bin_size_c
                
            return k
            
            #return np.sum(kl_div(self.pc, paux)*self.bin_size_c)

        # max of vectorial kullback leibler
        if self.model == 1:
            
            k = np.zeros(len(self.c))
            for i in range(len(self.c)):
                if paux[i] > 0. and self.pc[i] > 0.:
                    k[i] = self.pc[i] * np.log(self.pc[i]/paux[i])*self.bin_size_c
                elif paux[i] > 0. and self.pc[0] == 0.:
                    k[i] = 1e-10*np.log(1e-10/paux[i]) * self.bin_size_c
                elif paux[i] == 0. and self.pc[i] > 0.:
                    k[i] = self.pc[i] * np.log(self.pc[i]/1e-10)
                else:
                    k[i] = 0.
                    
            
            return np.amax(k)
            
            
            #return np.sum(kl_div(paux, self.pc)*self.bin_size_c)

        # Variance
        if self.model == 2:
            return np.sum((paux - self.pc)**2 * self.bin_size_c)
        
        # Infinite Norm
        if self.model == 3:
            return np.amax(np.abs(paux - self.pc))
    
    def PrintLoss(self, params):
        loss = self.loss(params)
        a, b, c = params
        print(f'[{self.temperature_index+1}][{self.attempt}] {loss:.4e}\t[{a:.2e} {b:.2e} {c:.2e}]')
    
    def Optimizer(self, range_a, range_b, range_c, initial_guess, tolerance):
        
        #print(self.f())
        bounds = (range_a, range_b, range_c)
        res = minimize(fun=self.loss,
                       #method='COBYLA',
                       x0=initial_guess,
                       bounds=bounds,
                       #tol=tolerance,
                       callback=self.PrintLoss)
        
        #print(f'\n\n{res.x}\n\n')
        a, b, c = res.x
        paux = self.paux(res.x)

        
        return a, b, c, paux, self.loss(res.x)
