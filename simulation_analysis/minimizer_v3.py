import numpy as np
from scipy.optimize import minimize
import sys


model = int(sys.argv[1])
bins = 1000
bin_size = 2./bins
def genGaussian(mean=0, std=1, edge=1., samples=1e6, bins=1000):
    g = np.random.normal(mean, std, int(samples))
    bin_size = 2*edge/bins
    h = np.histogram(g, bins=bins, range=[-edge, edge], density=True)[0]
    x = np.linspace(-edge + bin_size, edge-bin_size, bins)
    file = open('gaussian.dat', "w")
    np.savetxt(file, np.c_[x, h])
    file.close()
    h2 = np.histogram(g**2, bins=bins, range=[-edge, edge], density=True)[0]
    file = open('gaussian2.dat', "w")
    np.savetxt(file, np.c_[x, h2])
    file.close()
    return x, h, h2


x, pq, experimental_data = genGaussian(edge=1., samples=1e8, bins=bins)

def Mean(distribution, x):
    return np.sum(distribution * x) / np.sum(distribution)

def f(x):
    #return np.sqrt(1 - x**2)
    #return x**2 + np.abs(x) -1
    return 1.

def g3(q, a, b, c, vec_flag=False):
        
    if vec_flag:
        g = []
        for qi in x:
            g.append(1./(np.sqrt(2.)*c) * (qi - a * x**2 - b))
        return np.array(g, dtype=np.float64)
    else:
        return 1./(np.sqrt(2.)*c) * (x[q] - a * x**2 - b)

def pauxGeneral(a, b, c):
    p = []
    for q in range(len(x)):
        p.append(1./(c*np.sqrt(2.))*Mean(pq, np.exp(-g3(q, a, b, c)**2/f(x)**2)/f(x)))
    return np.array(p, dtype=np.float64)


def objective_function(params):
    a, b, c = params
    analytical_data = pauxGeneral(a, b, c)
    if model == 1:
        L = np.abs(analytical_data-experimental_data)
        return np.amax(L)
    elif model == 2:
        return np.square(analytical_data - experimental_data).sum() * bin_size
    elif model==0:
        k = 0.
        for i in range(len(x)):
            if analytical_data[i] > 0 and experimental_data[i] > 0:
                #k += analytical_data[i] * np.log(analytical_data[i]/experimental_data[i]) * bin_size
                k += experimental_data[i] * np.log(experimental_data[i]/analytical_data[i]) * bin_size
        return k

def print_values(params):
    print(f'{objective_function(params):.4e}\t{params}')

if model == 0:
    t = 'KL'
elif model == 1:
    t = 'NormInf'
elif model == 2:
    t = 'Norm2'
# Vincoli sui parametri
bounds = ((1e-8, None), (None, None), (1e-8, None))

# Valori iniziali dei parametri
initial_guess = (1, 0.05, 0.05)

# Minimizzazione della funzione obiettivo
result = minimize(objective_function, initial_guess, bounds=bounds, callback=print_values)
print("Valori dei parametri a, b, c:", result.x)

a, b, c = result.x
file = open(f'aux_scipy_{t}_{bins}bins.dat', "w")
paux = pauxGeneral(a, b, c)

np.savetxt(file, np.c_[x, paux])
file.close()