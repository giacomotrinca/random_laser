
import numpy as np
import os

def LoadOverlap(sample, npt):
    overlap = np.loadtxt(f'N18/DEV0/20230122_131730_PT/data/overlaps_PT_size18_sample{sample}.dat', dtype = np.float32)

    overlap = np.reshape(overlap, newshape=(npt, -1, 3))
    
    

    return overlap[0, :, 0], overlap[0, :, 1]


# Definizione di p(x) e q(f(x, a, b))
def p(x):
    return np.histogram(x, range=[-1, 1], bins=80, density=True)[0]

def q(x, a, b):
    return a* x**2 + b

def check(x, a, b):
    if np.max(q(x, a, b)) <= 1 and np.min(q(x, a, b)) >= -1:
        return True
    else:
        return False


# Definizione della funzione di perdita come la distanza di KL
def kl_divergence(ifo, aux):
    bin_size = 2./80
    
    k = 0.
    for i in range(len(ifo)):
        #print(ifo[i], aux[i])
        if ifo[i] != 0 and aux[i] != 0:
            k += bin_size*aux[i] * np.log(aux[i]/ifo[i])
    
    #print(k)
    return np.abs(k)
        
def grad(loss, step, ifo, x, a, b):
    
    #print(np.shape(ifo))
    aux_a = p(q(x, a+step, b))

    k_a = kl_divergence(p(ifo), aux_a)
    aux_b = p(q(x, a, b+step))
    k_b = kl_divergence(p(ifo), aux_b)
    #print(loss, k_a, k_b)
    da = loss - k_a
    db = loss - k_b

    return da/step, db/step



x, ifo  = LoadOverlap(1, 20)
# Definizione dei parametri a e b come array numpy
a = 2.0
b = 0.7

# Limite di a a un valore minimo positivo
a = np.maximum(a, 1e-6)

# Definizione del rate di apprendimento
learning_rate = 0.01
step = 0.1
# Ciclo di addestramento
loss_min = 1000.
a_min = 0.
b_min = 0.
file = open('space.dat', "w")
for i in range(1000000000):
    if check(x, a, b):
        file.write(f'{a:.4e}\t{b:.4e}')
        loss = kl_divergence(p(ifo), p(q(x, a, b)))
        if loss == 0:
            continue
        if loss < loss_min:
            a_min, b_min, loss_min = a, b, loss
        #print(p(q(x, a, b)))
        da, db = grad(loss, step, ifo, x, a, b)
        print(i, loss, a, b, check(x, a, b))
        
        a -= learning_rate * da
        b -= learning_rate * db
        a = np.maximum(a, 1)

    
    else:
        #print(a, b)
        while check(x, a, b):
            a = np.random.uniform(0, 2)
            b = np.random.uniform(0, 2)
        #a = np.maximum(a, 1e-6)
        
file.close()

bin_size = 2./80
x_of_histogram = np.linspace(bin_size, 1-bin_size, num = 80)
np.savetxt('aux_dist.dat', np.c_[x_of_histogram, p(q(x, a_min, b_min))])
# Stampa dei valori ottimizzati di a e b
print("\n\nLoss =", loss_min)
print("a =", a_min)
print("b =", b_min)


