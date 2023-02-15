
import numpy as np
import os


def LoadOverlap(sample, npt):
    overlap = np.loadtxt(f'N18/DEV0/20230122_131730_PT/data/\
                        overlaps_PT_size18_sample{sample}.dat',
                        dtype=np.float32)

    overlap = np.reshape(overlap, newshape=(npt, -1, 3))

    return overlap[0, :, 0], overlap[0, :, 1]


# Definizione di p(x) e q(f(x, a, b))
def p(x):
    return np.histogram(x, range=[-1, 1], bins=80, density=True)[0]


def q(x, a, b):
    return a * x**2 + b


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
    return k
        
def grad(loss, step, ifo, x, a, b):
    
    #print(np.shape(ifo))
    aux_a = p(q(x, a+step, b))

    k_a = kl_divergence(p(ifo), aux_a)
    aux_b = p(q(x, a, b+step))
    k_b = kl_divergence(p(ifo), aux_b)
    #print(loss, k_a, k_b)
    da = k_a - loss
    db = k_b - loss

    return da/step, db/step



x, ifo  = LoadOverlap(1, 20)
# Definizione dei parametri a e b come array numpy
a = 1.5
b = 0.7
bin_size = 2./80
x_of_histogram = np.linspace(bin_size-1, 1-bin_size, num = 80)

file = open("dist4.dat", "w")
par = p(x)
ifo_r = p(100 * x**2 + 0.1)

np.savetxt(file, np.c_[x_of_histogram, par, ifo_r])
file.close()
# Limite di a a un valore minimo positivo
a = np.maximum(a, 1e-6)

# Definizione del rate di apprendimento
learning_rate = 0.001
step = 0.02
# Ciclo di addestramento
loss_min = 1000.
a_min = 0.
b_min = 0.
file = open('space.dat', "w")
for i in range(100000):
    if check(x, a, b):
        
        loss = kl_divergence(p(ifo), p(q(x, a, b)))
        if loss == 0:
            continue
        if loss < loss_min:
            a_min, b_min, loss_min = a, b, loss
        #print(p(q(x, a, b)))
        da, db = grad(loss, step, ifo, x, a, b)
        print(f'[{i}]\t{loss:.4e}\t{a:.4f}\t{b:.4f}\t[{da:.4f}, {db:.4f}]\t{check(x, a, b)}')
        file.write(f'{i}\t{loss:.4e}\t{a:.4e}\t{b:.4e}\t{da:.4e}\t{db:.4e}\n')
        a -= learning_rate * da
        b -= learning_rate * db
        a = np.maximum(a, 1e-2)

    
    else:
        #print(a, b)
        r = np.random.uniform(0,1) * np.sqrt(a**2 + b**2)
        phi = np.random.uniform(0, 2*np.pi)
        a = r * np.cos(phi)
        b = r * np.sin(phi)
        a = np.maximum(a, 1e-2)
        
file.close()


np.savetxt('aux_dist.dat', np.c_[x_of_histogram, p(q(x, a_min, b_min))])
# Stampa dei valori ottimizzati di a e b
print("\n\nLoss =", loss_min)
print("a =", a_min)
print("b =", b_min)


