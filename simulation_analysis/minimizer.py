import KL2
import sys
import numpy as np
import os


size = int(sys.argv[1])
if size == 100:
    gaussian_flag=True
else:
    gaussian_flag=False

def LoadDist(type, size, pt_flag, block=None):
    if pt_flag:
        pt = 'PT'
    else:
        pt = f'block{block}'

    filename = f'N{size}/dis_ave/dis_ave_{type}_dist_{pt}_size{size}.dat'
    p = np.loadtxt(filename, dtype=np.float64, delimiter="\t", usecols=[0, 1])
    return p[:, 0], p[:, 1]

def GetTemperatures(tmin, tmax, npt):
    step = (tmax - tmin)/npt
    temperatures = np.array([tmin + (i+1) * step
                            for i in range(npt)],
                            dtype=np.float64)
    return temperatures

def genGaussian(mean=0, std=1, edge=1., samples=1e6, bins=1000):
    g = np.random.normal(mean, std, int(samples))
    bin_size = 2*edge/bins
    h = np.histogram(g, bins=bins, range=[-edge, edge], density=True)[0]
    x = np.linspace(-edge + bin_size, edge-bin_size, bins)
    file = open('gaussian.dat', "w")
    np.savetxt(file, np.c_[x, h])
    file.close()
    #g = np.random.normal(mean, std, int(samples))
    h2 = np.histogram(g**2, bins=bins, range=[-edge, edge], density=True)[0]
    
    file = open('gaussian2.dat', "w")
    np.savetxt(file, np.c_[x, h2])
    file.close()
    #print(np.sum(h2*bin_size))
    return x, h, h2

size = int(sys.argv[1])
pt = int(sys.argv[2])
npt = int(sys.argv[3])
tmin = float(sys.argv[4])
tmax = float(sys.argv[5])
model = int(sys.argv[6])
var_flag = int(sys.argv[7])
if var_flag==0:
    var_flag = False
elif var_flag == 1:
    var_flag = True
else:
    sys.exit(-1)
if model == 0:
    if var_flag:
        t = 'var_KL'
    else:
        t = 'KL'
elif model == 1:
    if var_flag:
        t = 'var_vecKL'
    else:
        t = 'vecKL'
elif model == 2:
    if var_flag:
        t = 'var_Norm2'
    else:
        t = 'Norm2'
elif model == 3:
    if var_flag:
        t = 'var_NormInf'
    else:
        t = 'NormInf'
else:
    sys.exit(-1)

temperatures = GetTemperatures(tmin, tmax, npt)
# os.system('clear && pfetch')

edge = 1.
if gaussian_flag:
    phi, pq, pc = genGaussian(bins=1000, samples=int(1e8), edge=edge)
    overlap = phi
    pc = [pc, 1] 
    pq = [pq, 1]
    #print(pq[0])
    npt = 1
# getting the overlap distributions

else:
    edge = 1.
    phi, pq = LoadDist('parisi', size, pt)
    c, pc = LoadDist('theo_ifo', size, pt)
    phi = np.reshape(phi, newshape=(npt, -1))
    c = np.reshape(c, newshape=(npt, -1))
    pq = np.reshape(pq, newshape=(npt, -1))
    pc = np.reshape(pc, newshape=(npt, -1))

    overlap = c[0]
    phi = phi[0]



bins_q = len(phi)
bins_c = len(overlap)

bin_size_q = 2.*edge/bins_q
bin_size_c = 2.*edge/bins_c

range_a = (1e-8, None)
range_b = (0., 1-bin_size_c/2)
range_c = (0., 1.)

attempts = int(sys.argv[8])

file_param = open(f'parameters_model_{t}_attempts{attempts}_size{size}.dat', "w")
file_aux = open(f'aux_dist_model_{t}_attempts{attempts}_size{size}.dat', "w")

for k in range(npt):
    norm_q = np.sum(pq[k]*bin_size_q)
    norm_c = np.sum(pc[k]*bin_size_c)
    pq[k] /= norm_q
    pc[k] /= norm_c
    
    K = 1000.
    aux_dist = np.zeros(len(pc))
    params = np.zeros(3)
    for j in range(attempts):
        initial_guess = (np.random.uniform(1e-8, 10), np.random.uniform(-1., 1.), np.random.uniform(0., 1.))
        #initial_guess = (2., -1., 0.001)
        optimizer = KL2.Minimizer(pq = pq[k],
                              pc = pc[k],
                              phi = phi,
                              c = overlap,
                              bins_c=bins_c,
                              bins_q=bins_q,
                              bin_size_c=bin_size_c,
                              bin_size_q=bin_size_q,
                              temperature_index=k,
                              temperature=temperatures,
                              model=model,
                              var_flag=var_flag,
                              edge=edge,
                              attempt = j+1)
        ai, bi, ci, paux, loss = optimizer.Optimizer(range_a=range_a,
                                                    range_b=range_b,
                                                    range_c=range_c,
                                                    tolerance=1e-10,
                                                    initial_guess=initial_guess)
        
        if loss == 0.:
            continue
        if loss < K:
            K = loss
            params[0] = ai
            params[1] = bi
            params[2] = ci
            aux_dist = paux

        
    print('\n')
    np.savetxt(file_aux, np.c_[overlap, aux_dist, np.full(shape = np.shape(overlap), fill_value=temperatures[k])], delimiter="\t", newline="\n", fmt="%4e")
    file_param.write(f'{temperatures[k]:.4e}\t{loss:.4e}\t{params[0]:.4e}\t{params[1]:.4e}\t{params[2]:.4e}\n')

file_aux.close()
file_param.close()


    


