import numpy as np
import os
import sys
import KL


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


size = int(sys.argv[1])
pt = int(sys.argv[2])
npt = int(sys.argv[3])
tmin = float(sys.argv[4])
tmax = float(sys.argv[5])
temperatures = GetTemperatures(tmin, tmax, npt)
os.system('clear && pfetch')
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

bin_size_q = 2./bins_q
bin_size_c = 2./bins_c

for k in range(npt):
    norm_q = np.sum(pq[k]*bin_size_q)
    norm_c = np.sum(pc[k]*bin_size_c)
    pq[k] /= norm_q
    pc[k] /= norm_c

momentum_flag = int(sys.argv[6])
model = int(sys.argv[7])
if momentum_flag == 0 or momentum_flag == 1:
    if momentum_flag == 0:
        momentum_bool = False
    elif momentum_flag == 1:
        momentum_bool = True

    if momentum_bool:
        if model == 0:
            filename_loss = f'loss_momentum_size{size}.dat'
            file_par = open(f'parameters_momentum_size{size}.dat', "w")
            file = open(f'aux_dist_momentum_size{size}.dat', "w")
        elif model == 1:
            filename_loss = f'loss_alt_momentum_size{size}.dat'
            file_par = open(f'parameters_alt_momentum_size{size}.dat', "w")
            file = open(f'aux_dist_alt_momentum_size{size}.dat', "w")
        elif model == 2:
            filename_loss = f'loss_normInf_momentum_size{size}.dat'
            file_par = open(f'parameters_normInf_momentum_size{size}.dat', "w")
            file = open(f'aux_dist_normInf_momentum_size{size}.dat', "w")
        elif model == 3:
            filename_loss = f'loss_norm2_momentum_size{size}.dat'
            file_par = open(f'parameters_norm2_momentum_size{size}.dat', "w")
            file = open(f'aux_dist_norm2_momentum_size{size}.dat', "w")
    else:
        if model == 0:
            filename_loss = f'loss_size{size}.dat'
            file_par = open(f'parameters_size{size}.dat', "w")
            file = open(f'aux_dist_size{size}.dat', "w")
        elif model == 1:
            filename_loss = f'loss_alt_size{size}.dat'
            file_par = open(f'parameters_alt_size{size}.dat', "w")
            file = open(f'aux_dist_alt_size{size}.dat', "w")
        elif model == 2:
            filename_loss = f'loss_normInf_size{size}.dat'
            file_par = open(f'parameters_normInf_size{size}.dat', "w")
            file = open(f'aux_dist_normInf_size{size}.dat', "w")
        elif model == 3:
            filename_loss = f'loss_norm2_size{size}.dat'
            file_par = open(f'parameters_norm2_size{size}.dat', "w")
            file = open(f'aux_dist_norm2_size{size}.dat', "w")
    file_loss = open(filename_loss, "w")
    #npt = 1
    for k in range(npt):
        minimizer = KL.Minimizer(pq=pq, pc=pc,
                                 phi=phi, c=overlap,
                                 bins_c=bins_c,
                                 bins_q=bins_q,
                                 bin_size_c=bin_size_c,
                                 bin_size_q=bin_size_q,
                                 temperature_index=k,
                                 temperature=temperatures[k])

        a, b, c, paux = minimizer.GD(epochs=100000,
                                     model=model,
                                     learning_rate=0.05,
                                     print_rate=1,
                                     momentum=momentum_bool,
                                     tolerance=1e-6,
                                     loss_file=file_loss)
        file_par.write(f'{temperatures[k]:.4e}\t{a:.4e}\t{b:.4e}\t{c:.4e}\n')
        file_loss.write("\n\n")
        # print(np.shape(paux), np.shape(c))
        np.savetxt(file, np.c_[overlap, paux, np.full(shape=np.shape(overlap),
                                                      fill_value=temperatures[k])])
        file.write("\n\n")
        # os.system("sleep 10s")
    file.close()
    file_par.close()
    file_loss.close()
else:
    file_map = open(f'map_size{size}.dat', "w")
    for k in range(npt):
        minimizer = KL.Minimizer(pq=pq, pc=pc,
                                 phi=phi, c=overlap,
                                 bins_c=bins_c,
                                 bins_q=bins_q,
                                 bin_size_c=bin_size_c,
                                 bin_size_q=bin_size_q,
                                 temperature_index=k,
                                 temperature=temperatures[k])
        minimizer.Mapper(range_a=2., range_b=2., range_c=2., points=100, model=model, file=file_map)
    file_map.close()
