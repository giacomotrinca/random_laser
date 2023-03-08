import numpy as np
import KL
import sys


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


bins = 80
x, pq, pc = genGaussian(edge=1., samples=1e8, bins=bins)
momentum_flag = int(sys.argv[1])
size = 100
model = int(sys.argv[2])
bin_size = 2./bins
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
    pc = [pc, 1]
    pq = [pq, 1]
    minimizer = KL.Minimizer(pc=pc,
                             pq=pq,
                             phi=x,
                             c=x,
                             bin_size_c=bin_size,
                             bin_size_q=bin_size,
                             bins_c=bins,
                             bins_q=bins,
                             temperature_index=0,
                             temperature=1.)
    a, b, c, paux = minimizer.GD(epochs=int(1e5),
                                 model=model,
                                 learning_rate=0.01,
                                 print_rate=1,
                                 momentum=momentum_bool,
                                 tolerance=1e-5,
                                 loss_file=file_loss)
    file_par.write(f'{a:.4e}\t{b:.4e}\t{c:.4e}\n')
    # file_loss.write("\n\n")
    np.savetxt(file, np.c_[x, paux])

