import os
import numpy as np

def list_directories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

def find_PT(string):
    if "PT" in string:
        return True
    else:
        return False

def print_settings(settings):
    N, nr, print_start, niter, first_iter, npt, pt_flag, print_every = settings

    print(f'SIZE: {N}\n NREP: {nr}\n Start: {print_start}\n NITER: {niter}\n first: {first_iter}\n NPT: {npt} \n pt_flag: {pt_flag}\n print_every: {print_every}\n')
    

def load_settings(path):
    settings = np.loadtxt(f'{path}/analysis_input.info', dtype=np.float64)
    return settings