import numpy as np
import tensorflow as tf

from loading import NPT
from loading import SIZE
from loading import NREP

def energy(energy_array, 
            number_of_sample, 
            real_replica_index, 
            temperature_array,
            print_flag = False):
    npt = NPT() - 1
    size = SIZE()
    max_iteration = len(energy_array)
    temp_min = 0
    temp_max = 2
    
    mean_block_energy = []
    std_block_energy = []
    blocks = 0
    while temp_max <= max_iteration:
        #print(f'{temp_min}\t{temp_max}\n')
        temp_energy = size*energy_array[temp_min: temp_max, :]
        #print(np.shape(temp_energy))
        mean_block_energy.append(np.mean(temp_energy, axis=0))
        std_block_energy.append(np.std(temp_energy, axis=0))
        temp_min = temp_max
        temp_max *= 2
        blocks += 1

    if print_flag:
        path = f'energy_nrep{real_replica_index}_sample{number_of_sample + 1}.dat'
        file = open(path, "w")
    
        for k in range(npt):
            for i in range(blocks):
                file.write(f'{64*2**i}\t{mean_block_energy[i][npt - 1 - k]:.4e}\t{std_block_energy[i][npt - 1 - k]:.4e}\t{temperature_array[k]}\n')
            file.write("\n\n")
        file.close()
    #print(np.shape(mean_block_energy))
    #print(mean_block_energy)
    return temp_energy


def specific_heat(energy, temperature_array, number_of_sample):
    #energy [replica, time, npt]
    nrep = NREP()
    size = SIZE()
    energy = tf.convert_to_tensor(energy)
    cv = tf.math.reduce_std(energy, axis = 1)**2
    
    npt = len(cv[0])
    path = f'specific_heat_sample{number_of_sample + 1}.dat'
    file = open(path, "w")
    
    for r in range(0, nrep):
        for k in range(npt):
            file.write(f'{temperature_array[k]:.4e}\t{cv[r, npt-k-1]/(size*temperature_array[k]**2):.4e}\t{r}\n')
        file.write("\n\n")
    file.close()
        

    
