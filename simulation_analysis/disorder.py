from re import X
import numpy as np
import tensorflow as tf
from tqdm import tqdm
from loading import NPT
from loading import SIZE
from loading import NOV_EQ
from loading import NREP
from loading import NFILE
from loading import BINS

def dist_disorder_average(file_name, number_of_samples, temperature_array):
    
    bins = BINS()
    npt = NPT()
    bin_size = 2 / bins
    min_ov = -1.
    max_ov = 1.
    '''
    if file_name == 'delta':
        min_ov = -10.
        max_ov = 10.
    '''
    
    x_of_histogram = tf.range(min_ov + bin_size/2, max_ov, bin_size)
    h = []
    for n in range(0, number_of_samples):
        path = file_name + f'_sample{n+1}.dat'
        h.append(tf.convert_to_tensor(np.loadtxt(path, usecols=[1, 2])))
    
    h = tf.convert_to_tensor(h)
    h = tf.reshape(h, shape=(number_of_samples, npt, bins, 2))
    #print(tf.shape(h))
    std = tf.math.reduce_std(h, axis=0)
    h   = tf.math.reduce_mean(h, axis=0)
    
    path = 'dis_ave_'+ file_name + '.dat'
    file = open(path, "w")
    file.write(f'#{file_name} distribution averaged over {number_of_samples} samples\n')
    file.write(f'#There are {npt} blocks of {bins} bins\n')
    file.write(f'# bin\thisto\tcounter\terror\ttemperature\n')
    for k in tqdm(range(npt), f'{file_name} disorder average'):
        for i in range(0, bins):
            file.write(f'{x_of_histogram[i]:.8e}\t{h[k][i][0]:.8e}\t{h[k][i][1]:.4f}\t{std[k][i][0]:.8e}\t{std[k][i][1]:.8e}\t{temperature_array[k]:.8e}\n')
        file.write("\n\n")
    
    file.close()
    
    return h;


def disorder_average_specific_heat(number_of_samples, temperature_array):
    
    cv = []
    real_replicas = NREP()
    npt = NPT() -1
    for n in tqdm(range(number_of_samples), 'specific_heat disorder average'):
        path = f'specific_heat_sample{n + 1}.dat'
        cv.append(np.loadtxt(path, usecols=[1]))
    
    cv = tf.convert_to_tensor(cv)
    cv = tf.reshape(cv, shape = (number_of_samples, real_replicas, npt))
    
    std = tf.math.reduce_std(cv, axis=0)
    cv = tf.reduce_mean(cv, axis = 0)
    
    std = tf.reduce_mean(std, axis=0)
    cv = tf.reduce_mean(cv, axis=0)
    
    path = f'dis_ave_specific_heat.dat'
    file = open(path, "w")
    
    for i in range(npt):
        file.write(f'{temperature_array[i]}\t{cv[i]}\t{std[i]}\n')
    file.close()
    
    
    
    
        
        
