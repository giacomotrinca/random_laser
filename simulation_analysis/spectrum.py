import numpy as np
import tensorflow as tf
from tqdm import tqdm

def spectrum(configuration):
    intensity = tf.math.multiply(configuration[:, :, :, 0], configuration[:, :, :, 0]) + tf.math.multiply(configuration[:, :, :, 1], configuration[:, :, :, 1])
    
    return intensity


def print_mean_spectrum(mean_spectrum, frequencies, number_of_sample, temperature_array):
    size = len(frequencies)
    n_blocks = int(len(mean_spectrum) / size)
    
    file = open(f'spectrum_sample{number_of_sample+1}.dat', "w")
    
    l = 0
    
    for k in tqdm(range(0, n_blocks), 'Print spectrum'):
        for j in range(0, size):
            file.write(f'{frequencies[j]}\t{mean_spectrum[l]/np.sqrt(temperature_array[k])}\t{temperature_array[k]}\n')
            l += 1
        file.write(f'\n\n')
    file.close()
            
    
    
    
    
    