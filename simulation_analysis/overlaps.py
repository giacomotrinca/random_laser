#from statistics import mean
from asyncio import FastChildWatcher
import tensorflow as tf
import numpy as np
from loading import NPT
from loading import SIZE
from loading import NOV_EQ
from loading import BIN_SIZE
from loading import NREP
from loading import NFILE
from loading import BINS
from tqdm import tqdm

def print_dist(overlap_tensor, file_name, number_of_sample, temperature_array, experimental_flag):
    
    #the overlap_tensor must be a rank 4 tensor [r1][r2][]
    npt = len(overlap_tensor[0])
    bins = BINS()
    bin_size = 2 / bins
    min_ov = -1.
    max_ov = 1.
    

    

def print_overlap(file_name, overlap, number_of_sample):
    
    #print(tf.shape(overlap))
    overlap = tf.reshape(overlap, (len(overlap), -1))
    file = f'{file_name}_overlaps_sample{number_of_sample+1}.dat'

    np.savetxt(file, overlap)
            
        
        
    
def compute_ifo_PT(intensity, 
                   size, 
                   npt,
                   temperature_array,
                   number_of_sample,
                   print_overlap_flag = True,
                   print_dist_flag = True,
                   experimental_flag = False,
                   return_type = None):
    
    nr  = len(intensity) #number of real replicas
    niter = len(intensity[0])
    
    if not (print_overlap_flag or print_dist_flag):
        print(f'Warning: one of print flag must be {True}.')

    #                                      shape =(replica,time,temp,mode)  
    intensity      = tf.reshape(intensity, shape = (nr, niter, npt, size))
    mean_intensity = tf.reduce_mean(intensity, axis=1, keepdims=(not experimental_flag))

    if experimental_flag:
        ifo_type = 'exp_ifo'
        mean_intensity_replica = tf.reduce_mean(mean_intensity, axis=0, keepdims=True)
        delta = mean_intensity-mean_intensity_replica
        #print(tf.shape(delta))00
    else:
        ifo_type = 'theo_ifo'
        delta = intensity - mean_intensity

    

    #computing sum_j delta[r1][i][k][j]*delta[r2][i][k][j]

    ifo = []
    for r1 in tqdm(range(nr), f'Computing {ifo_type}...'):
        ifo2 = []
        for r2 in range(nr):
            if experimental_flag:
                sum_element = delta[r1, :, :] * delta[r2, :, :]
                sum1 = tf.reduce_sum(tf.math.square(delta[r1, :, :]), axis = 1)
                sum2 = tf.reduce_sum(tf.math.square(delta[r2, :, :]), axis = 1)
                norm = tf.multiply(sum1, sum2)
                del sum1
                del sum2
                sum_element = tf.reduce_sum(sum_element, axis = 1)/tf.math.sqrt(norm)
                del norm
            else:
                sum_element = delta[r1, :, :, :] * delta[r2, :, :, :]
                sum1 = tf.reduce_sum(tf.math.square(delta[r1, :, :, :]), axis = 2)
                sum2 = tf.reduce_sum(tf.math.square(delta[r2, :, :, :]), axis = 2)
                norm = tf.multiply(sum1, sum2)
                del sum1
                del sum2
                sum_element = tf.reduce_sum(sum_element, axis = 2)/tf.math.sqrt(norm)
                del norm
            #print(sum_element)
            ifo2.append(sum_element)     
            
        ifo.append(ifo2)
        del ifo2
    ifo = tf.convert_to_tensor(ifo, dtype=tf.float64)
    #print(tf.shape(ifo))
    #temp time rep rep
    if experimental_flag:
        ifo = tf.transpose(ifo, perm = [2, 1, 0])
    else:
        ifo = tf.transpose(ifo, perm = [3, 2, 1, 0])
    if print_dist_flag or print_overlap_flag:
        new_ifo = []
        for k in range(npt):
            new_ifo_k = []
            for r1 in range(nr):
                for r2 in range(r1):
                    if experimental_flag:
                        new_ifo_k.append(tf.reshape(ifo[k, r1, r2], shape = (-1)))
                    else:
                        new_ifo_k.append(tf.reshape(ifo[k, :, r1, r2], shape = (-1)))
            new_ifo.append(new_ifo_k)
            del new_ifo_k
        new_ifo = tf.reshape(tf.convert_to_tensor(new_ifo, dtype=tf.float64), shape=(npt, -1))
    del ifo
    
    if print_dist_flag:
        bins = BINS()
        bin_size = 2./bins
        x_of_hist = np.arange(-1 + bin_size/2, 1, bin_size)
        file_handle = open(f'dist_{ifo_type}_sample{number_of_sample+1}.dat', "w")
        file_handle.write('#bins\tdist\ttemp\n')
        save_dist = []

        for k in tqdm(range(npt), f'Printing {ifo_type} dist'):
            dist = np.histogram(a=new_ifo[k].numpy(), range=[-1, 1], bins=bins, density=True)[0]
            temp = np.full(shape = np.shape(dist), fill_value=temperature_array[k], dtype=np.float64)
            ifo_print = np.transpose(np.vstack((x_of_hist, dist, temp), dtype=np.float64))
            np.savetxt(fname=file_handle, X=ifo_print, fmt='%.4e', delimiter='\t', newline='\n')
            del ifo_print
            file_handle.write("\n\n")
            save_dist.append(dist)
        file_handle.close()
        dist = tf.convert_to_tensor(save_dist, dtype=tf.float64)
        del save_dist
        #print(tf.shape(dist))

    if print_overlap_flag:
        file_handle = open(f'overlap_{ifo_type}_sample{number_of_sample+1}.dat', "w")
        file_handle.write('#ifo\t#temp\n')
        for k in tqdm(range(npt), f'Printing {ifo_type}'):
            temp = np.full(shape = np.shape(new_ifo[k].numpy()), fill_value=temperature_array[k], dtype=np.float64)
            ifo_print = np.transpose(np.vstack((new_ifo[k].numpy(), temp), dtype=np.float64))
            np.savetxt(fname=file_handle, X=ifo_print, fmt='%.4e', delimiter='\t', newline='\n')
            del ifo_print
            file_handle.write("\n\n")
        file_handle.close()
    

    if return_type == 'overlap':
        return new_ifo
    elif return_type == 'dist':
        return dist
    


def compute_parisi_PT(configuration, 
                    number_of_sample,
                    temperature_array,
                    npt,
                    size,
                    return_type = None,
                    print_overlap_flag = True,
                    print_dist_flag = True):
    
    #print(configuration)
    #configuration must be a 4-rank tensor with shape (nr, niter, npt*n, 2) where 2 is for real and imaginary part respectively
    nr    = len(configuration) #number of real replicas
    niter = len(configuration[0])
    sigma = configuration[:, :, :, 0]
    tau   = configuration[:, :, :, 1]
    
    sigma = tf.reshape(sigma, shape=(nr, niter, npt, size))
    tau   = tf.reshape(tau  , shape=(nr, niter, npt, size))

    #overlap^(1,2) = 1/N \sum (sigma^1 sigma^2 + tau^1 tau^2)
    overlap = []
    for r1 in range(nr):
        overlap2 = []
        for r2 in range(nr):
            overlap2.append(tf.reduce_mean(sigma[r1, :, :, :]*sigma[r2, :, :, :] + tau[r1, :, :, :]*tau[r2, :, :, :], axis = 2))
        overlap.append(overlap2)
        del overlap2
    del sigma
    del tau
    overlap = tf.convert_to_tensor(overlap, dtype = tf.float64)
    #print(tf.shape(overlap))
    overlap = tf.transpose(overlap, perm=[3, 2, 1, 0])
    if print_dist_flag or print_overlap_flag:
        new_overlap = []
        for k in range(npt):
            new_overlap_k = []
            for r1 in range(nr):
                for r2 in range(r1):          
                    new_overlap_k.append(tf.reshape(overlap[k, :, r1, r2], shape = (-1)))
            new_overlap.append(new_overlap_k)
            del new_overlap_k
        new_overlap = tf.reshape(tf.convert_to_tensor(new_overlap, dtype=tf.float64), shape=(npt, -1))
    del overlap
    overlap_type = 'parisi'
    if print_dist_flag:
        bins = BINS()
        bin_size = 2./bins
        x_of_hist = np.arange(-1 + bin_size/2, 1, bin_size)
        file_handle = open(f'dist_{overlap_type}_sample{number_of_sample+1}.dat', "w")
        file_handle.write('#bins\tdist\ttemp\n')
        save_dist = []

        for k in tqdm(range(npt), f'Printing {overlap_type} dist'):
            dist = np.histogram(a=new_overlap[k].numpy(), range=[-1, 1], bins=bins, density=True)[0]
            temp = np.full(shape = np.shape(dist), fill_value=temperature_array[k], dtype=np.float64)
            overlap_print = np.transpose(np.vstack((x_of_hist, dist, temp), dtype=np.float64))
            np.savetxt(fname=file_handle, X=overlap_print, fmt='%.4e', delimiter='\t', newline='\n')
            del overlap_print
            file_handle.write("\n\n")
            save_dist.append(dist)
        file_handle.close()
        dist = tf.convert_to_tensor(save_dist, dtype=tf.float64)
        del save_dist
        #print(tf.shape(dist))

    if print_overlap_flag:
        file_handle = open(f'overlap_{overlap_type}_sample{number_of_sample+1}.dat', "w")
        file_handle.write('#overlap\t#temp\n')
        for k in tqdm(range(npt), f'Printing {overlap_type}'):
            temp = np.full(shape = np.shape(new_overlap[k].numpy()), fill_value=temperature_array[k], dtype=np.float64)
            overlap_print = np.transpose(np.vstack((new_overlap[k].numpy(), temp), dtype=np.float64))
            np.savetxt(fname=file_handle, X=overlap_print, fmt='%.4e', delimiter='\t', newline='\n')
            del overlap_print
            file_handle.write("\n\n")
        file_handle.close()
    

    if return_type == 'overlap':
        return new_overlap
    elif return_type == 'dist':
        return dist
    
    
    


def printAuxDist(file, aux_dist, temperature, a, b, c, k, attempts_flag = False, attempts = None, kld = None):
    bins = BINS()
    bin_size = BIN_SIZE()
    x_of_histogram = tf.range(-1. + bin_size/2, 1., bin_size)
    parameters = [a, b, c]
    if attempts_flag and (attempts == None or kld == None):
        print('Warning: attempts_flag variable requires attempts and kld variable nequal to None.')

    if not attempts_flag:
        file.write(f'#AUX DISTRIBUTION for T = {temperature}, a = {a:.4e}, b = {b:.4e}, c = {c:.4e}, KLD = {k:.4e}\n')

    for i in range(bins):
        if attempts_flag:
            file.write(f'{x_of_histogram[i]:.4e}\t{aux_dist[i]:.4e}\t{temperature:.4e}\t{kld:.4e}\t{attempts}\t{parameters[0]:.4e}\t{parameters[1]:.4e}\t{parameters[2]:.4e}\n')
        else:
            file.write(f'{x_of_histogram[i]:.4e}\t{aux_dist[i]:.4e}\t{temperature:.4e}\n')
    file.write("\n\n")
    #file.close()
    del parameters
    

def load():
    def Overlap(file_name, npt, number_of_sample):
    
        file = f'overlap_{file_name}_sample{number_of_sample+1}.dat'
        overlap = np.loadtxt(file, usecols=[0])
        #print(np.shape(overlap))
        overlap = tf.reshape(tf.convert_to_tensor(overlap), shape=(npt, -1))
        #print(tf.shape(overlap))
        #print(tf.shape(overlap))
        return overlap


    def Dist(file_name, number_of_sample):
        file = f'{file_name}_sample{number_of_sample+1}.dat'
        dist = np.loadtxt(file, usecols=[0, 1])
    
        prob = dist[:,1]
        npt = NPT()
        prob = tf.convert_to_tensor(prob)
        #print(tf.shape(prob))
        prob = tf.reshape(prob, shape = (npt, -1))
    
        return prob

def Scatter(nsamples, temp):
    bins = BINS()
    bin_size = BIN_SIZE()
    x_of_histogram = tf.range(-1. + bin_size/2, 1., bin_size)
    overlap = []
    ifo = []
    for n in tqdm(range(nsamples), 'Loading Overlaps'):
        overlap.append(load.Overlap(file_name='parisi', npt=NPT, number_of_sample=n))
        ifo.append(load.Overlap(file_name='theo_ifo', npt = NPT, number_of_sample=n))
    overlap = tf.transpose(tf.convert_to_tensor(overlap, dtype=tf.float64), perm=[1, 0, 2])
    ifo = tf.transpose(tf.convert_to_tensor(ifo, dtype=tf.float64), perm = [1, 0, 2])
    file_dist = open('mod_dist_parisi.dat', "w")
    for k in tqdm(range(NPT), 'Printing Overlaps'):
        file_handle = open(f'mod_overlap_temp{k}.dat', "w")
        file_handle.write('#overlap\tifo\tsample\n')
        for n in range(nsamples):
            sample = np.full(shape = np.shape(overlap[k, n].numpy()), fill_value=n)
            print_mat = np.transpose(np.vstack((overlap[k, n].numpy(), ifo[k, n], sample)))
            np.savetxt(fname=file_handle, X=print_mat, fmt='%.4e', delimiter='\t', newline='\n')
            file_handle.write("\n\n")
        file_handle.close()
        del sample
        del print_mat
        h = np.histogram(a = overlap[k, :, :], range = [-1,1], bins=bins, density=True)[0]
        #print(np.shape(h))
        t = np.full(shape = np.shape(h), fill_value=temp[k])
        print_dist = np.transpose(np.vstack((x_of_histogram.numpy(), h, t)))
        np.savetxt(fname=file_dist, X=print_dist, fmt='%.4e', delimiter='\t', newline='\n')
        file_dist.write("\n\n")
    file_dist.close()