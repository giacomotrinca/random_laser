import numpy as np
import tensorflow as tf
import pandas as pd
import os
from tqdm import tqdm


def WelcomeScreen(N, eq, NPT, nRep, nsamples):
    print(f'°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°')
    print(f'°°°°°° SMrandomTETRADS Analysis °°°°°°°°°')
    print(f'°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°')
    print(f'Samples = {nsamples} N = {N} nREP = {nRep} NPT = {NPT}')
    if eq == 1:
        print(f'**** REPLICA EXCHANGE ****')
    else:
        print(f'**** DYNAMICS METROPOLIS ****')
    print(f'°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°')


def BINS():
    return 90

def BIN_SIZE():
    return 2./BINS()

    
def import_settings():
    path = f'SOURCE/SMrandomTetradsRUNCHECK.h'
    settings = pd.read_csv(path, sep = "\t", header = None, comment = "/", usecols = [2], skiprows = [21, 24])
    #print(settings)
    return settings

def NPT():
    data = np.array(import_settings())
    return int(data[6])

def SIZE():
	data = np.array(import_settings())
	return int(data[0])

def NREP():
	data = np.array(import_settings())
	return int(data[1])
	
def NOV_EQ():
	data = np.array(import_settings())
	NSTEP = int(data[2])
	NITERATIONS = int(data[3]/NSTEP)
	NITER_PRINT_CONF = int(data[4])	
	NITER_MIN_PRINT = int(NITERATIONS/2)
	dITER = NITER_PRINT_CONF
	IT_MIN = int(NITER_MIN_PRINT)
	IT_MAX = int((NITERATIONS-NITER_PRINT_CONF))
	iter = int((IT_MAX - IT_MIN)/dITER) + 1
	
	nr = NREP()
	
	return int((nr*(nr-1)/2)*iter)

def NFILE():
    data = np.array(import_settings())
    NSTEP = int(data[2])
    NITERATIONS = int(data[3]/NSTEP)
    NITER_PRINT_CONF = int(data[4])	
    NITER_MIN_PRINT = int(NITERATIONS/2)
    dITER = NITER_PRINT_CONF
    IT_MIN = int(NITER_MIN_PRINT)
    IT_MAX = int((NITERATIONS-NITER_PRINT_CONF))
    niter = int((IT_MAX - IT_MIN)/dITER) + 1
    
    return niter
	
 
def NOV_NEQ():
	data = np.array(import_settings())
	NSTEP = int(data[2])
	NITERATIONS = int(data[3]/NSTEP)
	NITER_PRINT_CONF = int(data[4])	
	NITER_MIN_PRINT = 0
	dITER = NITER_PRINT_CONF
	IT_MIN = int(NITER_MIN_PRINT)
	IT_MAX = int((NITERATIONS-NITER_PRINT_CONF))
	iter = int((IT_MAX - IT_MIN)/dITER) + 1
	
	nr = NREP()
	
	return int((nr*(nr-1)/2)*iter)	

def AQUIRE_DOMAIN():
    return 0

def getTemp(tmin, tmax, NPT):
    temp = np.zeros(NPT)
    dt = (tmax - tmin)/NPT

    for i in range(0, NPT):
        temp[NPT - i -1] = tmax - i*dt

    temp = tf.convert_to_tensor(temp)
    return temp

def getFrequencies(equilibrium_flag, number_of_sample):
    
    size = SIZE()
    if equilibrium_flag == 1:
        path = f'DEV_LOCAL/random-eq/N{size}/sample{number_of_sample+1}/frequencies.dat'
    else:
        path = f'DEV_LOCAL/random-neq/N{size}/sample{number_of_sample+1}/frequencies.dat'
        
    f = np.loadtxt(path, dtype=np.float64)
    #print(np.shape(f))
    return f[:, 1]

def loadconf(number_of_sample, real_replica_index, mcs_index, equilibrium_flag, size, npt):
    
    if equilibrium_flag == 1:
        path = f'./DEV_LOCAL/random-eq/N{size}/sample{number_of_sample+1}/config_nrep{real_replica_index}_iter_{mcs_index}.dat'
    else:
        path = f'./DEV_LOCAL/random-neq/N{size}/sample{number_of_sample + 1}/config_nrep{real_replica_index}_iter_{mcs_index}.dat'

    conf = tf.Variable(np.loadtxt(path, usecols=[2, 3]), dtype=tf.float64, shape=(size*npt, 2))

    return conf


def load_parallel(number_of_sample, real_replica_index, equilibrium_flag):
    
    size = SIZE()
    npt = NPT()
    
    if equilibrium_flag == 1:
        path = f'DEV_LOCAL/random-eq/N{size}/sample{number_of_sample+1}/parallel_tempering{real_replica_index}.dat'
    else:
        path = f'DEV_LOCAL/random-neq/N{size}/sample{number_of_sample+1}/parallel_tempering{real_replica_index}.dat'
        
        
    cols = np.arange(3, 3 * npt, 3)
    parallel = np.loadtxt(path, usecols=cols)
    
    #parallel = tf.convert_to_tensor(parallel)
    
    #parallel[time, temp]
    
    
    return parallel
    
def load_dist(file_name):
    
    bins = BINS()
    npt = NPT()
    
    path = f'dis_ave_{file_name}.dat'
    
    dist = np.loadtxt(path, usecols=[0, 1])
    
    dist = tf.convert_to_tensor(dist)
    
    dist = tf.reshape(dist, shape = (npt, bins, 2))
    
    return dist[:, :, 0], dist[:, :, 1]
    
    
    


def createstructure(eq, N):
    if os.path.exists(f'./data') == False:
        os.system('mkdir data')
    if eq == 1:
        os.system('mkdir data/eq')
        os.system(f'mkdir data/eq/N{N}')
        os.system(f'mkdir data/eq/N{N}/parisi')
        os.system(f'mkdir data/eq/N{N}/spectrum')
        #os.system(f'mkdir data/eq/N{N}/kld')
        os.system(f'mkdir data/eq/N{N}/ifo')
        os.system(f'mkdir data/eq/N{N}/delta')
        os.system(f'mkdir data/eq/N{N}/spec_heat')
        os.system(f'mkdir data/eq/N{N}/energy')
        #os.system(f'mkdir data/eq/N{N}/dist')

    else:
        os.system('mkdir data/neq')
        os.system(f'mkdir data/neq/N{N}')
        os.system(f'mkdir data/neq/N{N}/parisi')
        os.system(f'mkdir data/neq/N{N}/spectrum')
        #os.system(f'mkdir data/neq/N{N}/kld')
        os.system(f'mkdir data/neq/N{N}/ifo')
        os.system(f'mkdir data/neq/N{N}/delta')
        os.system(f'mkdir data/neq/N{N}/spec_heat')
        os.system(f'mkdir data/neq/N{N}/energy')
        #os.system(f'mkdir data/neq/N{N}/dist')

def organizeFiles(eq, N):
    if eq == 1:
        os.system(f'mv spectrum_*.dat data/eq/N{N}/spectrum/')
        os.system(f'mv parisi_*dat data/eq/N{N}/parisi/')
        os.system(f'mv *_ifo_*.dat data/eq/N{N}/ifo/')
        os.system(f'mv delta_*.dat data/eq/N{N}/delta/')
        os.system(f'mv energy_*.dat data/eq/N{N}/energy/')
        os.system(f'mv specific_*.dat data/eq/N{N}/spec_heat/')
        #os.system(f'mv dist_*.dat data/eq/N{N}/dist/')
        #os.system(f'mv min_kld_*.dat data/eq/N{N}/kld/')
        os.system(f'mv dis_ave_*.dat data/eq/N{N}/')

    else:
        os.system(f'mv spectrum_*.dat data/neq/N{N}/spectrum/')
        os.system(f'mv parisi_*dat data/neq/N{N}/overlaps/')
        os.system(f'mv *_ifo_*.dat data/neq/N{N}/ifo/')
        os.system(f'mv delta_*.dat data/neq/N{N}/delta/')
        os.system(f'mv energy_*.dat data/neq/N{N}/energy/')
        os.system(f'mv specific_*.dat data/neq/N{N}/spec_heat/')
        #os.system(f'mv min_kld_*.dat data/neq/N{N}/kld/')
        #os.system(f'mv dist_*.dat data/neq/N{N}/dist/')
        os.system(f'mv dis_ave_*.dat data/neq/N{N}/')

    #os.system(f'rm *.dat')
    


def fill_histogram(dist,
                   points_per_bin = 100,
                   bins = BINS(),
                   bin_size = BIN_SIZE(), 
                   print_to_file_flag = False, 
                   file_name = None,
                   min_ov = -1.,
                   max_ov = 1.):
  
    fill_x = tf.range(min_ov, max_ov, bin_size/points_per_bin)
    
    fill_y = []
    
    for i in range(bins):
        for j in range(i*points_per_bin, (i+1)*points_per_bin):
            fill_y.append(dist[i])
            
    
    fill_y = tf.convert_to_tensor(fill_y)
   
    
    if print_to_file_flag:
        path = f'{file_name}.dat'
        file = open(path, "w")
        for i in range(len(fill_y)):
            file.write(f'{fill_x[i]:.4e}\t{fill_y[i]:.4e}\n')
    
    return fill_x, fill_y
    