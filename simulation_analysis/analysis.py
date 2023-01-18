from tqdm import tqdm
import tensorflow as tf
import numpy as np
import os
import random

#custom libraries
import loading
import preliminary
import spectrum
import overlaps
import disorder
import kullback



gpu_devices = tf.config.experimental.list_physical_devices('GPU')

only_kld = 1
overlap_scatter = 0

#eq = 1 PT, eq = 0 NO PT
os.system("clear")
tmax = 1.3
tmin = 0.4
nsamples = 50
settings = np.array(loading.import_settings())
N = int(settings[0])
nRep = int(settings[1])
NSTEP = int(settings[2])
NITERATIONS = int(settings[3]/NSTEP)
NITER_PRINT_CONF = int(settings[4])
eq = int(settings[7])


if eq == 1:
    NITER_MIN_PRINT = int(NITERATIONS/2)
else:
    NITER_MIN_PRINT = 0

NPT = int(settings[6])

del settings

dITER = NITER_PRINT_CONF
IT_MIN = int(NITER_MIN_PRINT)
IT_MAX = int((NITERATIONS-NITER_PRINT_CONF))
loading.WelcomeScreen(N, eq, NPT, nRep, nsamples)
#IT_MAX = 16
iter = list(np.arange(IT_MIN, IT_MAX + dITER, dITER))
#print(iter)
temp = loading.getTemp(tmin, tmax, NPT)
#nsamples = 1
#nRep = 1


if eq == 1:
    if os.path.exists(f'./data/eq/N{N}') == False:
        loading.createstructure(eq, N)
else:
    if os.path.exists(f'./data/neq/N{N}') == False:
        loading.createstructure(eq, N)



if only_kld == 0:
    ifo = []
    q = []
    for n in range(0, nsamples):
        f = loading.getFrequencies(eq, n)   #frequency tensor
         
        equilibrium_energy = []
        a = []
        for r in tqdm(range(0, nRep), f'Loading Sample {n+1}'):
            ar = []
            for i in range(0, len(iter)):
                ar.append(loading.loadconf(n, r, iter[i], eq, N, NPT))  #aquiring configurations from files
            a.append(ar)
            energy = loading.load_parallel(n, r, eq) #loading parallel tempering files (only energies) [time, temp]
            equilibrium_energy.append(preliminary.energy(energy, n, r, temp)) #time block average of energy to check if the sample is thermalized
            del energy
        print(f'Sample {n+1} loaded in the GPU')
        a = tf.convert_to_tensor(a, dtype=tf.float64)
        #print(a)
        intensity = spectrum.spectrum(a)
        
        #print(intensity)
        preliminary.specific_heat(equilibrium_energy, temp, n)
        if eq == 1:
            
            mean_intensity = tf.math.reduce_mean(tf.math.reduce_mean(intensity, axis = 0), axis = 0) #compute mean intensity
            spectrum.print_mean_spectrum(mean_intensity, f, n, temp) #print spectrum doing the mean over real replicas as well
            #print(mean_intensity)
            if nRep > 1 :
                ifo.append(overlaps.compute_ifo_PT(intensity = intensity, 
                                        temperature_array=temp,
                                        number_of_sample=n,
                                        print_dist_flag=True,
                                        print_overlap_flag=True,
                                        experimental_flag=False,
                                        #return_type='overlap',
                                        size = N, 
                                        npt = NPT))  #computing IFO
                
                overlaps.compute_ifo_PT(intensity = intensity, 
                                        temperature_array=temp,
                                        number_of_sample=n,
                                        print_dist_flag=True,
                                        print_overlap_flag=True,
                                        experimental_flag=True,
                                        size = N, 
                                        npt = NPT)  #computing IFO
                
                q.append(overlaps.compute_parisi_PT(configuration=a,
                                                number_of_sample=n,
                                                temperature_array=temp,
                                                npt = NPT,
                                                size=N,
                                                print_overlap_flag=True,
                                                return_type=None,
                                                print_dist_flag=True)) #computing parisi overlaps
        print('------------------------------')

    disorder.disorder_average_specific_heat(nsamples, temp)
if nRep > 1:
    if only_kld == 0:
        disorder.dist_disorder_average('dist_parisi', nsamples, temp)  #disorder average for parisi overlaps
        disorder.dist_disorder_average('dist_theo_ifo', nsamples, temp)   #disorder average for IFOs
        disorder.dist_disorder_average('dist_exp_ifo', nsamples, temp)
        
if overlap_scatter == 1:
    overlaps.Scatter(nsamples, temp)



