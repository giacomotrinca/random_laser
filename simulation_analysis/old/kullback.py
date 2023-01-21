import numpy as np
import tensorflow as tf
import overlaps
from loading import BINS
from loading import NPT
from loading import BIN_SIZE
from loading import NFILE
from tqdm import tqdm
import math
import matplotlib
import matplotlib.pyplot as plt
import loading
import random

def p(overlap):
    h = np.histogram(overlap, bins=BINS(), density=True)
    
    h = tf.convert_to_tensor(h[0])
    
    return h
    

def q_aux(overlap_tensor,  
          a, 
          b, 
          c,
          gaussian_tensor):    
    return a * tf.multiply(overlap_tensor, overlap_tensor) + b*tf.ones(shape=tf.shape(overlap_tensor)) + c * tf.math.multiply(gaussian_tensor, tf.math.sqrt(tf.ones(shape=tf.shape(overlap_tensor)) - tf.multiply(overlap_tensor, overlap_tensor)))

def buildDistribution(parisi_binned,
                      x_of_histogram = tf.range(-1 + loading.BIN_SIZE()/2, 1, loading.BIN_SIZE()),
                      equilibrium_flag = True,
                      factor_random_variable = 1e3):

    norm = np.sum(parisi_binned.numpy())*BIN_SIZE()
    if equilibrium_flag == True:
        multinom = np.random.multinomial(factor_random_variable*loading.NOV_EQ(), pvals = parisi_binned.numpy()*BIN_SIZE()/norm, size=1)[0]
    else:
        multinom = np.random.multinomial(factor_random_variable*loading.NOV_NEQ(), pvals = parisi_binned.numpy()*BIN_SIZE()/norm, size=1)[0]
    
    q = np.zeros(shape=(0))
    for i in range(len(multinom)):
        if multinom[i] != 0:
            q = np.hstack((q, np.random.uniform(x_of_histogram[i]-loading.BIN_SIZE()/2, x_of_histogram[i] + loading.BIN_SIZE()/2, multinom[i])))
    del norm
    del multinom

    q = tf.convert_to_tensor(q)
    q = tf.cast(q, tf.float32)
    q = tf.sort(q)
    return q



def check_constraint(l):
    if tf.reduce_max(l, axis=None) > 1 or tf.reduce_min(l, axis=None) < -1:
        return False
    else:
        return True

def kld(dist_to_test, 
        analytical_dist):
    
    bin_size = BIN_SIZE()
    k = 0.
    q = tf.cast(analytical_dist, tf.float32)
    p = tf.cast(dist_to_test, tf.float32)
    #print(p)
    #print(q)
    for i in range(0, len(q)):
        #print(f'{i}\t{k}')
        if q[i] > 0 and p[i] > 0:
            k += bin_size * p[i] * np.log(p[i] / q[i]) 
             
    
    return k

def truthmatrix(overlap, 
                gaussian_numbers,
                range_of_a = 2., 
                step = 1e-3,
                print_to_file_flag = True,
                print_to_screen_flag = False,
                random_flag = True):
    
    #we observe that the domain is compact. So we need to find the minimum only for
    # -2.<= a <= 2.
    # -1.<= b <= 1
    # -0.5 <= c <= 0.5
    q = overlap
    z = gaussian_numbers
    #a = tf.range(-range_of_a, range_of_a, step)
    a = tf.range(0, range_of_a, step)
    b = tf.range(-range_of_a/2., range_of_a/2., step)
    tmat = []
    counter = 0.
    a_i = len(a)
    b_j = len(b)
    
    c = tf.range(-range_of_a/4., range_of_a/4., step)
    c_k = len(c)
    if print_to_file_flag:
        if random_flag:
            file = open(f'space_R{range_of_a}_step{step}_random.dat', "w")
        else:
            file = open(f'space_R{range_of_a}_step{step}_deterministic.dat', "w")
    
    
    counter_2 = 0
    
    for i in range(a_i):
        for j in range(b_j):
            for k in range(c_k):
                if a[i] != 0:
                    counter_2 += 1
                    aux = q_aux(overlap_tensor= q,
                                a = a[i],
                                b = b[j],
                                c = c[k],
                                gaussian_tensor=z)
                    if check_constraint(aux):
                        tmat.append([a[i], b[j], c[k]])
                        if print_to_screen_flag:
                            print(f'[{counter_2/(a_i * b_j * c_k):.2%}] {a[i]:.4f}\t{b[j]:.4f}\t{c[k]:.4f}\tOK!\t{counter/(a_i * b_j * c_k)}')
                        if print_to_file_flag:
                            file.write(f'{a[i]:.10e}\t{b[j]:.10e}\t{c[k]:.10e}\n')
                        counter += 1.
                    else:
                        if print_to_screen_flag:
                            print(f'[{counter_2/(a_i * b_j * c_k):.2%}] {a[i]:.4f}\t{b[j]:.4f}\t{c[k]:.4f}\tNO!\t{counter/(a_i * b_j * c_k)}')
         
                
                
    tmat = tf.convert_to_tensor(tmat)
    if print_to_file_flag:
        file.close()
    if print_to_screen_flag:
        print(f'Acceptance rate: {counter/(a_i * b_j * c_k)}')
    
    
    
    return tmat

def aquire_truth_matrix(range_of_a, step):
    path = f'space_R{range_of_a}_step{step}.dat'
    t = np.loadtxt(path)
    t = tf.convert_to_tensor(t)
    t = tf.cast(t, tf.float32)
    return t
    
def bubbleSort(arr, arr2):
    n = len(arr)
 
    # Traverse through all array elements
    for i in range(n):
 
        # Last i elements are already in place
        for j in range(0, n-i-1):
 
            # traverse the array from 0 to n-i-1
            # Swap if the element found is greater
            # than the next element
            if arr[j] > arr[j+1]:
                arr[j], arr[j+1] = arr[j+1], arr[j]
                arr2[j], arr2[j+1] = arr2[j+1], arr2[j]
                
    return arr, arr2
                
def InitializeParameters(range_of_a, step, c_flag = True):
    a = tf.range(step, range_of_a, step)
    a = tf.cast(a, tf.float32)
    #print(a)
    b = tf.range(0., range_of_a/2., step)
    b = tf.cast(b, tf.float32)
    if c_flag:
        c = tf.range(0., range_of_a/4., step)    
        c = tf.cast(c, tf.float32)
        return a, b, c
    else:
        return a, b

def DeltaF(overlap, a, b, c, z):
    
    bin_size = BIN_SIZE()
    
    s1 = tf.sqrt(1 - (overlap + bin_size/2)**2)
    s2 = tf.sqrt(1 - (overlap - bin_size/2)**2)
    
    return a * bin_size * overlap + c * z * (s1 - s2)


def getRandomFirstGuess(overlap_tensor, gaussian_tensor, range_a = 2.):
    check = False

    while check == False:
        a = random.uniform(-range_a, range_a)
        b = random.uniform(-range_a/2, range_a/2)
        c = random.uniform(-range_a/4, range_a/4)
        auxiliary = q_aux(overlap_tensor=overlap_tensor,
                        gaussian_tensor=gaussian_tensor,
                        a = a,
                        b = b, 
                        c = c)
        check = check_constraint(auxiliary)
    print(f'First guess: [{float(a):.4e}, {float(b):.4e}, {float(c):.4e}]')

    return a, b, c, check

def kldGradient(a, b, c, overlap_tensor, gaussian_tensor, ifo, derive_step = 1e-2):
    h = derive_step
    a_mod = a + h
    b_mod = b + h
    c_mod = c + h

    auxiliary = q_aux(overlap_tensor=overlap_tensor,
                      a = a,
                      b = b,
                      c = c,
                      gaussian_tensor=gaussian_tensor)
    
    auxiliary_a = q_aux(overlap_tensor=overlap_tensor,
                      a = a_mod,
                      b = b,
                      c = c,
                      gaussian_tensor=gaussian_tensor)

    auxiliary_b = q_aux(overlap_tensor=overlap_tensor,
                      a = a,
                      b = b_mod,
                      c = c,
                      gaussian_tensor=gaussian_tensor)

    auxiliary_c = q_aux(overlap_tensor=overlap_tensor,
                      a = a,
                      b = b,
                      c = c_mod,
                      gaussian_tensor=gaussian_tensor)

    check = check_constraint(auxiliary)
    dist = p(auxiliary)
    dist_a = p(auxiliary_a)
    dist_b = p(auxiliary_b)
    dist_c = p(auxiliary_c)
    k = kld(dist, ifo)
    k_a = (kld(dist_a, ifo) - k)/h
    k_b = (kld(dist_b, ifo) - k)/h
    k_c = (kld(dist_c, ifo) - k)/h

    return k, k_a, k_b, k_c, check


def minimize_ave(parisi_binned,
             ifo_binned,
             temperature,
             temperature_index,
             range_of_a = 2.,
             step_for_truthmatrix = 5e-2,
             randomize_z = True,
             mean_of_z = 0.,
             std_of_z = 1.,
             print_truthmatrix_to_screen_flag = True,
             print_to_file_space = False,
             print_to_file_aux = False,
             print_to_file_kld = True,
             method = 'Brutal',
             learning_rate = 1e-3,
             derive_step=1e-2,
             max_iterations=1e4,
             tolerance = 1e-3):
    
    x_of_histogram = tf.range(-1. + BIN_SIZE()/2, 1., BIN_SIZE())
    
    #transforming parisi histogram to a Heaviside composition function
    
    overlap_tensor = buildDistribution(parisi_binned=parisi_binned,
                                        x_of_histogram=x_of_histogram,
                                        factor_random_variable=1e2)
    
    if randomize_z:
        gaussian_tensor = tf.random.normal(shape = tf.shape(overlap_tensor), mean = mean_of_z, stddev=std_of_z)
    else:
        gaussian_tensor = tf.zeros(shape = (tf.shape(overlap_tensor)))
        
    if method == 'Brutal':
        print('Brutal method chosen')
        #in this method we compute the KL-Divergence for all points a,b,c such that q_aux(overlap_tensor, gaussian_tensor, a, b, c) \in [-1, 1]
        #and then we take the minimum
        
        parameters = truthmatrix(overlap=overlap_tensor,                                # tutte le possibili combinazioni delle a,b,c
                                 gaussian_numbers=gaussian_tensor,                      # vengono controllate. Prendiamo poi il minimo tra queste 
                                 range_of_a=range_of_a,                                 #
                                 step = step_for_truthmatrix,                           #
                                 print_to_file_flag=print_to_file_space,
                                 print_to_screen_flag=print_truthmatrix_to_screen_flag,
                                 random_flag=randomize_z)                
         
        k_min = 1000.
        a_min = 0.
        b_min = 0.
        c_min = 0.
        
        if print_to_file_aux:
            file_aux = open(f'aux_dist_attempts_temp{temperature_index}.dat', "w")
            file_aux.write(f'#Attempts for T = {temperature}\n')
            file_aux.write(f'overlap\tdistribution\tattempt\tkld\n')
        if print_to_file_kld:
            file_kld = open(f'kld_temp{temperature_index}.dat', "w")
            file_kld.write(f'#Kullback-Leibler in function of a, b, c\n')

        for i in range(len(parameters)):
            a, b, c                 = parameters[i]                                     # per comodità assegno i parametri ad altre tre variabili
            auxiliary_overlap       = q_aux(overlap_tensor=overlap_tensor,              # Dato l'array degli overlap espanso calcolo f(q| a,b,c)
                                        a = a,
                                        b = b,
                                        c = c,
                                        gaussian_tensor=gaussian_tensor)
            
            auxiliary_distribution  = p(auxiliary_overlap)
            norm = tf.reduce_sum(auxiliary_distribution*BIN_SIZE())
            auxiliary_distribution /= norm

            #print(f'{norm}\t{norm1}')

            k                       = kld(dist_to_test=auxiliary_distribution,          # Calcolo della KL-Divergence tra la distribuzione
                                          analytical_dist=ifo_binned)                 # degli overlap ausialiari e la ifo teorica
            if print_to_file_kld:
                file_kld.write(f'{float(a)}\t{float(b)}\t{float(c)}\t{k}\n')                             
            if print_to_file_aux:
                for j in range(len(x_of_histogram)):
                    file_aux.write(f'{x_of_histogram[j]}\t{auxiliary_distribution[j]}\t{i}\t{float(a)}\t{float(b)}\t{float(c)}\t{k}\n')
                file_aux.write("\n\n")
            if k <= k_min:                                                              # Trovo il minimo 
                k_min = k
                a_min = a
                b_min = b
                c_min = c
                aux_min = auxiliary_distribution
                print(f'[{float(i+1)/len(parameters):.2%}] KL:{k:.4e}*\t\t[{float(a):.4e}, {float(b):.4e}, {float(c):.4e}]\t\tK_min:\t{float(k_min):.4e}\t\t[{float(a_min):.4e}, {float(b_min):.4e}, {float(c_min):.4e}]')              
            else:
                print(f'[{float(i+1)/len(parameters):.2%}] KL:{k:.4e} \t\t[{float(a):.4e}, {float(b):.4e}, {float(c):.4e}]\t\tK_min:\t{float(k_min):.4e}\t\t[{float(a_min):.4e}, {float(b_min):.4e}, {float(c_min):.4e}]')
        if print_to_file_aux:
            file_aux.close()
        if print_to_file_kld:
            file_kld.close()
    elif method == 'GradientDescent':
        a, b, c, check = getRandomFirstGuess(overlap_tensor=overlap_tensor,
                                                            gaussian_tensor=gaussian_tensor,
                                                            range_a=range_of_a)
        a, b, c = tf.convert_to_tensor([a, b, c])
        lr = learning_rate
        for i in range(int(max_iterations)):
            while not check:
                kl, ka, kb, kc, check = kldGradient(a = a,
                               b = b,
                               c = c, 
                               overlap_tensor=overlap_tensor,
                               gaussian_tensor=gaussian_tensor,
                               ifo = ifo_binned,
                               derive_step=derive_step)
                grad = tf.convert_to_tensor([ka, kb, kc])
                diff = lr * grad
            
            print(f'[{i}] KLD:{float(kl):.3e}\t[{float(a):.2e} {float(b):.2e} {float(c):.2e}]\t[{float(ka):.4e} {float(kb):.4e} {float(kc):.4e}]\t{check}')
            a_old, b_old, c_old = a, b, c
            if tf.math.abs(tf.reduce_max(tf.math.abs(grad)))<tolerance:
                print('Critical Point reached!')
                break
            a, b, c = a-diff[0], b-diff[1], c-diff[2]

            if tf.math.reduce_all(tf.equal(a_old,a)) and tf.math.reduce_all(tf.equal(b_old,b)) and tf.math.reduce_all(tf.equal(c_old,c)):
                print('Try to tune parameters...')
                break

        a_min = a
        b_min = b
        c_min = c
        k_min = kl
        auxiliary = q_aux(overlap_tensor=overlap_tensor,
                        gaussian_tensor=gaussian_tensor,
                        a = a_min,
                        b = b_min,
                        c = c_min)

        aux_min = p(auxiliary)
        
    else:
        print(f'Method {method} is not valid. Exit')
        exit(-1)
    return aux_min, k_min, a_min, b_min, c_min
    

def singleSample3paramMinimize(overlap_tensor,
                         ifo_distribution,
                         file_attemps,
                         temperature_index,
                         temperature,
                         number_of_sample,
                         step = 5e-2,
                         range_of_a = 2.,
                         attempts_flag = False):
    
    overlap_tensor = tf.cast(overlap_tensor, tf.float32)

    gaussian_tensor = tf.random.normal(shape = tf.shape(overlap_tensor))
    gaussian_tensor = tf.cast(gaussian_tensor, tf.float32)
    a, b, c = InitializeParameters(range_of_a = range_of_a, 
                                    step = step)

    na = len(a)
    nb = len(b)
    nc = len(c)

    a_min = 0.
    b_min = 0.
    c_min = 0.
    k_min = 1000.
    count_attempts = 0 
    count = 0.
    tot = na*nb*nc
    file_space = open(f'space_step{step}_temp{temperature_index}_sample{number_of_sample+1}.dat', "w")
    file_space.write(f'#space for T = {temperature}\n')
    for i in range(na):
        for j in range(nb):
            for k in range(nc):
                count += 1.
                
                aux_overlap = q_aux(overlap_tensor=overlap_tensor,
                                    a = a[i],
                                    b = b[j],
                                    c = c[k],
                                    gaussian_tensor=gaussian_tensor)
                if check_constraint(aux_overlap):
                    aux_dist = p(aux_overlap)
                    K = kld(dist_to_test=aux_dist,
                            analytical_dist= ifo_distribution)
                    if attempts_flag:
                        overlaps.printAuxDist(file=file_attemps,
                                            aux_dist=aux_dist,
                                            temperature=k,
                                            a = a[i],
                                            b = b[j],
                                            c = c[k],
                                            k = temperature_index,
                                            attempts_flag=True,
                                            attempts = count_attempts,
                                            kld = K)
                    count_attempts += 1
                    
                    file_space.write(f'{a[i]:.4e}\t{b[j]:.4e}\t{c[k]:.4e}\t{K:.4e}\n')

                    if K <= k_min:
                        k_min = K
                        a_min = a[i]
                        b_min = b[j]
                        c_min = c[k]
                        aux_min = aux_dist
                        print(f'[{count/tot:.2%}]KLD: {K:.4e}*[{a[i]:.2e} {b[j]:.2e} {c[k]:.2e}]\t\tmin KLD: {k_min:.4e} for [{a_min:.2e} {b_min:.2e} {c_min:.2e}] PT:{temperature_index} Ns: {number_of_sample}')
                    else:
                        print(f'[{count/tot:.2%}]KLD: {K:.4e} [{a[i]:.2e} {b[j]:.2e} {c[k]:.2e}]\t\tmin KLD: {k_min:.4e} for [{a_min:.2e} {b_min:.2e} {c_min:.2e}] PT:{temperature_index} Ns: {number_of_sample}')
                
    file_space.close()
    print(f'minimum: {k_min:.4e} for [{a_min:.2e} {b_min:.2e} {c_min:.2e}]')

  

    return aux_min, k_min, a_min, b_min, c_min


    
    
    
    

