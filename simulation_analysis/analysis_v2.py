import loadingModule
import functionsModule
import sys
import numpy as np

# Getting parameters from bash script ######################
N = np.int32(sys.argv[1])             #Size of the system  #
t_min = np.float64(sys.argv[2])       #Minimal Temperature #
t_max = np.float64(sys.argv[3])       #Maximal Temperature #
nthreads=np.int32(sys.argv[4])        #Number of threads   #
#devices = np.int32(sys.argv[4])      #Number of devices   #
#-------------------------------------######################

devices = loadingModule.list_directories(f'N{N}')
for dev in devices:
    simulations = loadingModule.list_directories(f'N{N}/{dev}') 

    for sim in simulations:
        path =  f'N{N}/{dev}/{sim}'    
        samples = loadingModule.list_directories(path)
        options = loadingModule.Settings(path)
        options.print_settings()
        

        if len(samples) > 1 :
            print(f'Found {len(samples)} samples')
        elif len(samples) == 1:
            print(f'Found {len(samples)} sample')
        else :
            print(f'Warning: No samples found.', file=sys.stderr)
            sys.exit(-1)
        

        

        

