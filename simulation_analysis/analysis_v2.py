import loadingModule
import functionsModule
import sys
import numpy as np

#getting parameters from bash script######################
N = np.int32(sys.argv[1])           #Size of the system  #
t_min = np.float64(sys.argv[2])     #Minimal Temperature #
t_max = np.float64(sys.argv[3])     #Maximal Temperature #
#devices = np.int32(sys.argv[4])    #Number of devices   #
#-----------------------------------######################

devices = loadingModule.list_directories(f'N{N}')
for dev in devices:
    simulations = loadingModule.list_directories(f'N{N}/{dev}') 
    for sim in simulations:     
        samples = loadingModule.list_directories(f'N{N}/{dev}/{sim}')
        options = loadingModule.Settings(f'N{N}/{dev}/{sim}')
        options.print_settings()
        print(options.get_first())
        if len(samples) > 1 :
            print(f'Found {len(samples)} samples')
        elif len(samples) == 1:
            print(f'Found {len(samples)} sample')
        else :
            print(f'Warning: No samples found.', file=sys.stderr)
        for sample in samples:
            path = f'N{N}'

        

