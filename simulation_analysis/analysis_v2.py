import loadingModule
import functionsModule
import sys
import numpy as np

#getting parameters from bash script######################
N = np.int32(sys.argv[1])           #Size of the system  #
t_min = np.float64(sys.argv[2])     #Minimal Temperature #
t_max = np.float64(sys.argv[3])     #Maximal Temperature #
#devices = np.int32(sys.argv[4])     #Number of devices   #
#-----------------------------------######################

devices = loadingModule.list_directories(f'N{N}')
for dev in devices:
    simulations = loadingModule.list_directories(f'N{N}/{dev}') 
    for sim in simulations:     
        samples = loadingModule.list_directories(f'N{N}/{dev}/{sim}')
        for sample in samples:
            print(f'Directory ./N{N}/{dev}/{sim}/{sample} found!')

