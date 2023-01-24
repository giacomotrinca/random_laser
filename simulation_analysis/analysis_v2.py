import loadingModule
import functionsModule
import sys
import numpy as np

# Getting parameters from bash script ######################
N = np.int32(sys.argv[1])          # Size of the system
t_min = np.float64(sys.argv[2])    # Minimal Temperature
t_max = np.float64(sys.argv[3])    # Maximal Temperature
nthreads = np.int32(sys.argv[4])   # Number of threads
# ------------------------------------######################
size_path = f'N{N}'

devices = loadingModule.GetDirectories(path = size_path)

for dev in devices:
    dev_path = size_path + f'/{dev}'
    simulations = loadingModule.GetDirectories(path = dev_path)

    for simulation in simulations:
        simulation_path = dev_path + f'/{simulation}'
        
        samples = loadingModule.GetDirectories(path = simulation_path)
        functionsModule.checkSamples(samples=samples, path=simulation_path)
        options = loadingModule.Settings(simulation_path)
        
        full_samples_directories = [simulation_path + f'/{d}' for d in samples]
        for sample_path in full_samples_directories:
            analysis = functionsModule.Analysis(path=sample_path, param=options.get_all())        
            analysis.LoadWholeSample()
            print(f'Loaded {sample_path}')
            analysis.LoadFrequencies()
            #analysis.print_frequencies()
            #analysis.print_config()
        


        





