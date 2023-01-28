import loadingModule
import functionsModule
import sys
import numpy as np
import os
# Getting parameters from bash script ######################
N = np.int32(sys.argv[1])          # Size of the system
t_min = np.float64(sys.argv[2])    # Minimal Temperature
t_max = np.float64(sys.argv[3])    # Maximal Temperature
# ------------------------------------######################
size_path = f'N{N}'

devices = loadingModule.GetDirectories(path = size_path)

for dev in devices:
    dev_path = size_path + f'/{dev}'
    simulations = loadingModule.GetDirectories(path = dev_path)

    for simulation in simulations:
        
        simulation_path = dev_path + f'/{simulation}'
        
        samples = loadingModule.GetDirectories(path = simulation_path)
        number_of_sample = functionsModule.GetSampleIndex(samples)
        os.system(f'mkdir -p {simulation_path}/data')
        functionsModule.checkSamples(samples=samples, path=simulation_path)
        options = loadingModule.Settings(simulation_path, t_min=t_min, t_max=t_max)
        
        full_samples_directories = [simulation_path + f'/{d}' for d in samples]
        count = 0
        for sample_path in full_samples_directories:
            analysis = functionsModule.Analysis(path=sample_path, param=options.get_all(), sample = number_of_sample[count])        
            analysis.print_path()
            analysis.LoadWholeSampleFiles()
            print(f'Loaded {sample_path}')
            analysis.GetFrequenciesFile()
            analysis.GetParallelTemperingFile()
            analysis.GetTemperatures()
            analysis.DumpEnergy()
            analysis.ComputeSpecificHeat()
            analysis.MakeSpectrum(mean_flag=True)
            analysis.DumpSpectrum(print_instant=True)
            analysis.ComputeParisiOverlaps()
            analysis.ComputeTheoIFO()
            analysis.PrintDistributions(bins = 80)


            count += 1
        
        os.system(f'mv *.dat {simulation_path}/data')

            
        


        





