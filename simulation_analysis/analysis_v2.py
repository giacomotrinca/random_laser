import loadingModule
import functionsModule
import sys
import numpy as np
import os
from tqdm import tqdm
# Getting parameters from bash script ######################
N = np.int32(sys.argv[1])          # Size of the system
t_min = np.float64(sys.argv[2])    # Minimal Temperature
t_max = np.float64(sys.argv[3])    # Maximal Temperature
# ------------------------------------######################
size_path = f'N{N}'
only_disorder = 1

if only_disorder == 0:
    devices = loadingModule.GetDirectories(path = size_path)
    directories = []
    sim_directories = []
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
            directories.append(full_samples_directories)
            count = 0
            for sample_path in tqdm(full_samples_directories, f'Processing {simulation_path}'):
                analysis = functionsModule.Analysis(path=sample_path, param=options.get_all(), sample = number_of_sample[count])        
                analysis.LoadWholeSampleFiles()
                #print(f'Loaded {sample_path}')
                analysis.GetFrequenciesFile()
                analysis.GetParallelTemperingFile()
                analysis.GetTemperatures()
                analysis.DumpEnergy()
                analysis.ComputeSpecificHeat()
                analysis.MakeSpectrum(mean_flag=True)
                analysis.DumpSpectrum(print_instant=True)
                analysis.ComputeParisiOverlaps()
                analysis.ComputeTheoIFO()
                analysis.ComputeExpIFO()
                analysis.PrintDistributions(bins = 100)
                analysis.PrintOverlap()
                count += 1
        
            os.system(f'mv *.dat {simulation_path}/data')


os.system(f'mkdir -p {size_path}/dis_ave')
disorder = functionsModule.DisorderAverage(loadingModule.search_for_data_folder(size_path), t_min=t_min, t_max=t_max)
disorder.SpecificHeat()

os.system(f'mv *.dat {size_path}/dis_ave/')





        





