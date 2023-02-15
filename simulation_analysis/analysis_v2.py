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
bins = np.int32(sys.argv[4])
only_disorder = np.int32(sys.argv[5])
# ------------------------------------######################
size_path = f'N{N}'


if only_disorder == 0:
    devices = loadingModule.GetDirectories(path=size_path)
    directories = []
    sim_directories = []
    for dev in devices:
        dev_path = size_path + f'/{dev}'
        simulations = loadingModule.GetDirectories(path=dev_path)

        for simulation in simulations:

            simulation_path = dev_path + f'/{simulation}'

            samples = loadingModule.GetDirectories(path=simulation_path)
            number_of_sample = functionsModule.GetSampleIndex(samples)
            os.system(f'mkdir -p {simulation_path}/data')
            functionsModule.checkSamples(samples=samples, path=simulation_path)
            options = loadingModule.Settings(simulation_path,
                                             t_min=t_min,
                                             t_max=t_max)

            full_samples_directories = [simulation_path + f'/{d}'
                                        for d in samples]
            directories.append(full_samples_directories)
            c = 0
            for sample_path in tqdm(full_samples_directories,
                                    f'Processing {simulation_path}'):
                analysis = functionsModule.Analysis(path=sample_path,
                                                    param=options.get_all(),
                                                    sample=number_of_sample[c])
                analysis.LoadWholeSampleFiles()
                # print(f'Loaded {sample_path}')
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
                analysis.PrintDistributions(bins=bins)
                analysis.PrintOverlap()
                c += 1

            os.system(f'mv *.dat {simulation_path}/data')
else:
    print("Disorder average only mode.")

os.system(f'mkdir -p {size_path}/dis_ave')
disorder = functionsModule.DisorderAverage(loadingModule.search(size_path),
                                           t_min=t_min,
                                           t_max=t_max,
                                           bins=bins)
disorder.SpecificHeat()
disorder.Spectrum()
disorder.Distribution(type='parisi_dist')
disorder.Distribution(type='theo_ifo_dist')
disorder.Distribution(type='exp_ifo_dist')
os.system(f'mv *.dat {size_path}/dis_ave/')
