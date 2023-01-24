import numpy as np
import os
import sys
import loadingModule 
import concurrent.futures as multithreading



class Analysis:
    def __init__(self, param = None, path = None):
        if param and path:
            self.parameters = param
            self.path=path
            self.size = param[0]
            self.replicas = param[1]
            self.pt_rate = param[2]
            self.iters = param[3]
            self.first = param[4]
            self.npt = param[5]
            self.pt_flag = param[6]
            self.print_rate = param[7]
        else:
            print(f'You have to initialize the analyzer!')
            sys.exit(-1)
    
    def print_parameters(self):
        [print(p) for p in self.parameters]
    
    def print_path(self):
        [print(p) for p in self.path]
    
    def get_parameters(self):
        return [p for p in self.parameters]
    
    def get_path(self):
        return [p for p in self.path]
    
    def LoadWholeSample(self):

        configurations = []
        iters = [i for i in range(self.first, self.iters, self.print_rate)]
        
        for replica in range(0, self.replicas):
            configuration_times = []
            for i in iters:
                configuration_path = self.path + f'/config_nrep{replica}_iter_{i}.dat'
                configuration_times.append(loadingModule.GetConfig(configuration_path))
            configurations.append(configuration_times)
        configurations=np.array(configurations, dtype=np.float64)
        configurations = np.reshape(configurations, newshape=(self.replicas, len(iters), self.npt, self.size, 2))

        self.configurations = configurations
    def LoadFrequencies(self):
        self.frequencies = np.loadtxt(f'{self.path}/frequencies.dat', dtype=np.float64)[:, 1]

    def print_frequencies(self):
        print(self.frequencies)
        

    def print_config(self):
        print(self.configurations)
        

       

    
    

    
    

def checkSamples(samples, path):
    if len(samples) > 1 :
        print(f'Found {len(samples)} samples in {path}')
    elif len(samples) == 1:
        print(f'Found {len(samples)} sample in {path}')
    else :
        print(f'Warning: No samples found.', file=sys.stderr)
        sys.exit(-1)
        
    


        
