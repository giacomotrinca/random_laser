import numpy as np
import os
import sys
import loadingModule 
import concurrent.futures as multithreading



class Analysis:
    def __init__(self, param = None, path = None):
        # Paramaters Array
        # 0 -> Size
        # 1 -> Replicas
        # 2 -> pt_rate
        # 3 -> iter
        # 4 -> first
        # 5 -> npt
        # 6 -> pt_flag
        # 7 -> print_rate
        
        if param and path:
            self.parameters = param
            self.path=path
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
        iters = [i for i in range(self.parameters[4], self.parameters[3], self.parameters[7])]
        
        for replica in range(0, self.parameters[1]):
            configuration_times = []
            for i in iters:
                configuration_path = self.path + f'/config_nrep{replica}_iter_{i}.dat'
                configuration_times.append(loadingModule.GetConfig(configuration_path))
            configurations.append(configuration_times)
        configurations=np.array(configurations, dtype=np.float64)
        configurations = np.reshape(configurations, newshape=(self.parameters[1], len(iters), self.parameters[5], self.parameters[0], 2))

        self.configurations = configurations
    def LoadFrequencies(self):
        self.frequencies = np.loadtxt(f'{self.path}/frequencies.dat', dtype=np.float64)[:, 1]

    def print_frequencies(self):
        #if self.frequencies:
        print(self.frequencies)
        

    def print_config(self):
        #if self.configurations:
        print(self.configurations)
        

       

    
    

    
    

def checkSamples(samples, path):
    if len(samples) > 1 :
        print(f'Found {len(samples)} samples in {path}')
    elif len(samples) == 1:
        print(f'Found {len(samples)} sample in {path}')
    else :
        print(f'Warning: No samples found.', file=sys.stderr)
        sys.exit(-1)
        
    


        
