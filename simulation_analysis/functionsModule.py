import numpy
import os
import sys
import loadingModule

class Analysis:
    def __init__(self, param = None, paths = None):
        if param and paths:
            self.parameters = param
            self.path=paths
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
    

    
    

def checkSamples(samples, path):
    if len(samples) > 1 :
        print(f'Found {len(samples)} samples in {path}')
    elif len(samples) == 1:
        print(f'Found {len(samples)} sample in {path}')
    else :
        print(f'Warning: No samples found.', file=sys.stderr)
        sys.exit(-1)
        
    


        
