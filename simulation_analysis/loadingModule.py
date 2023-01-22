import os
import numpy as np

class Settings:

    def __init__(self, file_path= None):
        if file_path:
            self.load_settings(file_path)
        else:
            print(f'No file path. Exit')
            return -1

    def load_settings(self, file_path):
        self.size, self.replicas, self.pt_rate, self.iter, self.first, self.npt, self.pt,  self.print_rate = np.loadtxt(f'{file_path}/analysis_input.info', dtype=np.int32)
        
    def print_settings(self):
        print("Options")
        print(f' SIZE: {self.size}\n NREP: {self.replicas}\n PT_rate: {self.pt_rate}\n NITER: {self.iter}\n first: {self.first}\n NPT: {self.npt} \n pt_flag: {self.pt}\n print_every: {self.print_rate}\n')
    
    def get_size(self):
        return self.size
    
    def get_replicas(self):
        return self.replicas

    def get_pt_rate(self):
        return self.pt_rate

    def get_iter(self):
        return self.iter

    def get_first(self):
        return self.first
    
    def get_npt(self):
        return self.npt
    
    def get_pt(self):
        return self.pt
    
    def get_print_rate(self):
        return self.print_rate

    



def list_directories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

def find_PT(string):
    if "PT" in string:
        return True
    else:
        return False


