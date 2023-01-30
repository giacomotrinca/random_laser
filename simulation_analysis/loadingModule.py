import os
import numpy as np
import sys


class Settings:

    def __init__(self, file_path= None, t_min = None, t_max = None):
        if file_path and t_min and t_max:
            self.load_settings(file_path)
            self.tmin = t_min
            self.tmax = t_max
        else:
            print(f'No file path. Exit', file=sys.stderr)
            sys.exit(1)

    def load_settings(self, file_path):
        self.size, self.replicas, self.pt_rate, self.iter, self.first, self.npt, self.pt,  self.print_rate = np.loadtxt(f'{file_path}/analysis_input.info', dtype=np.int32)
        
    def print_settings(self):
        print("Options")
        print(f' size: {self.size}\n replicas: {self.replicas}\n pt_rate: {self.pt_rate}\n iter: {self.iter}\n first: {self.first}\n npt: {self.npt} \n pt_flag: {self.pt}\n print_every: {self.print_rate}\n')
    
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
    
    def get_all(self):
        return self.size, self.replicas, self.pt_rate, self.iter, self.first, self.npt, self.pt, self.print_rate, self.tmin, self.tmax 


def GetConfig(path):
    return np.loadtxt(path, usecols=[2, 3], dtype=np.float64)


def GetDirectories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and d != "data"]

def search_for_data_folder(path):
    found = False
    directories = []
    for root, dirs, files in os.walk(path):
        if "data" in dirs:
            directories.append(os.path.join(root, 'data'))
            found = True
    if not found:
        pass
    return directories

 