import numpy as np
import os
import sys
import loadingModule 
import concurrent.futures as multithreading



class Analysis:
    def __init__(self, param = None, path = None, sample = None):
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
            self.tmin = param[8]
            self.tmax = param[9]
            self.iterations = np.array([k for k in range(self.first, self.iters, self.print_rate)], dtype=np.int32)
            self.sample = sample
            
        else:
            print(f'You have to initialize the analyzer!')
            sys.exit(-1)
    
    def print_parameters(self):
        [print(p) for p in self.parameters]
    
    def GetTemperatures(self):
        step = (self.tmax - self.tmin)/self.npt
        self.temperatures = np.array([self.tmin + (i+1)* step for i in range(self.npt)], dtype = np.float64)
    
    def print_temperatures(self):
        print(self.temperatures)
    
    def print_path(self):
        print(self.path)
    
    def get_parameters(self):
        return np.array([p for p in self.parameters], dtype=np.float64)
    
    def get_path(self):
        return self.path
    
    def LoadWholeSampleFiles(self):

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
    def GetFrequenciesFile(self):
        self.frequencies = np.loadtxt(f'{self.path}/frequencies.dat', dtype=np.float64)[:, 1]

    def print_frequencies(self):
        print(self.frequencies)
        

    def print_config(self):
        print(self.configurations)

    def GetParallelTemperingFile(self):
        cols = np.arange(3, 3 * self.npt, 3)
        energy = []
        for replica in range(self.replicas):
            parallel_tempering_path = self.path + f'/parallel_tempering{replica}.dat'
            energy.append(np.loadtxt(parallel_tempering_path, usecols=cols))
        energy = np.array(energy, dtype=np.float64)
        self.energy = np.einsum("ijk->ikj", energy)
        #print(np.shape(self.energy))

    def print_energy(self):
        print(np.shape(self.energy))

    def find_PT(self):
        if "PT" in self.path:
            return True
        else:
            return False
    
    def MakeSpectrum(self, mean_flag = True):
        
        self.spectrum = np.einsum("ijkl -> lijk", self.configurations[:, :, :, :, 0]**2 + self.configurations[:, :, :, :, 1]**2)
        norm = np.sum(self.spectrum, axis=0)
        self.spectrum = self.spectrum/norm[np.newaxis, :,:,:]
        self.spectrum = self.spectrum/np.sqrt(self.temperatures[np.newaxis, np.newaxis, np.newaxis, :])
        #print(np.shape(self.spectrum))
        
        if self.find_PT():
            if mean_flag:
                self.mean_spectrum = np.mean(self.spectrum, axis = 2)
                self.std_spectrum = np.std(self.spectrum, axis=2)
                #print(np.shape(self.std_spectrum))
        else:
            if mean_flag:
                self.mean_spectrum=[]
                self.std_spectrum = []
                index = 0
                lenght=1
                for i in range(int(np.log2(len(self.spectrum[0, 0, :, 0])))):
                    self.mean_spectrum.append(np.mean(self.spectrum[:, :, index:index+lenght, :], axis=2))
                    self.std_spectrum.append(np.std(self.spectrum[:, :, index:index+lenght, :], axis=2))
                    index += lenght
                    lenght *= 2
                self.mean_spectrum = np.array(self.mean_spectrum, dtype=np.float64)
                self.std_spectrum = np.array(self.std_spectrum, dtype=np.float64)

    
    def DumpSpectrum(self, print_instant = False):
        
        if self.find_PT():
            self.print_mean_spectrum = np.mean(self.mean_spectrum, axis = 1)
            self.print_std_spectrum = np.mean(self.std_spectrum, axis = 1)
            self.print_mean_spectrum = np.einsum("ij -> ji", self.print_mean_spectrum)
            self.print_std_spectrum = np.einsum("ij -> ji", self.print_std_spectrum)
            #print(np.shape(self.print_mean_spectrum))
            filename = f'mean_spectrum_PT_size{self.size}_sample{self.sample}.dat'
            file_handle = open(filename, "w")
            #file_handle.close()

            for k in range(self.npt):
                np.savetxt(file_handle, np.c_[self.frequencies, self.print_mean_spectrum[k, :], self.print_std_spectrum[k, :], np.full(np.shape(self.print_mean_spectrum[k, :]), self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                file_handle.write("\n\n")
            file_handle.close()
            del self.print_mean_spectrum
            del self.print_std_spectrum
        
        else:
            for b in range(len(self.mean_spectrum)):
                self.print_mean_spectrum = np.mean(self.mean_spectrum[b], axis = 1)
                self.print_std_spectrum = np.mean(self.std_spectrum[b], axis = 1)
                self.print_mean_spectrum = np.einsum("ij -> ji", self.print_mean_spectrum)
                self.print_std_spectrum = np.einsum("ij -> ji", self.print_std_spectrum)
                #print(np.shape(self.print_mean_spectrum))
                filename = f'mean_spectrum_block{b}_size{self.size}_sample{self.sample}.dat'
                file_handle = open(filename, "w")
                #file_handle.close()

                for k in range(self.npt):
                    np.savetxt(file_handle, np.c_[self.frequencies, self.print_mean_spectrum[k, :], self.print_std_spectrum[k, :], np.full(np.shape(self.print_mean_spectrum[k, :]), self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                    file_handle.write("\n\n")
                file_handle.close()
                del self.print_mean_spectrum
                del self.print_std_spectrum
            
            if print_instant:
                print_spectrum = np.mean(self.spectrum, axis = 1)
                print_std = np.std(self.spectrum, axis=1)
                print(np.shape(print_spectrum))
            
                for i in range(len(print_spectrum[0])):
                    filename = f'spectrum_time{i}_size{self.size}_sample{self.sample}.dat'
                    file_handle = open(filename, "w")

                    for k in range(self.npt):
                        np.savetxt(file_handle, np.c_[self.frequencies, print_spectrum[:, i, k], print_std[:, i, k], np.full(np.shape(print_spectrum[:, i, k]), self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                        file_handle.write("\n\n")
                    file_handle.close()
                del print_spectrum
                del print_std
        
        
                
            





            
        



                
            
def GetSampleIndex(samples):
    k_list = []
    for string in samples:
        k = string.split("sample")[1]
        k_list.append(int(k))
    
    return np.array(k_list, dtype=np.int32)

       

       

    
    

    
    

def checkSamples(samples, path):
    if len(samples) > 1 :
        print(f'Found {len(samples)} samples in {path}')
    elif len(samples) == 1:
        print(f'Found {len(samples)} sample in {path}')
    else :
        print(f'Warning: No samples found.', file=sys.stderr)
        sys.exit(-1)
        
    


        
