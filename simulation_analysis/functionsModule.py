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
    
    def DumpEnergy(self):
        index = 0
        length = 1
        mean_energy = []
        std_energy = []
        times = []
        count = 0
        for i in range(int(np.log2(len(self.energy[0, 0])))):
            mean_energy.append(np.mean(self.energy[:, :, index:index+length], axis = 2))
            std_energy.append(np.std(self.energy[:, :, index:index+length], axis = 2)/np.sqrt(length))
            times.append((index+length) * self.pt_rate)
            index += length
            length *= 2
            count += 1
            if count == int(np.log2(len(self.energy[0, 0]))) -1 :
                self.raw_energy = self.size * self.energy[:, :, index:index+length]
                self.mean_energy = self.size * np.mean(self.energy[:, :, index:index+length], axis = 2)


        mean_energy = np.array(mean_energy, dtype=np.float64)
        std_energy = np.array(std_energy, dtype=np.float64)
        mean_energy = np.mean(mean_energy, axis = 1)
        times = np.array(times, dtype=np.int32)
        #print(times)
        std_energy = np.mean(std_energy, axis = 1)
        #print(np.shape(mean_energy))

        filename = f'mean_energy_size{self.size}_sample{self.sample}.dat'
        file_handle = open(filename, "w")
        for k in range(self.npt-1):
            np.savetxt(file_handle, np.c_[times, self.size * mean_energy[:, k], self.size* std_energy[:, k], np.full(np.shape(mean_energy[:, k]), self.temperatures[self.npt -1 - k])])
            file_handle.write("\n\n")
        file_handle.close()
        del mean_energy
        del std_energy
        del times
        del filename
        
        
    
    def MakeSpectrum(self, mean_flag = True):
        
        self.spectrum = np.einsum("ijkl -> lijk", self.configurations[:, :, :, :, 0]**2 + self.configurations[:, :, :, :, 1]**2)
        norm = np.sum(self.spectrum, axis=0)
        self.spectrum = self.spectrum/norm[np.newaxis, :,:,:]
        self.spectrum = self.spectrum/np.sqrt(self.temperatures[np.newaxis, np.newaxis, np.newaxis, :])
        #print(np.shape(self.spectrum))
        
        if self.find_PT():
            if mean_flag:
                self.mean_spectrum = np.mean(self.spectrum, axis = 2)
                self.std_spectrum = np.std(self.spectrum, axis=2)/np.sqrt(len(self.spectrum[0, 0]))
                #print(np.shape(self.std_spectrum))
        else:
            if mean_flag:
                self.mean_spectrum=[]
                self.std_spectrum = []
                index = 0
                length=1
                for i in range(int(np.log2(len(self.spectrum[0, 0, :, 0])))):
                    self.mean_spectrum.append(np.mean(self.spectrum[:, :, index:index+length, :], axis=2))
                    self.std_spectrum.append(np.std(self.spectrum[:, :, index:index+length, :], axis=2)/np.sqrt(length))
                    index += length
                    length *= 2
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
                #print(np.shape(print_spectrum))
            
                for i in range(len(print_spectrum[0])):
                    filename = f'spectrum_time{i}_size{self.size}_sample{self.sample}.dat'
                    file_handle = open(filename, "w")

                    for k in range(self.npt):
                        np.savetxt(file_handle, np.c_[self.frequencies, print_spectrum[:, i, k], print_std[:, i, k], np.full(np.shape(print_spectrum[:, i, k]), self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                        file_handle.write("\n\n")
                    file_handle.close()
                del print_spectrum
                del print_std

    def ComputeSpecificHeat(self):
        
        e = self.size * self.energy[:, :, int(self.iters/2):]
        
        e2 = e**2

        mean_e = np.mean(e, axis=2)
        sigma_e = np.std(e, axis=2)/np.sqrt(int(self.iters/2))
        mean_e2 = np.mean(e2, axis=2)
        #print(np.shape(mean_e))
        temperatures = np.delete(self.temperatures, -1)
        mean_e = np.fliplr(mean_e)
        mean_e2 = np.fliplr(mean_e2)
        cv = self.size * (mean_e2 - mean_e **2)/(self.size * temperatures[np.newaxis, :]**2)
        #print(np.shape(cv))
        sigma = np.sqrt(cv**2 * (2*mean_e2*sigma_e**2/(self.size*temperatures[np.newaxis, :]**4) + sigma_e**2/(self.size*temperatures[np.newaxis, :]**2)))
        sigma /= self.size
        #print(np.shape(sigma))
        if self.find_PT():
            filename = f'specific_heat_PT_sample{self.sample}.dat'
        else:
            filename = f'specific_heat_NOPT_sample{self.sample}.dat'

        file_handle = open(filename, "w")
        for r in range(self.replicas):
            np.savetxt(file_handle, np.c_[temperatures, cv[r, :], sigma[r, :], np.full(shape = np.shape(cv[r, :]), fill_value=r)], fmt="%4e", delimiter="\t", newline="\n")
            file_handle.write("\n\n")
        
        file_handle.close()
        del cv 
        del sigma 
        del mean_e
        del sigma_e
        del temperatures
        del e 
        del e2
    
    def ComputeParisiOverlaps(self):
        
        if self.replicas <= 1:
            print('I can not compute overlaps: replicas of the system are insufficients')
            return -1
        
        
        sigma = self.configurations[:, :, :, :, 0]
        tau = self.configurations[:, :, :, :, 1]
        q = (np.einsum('ajkl,bjkl->abjk', sigma, sigma) + np.einsum('ajkl,bjkl->abjk', tau, tau))/self.size
        i, j = np.triu_indices(self.replicas, k=1)
        self.parisi = q[i, j, :, :]
        del q
        del sigma
        del tau
        #print(np.shape(self.parisi))
    
    def ComputeTheoIFO(self):
        if self.find_PT():
            delta = self.spectrum - np.mean(self.spectrum, axis = 2)[:, :, np.newaxis, :]
            c = np.einsum('jaik, jbik -> abik', delta, delta)
            norm = np.einsum('jaik -> aik', delta*delta) 
            norm = np.sqrt(np.einsum('aij, bij -> abij', norm, norm))
            c /= norm
            i, j = np.triu_indices(self.replicas, k=1)
            self.theo_ifo = c[i, j, :, :]
        
        else:
            index = 0
            length = 2
            self.theo_ifo = []
            for b in range(int(np.log2(len(self.configurations[0])))):
                delta = self.spectrum - np.mean(self.spectrum[:, :, index:index+length, :], axis = 2)[:, :, np.newaxis, :]
                c = np.einsum('jaik, jbik -> abik', delta, delta)
                norm = np.einsum('jaik -> aik', delta*delta) 
                norm = np.sqrt(np.einsum('aij, bij -> abij', norm, norm))
                c /= norm
                i, j = np.triu_indices(self.replicas, k=1)
                self.theo_ifo.append(c[i, j, :, :])
                index += length
                length *= 2
            self.theo_ifo = np.array(self.theo_ifo, dtype=np.float64)
            #print(np.shape(self.theo_ifo))
        
            

            
        

    def ComputeExpIFO(self):
        pass

    def PrintDistributions(self, bins):
        self.bins = bins
        bin_size = 2./self.bins
        q = np.linspace(-1+bin_size, 1-bin_size, num=self.bins)
        if self.find_PT():
            #print(np.shape(self.parisi))
            #print(np.shape(self.theo_ifo))
            filename_pq = f'parisi_dist_PT_sample{self.sample}.dat'
            filename_pc = f'theo_ifo_dist_PT_sample{self.sample}.dat'
            file_handle_pc = open(filename_pc, "w")
            file_handle_pq = open(filename_pq, "w")

            for k in range(self.npt):
                pq = np.histogram(a = self.parisi[:, :, k], bins=self.bins, density=True)[0]
                pc = np.histogram(a = self.theo_ifo[:, :, k], bins=self.bins, density=True)[0]
                np.savetxt(file_handle_pq, np.c_[q, pq, np.full(shape=np.shape(pq),fill_value=self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                np.savetxt(file_handle_pc, np.c_[q, pc, np.full(shape=np.shape(pc),fill_value=self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                file_handle_pc.write("\n\n")
                file_handle_pq.write("\n\n")
            file_handle_pq.close()
            file_handle_pc.close()
        else:
            #print(np.shape(self.parisi))
            #print(np.shape(self.theo_ifo))
            index = 0
            length = 2
            for b in range(int(np.log2(len(self.configurations[0])))):
                filename_pq = f'parisi_dist_block{b}_size{self.size}_sample{self.sample}.dat'
                filename_pc = f'theo_ifo_dist_block{b}_size{self.size}_sample{self.sample}.dat'
                file_handle_pc = open(filename_pc, "w")
                file_handle_pq = open(filename_pq, "w")
                for k in range(self.npt):
                    pq = np.histogram(a = self.parisi[:, index:index+length, k], bins=self.bins, density=True)[0]
                    pc = np.histogram(a = self.theo_ifo[b, :, :, k], bins=self.bins, density=True)[0]
                    np.savetxt(file_handle_pq, np.c_[q, pq, np.full(shape=np.shape(pq),fill_value=self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                    np.savetxt(file_handle_pc, np.c_[q, pc, np.full(shape=np.shape(pc),fill_value=self.temperatures[k])], fmt="%4e", delimiter="\t", newline="\n")
                    file_handle_pc.write("\n\n")
                    file_handle_pq.write("\n\n")
                file_handle_pq.close()
                file_handle_pc.close()
                index += length
                length *= 2

            
    def PrintOverlap(self):
        pass


                
            
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
        
    


        
