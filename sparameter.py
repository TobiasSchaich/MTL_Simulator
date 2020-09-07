# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:37:04 2020

@author: Tobias
"""

import numpy as np
import matplotlib.pyplot as plt

class Sparameter:
    '''Class for Sparameter objects'''
    def __init__(self, freq=[], s=[[[]]], z0=[]):
        self.set_freq(freq)
        self.set_z0(z0)
        self.set_s(s)

    def set_freq(self,freq):
        self._freq=freq
        
    def get_freq(self):
        return self._freq
    
    def set_z0(self, z0):
        self._z0=z0
        
    def get_z0(self):
        return self._z0
    
    def set_s(self, s):
        if isinstance(s,list) or isinstance(s, np.ndarray):
                s=np.array(s)
                if np.size(s)==0:
                    self._num_freq=0
                    self._num_ports=0
                    self._s=np.array([[[]]])
                elif s.ndim==3:
                    shp=np.shape(s)
                    self._num_freq=shp[0]
                    self._num_ports=shp[1]
                    self._s=s
                else:
                    raise TypeError("s must be 3 dimensional list or ndarray")
        else:
            raise TypeError("s parameter matrix must be of type list or ndarray")
        
        
    def get_s(self):
        return self._s

    def get_num_freq(self):
        return self._num_freq
    
    def get_num_ports(self):
        return self._num_ports
    
    def create_from_abcd(self,freq, abcd, z0=50): 
        '''Takes ABCD Matrix and Transforms it into S-parameters with reference Impedance z0 Ohms at every port'''    
        num_freq, num_ports =np.shape(abcd)[0:2]
        sparams=np.zeros(np.shape(abcd))
        a, b = np.split(np.split(abcd, 2,axis=1)[0], 2,axis=2)
        c, d = np.split(np.split(abcd, 2,axis=1)[1], 2,axis=2)
        b = b/z0;
        c = c*z0;
        
       # Choose the optimal math depending on the size of the network parameter
        if num_ports == 2:
            delta = a+b+c+d
            s11= (a+b-c-d)/ delta
            s12 = 2*(a*d-b*c)/delta
            s21=2*np.ones(np.shape(b))/delta
            s22=(-a+b-c+d)/delta 
            sparams=np.concatenate((np.concatenate((s11,s12), axis=1), np.concatenate((s21, s22), axis=1)), axis=2)
        else:    
            identity_3d=np.zeros(np.shape(a))
            idx = np.arange(int(num_ports/2))
            identity_3d[:, idx, idx] = 1
            s21 = np.linalg.solve(a+b+c+d, 2.*identity_3d)
            s22 = np.linalg.solve(a+b+c+d, (-a+b-c+d))
            s11 = 1/2*(a+b-c-d)@s21
            s12 = 1/2*((a+b-c-d)@s22+(a-b-c+d))
            sparams=np.concatenate((np.concatenate((s11,s12), axis=1), np.concatenate((s21, s22), axis=1)), axis=2)
        self.set_freq(freq)
        self.set_z0(z0)
        self.set_s(sparams)
        
    def plot(self, ports=None):
        """Plots port pair specified as tuple in ports. If none are specified, will plot all sparameters."""
        if ports is None:
            self._plotall()
        elif isinstance(ports, tuple):
            f=self.get_freq()
            s=self.get_s()[:,ports[0]-1, ports[1]-1]# -1 due to 0 base definition in python not common for sparameters/ports
            plt.plot(f,20*np.log10(np.abs(s)))
            plt.xlabel('Frequency (Hz)')
            plt.ylabel(f'S{ports[0]}{ports[1]} (dB)')
            plt.grid(True)
        else:
            raise TypeError("ports must either be of type none or a tuple of ports to be plotted ")
	
    def _plotall(self):
            f=self.get_freq()
            s=self.get_s()
            for i in range(self.get_num_ports()):
                for j in range(self.get_num_ports()):
                    plt.plot(f,20*np.log10(np.abs(s[:,i,j])), label=f'S{i+1}{j+1}')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')		
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Magnitude (dB)')
            plt.grid(True)
            