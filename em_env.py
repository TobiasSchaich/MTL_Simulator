# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 13:16:30 2020

@author: Tobias
"""

import numpy as np
from scipy.constants import epsilon_0, pi

class EmEnv:
    """Environment in which simulation takes place including frequencies"""
    def __init__(self, frequencies: list, eps_rel: float = 1., loss_tan=0.):
        self.set_freq(frequencies)
        self.set_eps_rel(eps_rel) 
        self.set_loss_tan(loss_tan)
        
    def set_freq(self, frequencies):
        if isinstance(frequencies, float) or isinstance(frequencies, int):
            self._freq=np.array(float(frequencies))
        elif isinstance(frequencies, list):
            self._freq=np.array(frequencies)
        elif isinstance(frequencies, np.ndarray): #assume that these are numerical
            self._freq=frequencies
        else:
            raise TypeError("Frequencies must be int, float or list of floats")
            
    def get_freq(self):
        return self._freq
    
    def get_eps_rel(self):
        return self._eps_rel
    
    def set_eps_rel(self, eps_rel):
        if isinstance(eps_rel, float) or isinstance(eps_rel, int):
            if eps_rel>=1:
                self._eps_rel=float(eps_rel)
            else:
                raise ValueError("eps_rel must be greater or equal to one.")
        else:
            raise TypeError("relative permittivity must be of type float or integer")
    
    def get_loss_tan(self):
        return self._loss_tan
    
    def get_sigma_diel(self):
        return self.get_loss_tan()*self.get_eps_rel()*(2*pi*self.get_freq())*epsilon_0
    
    def set_loss_tan(self, loss_tan):
        if isinstance(loss_tan, float) or isinstance(loss_tan, int):
            if loss_tan>=0.:
                self._loss_tan=float(loss_tan)
            else:
                raise ValueError("loss_tan must be greater or equal to zero.")
        else:
            raise TypeError("loss tangent must be of type float or integer")