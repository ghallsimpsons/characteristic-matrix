# -*- coding: utf-8 -*-
"""
@author: Grantland Hall
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""

import numpy as np
from numpy import pi, sqrt, exp
import re

def unitize_f(freq):
    try:
        freq=str( '{0:f}'.format(freq) )
    except:#Can't format if already string, but cast to string just in case
        freq=str(freq)
    freqbase=float(re.findall(r"[0-9]*\.?[0-9]*", freq)[0])
    maxlen=0
    unit=""
    for u in re.findall(r"[a-zA-Z]*", freq):
        if(len(u)>maxlen):
            maxlen=len(u)
            unit=u
    if unit.lower() == "khz":
        freqbase*=1000
    elif unit.lower() == "mhz":
        freqbase*=1000000
    elif unit.lower() == "ghz":
        freqbase*=1000000000
    return freqbase
    
class Layer:
    """A single dialectric layer. These are compounded using the Interface class."""
    def __init__(self, eps=1, thickness=-1):
        self.thickness = thickness
        self.eps = eps
    def get_matrix(self, freq):
        if self.thickness == -1:
            thickness = 300000000.0/(4*sqrt(self.eps)*unitize_f(freq)) # Use quarter wavelength as default
        else:
            thickness = self.thickness
        return c_matrix(freq=freq, thickness=thickness, eps=self.eps)


class c_matrix:
    """Used for normal incidence TE wave on linear dialectric."""
    def __init__(self, freq, thickness, eps=1):
        self.wavelength=300000000.0/unitize_f(freq)
        self.permitivity=eps
        self.thickness=thickness
    def M(self, perm_next, z=0):
        if z==0:
            z=self.thickness
        r=(sqrt(self.permitivity)-sqrt(perm_next))/(sqrt(perm_next)+sqrt(self.permitivity))
        return np.matrix([[exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z),0],
                          [0, exp(-2j*pi/self.wavelength*sqrt(self.permitivity)*z)]])*np.matrix([[1,r],[r,1]])

#Construct interface of multiple dialectric slabs
class _InterfaceMatrix:
    '''Takes a tuple as an argument and returns an object with .t and .matrix attributes'''
    def __init__(self, layers):
        perm_list=[]
        perm_last=1
        for x in layers:
            perm_list.append(x.permitivity)
        self.t=1+(1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0]))
        perm_list.append(perm_last)
        self.matrix=np.matrix([[1,(1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0]))],[(1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0])),1]])
        i=0
        for x in layers:
            self.matrix=self.matrix*x.M(perm_list[i+1])
            self.t*=1+(sqrt(perm_list[i])-sqrt(perm_list[i+1]))/(sqrt(perm_list[i])+sqrt(perm_list[i+1]))
            i+=1

class Interface:
    """Wrapper around interface_matrix."""
    def __init__(self, *arg):
        self.layers = arg
    def construct(self, freq):
        matrices = [x.get_matrix(freq) for x in self.layers]
        return _InterfaceMatrix(matrices)
        
def trans(c_mat):
    """Calculate transmission ratio for medium/media."""
    return 1/abs(c_mat.matrix.item(0,0)/c_mat.t)**2