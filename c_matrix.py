# -*- coding: utf-8 -*-
"""
Created May 2013

@author: Grantland Hall
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
import re

#Expand SI units
def unitize_f(freq):
    try:
        freq=str( '{0:f}'.format(freq) )
    except:#Can't format if already string, but cast to string just in case
        freq=str(freq)
    freqbase=float(re.findall(r"[0-9]+\.?[0-9]*", freq)[0])
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
    
#Used for normal incidence TE wave on linear dialectric
class c_matrix:
    def __init__(self, freq, thickness, eps=1):
        self.wavelength=300000000.0/unitize_f(freq)
        self.permitivity=eps
        self.thickness=thickness
    def M(self, z=0):
        if z==0:
            z=self.thickness
        return np.array([[cos(2*pi/self.wavelength*sqrt(self.permitivity)*z),-1j/sqrt(self.permitivity)*sin(2*pi/self.wavelength*sqrt(self.permitivity)*z)],
                          [-1j*sqrt(self.permitivity)*sin(2*pi/self.wavelength*sqrt(self.permitivity)*z), cos(2*pi/self.wavelength*sqrt(self.permitivity)*z)]])

#Construct interface of multiple dialectric slabs
class interface:
    def __init__(self, *arg):
        self.matrix=np.array([[1,0],[0,1]])
        for x in arg:
            self.matrix=np.dot(self.matrix, x.M())
    def layer(self, new_m):
        self.matrix=np.dot(self.matrix, new_m)
        
#Calculate transmission ratio for medium/media.
def trans(c_mat, eps_1=1, eps_l=1):
    p_1=sqrt(eps_1)
    p_l=sqrt(eps_l)
    return abs(2*p_1/((c_mat.matrix.item(0,0)+c_mat.matrix.item(0,1)*p_l)*p_1+(c_mat.matrix.item(1,0)+c_mat.matrix.item(1,1)*p_l)))