# -*- coding: utf-8 -*-
"""
Created May 31, 2013

@author: Grantland Hall
grantlandhall@berkeley.edu
"""

import numpy as np
from numpy import pi, sqrt, exp
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
    def M(self, perm_next, z=0):
        if z==0:
            z=self.thickness
        r=(sqrt(self.permitivity)-sqrt(perm_next))/(sqrt(perm_next)+sqrt(self.permitivity))
#        print "M: "+str([[exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z),r*exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z)],
#                          [r*exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z), exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z)]])
        return np.matrix([[exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z),r*exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z)],
                          [r*exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z), exp(2j*pi/self.wavelength*sqrt(self.permitivity)*z)]])

#Construct interface of multiple dialectric slabs
class interface:
    def __init__(self, *arg):
        perm_list=[]
        perm_last=1
        self.t=1
        for x in arg:
            perm_list.append(x.permitivity)
        self.matrix=np.matrix([[1,(1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0]))],[(1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0])),1]])
        print self.matrix
        i=0
        for x in arg:
            try:
                self.matrix=self.matrix*x.M(perm_list[i+1])
                self.t*=sqrt(perm_list[i])*2/(sqrt(perm_list[i])+sqrt(perm_list[i+1]))
                print self.matrix
            except:
                self.matrix=self.matrix*x.M(perm_last)
                self.t*=sqrt(perm_list[i])*2/(sqrt(perm_list[i])+sqrt(perm_last))
                print "YIPEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE!"
            i+=1
        print self.matrix.item(0,0)
        
#Calculate transmission ratio for medium/media.
def trans(c_mat):
    return c_mat.t/abs(c_mat.matrix.item(0,0))**2