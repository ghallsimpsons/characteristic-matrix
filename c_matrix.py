# -*- coding: utf-8 -*-
"""
@author: Grantland Hall
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""

import numpy as np
from numpy import pi, sqrt, exp, sin, cos
import re
from operator import mul

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
    """
    A single dialectric layer. These are compounded using the
    Interface class. Birefringent layers can be formed by specifying
    eps2 and the orientation of the primary axis.
    """
    def __init__(self, eps=1, thickness=-1, eps2=-1, angle=0):
        if eps2 == -1:
            eps2 = eps
        self.eps2 = eps2
        self.angle = angle
        self.thickness = thickness
        self.eps = eps
    def get_matrix(self):
        thickness = self.thickness
        return c_matrix(thickness=thickness, eps=self.eps, eps2=self.eps2, angle=self.angle)

class c_matrix:
    """Used for normal incidence TE wave on linear dialectric."""
    def __init__(self, thickness, eps=1, eps2=1, angle=0):
        self.permitivity=eps
        self.permitivity2=eps2
        self.thickness=thickness
        self.angle=angle
    def M(self, layer_next):
        perm_next = layer_next.permitivity
        perm2_next = layer_next.permitivity2
        r_xx = (sqrt(self.permitivity)-sqrt(perm_next))/(sqrt(perm_next)+sqrt(self.permitivity))
        r_yy = (sqrt(self.permitivity2)-sqrt(perm2_next))/(sqrt(perm2_next)+sqrt(self.permitivity2))
        r_xy = (sqrt(self.permitivity)-sqrt(perm2_next))/(sqrt(perm2_next)+sqrt(self.permitivity))
        r_yx = (sqrt(self.permitivity2)-sqrt(perm_next))/(sqrt(perm_next)+sqrt(self.permitivity2))
        angle_diff = layer_next.angle - self.angle
        return np.matrix([[cos(angle_diff),r_xx*cos(angle_diff),sin(angle_diff),r_xy*sin(angle_diff)],
                        [r_xx*cos(angle_diff), cos(angle_diff),r_xy*sin(angle_diff),sin(angle_diff)],
                        [sin(angle_diff), r_yx*sin(angle_diff),cos(angle_diff),r_yy*cos(angle_diff)],
                        [r_yx*sin(angle_diff), sin(angle_diff),r_yy*cos(angle_diff),cos(angle_diff)],
                        ])
    def phase(self, freq):
        wavelength=300000000.0/unitize_f(freq)
        z=self.thickness
        return np.matrix([[exp(2j*pi/wavelength*sqrt(self.permitivity)*z),0,0,0],
                   [0, exp(-2j*pi/wavelength*sqrt(self.permitivity)*z),0,0],
                   [0,0, exp(2j*pi/wavelength*sqrt(self.permitivity2)*z),0],
                   [0,0,0, exp(-2j*pi/wavelength*sqrt(self.permitivity2)*z)]])

#Construct interface of multiple dialectric slabs
class _InterfaceMatrix:
    '''Takes a tuple as an argument and returns an object with .t and .matrix attributes'''
    def __init__(self, layers, freq):
        perm_list=[]
        perm_list2=[]
        perm_last=1
        for x in layers:
            perm_list.append(x.permitivity)
            perm_list2.append(x.permitivity2)
        self.t_x=1+(1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0]))
        self.t_y=1+(1-sqrt(perm_list2[0]))/(1+sqrt(perm_list2[0]))
        perm_list.append(perm_last)
        perm_list2.append(perm_last)
        r_xy = (1-sqrt(perm_list[0]))/(sqrt(perm_list[0])+1)
        r_yx = (1-sqrt(perm_list2[0]))/(sqrt(perm_list2[0])+1)
        theta_0 = layers[0].angle
        self.matrix=np.matrix([
            [cos(theta_0), cos(theta_0)*(self.t_x-1), sin(theta_0), sin(theta_0)*r_xy],
            [cos(theta_0)*(self.t_x-1),   cos(theta_0),   sin(theta_0)*r_xy, sin(theta_0)],
            [sin(theta_0), sin(theta_0)*r_yx, cos(theta_0), cos(theta_0)*(self.t_y-1)],
            [sin(theta_0)*r_yx,   sin(theta_0),   cos(theta_0)*(self.t_y-1), cos(theta_0)],
            ])
        i=0
        for x in layers:
            if i+1==len(layers):
                #last layer
                theta_1 = 0
                self.matrix=self.matrix*x.M(c_matrix(0))*x.phase(freq)
            else:
                theta_1 = layers[i+1].angle
                self.matrix=self.matrix*x.M(layers[i+1])*x.phase(freq)
            diff =  theta_1 - layers[i].angle
            self.t_x*=(1+(sqrt(perm_list[i])-sqrt(perm_list[i+1]))/(sqrt(perm_list[i])+sqrt(perm_list[i+1])))*cos(diff) \
                    + (1+(sqrt(perm_list2[i])-sqrt(perm_list[i+1]))/(sqrt(perm_list2[i])+sqrt(perm_list[i+1])))*sin(diff)
            self.t_y*=(1+(sqrt(perm_list2[i])-sqrt(perm_list2[i+1]))/(sqrt(perm_list2[i])+sqrt(perm_list2[i+1])))*cos(diff) \
                    + (1+(sqrt(perm_list[i])-sqrt(perm_list2[i+1]))/(sqrt(perm_list[i])+sqrt(perm_list2[i+1])))*sin(diff)
            i+=1

class Interface:
    """Wrapper around interface_matrix."""
    def __init__(self, *arg):
        self.layers = arg
    def build(self, freq):
        matrices = [x.get_matrix() for x in self.layers]
        self.c_mat = _InterfaceMatrix(matrices, freq)
    def trans(self, freq):
        """Calculate transmission ratio for medium/media."""
        self.build(freq)
        c_mat = self.c_mat
        return abs(c_mat.t_x/c_mat.matrix.item(0,0))**2/2 + abs(c_mat.t_y/c_mat.matrix.item(2,2))**2/2
