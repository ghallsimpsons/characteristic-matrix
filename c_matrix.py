# -*- coding: utf-8 -*-
"""
@author: Grantland Hall
grantlandhall@berkeley.edu
"""

import numpy as np
from numpy import pi, sqrt, exp, sin, cos
from unitutils import unitize_f, unitize_d, unitize_a
from vector_types import StokesVector, PolarizationTwoVector

CENTRAL_FREQ = 1.24E11
c = 3E8
mu_0 = 1.2566371E-6

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
        self.angle = unitize_a(angle)
        if (thickness==-1):
            # If no thickness is specified, choose a sane default
            self.thickness = c/(4*CENTRAL_FREQ*sqrt(eps))
        else:
            self.thickness = unitize_d(thickness)
        self.eps = eps  
    def __call__(self, angle):
        """
        Return a copy of self, but rotated by the given angle.
        """
        theta = unitize_a(angle)
        return Layer(self.eps, self.thickness, self.eps2, self.angle+theta)
    def _get_matrix(self):
        return c_matrix(thickness=self.thickness, eps=self.eps, eps2=self.eps2, angle=self.angle)

class c_matrix:
    """Used for normal incidence TE wave on linear dialectric or birefringent material."""
    def __init__(self, thickness, eps=1, eps2=1, angle=0):
        self.permitivity=eps
        self.permitivity2=eps2
        self.thickness=thickness
        self.angle=angle
    def M(self, freq, layer_next=None):
        """
        Returns the matrix operator which acts on the Poynting vector,
        (Ex,Ey,Hx,Hy). We do this by first projecting into E space, acting
        the phase propogator, and then projection back into S.
        We combine the last two steps, and a change of basis, into the 
        second matrix, because it lends itself to a compact form.
        """
        Cx = sqrt(self.permitivity)/(c*mu_0)
        Cy = sqrt(self.permitivity2)/(c*mu_0)
        E_to_S = np.matrix([ #Eq3.84
            [ 1 , 0 , 1 , 0 ],
            [ 0 , 1 , 0 , 1 ],
            [ 0 ,-Cy, 0 ,Cy ],
            [Cx , 0 ,-Cx, 0 ]
            ])
        wavelength = c/unitize_f(freq)
        z = self.thickness
        if layer_next is not None:
            theta = layer_next.angle - self.angle
        else:
            theta = -self.angle
        phase_x = exp(2j*pi/wavelength*sqrt(self.permitivity)*z)
        phase_y = exp(2j*pi/wavelength*sqrt(self.permitivity2)*z)
        inv_x = exp(-2j*pi/wavelength*sqrt(self.permitivity)*z)
        inv_y = exp(-2j*pi/wavelength*sqrt(self.permitivity2)*z)
        phase = 0.5*np.matrix([ #Eq3.85
            [phase_x*cos(theta),phase_x*sin(theta),-phase_x*sin(theta)/Cx,phase_x*cos(theta)/Cx],
            [-phase_y*sin(theta),phase_y*cos(theta),-phase_y*cos(theta)/Cy,-phase_y*sin(theta)/Cy],
            [inv_x*cos(theta),inv_x*sin(theta),inv_x*sin(theta)/Cx,-inv_x*cos(theta)/Cx],
            [-inv_y*sin(theta),inv_y*cos(theta),inv_y*cos(theta)/Cy,inv_y*sin(theta)/Cy]
            ])
        return E_to_S*phase #Eq3.86

class _InterfaceMatrix:
    """
    Turns a list of layers into an interaction matrix, representing the
    action on polarized light of a given frequency through the layer stack.
    There is no compelling reason for this to be a separate class from
    the Interface class.
    """
    def __init__(self, layers, freq):
        theta = layers[0].angle
        m = 1
        C = 1/(c*mu_0) #Eq3.73
        for i, this_layer in enumerate(layers[:-1]):
            m = m*this_layer.M(freq, layers[i+1]) #Eq3.87
        m = m*layers[-1].M(freq)
        self.matrix = 0.5*np.matrix([ #Eq3.94
            [cos(theta), sin(theta), -sin(theta)/C, cos(theta)/C],
            [-sin(theta), cos(theta), -cos(theta)/C, -sin(theta)/C],
            [cos(theta), sin(theta), sin(theta)/C, -cos(theta)/C],
            [-sin(theta), cos(theta), cos(theta)/C, sin(theta)/C]
                ]) * np.matrix([ #Eq3.93
            [m.item(0,0)+m.item(0,3)*C, m.item(0,1)-m.item(0,2)*C],
            [m.item(1,0)+m.item(1,3)*C, m.item(1,1)-m.item(1,2)*C],
            [m.item(2,0)+m.item(2,3)*C, m.item(2,1)-m.item(2,2)*C],
            [m.item(3,0)+m.item(3,3)*C, m.item(3,1)-m.item(3,2)*C]
            ])

class Interface:
    """Wrapper around interface_matrix."""
    def __init__(self, *arg):
        self.layers = arg
    def __call__(self, angle):
        """
        Returns a rotated representation of the entire interface by the
        given angle.
        """
        new_layers = [layer(angle) for layer in self.layers]
        return Interface(*new_layers)
    def build(self, freq):
        matrices = [x._get_matrix() for x in self.layers]
        self._built_freq = freq
        self.c_mat = _InterfaceMatrix(matrices, freq)
        return self
    def __mul__(self, vect):
        assert hasattr(self, "c_mat"), ("You must build the transmission matrix"
            " for some frequency before using the interface as an operator.")
        if isinstance(vect, StokesVector):
            # Convert to PolarizationTwoVector, then multiply
            v = vect.cartesian
            s_trans = self * v
            return s_trans.stokes
        elif isinstance(vect, PolarizationTwoVector):
            m = self.c_mat.matrix
            #Eq3.103-106
            t_xx = m.item(1,1)/(m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0))
            t_yy = m.item(0,0)/(m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0))
            t_xy = -m.item(0,1)/(m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0))
            t_yx = -m.item(1,0)/(m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0))
            v_x = t_xx*vect.v_x + t_xy*vect.v_y
            v_y = t_yy*vect.v_y + t_yx*vect.v_x
            return PolarizationTwoVector(v_x, v_y)
        else:   
            raise TypeError("Operand must be a StokesVector or PolarizationTwoVector")
    def __repr__(self):
        if hasattr(self, "c_mat"):
            return "Interface built at f={}\n{}".format(self._built_freq, str(self.c_mat.matrix))
        else:
            return "Interface has not been built yet."
