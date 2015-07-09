# -*- coding: utf-8 -*-
"""
@author: Grantland Hall
grantlandhall@berkeley.edu
"""

import numpy as np
from numpy import pi, sqrt, exp, sin, cos
from numpy.linalg import inv as invert
from unitutils import unitize_f, unitize_d, unitize_a
from vector_types import StokesVector, PolarizationTwoVector

CENTRAL_FREQ = 1.24E11
c = 3E8
mu_0 = 1.2566371E-6

def rot(matrix, theta):
    """
    Perform a tensor rotation of the operator matrix
    by an angle theta.
    """
    R = np.matrix([[cos(theta), 0, sin(theta), 0],
                      [0, cos(theta), 0, sin(theta)],
                      [-sin(theta), 0, cos(theta),0],
                      [0,-sin(theta), 0, cos(theta)]])
    # Hardcode Rinv, since the form is simple, and matrix inversion is hard
    Rinv = np.matrix([[cos(theta), 0, -sin(theta), 0],
                      [0, cos(theta), 0, -sin(theta)],
                      [sin(theta), 0, cos(theta),0],
                      [0,sin(theta), 0, cos(theta)]])
    return R*matrix*Rinv

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
        Return a copy of self, but rotated at a different angle.
        """
        return Layer(self.eps, self.thickness, self.eps2, angle)
    def get_matrix(self):
        return c_matrix(thickness=self.thickness, eps=self.eps, eps2=self.eps2, angle=self.angle)

class c_matrix:
    """Used for normal incidence TE wave on linear dialectric or birefringent material."""
    def __init__(self, thickness, eps=1, eps2=1, angle=0):
        self.permitivity=eps
        self.permitivity2=eps2
        self.thickness=thickness
        self.angle=angle
    def M(self, layer_next, freq):
        """
        Returns the matrix operator which acts on the Poynting vector,
        (Ex,Ey,Hx,Hy). We do this by first projecting into E space, acting
        the phase propogator, and then projection back into S.
        We combine the last two steps, and a change of basis, into the 
        second matrix, because it lends itself to a compact form.
        """
        Cx = sqrt(self.permitivity)/(c*mu_0)
        Cy = sqrt(self.permitivity2)/(c*mu_0)
        E_to_S = np.matrix([
            [ 1 , 0 , 1 , 0 ],
            [ 0 , 1 , 0 , 1 ],
            [ 0 ,-Cy, 0 ,Cy ],
            [Cx , 0 ,-Cx, 0 ]
            ])
        wavelength = c/unitize_f(freq)
        z = self.thickness
        theta = layer_next.angle - self.angle
        phase_x = exp(2j*pi/wavelength*sqrt(self.permitivity)*z)
        phase_y = exp(2j*pi/wavelength*sqrt(self.permitivity2)*z)
        phase = 0.5*np.matrix([
            [phase_x*cos(theta),phase_x*sin(theta),-phase_y*sin(theta)/Cx,phase_y*cos(theta)/Cx],
            [-phase_y*sin(theta),phase_y*cos(theta),-phase_x*cos(theta)/Cy,-phase_x*sin(theta)/Cy],
            [phase_x*cos(theta),phase_x*sin(theta),phase_y*sin(theta)/Cx,-phase_y*cos(theta)/Cx],
            [-phase_y*sin(theta),phase_y*cos(theta),phase_x*cos(theta)/Cy,phase_x*sin(theta)/Cy]
            ])
        return E_to_S*phase

class _InterfaceMatrix:
    """
    Turns a list of layers into an interaction matrix, representing the
    action on polarized light of a given frequency through the layer stack.
    There is no compelling reason for this to be a separate class from
    the Interface class.
    """
    def __init__(self, layers, freq):
        perm_list_x=[]
        perm_list_y=[]
        # What material the interface exits to (air)
        perm_last=1
        for x in layers:
            perm_list_x.append(x.permitivity)
            perm_list_y.append(x.permitivity2)
        theta = layers[0].angle
        m = 1
        C = 1/(c*mu_0)
        for i, this_layer in enumerate(layers[:-1]):
            m= m*this_layer.M(layers[i+1], freq)
        self.matrix = np.matrix([
            [cos(theta), sin(theta), -sin(theta)*C, cos(theta)*C],
            [-sin(theta), cos(theta), -cos(theta)*C, -sin(theta)*C],
            [cos(theta), sin(theta), sin(theta)*C, -cos(theta)*C],
            [-sin(theta), cos(theta), cos(theta)*C, sin(theta)*C]
                ]) * np.matrix([
            [m.item(0,0)+m.item(0,3)/C, m.item(0,1)-m.item(0,2)/C],
            [m.item(1,0)+m.item(1,3)/C, m.item(1,1)-m.item(1,2)/C],
            [m.item(2,0)+m.item(2,3)/C, m.item(2,1)-m.item(2,2)/C],
            [m.item(3,0)+m.item(3,3)/C, m.item(3,1)-m.item(3,2)/C]
            ])

class Interface:
    """Wrapper around interface_matrix."""
    def __init__(self, *arg):
        self.layers = arg
    def build(self, freq):
        matrices = [x.get_matrix() for x in self.layers]
        self._built_freq = freq
        self.c_mat = _InterfaceMatrix(matrices, freq)
    def trans(self, freq):
        """Calculate transmission ratio for medium/media."""
        self.build(freq)
        c_mat = invert(self.c_mat.matrix)
        return ( abs( 1/(c_mat.item(0,0)+c_mat.item(2,0)) )**2/2 +
                    abs( 1/(c_mat.item(2,2)+c_mat.item(0,2)) )**2/2 )
        #return abs(1/(c_mat.matrix.item(0,0)))**2/2 + abs(1/(c_mat.matrix.item(2,2)))**2/2
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
            t_xx = ((m.item(1,1)*m.item(2,0)-m.item(2,1)*m.item(1,0))/
                    (m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0)))
            t_yy = ((m.item(0,0)*m.item(2,1)-m.item(2,0)*m.item(0,1))/
                    (m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0)))
            t_xy = ((m.item(3,0)*m.item(1,1)-m.item(3,1)*m.item(1,0))/
                    (m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0)))
            t_yx = ((m.item(3,1)*m.item(0,0)-m.item(3,0)*m.item(0,1))/
                    (m.item(1,1)*m.item(0,0)-m.item(0,1)*m.item(1,0)))
            v_x = t_xx*vect.v_x - t_xy*vect.v_y
            v_y = t_yy*vect.v_y - t_yx*vect.v_x
            return PolarizationTwoVector(v_x, v_y)
        else:
            raise TypeError("Operand must be a StokesVector or PolarizationVector")
    def __repr__(self):
        if hasattr(self, "r_mat"):
            return "Interface built at f={}\n{}".format(self._built_freq, str(self.r_mat))
        else:
            return "Interface has not been built yet."
