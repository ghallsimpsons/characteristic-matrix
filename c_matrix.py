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
    def M(self, layer_next):
        """
        Builds the matrix operator representing the interface between
        self and layer_next.
        """
        perm_next = layer_next.permitivity
        perm2_next = layer_next.permitivity2
        r_xx = (sqrt(self.permitivity)-sqrt(perm_next))/(sqrt(perm_next)+sqrt(self.permitivity))
        r_yy = (sqrt(self.permitivity2)-sqrt(perm2_next))/(sqrt(perm2_next)+sqrt(self.permitivity2))
        r_xy = (sqrt(self.permitivity)-sqrt(perm2_next))/(sqrt(perm2_next)+sqrt(self.permitivity))
        r_yx = (sqrt(self.permitivity2)-sqrt(perm_next))/(sqrt(perm_next)+sqrt(self.permitivity2))
        angle_diff = layer_next.angle - self.angle
        return np.matrix([
             [cos(angle_diff)/(1+r_xx),  r_xx*cos(angle_diff)/(1+r_xx),
                 sin(angle_diff)/(1+r_xy),  r_xy*sin(angle_diff)/(1+r_xy)],
             [r_xx*cos(angle_diff)/(1+r_xx),  cos(angle_diff)/(1+r_xx),
                 r_xy*sin(angle_diff)/(1+r_xy),  sin(angle_diff)/(1+r_xy)],
             [-sin(angle_diff)/(1+r_yx),  -r_yx*sin(angle_diff)/(1+r_yx),
                 cos(angle_diff)/(1+r_yy),  r_yy*cos(angle_diff)/(1+r_yy)],
             [-r_yx*sin(angle_diff)/(1+r_yx),  -sin(angle_diff)/(1+r_yx),
                 r_yy*cos(angle_diff)/(1+r_yy),  cos(angle_diff)/(1+r_yy)]
             ])
    def phase(self, freq):
        """
        Returns the phase propogation operator through the layer.
        """
        wavelength=c/unitize_f(freq)
        z=self.thickness
        return np.matrix([
                   [exp(2j*pi/wavelength*sqrt(self.permitivity)*z),0,0,0],
                   [0, exp(-2j*pi/wavelength*sqrt(self.permitivity)*z),0,0],
                   [0,0, exp(2j*pi/wavelength*sqrt(self.permitivity2)*z),0],
                   [0,0,0, exp(-2j*pi/wavelength*sqrt(self.permitivity2)*z)]
                   ])

class _InterfaceMatrix:
    """
    Turns a list of layers into an interaction matrix, representing the
    action on polarized light of a given frequency through the layer stack.
    There is no compelling reason for this to be a separate class from
    the Interface class.
    """
    def __init__(self, layers, freq):
        perm_list=[]
        perm_list2=[]
        # What material the interface exits to (air)
        perm_last=1
        for x in layers:
            perm_list.append(x.permitivity)
            perm_list2.append(x.permitivity2)
        perm_list.append(perm_last)
        perm_list2.append(perm_last)
        # Calculate the reflection coefficients for the first layer from air
        r_xx = (1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0]))
        r_yy = (1-sqrt(perm_list2[0]))/(1+sqrt(perm_list2[0]))
        r_xy = (1-sqrt(perm_list[0]))/(1+sqrt(perm_list[0]))
        r_yx = (1-sqrt(perm_list2[0]))/(1+sqrt(perm_list2[0]))
        theta_0 = layers[0].angle
        # 4x4 analog of the Hou reflection matrix
        self.matrix=np.matrix([
             [cos(theta_0)/(1+r_xx),  r_xx*cos(theta_0)/(1+r_xx),
                 sin(theta_0)/(1+r_xy),  r_xy*sin(theta_0)/(1+r_xy)],
             [r_xx*cos(theta_0)/(1+r_xx),  cos(theta_0)/(1+r_xx),
                 r_xy*sin(theta_0)/(1+r_xy),  sin(theta_0)/(1+r_xy)],
             [-sin(theta_0)/(1+r_yx),  -r_yx*sin(theta_0)/(1+r_yx),
                 cos(theta_0)/(1+r_yy),  r_yy*cos(theta_0)/(1+r_yy)],
             [-r_yx*sin(theta_0)/(1+r_yx),  -sin(theta_0)/(1+r_yx),
                 r_yy*cos(theta_0)/(1+r_yy),  cos(theta_0)/(1+r_yy)]
            ])
        for i, this_layer in enumerate(layers):
            if i+1==len(layers):
                # Last layer
                angle_diff = this_layer.angle
                rot_mat = rot( this_layer.M(c_matrix(0)), angle_diff )
            else:
                angle_diff = this_layer.angle - layers[i+1].angle
                rot_mat = rot( this_layer.M(layers[i+1]), angle_diff )
            self.matrix=self.matrix * this_layer.phase(freq) * rot_mat

class Interface:
    """Wrapper around interface_matrix."""
    def __init__(self, *arg):
        self.layers = arg
    def build(self, freq):
        matrices = [x.get_matrix() for x in self.layers]
        self._built_freq = freq
        self.c_mat = _InterfaceMatrix(matrices, freq)
        self.r_mat = invert(self.c_mat.matrix)
    def trans(self, freq):
        """Calculate transmission ratio for medium/media."""
        self.build(freq)
        c_mat = invert(self.c_mat.matrix)
        return ( abs( 1/(c_mat.item(0,0)+c_mat.item(2,0)) )**2/2 +
                    abs( 1/(c_mat.item(2,2)+c_mat.item(0,2)) )**2/2 )
        #return abs(1/(c_mat.matrix.item(0,0)))**2/2 + abs(1/(c_mat.matrix.item(2,2)))**2/2
    def cross(self, freq):
        """
        Returns the total transmitted power along each axis,
        the cross x-y power, and cross y-x power. The i-j cross
        power is the output power along the j axis with an input
        polarization along the i axis.
        """
        self.build(freq)
        r_mat = self.r_mat
        return (
                abs(1/(r_mat.item(0,0)-r_mat.item(0,2)*r_mat.item(2,0)/r_mat.item(2,2))),#t_xx
                abs(1/(r_mat.item(2,2)-r_mat.item(0,2)*r_mat.item(2,0)/r_mat.item(0,0))),#t_yy
                abs(r_mat.item(0,2)/(r_mat.item(0,0)*r_mat.item(2,2)-r_mat.item(2,0)*r_mat.item(0,2))),#t_xy
                abs(r_mat.item(2,0)/(r_mat.item(0,0)*r_mat.item(2,2)-r_mat.item(2,0)*r_mat.item(0,2))) #t_yx
                )
    def __mul__(self, vect):
        assert hasattr(self, "c_mat"), ("You must build the transmission matrix"
            " for a given frequency before using the interface as an operator.")
        if isinstance(vect, StokesVector):
            # Convert to PolarizationTwoVector, then multiply
            v = vect.cartesian
            s_trans = self * v
            return s_trans.stokes
        elif isinstance(vect, PolarizationTwoVector):
            r_mat = self.r_mat
            t_xx = 1/(r_mat.item(0,0)-r_mat.item(0,2)*r_mat.item(2,0)/r_mat.item(2,2))
            t_yy = 1/(r_mat.item(2,2)-r_mat.item(0,2)*r_mat.item(2,0)/r_mat.item(0,0))
            t_xy = r_mat.item(0,2)/(r_mat.item(0,0)*r_mat.item(2,2)-r_mat.item(2,0)*r_mat.item(0,2))
            t_yx = r_mat.item(2,0)/(r_mat.item(0,0)*r_mat.item(2,2)-r_mat.item(2,0)*r_mat.item(0,2))
            v_x = t_xx*vect.v_x + t_yx*vect.v_y
            v_y = t_yy*vect.v_y + t_xy*vect.v_x
            return PolarizationTwoVector(v_x, v_y)
        else:
            raise TypeError("Operand must be a StokesVector or PolarizationVector")
    def __repr__(self):
        if hasattr(self, "r_mat"):
            return "Interface built at f={}\n{}".format(self._built_freq, str(self.r_mat))
        else:
            return "Interface has not been built yet."