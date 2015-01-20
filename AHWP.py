# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 11:05:30 2015

@author: Grantland
"""

"""
Calculation of cross-polarozation in AHWP.
Assume that the AR-saphire interface has ~0 reflecion.
Furthermore, if the AR layer is matched to the mean of the
ordinary and extroardinary axes, then mismatch on either
axis is roughly equivalent. Thus, we can assume that total
reflection = differential reflection = 0 at the AR interface.
Now take the intermediate HWP layers, where the extroardinary
axis of HWP_n is offset from HWP_(n-1) by an angle \phi.
Then the reflection can be broken into 4 components - 
the component along each of x, y in layer n-1 projected on the
axes x', y' in n. Assuming each HWP layer is eqivalent, 
R_eo given by R_eo = [(n_e-n_o)/(n_e+n_o)]^2 is equal to R_oe.
Then the reflection operator for layer n acting on layer n-1
is given by:
(x; y)(R_oo cos(\phi), R_eo sin(\phi); R_eo sin(\phi), R_ee cos(\phi))
where the y axis is the extroardinary axis in layer n. We can
see that if the layers are aligned, sin(\phi) = 0 and the form
reduces to R = (0, 0; 0, 0), as expected.

The steps in calculating the cross polarization is as follows:
1. Calculate the frequency dependent reflection, assuming a uniform
single-index layer for the HWP layers.
2. Starting at HWP_2, calculate the reflection matrix with HWP_(n-1).
3. Combine these matrices with Rot(\phi) projection at each layer.
4. Act the resulting operator on [cos(\theta); sin(\theta)] where
\theta is the polarization angle relative to the ordinary axis of
the first layer. Since the first and last HWP layers are aligned,
the resulting vector is the polarization state in the initial
coordinate system.

Phase consideration:
    Let R_0 be the reflection matrix at the top of a HWP interface,
    and R_1 the bottom. The phase offset from a second-order reflection
    is p=2pi*f*n*(d*2)/c and the resulting power along the x axis is
    a_x^2+b_x^2+2a_x*b_x*cos(p), where b_x=(R_1*R_0*(x,y))_x, the
    amplitude of the second-order reflected wave.
    
Things that need to be fixed:
    For n_x>n_y, R_yx should cause a pi phase shift.
    When creating stokes vector, it assumes axis offset = 0
    
Current Examples:
    Add two parallel vectors - 
        result = PolarizationVector(1, pi/2) + PolarizationVector(0.1, pi)
    Add two perpendicular vectors - 
        result = PolarizationTwoVector( \
                PolarizationVector(1,0), PolarizationVector(0.2, pi/3))
    Add two vectors rotated pi/6 -
        result = PolarizationTwoVector( .1, 0, 0 ) +
                PolarizationTwoVector( .4, pi/2, pi/6 )

"""

from numpy import cos, sin, pi, arctan2
import numpy as n

AHWP_e = 9 # extroardinary index
AHWP_o = 8 # ordinary index
AHWP_phi = [0, 15, 30, 15, 0] # orientation of each layer

c = 300000000

class AHWP:
    def __init__(self, index_o, index_e, d, angles):
        """
        @param index_o: index of refraction along ordinary axis
        @param index_e: index of refraction along extraordinary axis
        @param d: thickness of each HWP layer
        @param angles: list of orientation of each layer, in radians
        """
        self.n_o = index_o
        self.n_e = index_e
        self.R_oe = ((index_o-index_e)/(index_o+index_e))**2
        self.angles = angles
        self.matrix = np.matrix([[1, 0], [0, 1]])
    def build(self, freq):
        """
        TODO: Refactor to build as many matrices as possible in init,
        so they only need to be computed once, or memoize in build.
        """
        angles = self.angles
        R_oe = self.R_oe
        n_e = self.n_e
        n_o = self.n_o
        for i in xrange(len(angles)):
            # Compute reflection at interface, plus phase offset
            # from second order reflection.
            
            if (i == len(angles)):
                # Reflection at last AR interface
                R_1 = np.matrix([[0, 0], [0, 0]])
                R_0 = np.matrix([0, R_oe*cos(angles[i]-angles[i-1])],
                                [R_oe*cos(angles[i]-angles[i-1]), 0])
            elif (i == 0):
                # Reflection at first AR interface
                R_0 = np.matrix([[0, 0], [0, 0]])
                R_1 = np.matrix([0, R_oe*cos(angles[i+1]-angles[i])],
                                [R_oe*cos(angles[i+1]-angles[i]), 0])
            else:
                R_0 = np.matrix([0, R_oe*cos(angles[i]-angles[i-1])],
                                [R_oe*cos(angles[i]-angles[i-1]), 0])
                R_1 = np.matrix([0, R_oe*cos(angles[i+1]-angles[i])],
                                [R_oe*cos(angles[i+1]-angles[i]), 0])
            second_r = R_1*R_0
            p = 4*pi*freq*(n_e-n_o)*d/c
            transmit = R_1 + second_r

class PolarizationVector:
    def __init__(self, amplitude, phase):
        """
        @param amplitude: Amplitude of wave
        @param phase: Phase of wave in radians
        """
        self.amp = amplitude
        self.phase = phase
    def __add__(self, other):
        a1 = self.amp
        a2 = other.amp
        p1 = self.phase
        p2 = other.phase
        amp_new = (a1**2) + (a2**2) + (2*a1*a2*cos(p1-p2))
        phase_new = arctan2( a1*cos(p1) + a2*cos(p2), a1*sin(p1) + a2*sin(p2) )
        return PolarizationVector(amp_new, phase_new)
    def power(self):
        return self.amp**2

class StokesVector:
    def __init__(self, *args):
        if len(args) == 4:
            # stokes repr
            self.I = args[0]
            self.Q = args[1]
            self.U = args[2]
            self.V = args[3]
        elif len(args) == 2:
            x = args[0]
            y = args[1]
            I = x.power() + y.power()
            Q = x.power() - y.power()
            U = 2*x.amp*y.amp*(cos(x.phase)*cos(y.phase)+sin(x.phase)*sin(y.phase))
            V = 2*x.amp*y.amp*(cos(x.phase)*sin(y.phase)+sin(x.phase)*cos(y.phase))

class PolarizationTwoVector:
    """
    This class is similar to a Stokes Vector,
    but designed to allow phase-offset vector addition.
    """
    def __init__(self, vector_x, vector_y, axis_orientation = 0):
        self.v_x = vector_x
        self.v_y = vector_y
        self.orientation = axis_orientation

    def rot(self, rad):
        """Return rotated representation of self"""
        x_to_rot_x = vector_x.amp*cos(rad)
        x_to_rot_y = -1*vecrot_x.amp*sin(rad)
        y_to_rot_x = vector_y.amp*cos(sin)
        y_to_rot_y = vector_y.amp*cos(rad)
        new_x = PolarizationVector(x_to_rot_x, vector_x.phase) + \
            PolarizationVector(y_to_rot_x, vector_y.phase)
        new_y = PolarizationVector(x_to_rot_y, vector_x.phase) + \
            PolarizationVector(y_to_rot_y, vector_y.phase)
        return PolarizationTwoVector(new_x, new_y, self.orientation + rad)
    
    def _rot(self, rad):
        """Rotate self"""
        self = self.rot(rad)
    
    def __add__(self, other):
        """
        Returns two-vector in the orientation of the first two-vector
        """
        offset = other.orientation - self.orientation
        other_r = other.rot(-offset)
        return PolarizationTwoVector( self.v_x + other_r.v_x, 
                            self.v_y + other_r.v_y, self.orientation )

    def stokes(self):
        unrotated = self.rot( -self.orientation )
        return StokesVector(self.v_x, self.v_y)
