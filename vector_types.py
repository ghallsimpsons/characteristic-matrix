# -*- coding: utf-8 -*-
"""
Created on Wed Jun 03 14:32:21 2015

@author: Grantland
"""

from numpy import arctan2, sin, cos, arcsin, angle
from unitutils import unitize_f, unitize_d, unitize_a

class PolarizationVector:
    def __init__(self, *args):
        """
        Essentially just an abstraction of complex numbers, representing the
        electric field at a specific point in space-time.
        Takes amplitude and wave phase in radians, or a complex number.
        """
        # TODO: Just make the internal repr a complex number, and make
        #       amp and phase computed properties
        if len(args) == 2:
            self.amp = args[0]
            self.phase = args[1]
        elif len(args) ==1:
            self.amp = abs(args[0])
            self.phase = angle(args[0])
    def __add__(self, other):
        a1 = self.amp
        a2 = other.amp
        p1 = self.phase
        p2 = other.phase
        amp_new = (a1**2) + (a2**2) + (2*a1*a2*cos(p1-p2))
        phase_new = arctan2( a1*cos(p1) + a2*cos(p2), a1*sin(p1) + a2*sin(p2) )
        return PolarizationVector(amp_new, phase_new)
    def __mul__(self, num):
        assert isinstance(num, (int, long, float, complex))
        amp_new = self.amp * abs(num)
        phase_new = self.phase + angle(num)
        return PolarizationVector(amp_new, phase_new)
    def __rmul__(self, num):
        assert isinstance(num, (int, long, float, complex))
        amp_new = self.amp * abs(num)
        phase_new = self.phase + angle(num)
        return PolarizationVector(amp_new, phase_new)
    def power(self):
        return self.amp**2
    def __repr__(self):
        return "Amplitude: {}\nRelative phase: {}".format(self.amp, self.phase)

class StokesVector:
    def __init__(self, *args):
        """
        Takes either I, Q, U, V or two perpendicular PolarizationVectors,
        and returns the Stokes representation.
        """
        if len(args) == 4:
            # stokes repr
            self.I = args[0]
            self.Q = args[1]
            self.U = args[2]
            self.V = args[3]
        elif len(args) == 2:
            x = args[0]
            y = args[1]
            self.I = x.power() + y.power()
            self.Q = x.power() - y.power()
            self.U = 2*x.amp*y.amp*(cos(x.phase)*cos(y.phase)+sin(x.phase)*sin(y.phase))
            self.V = 2*x.amp*y.amp*(cos(x.phase)*sin(y.phase)+sin(x.phase)*cos(y.phase))
        else:
            raise TypeError("Arguments must either be I, Q, U, V or two " \
                            "perpendicular PolarizationVectors.")
    @property
    def cartesian(self):
        """
        Returns the PolarizationTwoVector corresponding to the Stokes vector.
        """
        tilt = arctan2(self.U, self.Q)/2
        elipticity = arcsin(self.V/self.I)/2
        E_x = cos(tilt)*cos(elipticity)-1j*sin(tilt)*sin(elipticity)
        E_y = sin(tilt)*cos(elipticity)+1j*cos(tilt)*sin(elipticity)
        v_x = PolarizationVector(E_x)
        v_y = PolarizationVector(E_y)
        return PolarizationTwoVector(v_x, v_y)
    def __repr__(self):
        return "I: {}, Q: {}, U: {}, V: {}".format(self.I, self.Q, self.U, self.V) 
    
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
        x_to_rot_x = self.vector_x.amp*cos(rad)
        x_to_rot_y = -self.vecrot_x.amp*sin(rad)
        y_to_rot_x = self.vector_y.amp*cos(rad)
        y_to_rot_y = self.vector_y.amp*cos(rad)
        new_x = PolarizationVector(x_to_rot_x, self.vector_x.phase) + \
            PolarizationVector(y_to_rot_x, self.vector_y.phase)
        new_y = PolarizationVector(x_to_rot_y, self.vector_x.phase) + \
            PolarizationVector(y_to_rot_y, self.vector_y.phase)
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
    @property
    def stokes(self):
        return StokesVector(self.v_x, self.v_y)
    def __repr__(self):
        return "X:\n{}\nY:\n{}".format(self.v_x, self.v_y)