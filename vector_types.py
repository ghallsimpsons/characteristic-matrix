# -*- coding: utf-8 -*-
"""
Created on Wed Jun 03 14:32:21 2015

@author: Grantland
"""

from numpy import arctan2, sin, cos, arcsin, angle, sqrt, pi, exp
from numpy import isclose, allclose
from unitutils import unitize_a

def _isclosemod(a, b, atol=1E-5, mod=2*pi):
    """
    Return whether two numbers (or arrays) are within atol of each other
    in the modulo space determined by mod.
    """
    return (isclose(a%mod, b%mod, atol=atol) 
            or isclose((a+atol)%mod, (b+atol)%mod, atol=atol))

class PolarizationVector:
    def __init__(self, *args):
        """
        Essentially just an abstraction of complex numbers, representing the
        electric field at a specific point in space-time.
        Takes amplitude and wave phase in radians, or a complex number.
        """
        if len(args) == 2:
            # Args are amplitude, phase
            self.pol = args[0]*(exp(1j*args[1]))
        elif len(args) == 1:
            # Arg is complex number
            self.pol = args[0]
    def __add__(self, other):
        a1 = self.amp
        a2 = other.amp
        p1 = self.phase
        p2 = other.phase
        amp_new = sqrt(a1**2 + a2**2 + 2*a1*a2*cos(p1-p2))
        phase_new = arctan2( a1*sin(p1) + a2*sin(p2), a1*cos(p1) + a2*cos(p2) )
        return PolarizationVector(amp_new, phase_new)
    def __sub__(self, other):
        a1 = self.amp
        a2 = -other.amp
        p1 = self.phase
        p2 = other.phase
        amp_new = sqrt(a1**2 + a2**2 + 2*a1*a2*cos(p1-p2))
        phase_new = arctan2( a1*cos(p1) + a2*cos(p2), a1*sin(p1) + a2*sin(p2) )
        return PolarizationVector(amp_new, phase_new)
    def __mul__(self, num):
        """Scalar multiplication."""
        assert isinstance(num, (int, long, float, complex))
        amp_new = self.amp * abs(num)
        phase_new = self.phase + angle(num)
        return PolarizationVector(amp_new, phase_new)
    def __rmul__(self, num):
        """Scalar multiplication."""
        assert isinstance(num, (int, long, float, complex))
        amp_new = self.amp * abs(num)
        phase_new = self.phase + angle(num)
        return PolarizationVector(amp_new, phase_new)
    @property
    def power(self):
        """
        Power is the square of the amplitude of the electric field.
        """
        return self.amp**2
    @property
    def amp(self):
        return abs(self.pol)
    @property
    def phase(self):
        return angle(self.pol)
    def __eq__(self, other):
        """
        Check if vector is essentially equal to another. This makes it
        easy to confirm that vector transformations are behaving as they
        should.
        """
        if not _isclosemod(self.phase, other.phase):
            if _isclosemod(self.phase, other.phase+pi):
                # If they're pi out of phase, flip one of the amplitudes
                return isclose(self.amp, -other.amp, atol=1E-5)
            else:
                return False
        return isclose(self.amp, other.amp)
    def __ne__(self, other):
        #This isn't the default behavior because <insert bogus explanation here>
        return not self==other
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
            self.phase = 0
        elif len(args) == 5:
            # stokes repr with phase
            self.I = args[0]
            self.Q = args[1]
            self.U = args[2]
            self.V = args[3]
            self.phase = args[4]
        elif len(args) == 2:
            x = args[0]
            y = args[1]
            self.I = x.power + y.power
            self.Q = x.power - y.power
            self.U = 2*x.amp*y.amp*(cos(x.phase-y.phase))
            self.V = -2*x.amp*y.amp*(sin(x.phase-y.phase))
            if isclose(x.amp, 0, atol=1E-5):
                if isclose(y.amp, 0, atol=1E-5):
                    self.phase = 0
                else:
                    """
                    If x is zero, use y phase. Since stokes vectors can't be
                    modified except by casting to cartesian coordinates, the
                    result of this check will be preserved.
                    """
                    self.phase = y.phase
            else:
                self.phase = x.phase
        else:
            raise TypeError("Arguments must either be I, Q, U, V or two " \
                            "perpendicular PolarizationVectors.")
    @property
    def cartesian(self):
        """
        Returns the PolarizationTwoVector corresponding to the Stokes vector.
        """
        tilt = arctan2(self.U, self.Q)/2
        if self.I == 0:
            elipticity = 0
        else:
            elipticity = arcsin(self.V/self.I)/2
        E_x = sqrt(self.I)*(cos(tilt)*cos(elipticity)-1j*sin(tilt)*sin(elipticity))
        E_y = sqrt(self.I)*(sin(tilt)*cos(elipticity)+1j*cos(tilt)*sin(elipticity))
        # Since pure Stokes vectors don't preserve phase, reconstruct from saved phase
        if isclose(E_x, 0, atol=1E-5):
            offset = self.phase-angle(E_y)
        else:
            offset = self.phase-angle(E_x)
        E_x = E_x*exp(1j*offset)
        E_y = E_y*exp(1j*offset)
        v_x = PolarizationVector(E_x)
        v_y = PolarizationVector(E_y)
        return PolarizationTwoVector(v_x, v_y)
    @property
    def vect(self):
        return [self.I, self.Q, self.U, self.V]
    def rot(self, angle):
        """Returns rotated copy of self."""
        return self.cartesian.rot(angle).stokes
    def _rot(self, angle):
        """Rotates self by angle."""
        self = self.rot(angle)
    def __add__(self, other):
        if not isinstance(other, PolarizationTwoVector):
            other = other.cartesian
        return (self.cartesian+other).stokes
    def __sub__(self, other):
        if not isinstance(other, PolarizationTwoVector):
            other = other.cartesian
        return (self.cartesian-other).stokes
    def __eq__(self, other):
        """
        Returns whether or not two stokes vectors are essentially equal.
        """
        if isinstance(other, StokesVector):
            return allclose([self.I, self.Q, self.U, self.V], 
                [other.I, other.Q, other.U, other.V], atol=1E-5) and (
                _isclosemod(self.phase, other.phase) or isclose(self.I, 0, atol=1E-5))
        elif isinstance(other, PolarizationTwoVector):
            return self == other.stokes
    def __ne__(self, other):
        #This isn't the default behavior because <insert bogus explanation here>
        return not self==other
    def __repr__(self):
        return "I: {}, Q: {}, U: {}, V: {}, Phase: {}".format(
            self.I, self.Q, self.U, self.V, self.phase) 
    
class PolarizationTwoVector:
    """
    This class is similar to a Stokes Vector,
    but designed to allow phase-offset vector addition.
    """
    def __init__(self, vector_x, vector_y):
        self.v_x = vector_x
        self.v_y = vector_y
    def rot(self, angle):
        """Return rotated copy of self."""
        rad = unitize_a(angle)
        x_to_rot_x = self.v_x.amp*cos(rad)
        x_to_rot_y = -self.v_x.amp*sin(rad)
        y_to_rot_x = self.v_y.amp*sin(rad)
        y_to_rot_y = self.v_y.amp*cos(rad)
        new_x = PolarizationVector(x_to_rot_x, self.v_x.phase) + \
            PolarizationVector(y_to_rot_x, self.v_y.phase)
        new_y = PolarizationVector(x_to_rot_y, self.v_x.phase) + \
            PolarizationVector(y_to_rot_y, self.v_y.phase)
        return PolarizationTwoVector(new_x, new_y)   
    def _rot(self, angle):
        """Rotate self by angle."""
        self = self.rot(angle)
    def __add__(self, other):
        if not isinstance(other, PolarizationTwoVector):
            other = other.cartesian
        return PolarizationTwoVector( self.v_x + other.v_x, 
                                      self.v_y + other.v_y )
    def __sub__(self, other):
        if not isinstance(other, PolarizationTwoVector):
            other = other.cartesian
        return PolarizationTwoVector( self.v_x - other.v_x,
                                      self.v_y - other.v_y)
    @property
    def stokes(self):
        return StokesVector(self.v_x, self.v_y)
    @property
    def I(self):
        return self.stokes.I
    @property
    def Q(self):
        return self.stokes.Q
    @property
    def U(self):
        return self.stokes.U
    @property
    def V(self):
        return self.stokes.V
    def __eq__(self, other):
        if isinstance(other, PolarizationTwoVector):
            return self.v_x == other.v_x and self.v_y == other.v_y
        elif isinstance(other, StokesVector):
            return self.stokes == other
    def __ne__(self, other):
        #This isn't the default behavior because <insert bogus explanation here>
        return not self==other
    def __repr__(self):
        return "X:\n{}\nY:\n{}".format(self.v_x, self.v_y)
