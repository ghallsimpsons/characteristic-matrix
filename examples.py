# -*- coding: utf-8 -*-
"""
@author: Grantland
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""

import numpy as np
from numpy import pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from unitutils import unitize_f, unitize_a
from vector_types import (PolarizationVector, PolarizationTwoVector,
                          StokesVector)

central_freq="150GHz"

base_thick=300000000/(4*unitize_f("150ghz"))
default_input = StokesVector(1, 1, 0, 0) # I=Q=1, U=V=0

def rot_sweep(interface, freq, step='1 deg', pol_in=default_input):
    """
    Plots the output Q and U vs rotation angle given an input polarization
    vector. The output is normalized by the input intensity. The interface
    must be built before calling this sweep.
    """
    Q=np.array([])
    U=np.array([])

    f = unitize_f (freq)
    radstep = unitize_a(step)
    angle_range = np.arange(0, 2*pi, radstep)
    
    for angle in angle_range:
        rot_int = interface(angle).build(f)
        pol_out = rot_int*pol_in
        Q = np.append(Q, pol_out.Q)
        U = np.append(U, pol_out.U)

    Q = Q/pol_in.I
    U = U/pol_in.I

    sin_to_fit = lambda x, phase: np.sin(phase + 4*x)

    phase_fit = curve_fit(sin_to_fit, angle_range, Q)
    sin_fit = np.sin(4*angle_range+phase_fit[0])
    
    diff = np.sum((Q-sin_fit)**2)*radstep/(2*pi)
    print "Power difference: {}".format(diff)
    
    plt.plot(angle_range, Q, label="Q")
    plt.plot(angle_range, U, label="U")
    plt.plot(angle_range, sin_fit, label="Fit")
    plt.xlabel("Angle (radians)")
    plt.ylabel("Normalized Transmission")
    plt.ylim(ymax=1)
    plt.legend()
    plt.show()
    


def graph(interface, start="1ghz", stop="300ghz", step="0.1ghz"):
    """Show frequency spectrum of interface."""
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    x_vals=[]
    # Choose a sane base for the plot
    divisor=1
    frq_range="Hz"
    if (stop+start)/2>1000000000:
        divisor=1000000000
        frq_range="GHz"
    elif (stop+start)/2>1000000:
        divisor=1000000
        frq_range="MHz"
    elif (stop+start)/2>1000:
        divisor=1000
        frq_range="kHz"
    for x in xrange(int((stop-start)/step)):
        x_vals.append((x*step+start)/divisor)
    (trans, xy_cross)=get_data(interface, start, stop, step)
    plt.plot(x_vals, trans, label="Total transmission")
    plt.plot(x_vals, xy_cross, label="Cross polarization")
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Transmission Ratio")
    plt.ylim(ymax=1)
    plt.legend()
    plt.show()

def get_data(interface, start, stop, step):
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    trans=[]
    xy_cross=[]
    yx_cross=[]
    unit = PolarizationVector(1)
    zero = PolarizationVector(0)
    vect_x = PolarizationTwoVector(unit,zero)
    vect_y = PolarizationTwoVector(zero,unit)
    for f in xrange(int((stop-start)/step)):
        interface.build(start+f*step)
        out_x = interface*vect_x
        trans.append(out_x.v_x.power)
        xy_cross.append(out_x.v_y.power)
    return (trans, xy_cross)
