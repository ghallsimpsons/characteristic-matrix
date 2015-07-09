# -*- coding: utf-8 -*-
"""
@author: Grantland
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""

import numpy as np
import matplotlib.pyplot as plt
from unitutils import unitize_f
from vector_types import PolarizationVector, PolarizationTwoVector

central_freq="150GHz"

base_thick=300000000/(4*unitize_f("150ghz"))    

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
