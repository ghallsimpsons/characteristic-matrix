# -*- coding: utf-8 -*-
"""
@author: Grantland
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""

import numpy as np
#from scipy.signal import argrelextrema
from c_matrix import unitize_f
import matplotlib.pyplot as plt
from unitutils import *
from numpy import sqrt

central_freq="150GHz"


base_thick=300000000/(4*unitize_f("150ghz"))#central_freq))
    
def graph(interface, start="1ghz", stop="300ghz", step="0.1ghz"):
    """Show frequency spectrum of interface."""
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    x_vals=[]
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
    (trans, xy_cross, yx_cross)=get_data(interface, start, stop, step)
    plt.plot(x_vals, trans, label="Total transmission")
    plt.plot(x_vals, xy_cross, label="Cross polarization")
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Transmission Ratio")
    #plt.legend()
    plt.show()

def get_data(interface, start, stop, step):
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    trans=[]
    xy_cross=[]
    yx_cross=[]
    for x in xrange(int((stop-start)/step)):
        t_vals = interface.cross(x*step+start)
        trans.append((t_vals[0])**2/2+(t_vals[1])**2/2+(t_vals[2])**2/2+(t_vals[3])**2/2)
        xy_cross.append((t_vals[2])**2)
        yx_cross.append((t_vals[3])**2)
    return (trans, xy_cross, yx_cross)
    
#If I recall, this isn't working and I don't have time to figure out why.
"""    
def statistics_of(layer, layer_thick, start="1ghz", stop="300ghz", step="0.1ghz"):
    y_vals=get_data(layer, layer_thick, start, stop, step)
    local_mins=argrelextrema(np.array(y_vals), np.less)[0]
    envelope_maxima=[local_mins[loc] for loc in argrelextrema(np.array([y_vals[i] for i in local_mins]),np.greater)]
    start=unitize_f(start)
    stop=unitize_f(stop)
    step=unitize_f(step)
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
    bandwidth=str((envelope_maxima[0][-1]-envelope_maxima[0][0])*step/1000000000)
    print "Bandwidth: "+bandwidth+"GHz"
    print "Area: "+str((envelope_maxima[0][-1]-envelope_maxima[0][0]-np.sum(sqrt(y_vals[envelope_maxima[0][0]:envelope_maxima[0][-1]])))*step/10000000/float(bandwidth))

def run_trial(layer, layer_start, layer_stop, layer_step, start="1ghz", stop="300ghz", step="0.1ghz"):
    print "Data for thickness of layer "+str(layer)
    for x in xrange(int((layer_stop-layer_start)/layer_step)):
        print "Thickness: "+str(x*layer_step+layer_start)+"m"
        statistics_of(layer, x*layer_step+layer_start, start, stop, step)
    """
    
def sma(interval, window_size):
    """Smoothes a data set of 'window_size' length by averaging over 'interval' data points."""
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
