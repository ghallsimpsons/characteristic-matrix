# -*- coding: utf-8 -*-
"""
@author: Grantland
May 31, 2013

Created with Spyder IDE
"""

import numpy as np
from scipy.signal import argrelextrema
import c_matrix as cm
from c_matrix import unitize_f, trans
import matplotlib.pyplot as plt
from numpy import sqrt

central_freq="150GHz"


base_thick=300000000/(4*unitize_f("150ghz"))#central_freq))

#default for optimal results
def interface_default(freq, layer=0, thick=base_thick):
    layer_1=cm.c_matrix(freq, thick/sqrt(2), 2)
    layer_2=cm.c_matrix(freq, thick/sqrt(4), 4)
    layer_3=cm.c_matrix(freq, thick/sqrt(7), 7)
    layer_4=cm.c_matrix(freq, .045, 9.6)
    return cm.interface(layer_1, layer_2, layer_3, layer_4, layer_3, layer_2, layer_1)

#default for current materials
def interface_current(freq, layer=0, thick=base_thick):
    layer_perms=[2, 4, 7, 9.6]
    air=cm.c_matrix(freq, 1000, 1)
    layer_1=cm.c_matrix(freq, thick/sqrt(layer_perms[0]), layer_perms[0])
    layer_2=cm.c_matrix(freq, thick/sqrt(layer_perms[1]), layer_perms[1])
    layer_3=cm.c_matrix(freq, thick/sqrt(layer_perms[2]), layer_perms[2])
    layer_4=cm.c_matrix(freq, .006, layer_perms[3])
    if layer>0:
        exec("layer_"+str(layer)+"=cm.c_matrix(freq, thick, layer_perms[layer-1] )")
    return cm.interface(layer_1, layer_2, layer_3, layer_4, layer_3, layer_2, layer_1)
    
#Show frequency spectrum of medium defined in trans_y
def graph(layer=0, layer_thick=base_thick, start="1ghz", stop="300ghz", step="0.1ghz"):
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
    y_vals=get_data(layer, layer_thick, start, stop, step)
    print y_vals
    #plt.plot(x_vals, y_vals)
    #plt.xlabel("Frequency ("+frq_range+")")
    #plt.ylabel("Transmission Ratio")
    #plt.show()

def get_data(layer, layer_thick, start, stop, step):
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    y_vals=[]
    for x in xrange(int((stop-start)/step)):
        y_vals.append(trans(interface_current(x*step+start, layer, layer_thick)))
    return y_vals
    
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
        
def sma(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')