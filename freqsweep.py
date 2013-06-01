# -*- coding: utf-8 -*-
"""
@author: Grantland
May 31, 2013

Created with Spyder IDE
"""

import c_matrix as cm
from c_matrix import unitize_f, trans
import matplotlib.pyplot as plt

thick=.0005 #thickness of outer layers, in m

def trans_y(freq):
    layer_1=cm.c_matrix(freq, thick, 2)
    layer_2=cm.c_matrix(freq, thick, 4)
    layer_3=cm.c_matrix(freq, thick, 7)
    layer_4=cm.c_matrix(freq, .006, 9.6)
    ar = cm.interface(layer_1, layer_2, layer_3, layer_4, layer_3, layer_2, layer_1)
    return trans(ar, 2, 2)
    
def graph(start, stop, step):
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    x_vals=[]
    y_vals=[]
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
        y_vals.append(trans_y(x*step+start))
    plt.plot(x_vals, y_vals)
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Transmission Ratio")
    plt.show()