# -*- coding: utf-8 -*-
"""
@author: Grantland
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""
import time
from multiprocessing import Process, Pool, cpu_count

import numpy as np
from numpy import pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from unitutils import unitize_f, unitize_a, unitize_d, frange
from vector_types import (PolarizationVector, PolarizationTwoVector,
                          StokesVector)

central_freq="150GHz"
c=300000000.0

base_thick=c/(4*unitize_f("150ghz"))
default_input = StokesVector(1, 1, 0, 0) # I=Q=1, U=V=0

rms = lambda arr: np.sqrt(np.mean(np.square(arr)))

def contrived_ahwp(f_range,n1,n2,d):
    a_range = np.arange(0, pi/2, unitize_a("1deg"))
    d=unitize_d(d)
    center_angle = unitize_a("58deg")
    mods = []
    for f in f_range:
        mueller = mueller_wp((n2-n1)*2*pi*f*d/c)
        p = []
        for a in a_range:
            s = StokesVector(1,1,0,0)
            v = mueller.dot(np.array([s.rot(a).vect]).reshape([4,1])).reshape([4,])
            s = StokesVector(*v)
            v = mueller.dot(np.array([s.rot(a+center_angle).vect]).reshape([4,1])).reshape([4,])
            s = StokesVector(*v)
            v = mueller.dot(np.array([s.rot(a-center_angle).vect]).reshape([4,1])).reshape([4,])
            p.append(np.sqrt(v[1]**2+v[2]**2))
        mods.append(rms(p))
    return mods
        
def mueller_wp(phase):
    """
    Returns the Mueller matrix in the standard basis for a waveplate which
    induces the given phase offset.
    """
    return np.array([[1,0,0,0],
                      [0,1,0,0],
                      [0,0,np.cos(phase),-np.sin(phase)],
                      [0,0,np.sin(phase),np.cos(phase)]])

def rot_sweep_data(interface, freq, step='1 deg', a_range=2*pi, pol_in=default_input):
    """
    Returns the I, Q, U and V over a rotation of 2pi.
    """
    I=np.array([])
    Q=np.array([])
    U=np.array([])
    V=np.array([])

    a_range=unitize_a(a_range)
    f = unitize_f(freq)
    interface.build(f)
    
    radstep = unitize_a(step)
    angle_range = np.arange(0, a_range, radstep)
    
    for angle in angle_range:
        # Rotating the input polarization rather than rebuilding the interface
        # at each rotation gives a factor of 2 improvement in performance.
        rot_pol = pol_in.rot(angle)
        pol_out = (interface*rot_pol).rot(-angle)
        I = np.append(I, pol_out.I)
        Q = np.append(Q, pol_out.Q)
        U = np.append(U, pol_out.U)
        V = np.append(V, pol_out.V)  

    I = I/pol_in.I
    Q = Q/pol_in.I
    U = U/pol_in.I
    V = V/pol_in.I
    
    return (I, Q, U, V)
    
def mp_rot_sweep_data(interface, freq, step='1 deg', a_range=2*pi, pol_in=default_input):
    """
    Returns the I, Q, U and V over a rotation of 2pi.
    """
    Q=np.array([])
    U=np.array([])

    a_range=unitize_a(a_range)
    f = unitize_f(freq)
    interface.build(f)
    
    radstep = unitize_a(step)
    angle_range = np.arange(0, a_range, radstep)
    
    for angle in angle_range:
        # Rotating the input polarization rather than rebuilding the interface
        # at each rotation gives a factor of 2 improvement in performance.
        rot_pol = pol_in.rot(angle)
        pol_out = (interface*rot_pol).rot(-angle)
        Q = np.append(Q, pol_out.Q)
        U = np.append(U, pol_out.U)

    Q = Q/pol_in.I
    U = U/pol_in.I
    
    eff = rms(np.sqrt(Q**2+U**2)) 
    
    return eff
    
def rot_sweep(interface, freq, step='1 deg', pol_in=default_input):
    """
    Plots the output Q and U vs rotation angle given an input polarization
    vector. The output is normalized by the input intensity. The interface
    must be built before calling this sweep.
    """
    I, Q, U, V = rot_sweep_data(interface, freq, step=step, pol_in=pol_in)
    
    sin_to_fit = lambda x, phase: np.sin(phase + 4*x)

    radstep = unitize_a(step)
    angle_range = np.arange(0, 2*pi, radstep)
    phase_fit = curve_fit(sin_to_fit, angle_range, Q)
    sin_fit = np.sin(4*angle_range+phase_fit[0])
    
    diff = np.sum((Q-sin_fit)**2)*radstep/(2*pi)
    print "Power difference: {}".format(diff)
    print "Minimum transmission: {}".format(min(I))
    
    plt.plot(angle_range, Q, label="Q")
    plt.plot(angle_range, U, label="U")
    plt.plot(angle_range, I, label="I")
    plt.plot(angle_range, V, label="V")
    plt.plot(angle_range, sin_fit, label="Fit")
    plt.xlabel("Angle (radians)")
    plt.ylabel("Normalized Transmission")
    plt.ylim(ymax=1)
    plt.legend()
    plt.show()
    
def mod_vs_freq(interface, fstart, fstop, fstep, astep="1deg"):
    """
    Calculates and plots the modulation efficiency vs frequency.
    """
    start_time = time.time()
    fstart = unitize_f(fstart)
    fstop = unitize_f(fstop)
    fstep = unitize_f(fstep)
    f_vals = []
    efficiency = []
    divisor, frq_range = frange(fstart, fstop)
    #contrived_data = contrived_ahwp(
    #    np.arange(fstart,fstop,fstep),np.sqrt(9.26),np.sqrt(11.56),"3.6mm")
    for f in xrange(int((fstop-fstart)/fstep)):
        f_vals.append((f*fstep+fstart)/divisor)
        I,Q,U,V=rot_sweep_data(interface, fstart+f*fstep, step=astep, a_range=pi/2)
        efficiency.append(rms(np.sqrt(Q**2+U**2)))
    plt.plot(f_vals, efficiency, label="Modulation Efficiency")
    #plt.plot(f_vals, contrived_data, label="Contrived AHWP")
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Modulation Efficiency")
    plt.ylim(ymin=0,ymax=1)
    plt.legend()
    plt.show()   
    end_time = time.time()
    print "elapsed: {}s".format((end_time-start_time))

def mp_mod_vs_freq(interface, fstart, fstop, fstep, astep="1deg"):
    """
    Calculates and plots the modulation efficiency vs frequency.
    """
    nprocs = cpu_count()-1
    pool = Pool(nprocs)
    print "Pooling {} procs.".format(nprocs)
        
    start_time = time.time()
    fstart = unitize_f(fstart)
    fstop = unitize_f(fstop)
    fstep = unitize_f(fstep)
    divisor, frq_range = frange(fstart, fstop)
    
    efficiency = [pool.apply(mp_rot_sweep_data,(interface,f*fstep+fstart,astep,pi/2))
                            for f in range(int((fstop-fstart)/fstep))]
    f_vals = [(f*fstep+fstart)/divisor for f in range(int((fstop-fstart)/fstep))]
    
    plt.plot(f_vals, efficiency, label="Modulation Efficiency")
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Modulation Efficiency")
    plt.ylim(ymin=0,ymax=1)
    plt.legend()
    plt.show()   
    end_time = time.time()
    print "elapsed: {}s".format((end_time-start_time))

def graph(interface, start="1ghz", stop="300ghz", step="0.1ghz"):
    """Show frequency spectrum of interface."""
    start = unitize_f(start)
    stop = unitize_f(stop)
    step = unitize_f(step)
    x_vals=[]
    # Choose a sane base for the plot
    divisor, frq_range = frange(start,stop)
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
    """Get polarization transmittion data over a range of frequencies."""
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
