# -*- coding: utf-8 -*-
"""
@author: Grantland
grantlandhall@berkeley.edu

Last Modified: Nov 2, 2013
"""
import time
# Import unqualified module for pp
import numpy
import numpy as np
from numpy import pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from unitutils import unitize_f, unitize_a, unitize_d, frange, deg
from vector_types import (PolarizationVector, PolarizationTwoVector,
                          StokesVector)
from c_matrix import Layer, Interface

from mp import mpexec

hwp_axis_1 = 11.56
hwp_axis_2 = 9.36
hwp = Layer(eps=hwp_axis_1, eps2=hwp_axis_2, thickness="3.6mm")
coating_2a = 2.32
coating_2b = 6.79
ar2a = Layer(coating_2a, "400um")
ar2b = Layer(coating_2b, "230um")
dual_band_ahwp = Interface(ar2a, ar2b, hwp, hwp("58deg"), hwp, ar2b, ar2a)

central_freq="150GHz"
c=300000000.0

base_thick=c/(4*unitize_f("150ghz"))
default_input = StokesVector(1, 1, 0, 0) # I=Q=1, U=V=0

rms = lambda arr: np.sqrt(np.mean(np.square(arr)))

def contrived_ahwp_mod(f_range,n1,n2,d):
    a_range = np.arange(0, pi/2, unitize_a("1deg"))
    d=unitize_d(d)
    center_angle = unitize_a("58deg")
    mods = []
    for f in f_range:
        mueller = mueller_wp((n2-n1)*2*pi*f*d/c)
        I = np.array([])
        Q = np.array([])
        for a in a_range:
            s = StokesVector(1,1,0,0)
            v = mueller.dot(np.array([s.rot(a).vect]).reshape([4,1])).reshape([4,])
            s = StokesVector(*v)
            v = mueller.dot(np.array([s.rot(center_angle).vect]).reshape([4,1])).reshape([4,])
            s = StokesVector(*v)
            v = mueller.dot(np.array([s.rot(-center_angle).vect]).reshape([4,1])).reshape([4,])
            s = StokesVector(*v).rot(-a)
            Q = np.append(Q,s.Q)
        mods.append((max(1+Q)-min(1+Q))/2)
    return mods

def contrived_ahwp_cross(f_range, n1, n2, d, angle):
    d=unitize_d(d)
    center_angle = unitize_a("58deg")
    out_angles = []
    for f in f_range:
        mueller = mueller_wp((n2-n1)*2*pi*f*d/c)
        p = []
        s = StokesVector(1,1,0,0)
        v = mueller.dot(np.array([s.rot(angle).vect]).reshape([4,1])).reshape([4,])
        s = StokesVector(*v)
        v = mueller.dot(np.array([s.rot(center_angle).vect]).reshape([4,1])).reshape([4,])
        s = StokesVector(*v)
        v = mueller.dot(np.array([s.rot(-center_angle).vect]).reshape([4,1])).reshape([4,])
        s = StokesVector(*v)
        a = s.rot(-angle).pol_angle
        a = (a+pi)%(pi)
        out_angles.append(a)
    return out_angles

def contrived_cross_sweep(fstart, fstop, fstep, n1, n2, d):
    d = unitize_d(d)
    fstart = unitize_f(fstart)
    fstop = unitize_f(fstop)
    fstep = unitize_f(fstep)
    divisor, frq_range = frange(fstart,fstop)
    arange = np.arange(0, pi/2, unitize_a("15deg"))
    f_range = np.arange(fstart, fstop, fstep)
    for angle in arange:
        angles = contrived_ahwp_cross(f_range, n1, n2, d, angle)
        label = "AHWP Angle: {0:.2f}".format(angle)
        plt.plot(f_range, angles, label=label)
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Polarization Angle")
    plt.ylim(ymin=0, ymax=pi)
    plt.legend()
    plt.show()
        
def mueller_wp(phase):
    """
    Returns the Mueller matrix in the standard basis for a waveplate which
    induces the given phase offset.
    """
    return np.array([[1,0,0,0],
                      [0,1,0,0],
                      [0,0,np.cos(phase),-np.sin(phase)],
                      [0,0,np.sin(phase),np.cos(phase)]])

# We can't use default args because pp was poorly designed
def rot_sweep_data(interface, freq, step, a_range, pol_in):
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
    
    return (I, Q, U, V)
    
def phase_vs_freq(interface, fstart, fstop, fstep="1GHz"):
    """
    Plot polarization rotation phase vs frequency for a given interface.
    """
    sin_to_fit = lambda x, a, phase: a*np.sin(phase + 4*x)
    pol_in = StokesVector(1,1,0,0)
    fstart = unitize_f(fstart)
    fstop = unitize_f(fstop)
    fstep = unitize_f(fstep)
    f_range = np.arange(fstart, fstop, fstep)
    divisor, frq_range = frange(fstart, fstop)

    radstep = unitize_a("1deg")
    angle_range = np.arange(0, pi/2, radstep)

    data = np.array(mpexec(rot_sweep_data, "freq", f_range,
        interface=interface, step="1deg", a_range=pi/2, pol_in=pol_in,
        globals=globals()))
    pol = []
    phase = []
    print "starting fit"
    for i, f in enumerate(f_range):
        # Get Q vs angle, indexed by frequency
        pol.append(data[i,1])
        phase_fit = curve_fit(sin_to_fit, angle_range, pol[i])
        phase.append(phase_fit[0][1])
    min_phase = min(phase)
    phase -= min_phase
    dphase = map(deg, phase)
    plt.plot(f_range/divisor, dphase)
    plt.xlabel("Frequency ({})".format(frq_range))
    plt.ylabel("Phase angle (deg)")
    plt.show()

def rot_sweep(interface, freq, pol_in, step='1 deg'):
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
    contrived_data = contrived_ahwp_mod(
        np.arange(fstart,fstop,fstep),np.sqrt(9.26),np.sqrt(11.56),"3.6mm")
    for f in xrange(int((fstop-fstart)/fstep)):
        f_vals.append((f*fstep+fstart)/divisor)
        #I,Q,U,V=rot_sweep_data(interface, fstart+f*fstep, step=astep, a_range=pi/2)
        #efficiency.append(rms(np.sqrt(Q**2+U**2)))
    end_time = time.time()
    print "elapsed: {}s".format((end_time-start_time))
    #plt.plot(f_vals, efficiency, label="Modulation Efficiency")
    plt.plot(f_vals, contrived_data, label="Contrived AHWP")
    plt.xlabel("Frequency ("+frq_range+")")
    plt.ylabel("Modulation Efficiency")
    plt.ylim(ymin=0,ymax=1.1)
    plt.legend()
    plt.show()   

def cross_pol(interface, fvals, astep="1deg"):
    """
    Calculates and plots the cross polarization vs angle for a number of freqs.
    """
    f_vals = map(unitize_f,fvals)
    divisor, frq_range = frange(min(f_vals), max(f_vals))

    radstep = unitize_a(astep)
    angle_range = np.arange(0, pi/2, radstep)
    deg_range = map(deg, angle_range)

    data = mpexec(rot_sweep_data, "freq", f_vals, interface=interface,
            a_range=pi/2, step=astep, pol_in=StokesVector(1,1,0,0),
            globals=globals())
    for i, f in enumerate(f_vals):
        I,Q,U,V = data[i]
        plt.plot(deg_range, U, label="{}{}".format(int(f/divisor), frq_range))
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Cross Polarization")
    plt.ylim(ymin=-1,ymax=1)
    plt.legend()
    plt.show()   

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
