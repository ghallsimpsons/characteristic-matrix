# -*- coding: utf-8 -*-
"""
Created on Tue Jun 02 23:52:17 2015

@author: Grantland
"""

import re
import numpy as np

def unitize_a(angle):
    """
    Takes a string such as "58 deg" or a number and returns the
    angle as a float, in radians.
    """
    try:
        s_angle=str( '{0:f}'.format(angle) )
    except:#Can't format if already string, but cast to string just in case
        s_angle=str(angle)
    ret_angle=float(re.findall(r"-?[0-9]*\.?[0-9]*", s_angle)[0])
    if re.match(r".*deg.*", s_angle, re.I):
        return np.deg2rad(ret_angle)
    return ret_angle
    
def unitize_f(freq):
    """
    Takes a string such as "150 GHz" or a number and returns the
    frequency as a float, in Hz.
    """
    try:
        s_freq=str( '{0:f}'.format(freq) )
    except:#Can't format if already string, but cast to string just in case
        s_freq=str(freq)
    ret_freq=float("0"+re.findall(r"[0-9]*\.?[0-9]*", s_freq)[0])
    unit = max(re.findall(r"[a-zA-Z]*", s_freq), key=len).lower()
    if unit == "khz":
        ret_freq*=1E3
    elif unit == "mhz":
        ret_freq*=1E6
    elif unit == "ghz":
        ret_freq*=1E9
    return ret_freq
    
def unitize_d(dist):
    """
    Takes a string such as "3 um" or a number and returns the
    distance as a float, in meters.
    """
    try:
        s_dist=str( '{0:f}'.format(dist) )
    except:#Can't format if already string, but cast to string just in case
        s_dist=str(dist)
    ret_dist=float("0"+re.findall(r"[0-9]*\.?[0-9]*", s_dist)[0])
    unit = max(re.findall(r"[a-zA-Z]*", s_dist), key=len).lower()
    if unit == "mm":
        ret_dist*=1E-3
    elif unit == "cm":
        ret_dist*=1E-2
    elif unit == "um":
        ret_dist*=1E-6
    return ret_dist