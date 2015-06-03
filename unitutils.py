# -*- coding: utf-8 -*-
"""
Created on Tue Jun 02 23:52:17 2015

@author: Grantland
"""

import re

def unitize_f(freq):
    """
    Takes a string such as "150 GHz" and returns the
    frequency as a float, in Hz.
    """
    try:
        s_freq=str( '{0:f}'.format(freq) )
    except:#Can't format if already string, but cast to string just in case
        s_freq=str(freq)
    ret_freq=float(re.findall(r"[0-9]*\.?[0-9]*", s_freq)[0])
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
    Takes a string such as "3 um" and returns the
    distance as a float, in meters.
    """
    try:
        s_dist=str( '{0:f}'.format(dist) )
    except:#Can't format if already string, but cast to string just in case
        s_dist=str(dist)
    ret_dist=float(re.findall(r"[0-9]*\.?[0-9]*", s_dist)[0])
    unit = max(re.findall(r"[a-zA-Z]*", s_dist), key=len).lower()
    if unit == "mm":
        ret_dist*=1E-3
    elif unit == "cm":
        ret_dist*=1E-2
    elif unit == "um":
        ret_dist*=1E-6
    return ret_dist