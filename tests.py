# -*- coding: utf-8 -*-
"""
Created on Thu Jun 04 20:13:13 2015

@author: grantlandhall
"""

from c_matrix import Layer, Interface
from numpy import pi
from vector_types import StokesVector, PolarizationTwoVector, PolarizationVector

def warn(msg):
    print "\033[91m{}\033[0m".format(msg)

hwp_axis_1 = 11.56
hwp_axis_2 = 9.36
optimal_coating = 3.23

ar1 = Layer(optimal_coating)
hwp = Layer(eps=hwp_axis_1, eps2=hwp_axis_2, thickness="3.6mm")

ideal_ahwp = Interface(ar1, hwp, hwp("58deg"), hwp, ar1)

pol_in = StokesVector(1, -1, 0, 0)
ideal_ahwp.build("124GHz")
pol_out = ideal_ahwp * pol_in

cart_in = PolarizationTwoVector(PolarizationVector(1),PolarizationVector(0)).rot("-90 deg")
cart_out = ideal_ahwp * cart_in

pol_in_offset = StokesVector(1, -1, 0, 0, pi)

zero = StokesVector(0,0,0,0)

failed_tests = []
if pol_in != cart_in:
    failed_tests.append(("Cartesian representation not equal to Stokes:\n"
                         "Cartesian: {}\nStokes: {}").format(pol_in, cart_in))
if pol_out != cart_out:
    failed_tests.append(("Cartesian representation not equal to Stokes:\n"
                         "Cartesian: {}\nStokes: {}").format(pol_out, cart_out))
if pol_in == pol_out:
    failed_tests.append(("Input vector equal to output vector:\n"
                         "Input: {}\nOutput: {}").format(pol_in, pol_out))
if pol_out.I > pol_in.I:
    failed_tests.append(("Output intensity greater than input intensity:\n"
                         "Input: {}\nOutput: {}").format(pol_in, pol_out))
if pol_in != pol_in.cartesian.stokes:
    failed_tests.append(("Doubly transformed Stokes vector not equal to "
                         "itself:\nBefore: {}\nAfter: {}").format(
                          pol_in, pol_in.cartesian.stokes))
if cart_in != cart_in.rot(2*pi):
    failed_tests.append(("Cartesian vector not equal to itself rotated 2pi:\n"
                         "Before:\n{}\nRotated:\n{}").format(cart_in, cart_in.rot(2*pi)))
if pol_in != pol_in.rot(2*pi):
    failed_tests.append(("Stokes vector not equal to itself rotated 2pi:\n"
                         "Before:\n{}\nRotated:\n{}").format(pol_in, pol_in.rot(2*pi)))
if pol_in - pol_in != zero:
    failed_tests.append(("StokesVector subtracted from itself not zero:\n"
                         "Vector:\n{}\nResult:\n{}").format(pol_in, pol_in-pol_in))
if cart_in - cart_in != zero:
    failed_tests.append(("PolarizationTwoVector subtracted from itself not zero:\n"
                         "Vector:\n{}\nResult:\n{}").format(cart_in, cart_in-cart_in))            
if pol_in - cart_in != zero:
    failed_tests.append(("Supposedly equal vectors subtracted from each other not zero:\n"
                         "Stokes:\n{}\nCartesian:\n{}\nResult:\n{}").format(
                         pol_in, cart_in, pol_in-cart_in))
if cart_in + zero != cart_in:
    failed_tests.append(("Adding zero vector modifies value:\n"
                         "Before:\n{}\nAfter:\n{}".format(cart_in, cart_in+zero)))
if pol_in + pol_in_offset != zero:
    failed_tests.append(("Phase destructive Stokes vectors not zero:\n"
                         "Vector:\n{}\nPhase offset vector:\n{}\nResult:\n{}".format(
                         pol_in, pol_in_offset, pol_in+pol_in_offset)))

if len(failed_tests):
    warn("Unit Tests Failed!")
    for failure in failed_tests:
        warn(failure)
else:
    if __name__ == "__main__":
        # If running deliberately, let user know all is good.
        # This way, tests can be included unobtrusively in other modules.
        print "Unit Tests Passed!"
