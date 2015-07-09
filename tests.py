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

if len(failed_tests):
    warn("Unit Tests Failed!")
    for failure in failed_tests:
        warn(failure)
else:
    if __name__ == "__main__":
        # If running deliberately, let user know all is good.
        # This way, tests can be included unobtrusively in other modules.
        print "Unit Tests Passed!"
