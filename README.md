Characteristic Matrix Solver for Dialectric Media
====================
Methods for finding analytic solutions for waves in dialectric media.



This is not meant to be an extensive library for Characteristic matrices, but an emphasis is placed on modularity, so that these modules may be reused.
There are several missing features which are not difficult to implement, they just require more arguments than were required for my purposes, and so were omitted. If somebody would like to prepare a fully-fledged scipy extension, they would likely need to be added. Examples include:

-Incidence angles: The mathematical methods assume a TE wave incident to the normal of the dialectric.

-Integral methods for continuous strata: Currently, approximations can be formed by finite layering of gradient permittivities, but integral solutions are preferable, both mathematically and computationally.

-Inclusion of arbitrary dialectrics: Current library uses analytic solution for dialectric slab as a starting point, and so is limited in breadth.


The current example:

freqsweep.graph creates a graph of the transmission (neglecting absorption) versus a range of frequencies. 
freqsweep.run_trial optimizes dialectric thickness by displaying bandwidth and relative reflectivity (a dimensionless constant).

EXAMPLE
====================
layer1=layer(2, .00035)
layer2=layer(4, .00025)
layer3=layer(7, .00018)
layer4=layer(9.6, .006)
iface=interface(layer1,layer2,layer3,layer4,layer3,layer2,layer1)
graph(iface)

Acknowledgements
====================
Created for POLARBEAR II project at UC Berkeley under the advisement of Aritoki Suzuki
