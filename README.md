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


Known bugs:

Frequencies less than one must include a leading 0 due to inproper parsing.



TODO:

Fix the units on axis (currently it only graphs in Hz)

Make more modular so it can do more different number of layers with a single command.

Subroutines to focus permittivity and thickness parameters to return the desired central frequency, maximize the the spread of the +- peak nodes, and/or minimize rms within a given range.
