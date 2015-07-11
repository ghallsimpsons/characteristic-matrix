#Characteristic Matrix Solver for Dialectric Media
This package provides methods for finding analytic solutions for electromagnetic
waves in planar dialectric and birefringent media.

Before using this package, make sure to run the test suite to ensure the
current build is functioning properly. You can run the tests from the
command line with `python tests.py`. If you plan on including this package
in your own project, you can include the tests module in your own tests. When
imported in another module, the tests will silently succeed, but will still
complain (loudly) if anything is awry.

##Dependencies
The core modules currently depend upon [NumPy](http://numpy.org) and have
only been tested in Python 2.7.
The included examples also depend upon [SciPy](http://scipy.org) and 
[matplotlib](http://matplotlib.org/).

#Contents
[Design](#design)<br />
[Examples](#examples)<br />
[Classes](#classes)<br />
[Unit Types](#unittypes)<br />

<a name="design" />
#Design
This package is broken into two main modules. The [c_matrix](#c_matrix) module
contains the classes related to dielectric interfaces. The [vector_types](#vector_types)
module provides representations of various vectors that are useful for
polarization sensitive simulations.

In general, an interface matrix is formed by stacking one or more dielectric
(birefringent) layers together. This matrix then acts on a polarization vector
(via multiplication) to produce an output polarization after passing through
the interface.

<a name="examples" />
#Examples

##Anti-Reflection Coating
```python
from c_matrix import Layer, Interface
from examples import graph
layer1=Layer(2, .00035)
layer2=Layer(4, .00025)
layer3=Layer(7, .00018)
layer4=Layer(9.6, .006)
iface=Interface(layer1,layer2,layer3,layer4,layer3,layer2,layer1)
graph(iface)
```

##3-Layer Achromatic Half Wave Plate
```python
from c_matrix import Layer, Interface
from examples import graph, rot_sweep
ar=Layer(3.23, "0.34mm")
hwp=Layer(eps=9.36, eps2=11.56, thickness="3.6mm")
ahwp=Interface(ar,hwp,hwp("58deg"),hwp,ar)
graph(ahwp)
rot_sweep(ahwp)
```

<a name="classes" />
#Classes
<a name="c_matrix" />
##c_matrix
###Layer
The Layer class represents a single dielectric/birefringent layer.
######Arguments
Although the Layer class can be instantiated with one or no arguments, this is
strongly discouraged. Layers can be initialized with the following properties:<br />
`eps`: The relative permittivity along the ordinary axis of the layer.<br />
`thickness (strongly encouraged)`: The thickness of the layer in [Distance Units](#distanceunits).
While the thickness will defaut to a sane value (1/4 wavelength at 124GHz) if not specified,
it is strongly encouraged that you specify the thickness to avoid unexpected behavior.<br />
`eps2 (optional)`: The relative permittivity along the extraordinary axis of the layer. If
unspecified, the layer will have only a single index.<br />
`angle (optional)`: The angle, in [Angle Units](#angleunits), at which to rotate the
ordinary axis of the layer relative to the x-axis.<br />
######Public Methods
<a name="layercall" />
`__call__(angle)`: Because it is undesirable to create multiple identical layers
differing only by the orientation axis, a rotated copy of a layer can be created
simple by calling the layer with the angle at which you would like to rotate.
This method does not modify the existing layer, unless of course you set the
existing layer equal to the rotated representation.

###Interface
The Interface class represents a stack of coincident Layers.
######Arguments
`layers`: An Interface is initialized by passing in an arbitrary number or Layer
objects. The Interface is generated with the Layers in the order by which they
are passed, but a single layer can be passed multiple times, e.g.
`Interface(cracker, chocolate, marshmallow, cracker)`
######Public Methods
`build(frequency)`: The optical properties of an interface in general depend upon
the frequency of the incoming light. You must therefore construct the operator
at a given frequency before acting it on an incoming light vector by calling
the build method.<br />
`__mul__(vect)`: The primary usage of an Interface is by acting it on a light
vector by multiplication. As with all well-behaved multiplication operators,
this returns the new polarization vector, but does not modify the original
incoming vector.<br />
`__call__(angle)`: See [Layer.\__call__](#layercall)

<a name="vector_types" />
##vector_types
The vector_types module provides a few conventient vector representations of
(partially) polarized light. 
<a name="polarizationvector" />
###PolarizationVector
Essentially represents a sine wave with a phase offset. PolarizationVectors
form a complete algebra (addition and scalar multiplication).
######Arguments
Either:<br />
`amplitude`: Unitless amplitude of electric field.<br />
`phase`: Phase angle in [Angle Units](#angleunits)<br />
OR:<br />
`phasor`: Complex number represention the amplitude and phase of the vector.
######Properties
`power`: The square of the electric field.<br />
`amp`: The amplitude of the electric field.<br />
`phase`: The absolute phase of the electric field.<br />
<a name="polarizationtwovector" />
###PolarizationTwoVector
Represents an arbitrary normal-incidence plane wave. PolarizationTwoVectors can
be added together to generate composite waves.
######Arguments
`vector_x`: A [PolarizationVector](#polarizationvector) representing the x
component of field.<br />
`vector_y`: A [PolarizationVector](#polarizationvector) representing the y
component of field.<br />
######Properties
`I`: The intensity of the field.<br />
`Q`: The Q polarized intensity of the field.<br />
`U`: The U polarized intensity of the field.<br />
`V`: The circularly polarized intensity of the field.<br />
`stokes`: The [StokesVector](#stokesvector) representation of the field.<br />
######Public Methods
`rot(angle)`: Returns the vector rotated by the given angle, specified in
[Angle Units](#angleunits)<br />
<a name="stokesvector" />
###StokesVector
The Stokes representation of a normal-incidence plane wave.
######Arguments
Either:<br />
`I`: The intensity of the field.<br />
`Q`: The Q polarized intensity of the field.<br />
`U`: The U polarized intensity of the field.<br />
`V`: The circularly polarized intensity of the field.<br />
`Phase (optional)`: The absolute phase of the [StokesVector](#stokesvector), 
in [Angle Units](#angleunits).<br />
The phase offset is relative to the x component if the x component is non-zero,
or the y component otherwise. It is zero by default.<br />
OR:<br />
`vector_x`: A [PolarizationVector](#polarizationvector) representing the x
component of field.<br />
`vector_y`: A [PolarizationVector](#polarizationvector) representing the y
component of field.<br />
######Properties
`I`: The intensity of the field.<br />
`Q`: The Q polarized intensity of the field.<br />
`U`: The U polarized intensity of the field.<br />
`V`: The circularly polarized intensity of the field.<br />
`cartesian`: The [PolarizationTwoVector](#polarizationtwovector) corresponding
to the [StokesVector](#stokesvector)<br />
`phase`: The phase offset in radians relative to the x component if the
x component is non-zero, or the y component otherwise.<br />
######Public Methods
`rot(angle)`: Returns a [StokesVector](#stokesvector) rotated by the angle
given in [Angle Units](#angleunits)<br />

<a name="unittypes" />
#Unit Types
<a name="angleunits" />
##Angle Units
Any angle can either be specified as an int/float in radians, or in degrees by
passing a string, e.g. "90 deg".
<a name="frequencyunits" />
##Frequency Units
Any frequency can either be specified as an int/float (optionally in scientific
notation, e.g. 120E9) in Hz, or in higher frequencies by passing a string, e.g.
"120 GHz".
<a name="distanceunits" />
##Distance Units
Any distance can either be specified as an int/float (optionally in scientific
notation, e.g. 10E-6) in meters, or in smaller unity by passing a string, e.g.
"30um".

#Acknowledgements
These modules were written with the help of Aritoki Suzuki and Charles Hill at UC Berkeley.
The mathematics draws extensively from Chapter 3 of Tomotake Matsumura's *A Cosmic Microwave Background Radiation Polarimeter Using Superconducting Magnetic Bearings*
