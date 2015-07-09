#Characteristic Matrix Solver for Dialectric Media
Methods for finding analytic solutions for waves in dialectric media.

##Contents
[Design](#design)
[Examples](#examples)
[Classes](#classes)
[Unit Types](#unittypes)

<a name="design" />
##Design
This package is broken into two main modules. The c_matrix module contains
the classes related to dielectric interfaces. The vector_types module provides
representations of various vectors that are useful for polarization sensitive
simulations.

In general, an interface matrix is formed by stacking one or more dielectric
(birefringent) layers together. This matrix then acts on a polarization vector
to produce an output polarization after passing through the interface.

<a name="examples" />
##Examples

###AR COATING EXAMPLE
```python
layer1=Layer(2, .00035)
layer2=Layer(4, .00025)
layer3=Layer(7, .00018)
layer4=Layer(9.6, .006)
iface=Interface(layer1,layer2,layer3,layer4,layer3,layer2,layer1)
graph(iface)
```

###3-Layer Achromatic Half Wave Plate Example
```python
AR1=Layer(2, "350um")
AR2=Layer(5, "250um")
HWP=Layer(eps=9.36, eps2=11.56, thickness="3.6mm")
iface=Interface(AR1,AR2,HWP,HWP("58deg"),HWP,AR2,AR1)
graph(iface)
```

<a name="classes" />
##Classes
###c_matrix
####Layer
The Layer class represents a single dielectric/birefringent layer.
#####Arguments
Although the Layer class can be instantiated with one or no arguments, this is
strongly discouraged. Layers can be initialized with the following properties:
`eps`: The relative permittivity along the ordinary axis of the layer.
`thickness`: The thickness of the layer in [Distance Units](#distanceunits)
`eps2`: The relative permittivity along the extraordinary axis of the layer. If
unspecified, the layer will has a single index.
`angle`: The angle, in [Angle Units](#angleunits), at which to rotate the ordinary axis of the layer
from the x-axis.
for a 90 degree rotation.
#####Public Methods
<a name="layercall" />
`__call__(angle)`: Because it is undesirable to create multiple identical layers
differing only by the orientation axis, a rotated copy of a layer can be created
simple by calling the layer with the angle at which you would like to rotate.
This method does not modify the existing layer, unless of course you set the
existing layer equal to the rotated representation.

####Interface
The Interface class represents a stack of coincident Layers.
#####Arguments
`layers`: An Interface is initialized by passing in an arbitrary number or Layer
objects. The Interface is generated with the Layers in the order by which they
are passed, but a single layer can be passed multiple times, e.g. 
`Interface(cracker, chocolate, marshmallow, cracker)`
#####Public Methods
`build(frequency)`: The optical properties of an interface in general depend upon
the frequency of the incoming light. You must therefore construct the operator
at a given frequency before acting it on an incoming light vector by calling
the build method.
`__mul__(vect)`: The primary usage of an Interface is by acting it on a light
vector by multiplication. As with all well-behaved multiplication operators,
this returns the new polarization vector, but does not modify the original
incoming vector.
`__call__(angle)`: See [Layer.__call__](#layercall)

###vector_types
The vector_types module provides a few conventient vector representations of
(partially) polarized light. 
####PolarizationVector
Essentially represents a sine wave with a phase offset. PolarizationVectors
form a complete algebra (addition and scalar multiplication).
#####Arguments
Either:
`amplitude`: Unitless amplitude of electric field.
`phase`: Phase angle in [angle units](#angleunits)
OR:
`phasor`: Complex number represention the amplitude and phase of the vector.

<a name="unittypes" />
##Unit Types
<a name="angleunits" />
###Angle Units
Any angle can either be specified as an int/float in radians, or in degrees by
passing a string, e.g. "90 deg".
<a name="frequencyunits" />
###Frequency Units
Any frequency can either be specified as an int/float (optionally in scientific
notation, e.g. 120E9) in Hz, or in higher frequencies by passing a string, e.g.
"120 GHz".
<a name="distanceunits" />
###Distance Units
Any distance can either be specified as an int/float (optionally in scientific
notation, e.g. 10E-6) in meters, or in smaller unity by passing a string, e.g.
"30um".

#Acknowledgements
These modules were written with the help of Aritoki Suzuki and Charles Hill at UC Berkeley.
The mathematics draws extensively from Chapter 3 of Tomotake Matsumura's *A Cosmic Microwave Background Radiation Polarimeter Using Superconducting Magnetic Bearings*
