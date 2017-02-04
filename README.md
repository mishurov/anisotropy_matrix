Anysotropy Matrix SOP node
=====
Anysotropy Matrix SOP node (multithreaded). SideFX® Houdini® Plugin. A base implementation of an anisotropy matrix made for research purposes by paper <b>"Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels", Jihun Yu (Industrial Light and Magic) and Greg Turk (Georgia Institute of Technology)</b>. This plugin doesn't build scalar density field or liquid surface, it just gets points and a geometry, calculates an anisotropic transformation matrix for the geometry, and copies it to an every point position of the first input.

![ScreenShot1](http://mishurov.5gbfree.com/github/anisotropy_matrix/kernels1.png)

## Dependencies
In order to compile plugin you need to link it with <a href="http://www.netlib.org/lapack/">LAPACK</a> library and install Boost LAPACK bindings from <a href="http://mathema.tician.de/software/boost-numeric-bindings/">Numeric Library Bindings for Boost UBlas</a> into an include directory of Boost in Houdini's toolkit include folder. That's neccessary to compute Singular Value Decomposition to get scale and rotation matrices. Also OpenMP library is neccessary for multithreading.

## Building
To compile on Linux, use Makefile. OSX and Windows are not tested.

![ScreenShot2](http://mishurov.5gbfree.com/github/anisotropy_matrix/kernels2.png)
