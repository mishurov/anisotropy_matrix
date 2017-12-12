Anysotropy Matrix SOP node
=====
Anysotropy Matrix SOP node (multithreaded). SideFX® Houdini® Plugin. A base implementation of an anisotropy matrix by paper <b>"Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels", Jihun Yu (Industrial Light and Magic) and Greg Turk (Georgia Institute of Technology)</b>. The plugin doesn't build a scalar density field or a liquid surface, it merely gets points and geometry, calculates an anisotropic transformation matrix for the geometry, and copies it on every point position of the first input.

![ScreenShot1](http://mishurov.co.uk/images/github/anisotropy_matrix/kernels1.png)

## Dependencies
I removed all the external dependencies and use HDK's methods for SVD and multithreading

## Building
There's no more a Makefile and Unux limitations, hcustom works out the box

![ScreenShot2](http://mishurov.co.uk/images/github/anisotropy_matrix/kernels2.png)
