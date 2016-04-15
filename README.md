Anisotropy Matrix SOP node
=====
Anysotropy Matrix SOP node (multithreaded). SideFX® Houdini® Plugin. A base implementation of an anisotropy matrix made for research purposes by paper <b>"Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels", Jihun Yu (Industrial Light and Magic) and Greg Turk (Georgia Institute of Technology)</b>. This plugin doesn't build scalar density field or liquid surface, it just gets points and a geometry, calculates an anisotropic transformation matrix for the geometry, and copies it to an every point position of the first input.

Based on original Anisotropy Matrix SOP by Alexander Mishurov.

## Dependencies
No library dependencies apart from the HDK. Used GCC 4.8 under Linux.

## Building
To compile on Linux, use hcustom. OSX and Windows haven't been tested, but as long as you can build the SOP Star example, this should compile.
