Anysotropy Matrix SOP node
=====
Anysotropy Matrix SOP node. SideFX® Houdini® Plugin. A base implementation of an anisotropy matrix made for research purposes by paper "Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels" from Jihun Yu (Industrial Light and Magic) and Greg Turk (Georgia Institute of Technology). This plugin doesn't build scalar density field or liquid surface, it's just get points and geometry then transforms geometry by the anisotropy matrix and copy it to an every point position.

![ScreenShot1](https://dl.dropboxusercontent.com/u/20988720/CG/anisotropy_matrix/kernels1.png)

## Dependencies
In order to compile plugin you need to link it with static <a href="http://www.netlib.org/lapack/">LAPACK</a> library and install Boost LAPACK bindings form <a href="http://mathema.tician.de/software/boost-numeric-bindings/">Numeric Library Bindings for Boost UBlas</a> into an include directory of Boost in Houdini's toolkit include folder. That's neccessary to compute Singular Value Decomposition to get scale and rotation matrices.

## Building
I'd compiled the plugin with hcustom.
```bash
hcustom SOP_anisotropy_matrix.cc -llapack
```

![ScreenShot2](https://dl.dropboxusercontent.com/u/20988720/CG/anisotropy_matrix/kernels2.png)
