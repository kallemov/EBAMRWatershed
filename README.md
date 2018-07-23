# EWBAMRWatershed
## Author: Baky Kallemov, Lawrence Berkeley National Laboratory

This code implements a model for the simulation of coupled surface and subsurface flow
using the multi-dimensional library CHOMBO.

- Subsurface flow is solved with Richards equation in an 3D embedded boundary geometry with AMR

- Surface flow is represented by the kinematic wave approximation and Manning equations on a 2D
  regular grid with AMR, which follows the underlying surface grid
