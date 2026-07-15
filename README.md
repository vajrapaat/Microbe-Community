# Spacial Dynamics in Microbe Community
Simulation codes for some of the papers I am going to present for my proposition work.
I'll put in the codes for making the simulations for each paper, mainly for some visuals.

`simplesim.py`
This file contains a very simple simulation of "cells" which are points on the lattice. It is of open boundary. They move randomly in discrete time steps and detect cells in proximity of a certain radius of detection (here 1.5). This is just to check if the neighbour detection is working, and the plotting is for confirming the neighbourhood detection and the timestep updates are working fine.

Next steps to this would be to build a periodic boundary (moving out of one edge brings the particle to the parallel edge), add growth rate and a proper shape to the bacteria, add hertzian forces and interactions finally interactions between two kinds.
