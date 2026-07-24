# Spacial Dynamics in Microbe Community
Simulation codes for some of the papers I am going to present for my proposition work.
I'll put in the codes for making the simulations for each paper, mainly for some visuals.

`cell_sim.py`
This file contains a partially recreated simulation of "cells" which have an area (are rod shaped). Lattice is of open boundary. They move randomly in discrete time steps and detect cells in proximity of a certain radius of detection (here 1.5). This implements the hertzian forces, cell movement, cell growth, cell division in the simulation.
Next steps would be to implement two kinds of bacteria present, develop interactions, growth rate influenced by the same and finally a nutrient gradient which can be implemented.

Note: I should probably fix the overlap between the bacteria. In my current simulation I see a lot more of bacteria overlapping, I will find a way to plot the actual shapes and observe if they overlap or not.

