# Spacial Dynamics in Microbe Community
Simulation codes for some of the papers I am going to present for my proposition work.
I'll put in the codes for making the simulations for each paper, mainly for some visuals.

`cell_sim.py`
This file contains a partially recreated simulation of "cells" which have an area (are rod shaped). Lattice is of open boundary, starts from a 20x20 grid and cells can move elsewhere. They move in discrete time steps and detect cells in proximity of a certain radius of detection. Eligible cells are then processed for the type of cell and the kind of interaction that occurs with the cell are defined, which then creates a varied growth rate/growth arrest. This implements the hertzian forces, cell movement, cell growth, cell division in the simulation. All 6 interactions types are accounted for currently, now I will build the system to do statistical analyses on this system. After that is done, I will move on to modify the model to account for chemotaxis, which is accounting for the follow up paper.

Note: I should probably fix the overlap between the bacteria. In my current simulation I see a lot more of bacteria overlapping, I will find a way to plot the actual shapes and observe if they overlap or not.
Note 2: Apparently it was a bad visual simulation. I learnt that matplotlib can scale shapes according to actual interaction boundaries, so I have implemented that instead of line thickness.

