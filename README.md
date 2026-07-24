# Spacial Dynamics in Microbe Community
Simulation codes for some of the papers I am going to present for my proposition work.
I'll put in the codes for making the simulations for each paper, mainly for some visuals.

`cell_sim.py`
This file contains a partially recreated simulation of "cells" which have an area (are rod shaped). Lattice is of open boundary. They move randomly in discrete time steps and detect cells in proximity of a certain radius of detection (here 1.5). This implements the hertzian forces, cell movement, cell growth, cell division in the simulation.

UPDATE: I added two kinds of species with neutral interactions. basically their growth rates do not change anything about the other species growth rates at all. It seems to be working fine, shown below:
Initial State:
<img width="607" height="754" alt="initial_state" src="https://github.com/user-attachments/assets/4d1c6f19-437d-446a-a80a-0cd1665158c4" />

Final State:
<img width="685" height="754" alt="savedsim" src="https://github.com/user-attachments/assets/5b0eb6d0-947d-46a1-966c-8c0a2ff8e158" />


Next steps would be to develop different types of interactions, growth rate influenced by the same, a nutrient gradient, then to chemotactic movement.

Note: I should probably fix the overlap between the bacteria. In my current simulation I see a lot more of bacteria overlapping, I will find a way to plot the actual shapes and observe if they overlap or not.
Note 2: Apparently it was a bad visual simulation. I learnt that matplotlib can scale shapes according to actual interaction boundaries, so I have implemented that instead of line thickness.

