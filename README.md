# Spacial Dynamics in Microbe Community
Simulation codes for some of the papers used as a basis in my proposition.
I'll put in the codes for making the simulations for each paper, mainly for some visuals.

`cell_sim.py`
This file contains a partially recreated simulation of "cells" which have an area (are rod shaped). Lattice is of open boundary, starts from a 20x20 grid and cells can move elsewhere. They move in discrete time steps and detect cells in proximity of a certain radius of detection. Eligible cells are then processed for the type of cell and the kind of interaction that occurs with the cell are defined, which then creates a varied growth rate/growth arrest. This implements the hertzian forces, cell movement, cell growth, cell division in the simulation. All 6 interactions types are accounted for currently, now I will build the system to do statistical analyses on this system.<br>
Note: I should probably fix the overlap between the bacteria. In my current simulation I see a lot more of bacteria overlapping, I will find a way to plot the actual shapes and observe if they overlap or not.<br>
Note 2: Apparently it was a bad visual simulation. I learnt that matplotlib can scale shapes according to actual interaction boundaries, so I have implemented that instead of line thickness.

I simulated all of them for 300 steps, except competition, where growth was very little, so I did it for 800 steps. The aim is to visualise how the interactions are occuring, shouldn't be used for comparison right now because all of them aren't on the same base/have a common line.
### Mutualism:<br> 
Both Red and Blue benefit each other.<br>
<img width="300" height="300" alt="savedsim" src="https://github.com/user-attachments/assets/7737c551-34a6-4045-a823-a9fdcd3e8436" /><br>
### Parasitism:<br>
Here, Red is the parasite, Blue is the victim.<br>
<img width="300" height="300" alt="savedsim" src="https://github.com/user-attachments/assets/bce1b6c9-354e-4992-8276-68230c84f496" />
### Amensalism:<br>
Red doesn't affect/is neutral to Blue, Blue negatively affects the growth of Red.<br>
<img width="300" height="300" alt="savedsim" src="https://github.com/user-attachments/assets/d537fbd1-bd2f-4fa2-9a81-d94c10e8e21b" /><br>
### Commensalism:<br>
Red is neutral/doesn't affect Blue, Blue positively affects growth of Red.<br>
<img width="300" height="300" alt="savedsim" src="https://github.com/user-attachments/assets/1fe19420-c0b8-4efe-890e-4699a3091edd" /><br>
### Competition:<br>
Both Red and Blue hinder each other's growth.<br>
<img width="300" height="300" alt="savedsim" src="https://github.com/user-attachments/assets/72a8a8d9-0186-4e55-abcd-ffbff2d73f8e" /><br>
### Neutralism (Control):<br> 
Both Red and Blue don't affect each other's growth.<br>
<img width="300" height="300" alt="savedsim" src="https://github.com/user-attachments/assets/5e25dafc-54f4-4200-9147-3185d80195db" /><br>


# Next Steps
What I will do is to figure out how to stop a simulation for a given number of cells instead of a given timestep. Afterwards I will update the use of `plot_history` to account for lineages, fractions and more stuff to properly assess how the model is working. After that is done, I will move on to modify the model to account for chemotaxis, which is accounting for the follow up paper.
