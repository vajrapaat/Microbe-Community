# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 17:15:02 2026

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def lv_dyn_allee(t, y, r1, r2, K1, K2, alpha12, alpha21, a1, d):
    X1, X2 = y
    dX1_dt = r1 * X1 * ((K1 - X1) * (X1 - a1) - alpha12 * X2) / K1 - d * X1
    dX2_dt = r2 * X2 * (K2 - X2 - alpha21 * X1) / K2 - d * X2
    return [dX1_dt, dX2_dt]

t_span = (0, 1000)
t_eval = np.linspace(0, 1000, 10000)

initial_conditions = []
grid_points = np.linspace(0.05, 0.95, 8)
for x0 in grid_points:
    for y0 in grid_points:
        initial_conditions.append([x0, y0])

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

pars_sim1 = {
    'r1': 0.3, 'r2': 0.3,
    'K1': 1.0, 'K2': 1.0,
    'alpha12': 0.8, 'alpha21': 1.1,
    'a1': -0.3, 'd': 0.0
}

axes[0].set_title("Simulation 1: Strong Bistability\n(Allee Effect present)", fontsize=12, fontweight='bold')

for yini in initial_conditions:
    sol = solve_ivp(lv_dyn_allee, t_span, yini, t_eval=t_eval, args=tuple(pars_sim1.values()), method='RK45')
    axes[0].plot(sol.y[0], sol.y[1], color='black', alpha=0.3, linewidth=1)
    axes[0].plot(yini[0], yini[1], 'o', color='#3C5CC2', alpha=0.5, markersize=3)

x_vals = np.linspace(0.001, 1.2, 200)
null1_sim1 = ((pars_sim1['K1'] - x_vals) * (x_vals - pars_sim1['a1'])) / pars_sim1['alpha12']
null2_sim1 = pars_sim1['K2'] - pars_sim1['alpha21'] * x_vals

axes[0].plot(x_vals, null1_sim1, '--', color='#a0ced9', linewidth=2.5, label='Species 1 Nullcline')
axes[0].plot(x_vals, null2_sim1, '--', color='#f9ac95', linewidth=2.5, label='Species 2 Nullcline')

pars_sim2 = pars_sim1.copy()
pars_sim2['a1'] = 0.1  

axes[1].set_title("Simulation 2: Glutamate Effect\n(Weakened/Eliminated Bistability)", fontsize=12, fontweight='bold')

for yini in initial_conditions:
    sol = solve_ivp(lv_dyn_allee, t_span, yini, t_eval=t_eval, args=tuple(pars_sim2.values()), method='RK45')
    axes[1].plot(sol.y[0], sol.y[1], color='black', alpha=0.3, linewidth=1)
    axes[1].plot(yini[0], yini[1], 'o', color='#3C5CC2', alpha=0.5, markersize=3)

null1_sim2 = ((pars_sim2['K1'] - x_vals) * (x_vals - pars_sim2['a1'])) / pars_sim2['alpha12']
null2_sim2 = pars_sim2['K2'] - pars_sim2['alpha21'] * x_vals

axes[1].plot(x_vals, null1_sim2, '--', color='#a0ced9', linewidth=2.5)
axes[1].plot(x_vals, null2_sim2, '--', color='#f9ac95', linewidth=2.5)

for ax in axes:
    ax.set_xlabel("Abundance Species 1 ($X_1$ - Se / Turquoise)", fontsize=11)
    ax.set_ylabel("Abundance Species 2 ($X_2$ - Cn / Orange)", fontsize=11)
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 1.2)
    ax.grid(True, linestyle=':', alpha=0.6)

axes[0].legend(loc='upper right')
plt.tight_layout()
plt.show()