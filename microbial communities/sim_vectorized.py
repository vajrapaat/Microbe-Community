import random
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
import diffusion


DIAMETER = 1.0
G = 1.0
growthfact = 2.0
KAPPA = 0.333
INTERACTION_STRENGTH = 1300.0
DT = 0.001
Color = {0: 'red', 1: 'seagreen'}

interactions = {
    'neutral': (0.0, 0.0),
    'commensalism': (0.0, -10),
    'amensalism': (0.0, 10),
    'mutualism': (-10, -10),
    'competition': (10, 10),
    'parasitism': (10, -10)
}


def division_len(n=1):
    l = np.random.normal(4.0, 0.3, size=n)
    return np.clip(l, 3.1, 4.9)


def part_init(n, box_size):
    half = box_size / 2.0
    x = np.random.uniform(-half, half, n)
    y = np.random.uniform(-half, half, n)
    angle = np.random.uniform(0, 2 * math.pi, n)
    length = np.random.uniform(2.0, 3.0, n)
    div_len = division_len(n)
    species = (np.arange(n) % 2).astype(np.int64)
    return dict(x=x, y=y, angle=angle, length=length, div_len=div_len, species=species)


def n_cells(cells):
    return len(cells['x'])


def capsule(ax, cells, i, alpha=1):
    x, y, angle, length, species = (
        cells['x'][i], cells['y'][i], cells['angle'][i],
        cells['length'][i], cells['species'][i]
    )
    color = Color[species]
    r = DIAMETER / 2.0
    half = length / 2.0
    dirx, diry = math.cos(angle), math.sin(angle)
    perp_x, perp_y = -diry, dirx

    p1 = (x - half * dirx + r * perp_x, y - half * diry + r * perp_y)
    p2 = (x + half * dirx + r * perp_x, y + half * diry + r * perp_y)
    p3 = (x + half * dirx - r * perp_x, y + half * diry - r * perp_y)
    p4 = (x - half * dirx - r * perp_x, y - half * diry - r * perp_y)
    ax.add_patch(Polygon([p1, p2, p3, p4], closed=True, color=color, alpha=alpha, linewidth=0))

    ax.add_patch(Circle((x - half * dirx, y - half * diry), r, color=color, alpha=alpha, linewidth=0))
    ax.add_patch(Circle((x + half * dirx, y + half * diry), r, color=color, alpha=alpha, linewidth=0))


def endpoints(cells):
    half = cells['length'] / 2.0
    dx = np.cos(cells['angle']) * half
    dy = np.sin(cells['angle']) * half
    p0 = np.stack([cells['x'] - dx, cells['y'] - dy], axis=1)
    p1 = np.stack([cells['x'] + dx, cells['y'] + dy], axis=1)
    return p0, p1





def cell_area(cells):
    r = DIAMETER / 2.0
    return DIAMETER * cells['length'] + math.pi * r ** 2


def find_pairs(cells, max_dist):
    x, y = cells['x'], cells['y']
    n = len(x)
    ii, jj = np.triu_indices(n, k=1)
    if len(ii) == 0:
        return ii, jj
    dx = x[ii] - x[jj]
    dy = y[ii] - y[jj]
    dist = np.hypot(dx, dy)
    keep = dist < max_dist
    return ii[keep], jj[keep]


def closest_batch(p1, q1, p2, q2, eps=1e-9):
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2

    a = np.einsum('ij,ij->i', d1, d1)
    e = np.einsum('ij,ij->i', d2, d2)
    f = np.einsum('ij,ij->i', d2, r)
    c = np.einsum('ij,ij->i', d1, r)
    b = np.einsum('ij,ij->i', d1, d2)
    denom = a * e - b * b

    s = np.where(denom > eps, np.clip((b * f - c * e) / np.where(denom > eps, denom, 1.0), 0.0, 1.0), 0.0)
    t = (b * s + f) / e

    below = t < 0.0
    above = t > 1.0
    t = np.clip(t, 0.0, 1.0)

    s_below = np.clip(-c / a, 0.0, 1.0)
    s_above = np.clip((b - c) / a, 0.0, 1.0)
    s = np.where(below, s_below, np.where(above, s_above, s))

    c1 = p1 + d1 * s[:, None]
    c2 = p2 + d2 * t[:, None]
    dist = np.linalg.norm(c1 - c2, axis=1)
    return c1, c2, dist


def force(cells):
    n = n_cells(cells)
    forces = np.zeros((n, 2))
    torques = np.zeros(n)

    if n == 0:
        return forces, torques, []

    max_len = cells['length'].max()
    broad_dist = max_len + DIAMETER
    ii, jj = find_pairs(cells, broad_dist)
    if len(ii) == 0:
        return forces, torques, []

    p0, p1 = endpoints(cells)
    pa, pb, dist = closest_batch(p0[ii], p1[ii], p0[jj], p1[jj])

    h = DIAMETER - dist
    contact_mask = h > 0
    if not contact_mask.any():
        return forces, torques, []

    ii = ii[contact_mask]
    jj = jj[contact_mask]
    pa = pa[contact_mask]
    pb = pb[contact_mask]
    h = h[contact_mask]

    sep = pa - pb
    norm = np.linalg.norm(sep, axis=1)
    zero = norm < 1e-9
    if zero.any():
        rand_dir = np.random.uniform(-1, 1, size=(zero.sum(), 2))
        sep = sep.copy()
        sep[zero] = rand_dir
        norm = norm.copy()
        norm[zero] = np.linalg.norm(rand_dir, axis=1)
    direction = sep / norm[:, None]

    mag = INTERACTION_STRENGTH * (h ** 1.5)
    li = cells['length'][ii]
    lj = cells['length'][jj]
    fi = mag / li
    fj = mag / lj

    force_i = fi[:, None] * direction
    force_j = -fj[:, None] * direction

    np.add.at(forces, ii, force_i)
    np.add.at(forces, jj, force_j)

    lever_i = pa - np.stack([cells['x'][ii], cells['y'][ii]], axis=1)
    lever_j = pb - np.stack([cells['x'][jj], cells['y'][jj]], axis=1)

    torque_i = (12 * fi / li ** 2) * (lever_i[:, 0] * direction[:, 1] - lever_i[:, 1] * direction[:, 0])
    torque_j = (12 * fj / lj ** 2) * (lever_j[:, 0] * (-direction[:, 1]) - lever_j[:, 1] * (-direction[:, 0]))

    np.add.at(torques, ii, torque_i)
    np.add.at(torques, jj, torque_j)

    contacts = list(zip(ii.tolist(), jj.tolist()))
    return forces, torques, contacts


def integrate(cells, forces, torques, dt):
    cells['x'] = cells['x'] + forces[:, 0] * dt
    cells['y'] = cells['y'] + forces[:, 1] * dt
    cells['angle'] = cells['angle'] + torques * dt


def grow_cells(cells, dt, xi1, xi2, n_field, chem_from_red, chem_from_green):
    area = cell_area(cells)
    x, y, species = cells['x'], cells['y'], cells['species']

    n_local = diffusion.lookup(n_field, x, y)
    basedl = G * area * n_local / (KAPPA + n_local) * dt

    is_red = species == 0
    T_if_red = diffusion.lookup(chem_from_green, x, y)
    T_if_green = diffusion.lookup(chem_from_red, x, y)
    T = np.where(is_red, T_if_red, T_if_green)
    xi = np.where(is_red, xi1, xi2)

    factor = np.clip(1 - xi * T, 0.0, growthfact)
    cells['length'] = cells['length'] + basedl * factor


def divide_cells(cells):
    length = cells['length']
    div_len = cells['div_len']
    dividing = length >= div_len
    n_div = int(dividing.sum())

    survivors = {k: v[~dividing] for k, v in cells.items()}

    if n_div == 0:
        return survivors

    x = cells['x'][dividing]
    y = cells['y'][dividing]
    angle = cells['angle'][dividing]
    length_d = length[dividing]
    species = cells['species'][dividing]

    frac = np.random.uniform(0.4, 0.6, size=n_div)
    l1 = length_d * frac
    l2 = length_d - l1

    dirx, diry = np.cos(angle), np.sin(angle)
    offset1 = length_d / 2.0 - l1 / 2.0
    offset2 = length_d / 2.0 - l2 / 2.0

    daughters = {
        'x': np.concatenate([x - offset1 * dirx, x + offset2 * dirx]),
        'y': np.concatenate([y - offset1 * diry, y + offset2 * diry]),
        'angle': np.concatenate([angle, angle]),
        'length': np.concatenate([l1, l2]),
        'div_len': division_len(2 * n_div),
        'species': np.concatenate([species, species]),
    }

    return {k: np.concatenate([survivors[k], daughters[k]]) for k in cells}


def plot_history(history, filename="fraction_green_history.png"):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(history, color="seagreen")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("step")
    ax.set_ylabel("fraction green cells")
    ax.set_ylim(0, 1)
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved visualization to", filename)


def plotting(cells, filename="savedsim.png"):
    fig, ax = plt.subplots(figsize=(6, 6))
    n = n_cells(cells)

    for i in range(n):
        capsule(ax, cells, i)

    if n:
        xs, ys = cells['x'], cells['y']
        max_len = cells['length'].max()
        pad = max_len / 2.0 + DIAMETER
        ax.set_xlim(xs.min() - pad, xs.max() + pad)
        ax.set_ylim(ys.min() - pad, ys.max() + pad)

    ax.set_aspect("equal")
    ax.set_title(f"{n} cells")
    ax.axis("off")
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved visualization to", filename)


def checkpoint_path(outdir, interaction_type):
    return os.path.join(outdir, f"checkpoint_{interaction_type}.npz")


def save_checkpoint(path, cells, step, history, n_field, chem_from_red, chem_from_green):
    base, ext = os.path.splitext(path)
    tmp_path = base + ".tmp" + ext
    np.savez(
        tmp_path,
        step=step,
        history=np.array(history),
        x=cells['x'], y=cells['y'], angle=cells['angle'],
        length=cells['length'], div_len=cells['div_len'], species=cells['species'],
        n_field=n_field, chem_from_red=chem_from_red, chem_from_green=chem_from_green,
    )
    os.replace(tmp_path, path)


def load_checkpoint(path):
    data = np.load(path)
    cells = dict(
        x=data['x'], y=data['y'], angle=data['angle'],
        length=data['length'], div_len=data['div_len'], species=data['species'],
    )
    step = int(data['step'])
    history = list(data['history'])
    n_field = data['n_field']
    chem_from_red = data['chem_from_red']
    chem_from_green = data['chem_from_green']
    return cells, step, history, n_field, chem_from_red, chem_from_green


def run_sim(n_particles=64, box_size=12, maxind=10000, n_steps=1000000, interaction_type="control",
            outdir=".", checkpoint_every=None, resume=False):
    xi1, xi2 = interactions[interaction_type]
    ckpt_path = checkpoint_path(outdir, interaction_type)

    if resume and os.path.exists(ckpt_path):
        cells, step, history, n_field, chem_from_red, chem_from_green = load_checkpoint(ckpt_path)
        print(f"[{interaction_type}] resumed from step {step} ({n_cells(cells)} cells)")
    else:
        cells = part_init(n_particles, box_size)
        plotting(cells, filename=os.path.join(outdir, f"initial_state_{interaction_type}.png"))
        history = []
        n_field = diffusion.nutrientfield()
        chem_from_red = diffusion.chemmfield()
        chem_from_green = diffusion.chemmfield()
        step = 0

    contacts = []
    while n_cells(cells) < maxind and step < n_steps:
        forces, torques, contacts = force(cells)
        integrate(cells, forces, torques, DT)

        areas = cell_area(cells)
        rho_red = diffusion.density(cells, areas, specfilt=0)
        rho_green = diffusion.density(cells, areas, specfilt=1)
        rho_total = rho_red + rho_green
        n_field = diffusion.nutrientstep(n_field, rho_total, DT)
        chem_from_red = diffusion.chemmstep(chem_from_red, rho_red, DT)
        chem_from_green = diffusion.chemmstep(chem_from_green, rho_green, DT)

        grow_cells(cells, DT, xi1, xi2, n_field, chem_from_red, chem_from_green)
        cells = divide_cells(cells)

        n = n_cells(cells)
        n_green = int(np.sum(cells['species'] == 1)) if n else 0
        history.append(n_green / n if n else 0.0)

        step += 1

        if step % 10 == 0:
            print(f"[{interaction_type}] step {step}: {n} cells")

        if checkpoint_every and step % checkpoint_every == 0:
            save_checkpoint(ckpt_path, cells, step, history, n_field, chem_from_red, chem_from_green)

    if checkpoint_every:
        save_checkpoint(ckpt_path, cells, step, history, n_field, chem_from_red, chem_from_green)

    print()
    return cells, contacts, history


if __name__ == "__main__":
    cells, contacts, history = run_sim(interaction_type="competition")
    plotting(cells)
    plot_history(history)
