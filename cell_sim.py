import random
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon


DIAMETER = 1.0         
G = 1.0                 
diffradius = 10.0
growthfact = 2.0
KAPPA = 0.333           
INTERACTION_STRENGTH = 1800.0  
DT = 0.01               
Color = {0:'red', 1:'steelblue'}
tmax = 0.1
nnorm = 20.0

interactions = {
    'neutral': (0.0,0.0),
    'commensalism':(-10,0),
    'amensalism': (10,0.0),
    'mutualism':(-10,-10),
    'competition': (10,10),
    'parasitism':(-10,10)
    }

def division_len():
    l = random.gauss(4.0, 0.3)
    return min(max(l, 3.1), 4.9)

def part_init(n, box_size):
    cells = []
    for i in range(n):
        x = random.uniform(0, box_size)
        y = random.uniform(0, box_size)
        angle = random.uniform(0, 2 * math.pi)
        length = random.uniform(2.0, 3.0)
        species = i % 2
        cells.append([x, y, angle, length, division_len(), species])
    return cells

def capsule(ax, cell, alpha=1):
    x, y, angle, length, species = cell[0], cell[1], cell[2], cell[3], cell[5]
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


def distance(p1, p2):
    dx2 = (p1[0] - p2[0])**2
    dy2 = (p1[1] - p2[1])**2
    return math.sqrt(dx2 + dy2)


def find_pairs(particles, max_dist):
    pair = []
    n = len(particles)
    for i in range(n):
        for j in range(i + 1, n):
            if distance(particles[i], particles[j]) < max_dist:
                pair.append((i, j))
    return pair

def endpoints(cell):
    x, y, angle, length = cell[0], cell[1], cell[2], cell[3]
    half = length / 2.0
    dx = math.cos(angle) * half
    dy = math.sin(angle) * half
    return (x - dx, y - dy), (x + dx, y + dy)

def cell_area(cell):
    length = cell[3]
    r = DIAMETER / 2.0
    return DIAMETER * length + math.pi * r ** 2

def closest(p1, q1, p2, q2, eps=1e-9):
    p1x, p1y = p1
    q1x, q1y = q1
    p2x, p2y = p2
    q2x, q2y = q2

    d1x, d1y = q1x - p1x, q1y - p1y
    d2x, d2y = q2x - p2x, q2y - p2y
    rx, ry = p1x - p2x, p1y - p2y

    a = d1x * d1x + d1y * d1y
    e = d2x * d2x + d2y * d2y
    f = d2x * rx + d2y * ry

    if a <= eps and e <= eps:
        return (p1x, p1y), (p2x, p2y), math.hypot(p1x - p2x, p1y - p2y)

    if a <= eps:
        s = 0.0
        t = min(max(f / e, 0.0), 1.0)
    else:
        c = d1x * rx + d1y * ry
        if e <= eps:
            t = 0.0
            s = min(max(-c / a, 0.0), 1.0)
        else:
            b = d1x * d2x + d1y * d2y
            denom = a * e - b * b
            s = min(max((b * f - c * e) / denom, 0.0), 1.0) if denom > eps else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = min(max(-c / a, 0.0), 1.0)
            elif t > 1.0:
                t = 1.0
                s = min(max((b - c) / a, 0.0), 1.0)

    c1x, c1y = p1x + d1x * s, p1y + d1y * s
    c2x, c2y = p2x + d2x * t, p2y + d2y * t
    return (c1x, c1y), (c2x, c2y), math.hypot(c1x - c2x, c1y - c2y)


def others(cells,radius):
    centers = [[c[0],c[1]] for c in cells]
    pairs = find_pairs(centers, radius)
    counts = [0]*len(cells)
    for i, j in pairs:
        if cells[i][5] != cells[j][5]:
            counts[i] += 1
            counts[j] += 1
    return counts
            
def force(cells):
    n = len(cells)
    forces = [[0.0, 0.0] for _ in range(n)]
    torques = [0.0 for _ in range(n)]
    centers = [[cell[0], cell[1]] for cell in cells]
    max_len = max((cell[3] for cell in cells), default=0.0)
    broad_dist = max_len + DIAMETER
    candidate_pairs = find_pairs(centers, broad_dist)

    contacts = []
    for (i, j) in candidate_pairs:
        a0, a1 = endpoints(cells[i])
        b0, b1 = endpoints(cells[j])
        pa, pb, dist = closest(a0, a1, b0, b1)

        h = DIAMETER - dist
        if h <= 0:
            continue
        contacts.append((i, j))

        sepx, sepy = pa[0] - pb[0], pa[1] - pb[1]
        norm = math.hypot(sepx, sepy)
        if norm < 1e-9:
            sepx, sepy = random.uniform(-1, 1), random.uniform(-1, 1)
            norm = math.hypot(sepx, sepy)
        dirx, diry = sepx / norm, sepy / norm

        mag = INTERACTION_STRENGTH * (h ** 1.5)
        li, lj = cells[i][3], cells[j][3]
        fi, fj = mag / li, mag / lj

        forces[i][0] += fi * dirx
        forces[i][1] += fi * diry
        forces[j][0] -= fj * dirx
        forces[j][1] -= fj * diry
        lever_ix, lever_iy = pa[0] - cells[i][0], pa[1] - cells[i][1]
        lever_jx, lever_jy = pb[0] - cells[j][0], pb[1] - cells[j][1]
        torques[i] += (12 * fi / li ** 2) * (lever_ix * diry - lever_iy * dirx)
        torques[j] += (12 * fj / lj ** 2) * (lever_jx * (-diry) - lever_jy * (-dirx))

    return forces, torques, contacts


def integrate(cells, forces, torques, dt):
    for i, cell in enumerate(cells):
        cell[0] += forces[i][0] * dt
        cell[1] += forces[i][1] * dt
        cell[2] += torques[i] * dt

def grow_cells(cells, dt, xi1 = 0.0, xi2 = 0.0, n_nutrient=1.0):
    othercn = others(cells, diffradius)
    for i, cell in enumerate(cells):
        area = cell_area(cell)
        basedl = G * area * n_nutrient / (KAPPA + n_nutrient) * dt
        xi = xi1 if cell[5] == 0 else xi2
        T = min((othercn[i]/nnorm), tmax)
        factor = min(max(1 - xi*T,0.0),growthfact)
        cell[3] += basedl*factor*dt

def divide_cells(cells):
    survivors = []
    new_cells = []

    for cell in cells:
        x, y, angle, length, div_len, species = cell
        if length < div_len:
            survivors.append(cell)
            continue

        frac = random.uniform(0.4, 0.6)
        l1 = length * frac
        l2 = length - l1

        dirx, diry = math.cos(angle), math.sin(angle)
        offset1 = length / 2.0 - l1 / 2.0
        offset2 = length / 2.0 - l2 / 2.0

        daughter1 = [x - offset1 * dirx, y - offset1 * diry, angle, l1, division_len(), species]
        daughter2 = [x + offset2 * dirx, y + offset2 * diry, angle, l2, division_len(), species]
        new_cells.append(daughter1)
        new_cells.append(daughter2)

    return survivors + new_cells


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


def run_sim(n_particles=20, box_size=20, n_steps=300, interaction_type="control"):
    cells = part_init(n_particles, box_size)
    plotting(cells, filename="initial_state.png")

    xi1, xi2 = interactions[interaction_type]
    history = []

    for step in range(n_steps):
        forces, torques, contacts = force(cells)
        integrate(cells, forces, torques, DT)
        grow_cells(cells, DT, xi1, xi2)
        cells = divide_cells(cells)

        n_green = sum(1 for c in cells if c[5] == 1)
        history.append(n_green / len(cells) if cells else 0.0)

        if step % 10 == 0:
            print(f"  step {step:3d}: {len(cells)} cells, {len(contacts)} in contact")

    print()
    return cells, contacts, history

def plotting(cells, filename="savedsim.png"):
    fig, ax = plt.subplots(figsize=(6, 6))

    for cell in cells:
        capsule(ax, cell)

    if cells:
        xs = [c[0] for c in cells]
        ys = [c[1] for c in cells]
        max_len = max(c[3] for c in cells)
        pad = max_len / 2.0 + DIAMETER
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
    
    #plt.grid(True, 'minor', 'both')
    ax.set_aspect("equal")
    ax.set_title(f"{len(cells)} cells")
    ax.axis("off")
    #plt.grid(True, 'minor', 'both')
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved visualization to", filename)

if __name__ == "__main__":
    cells, contacts, history = run_sim(interaction_type="parasitism")
    plotting(cells)
    plot_history(history)