import random
import math
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt

def part_init(n, box_size):
    particles = []
    for i in range(n):
        x = random.uniform(0, box_size)
        y = random.uniform(0, box_size)
        particles.append([x, y])
    return particles

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


def run_sim(n_particles=100, box_size=20, max_dist=1.5, n_steps=30, jitter=0.3):
    particles = part_init(n_particles, box_size)

    print("jitter")
    for step in range(n_steps):
        for p in particles:
            p[0] += random.uniform(-jitter, jitter)
            p[1] += random.uniform(-jitter, jitter)
        pairs = find_pairs(particles, max_dist)

        if step % 5 == 0:
            print(f"  step {step:2d}: {len(pairs)} pairs within {max_dist}")

    print()
    return particles, pairs

def plotting(particles, pairs, filename="savedsim.png"):
    fig, ax = plt.subplots(figsize=(5, 5))

    xs = [p[0] for p in particles]
    ys = [p[1] for p in particles]
    ax.scatter(xs, ys, s=15)

    for (i, j) in pairs:
        x_vals = [particles[i][0], particles[j][0]]
        y_vals = [particles[i][1], particles[j][1]]
        ax.plot(x_vals, y_vals, color="gray", linewidth=0.5)

    ax.set_aspect("equal")
    ax.set_title(f"{len(particles)} particles, {len(pairs)} pairs within range")
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved visualization to", filename)

if __name__ == "__main__":
    particles, pairs = run_sim()
    plotting(particles, pairs)