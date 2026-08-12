import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import sim_vectorized as sv

OUTDIR = "sim_output"
N_PARTICLES = 64
BOX_SIZE = 12
MAXIND = 10000
N_STEPS = 1000000
CHECKPOINT_EVERY = 500
RESUME = True


def run_one(interaction_type, seed):
    np.random.seed(seed)
    os.makedirs(OUTDIR, exist_ok=True)

    cells, contacts, history = sv.run_sim(
        n_particles=N_PARTICLES,
        box_size=BOX_SIZE,
        maxind=MAXIND,
        n_steps=N_STEPS,
        interaction_type=interaction_type,
        outdir=OUTDIR,
        checkpoint_every=CHECKPOINT_EVERY,
        resume=RESUME,
    )

    sv.plotting(cells, filename=os.path.join(OUTDIR, f"savedsim_{interaction_type}.png"))
    sv.plot_history(history, filename=os.path.join(OUTDIR, f"fracgreen{interaction_type}.png"))

    return interaction_type, sv.n_cells(cells)


if __name__ == "__main__":
    interaction_types = list(sv.interactions.keys())

    with ProcessPoolExecutor(max_workers=6) as executor:
        futures = {
            executor.submit(run_one, itype, seed): itype
            for seed, itype in enumerate(interaction_types)
        }
        for future in as_completed(futures):
            itype = futures[future]
            try:
                name, final_n = future.result()
                print(f"[done] {name}: {final_n} cells")
            except Exception as exc:
                print(f"[FAILED] {itype}: {exc!r}")
