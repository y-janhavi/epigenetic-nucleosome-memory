import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------
# Simulation Parameters
# --------------------------------------------
N = 60                             # Number of nucleosomes
F_values = np.logspace(-1, 2, 30)  # Feedback-to-noise ratios from 0.1 to 100 (log scale)
Tmax = 30000*N                     # Number of simulation steps
Teq = 10 * N                       # Equilibration time before measurement
num_runs = 10                      # Number of independent simulation runs

# Nucleosome states: A = Acetylated, U = Unmodified, M = Methylated
A, U, M = 0, 1, 2

# --------------------------------------------
# Model A: Standard model (no spatial constraint)
# --------------------------------------------
def standard(state, F):
    """Standard model: Random nucleosomes can recruit from anywhere in the array."""
    α = F / (F + 1)
    N = len(state)
    n1 = np.random.randint(N)  # Target nucleosome

    if np.random.rand() < α:
        n2 = np.random.randint(N)  # Recruiter chosen randomly from entire region
        if state[n2] != U:
            if state[n2] == M:
                if state[n1] == A:
                    state[n1] = U
                elif state[n1] == U:
                    state[n1] = M
            elif state[n2] == A:
                if state[n1] == M:
                    state[n1] = U
                elif state[n1] == U:
                    state[n1] = A
    else:
        # Noisy conversion (equal chance between accessible transitions)
        r = np.random.rand()
        if state[n1] == A and r < 1/3: state[n1] = U
        elif state[n1] == U:
            if r < 1/3: state[n1] = A
            elif r < 2/3: state[n1] = M
        elif state[n1] == M and r < 1/3: state[n1] = U
    return state

# --------------------------------------------
# Model B: Neighbor-limited recruitment (local only)
# --------------------------------------------
def neighbour_limited(state, F):
    """Only neighboring nucleosomes (left/right) can recruit."""
    α = F / (F + 1)
    N = len(state)
    n1 = np.random.randint(N)

    if np.random.rand() < α:
        # Recruit from either left or right neighbor
        n2 = (n1 - 1) % N if np.random.rand() < 0.5 else (n1 + 1) % N

        if state[n2] != U:
            if state[n2] == M:
                if state[n1] == A:
                    state[n1] = U
                elif state[n1] == U:
                    state[n1] = M
            elif state[n2] == A:
                if state[n1] == M:
                    state[n1] = U
                elif state[n1] == U:
                    state[n1] = A
    else:
        return standard(state, 0)  # Use noise update from standard model
    return state

# --------------------------------------------
# Model C: Power-law contact probability
# Simulates 3D chromatin folding where distant nucleosomes can contact
# --------------------------------------------
def power_law_contact(state, F):
    """Recruitment with probability ~1/distance^1.5 to mimic 3D chromatin contacts."""
    α = F / (F + 1)
    N = len(state)
    n1 = np.random.randint(N)

    if np.random.rand() < α:
        distances = np.arange(1, N)        # All distances from 1 to N−1
        probs = 1 / distances**1.5         # Power-law decay
        probs /= probs.sum()              # Normalize to make a probability distribution

        d = np.random.choice(distances, p=probs)  # Pick distance using power-law
        # Decide left or right direction with equal chance
        n2 = (n1 + d) % N if np.random.rand() < 0.5 else (n1 - d) % N

        # Recruited conversion
        if state[n2] != U:
            if state[n2] == M:
                if state[n1] == A:
                    state[n1] = U
                elif state[n1] == U:
                    state[n1] = M
            elif state[n2] == A:
                if state[n1] == M:
                    state[n1] = U
                elif state[n1] == U:
                    state[n1] = A
    else:
        return standard(state, 0)  # Noisy part same as in standard model
    return state

# --------------------------------------------
# Simulation Loop with Averaging
# --------------------------------------------
def run_simulation_avg(step_fn, F, Tmax, num_runs=10):
    """
    Runs the simulation with a given update function (step_fn) and F value.
    Computes the average gap score G over all runs.
    """
    G_vals_all_runs = []

    for _ in range(num_runs):
        state = np.random.choice([A, U, M], size=N)
        
        # Equilibration
        for _ in range(Teq):
            state = step_fn(state, F)

        G_vals = []

        for i in range(Tmax):
            state = step_fn(state, F)
            if i % N == 0:
                m = np.sum(state == M)
                a = np.sum(state == A)
                if (m + a) > 0:
                    G_vals.append(abs(m - a) / (m + a))  # Gap score at this timepoint

        if G_vals:
            G_vals_all_runs.append(np.mean(G_vals))

    return np.mean(G_vals_all_runs)

# --------------------------------------------
# Run simulations for each model and F value
# --------------------------------------------
Gc_A, Gc_B, Gc_C = [], [], []

for F in F_values:
    print(f"Running F={F:.2f}")
    Gc_A.append(run_simulation_avg(standard, F, num_runs))
    Gc_B.append(run_simulation_avg(neighbour_limited, F, num_runs))
    Gc_C.append(run_simulation_avg(power_law_contact, F, num_runs))

# --------------------------------------------
# Plot Results: G vs F for each spatial model
# --------------------------------------------
plt.figure(figsize=(12, 4))
x_min, x_max = 0.1, 100

for idx, (Gc, case) in enumerate(zip([Gc_A, Gc_B, Gc_C], ['A', 'B', 'C']), start=1):
    ax = plt.subplot(1, 3, idx)
    ax.plot(F_values, Gc, '-', color='darkblue')
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 1)
    ax.set_xlabel('F')
    ax.set_ylabel('Gap Score (G)')
    ax.set_title(f'Case {case}')
    ax.grid(True)

plt.tight_layout()
plt.show()
