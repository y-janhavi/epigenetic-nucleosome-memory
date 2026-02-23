import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# ---------------------------
# Parameters
# ---------------------------
N = 60                           # Number of nucleosomes
F_values = [1, 2.6, 6, 26, 77]   # Feedback-to-noise ratios to test
Tmax = 50000 * N                 # Simulation time
Teq = 10 * N                     # Equilibration time
runs = 10                        # Independent runs per F

# Nucleosome States
A, U, M = 0, 1, 2 # Acetylated, Unmodified, Methylated


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
def neighbor_limited(state, F):
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



# ---------------------------
# Simulation Function
# ---------------------------
def simulate_distribution(step_fn, F, runs):
    hist_total = Counter()
    for _ in range(runs):
        state = np.random.choice([A, U, M], size=N)
        for _ in range(Teq): # Equilibration
            state = step_fn(state.copy(), F)

        for t in range(Tmax):
            state = step_fn(state, F)
            if t % N == 0:
                m = np.sum(state == M)
                a = np.sum(state == A)
                delta = m - a
                hist_total[delta] = hist_total.get(delta, 0) + 1
    return hist_total


# ---------------------------
# Plotting Function
# ---------------------------
def plot_distributions(model_name, step_fn, runs):
    fig, ax = plt.subplots(figsize=(8, 5))
    for F in F_values:
        hist = simulate_distribution(step_fn, F, runs)
        x = sorted(hist.keys())
        y = np.array([hist[k] for k in x], dtype=float)
        y /= y.sum()
        ax.plot(x, y, label=f'F = {F}')
        print(f"Simulation done for {model_name}, F = {F}")
    ax.set_yscale('log')
    ax.set_xlim(-60, 60)
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel('M - A')
    ax.set_ylabel('P(M - A)')
    ax.set_title(f'{model_name} Model: Probability Distribution')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()


# ---------------------------
# Run All Models
# ---------------------------
plot_distributions("Standard (Case A)", standard, runs)
plot_distributions("Neighbor-limited (Case B)", neighbor_limited, runs)
plot_distributions("Power-law Contact (Case C)", power_law_contact, runs)
