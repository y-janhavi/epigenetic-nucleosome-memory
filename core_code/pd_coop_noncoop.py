import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# -----------------------
# Model Parameters
# -----------------------
N = 60             # Number of nucleosomes
F = 77             # Feedback-to-noise ratio
Tmax = 200000 * N  # Total simulation steps after equilibration
Teq = 10 * N       # Equilibration steps
runs = 10          # Number of independent simulations

# Nucleosome States
A, U, M = 0, 1, 2 # Acetylated, Unmodified, Methylated

# ----------------------------------------------------------
# STEP FUNCTIONS FOR COOPERATIVE AND NON-COOPERATIVE CASES
# ----------------------------------------------------------

# --------- Case A: Both modifying and demodifying enzymes ---------

def step_coop_A(state, F):
    """ Cooperative recruitment for both modification and demodification """
    α = F / (F + 1)
    n1 = np.random.randint(N)
    if np.random.rand() < α:
        n2, n3 = np.random.randint(N, size=2)
        if state[n2] == state[n3] and state[n2] != U:
            if state[n2] == M:
                if state[n1] == A: state[n1] = U
                elif state[n1] == U: state[n1] = M
            elif state[n2] == A:
                if state[n1] == M: state[n1] = U
                elif state[n1] == U: state[n1] = A
    else:
        # Random (noise) conversion
        r = np.random.rand()
        if state[n1] == A and r < 1/3: state[n1] = U
        elif state[n1] == U:
            if r < 1/3: state[n1] = A
            elif r < 2/3: state[n1] = M
        elif state[n1] == M and r < 1/3: state[n1] = U
    return state

def step_noncoop_A(state, F):
    """ Non-cooperative version of Case A """
    α = F / (F + 1)
    n1 = np.random.randint(N)
    if np.random.rand() < α:
        n2 = np.random.randint(N)
        if state[n2] != U:
            if state[n2] == M:
                if state[n1] == A: state[n1] = U
                elif state[n1] == U: state[n1] = M
            elif state[n2] == A:
                if state[n1] == M: state[n1] = U
                elif state[n1] == U: state[n1] = A
    else:
        return step_coop_A(state, 0)
    return state

# --------- Case B: Only modifying enzymes recruited ---------

def step_coop_B(state, F):
    """ Cooperative recruitment for modification only """
    α = F / (F + 1)
    n1 = np.random.randint(N)
    if np.random.rand() < α:
        n2, n3 = np.random.randint(N, size=2)
        if state[n2] == state[n3] and state[n2] != U:
            if state[n2] == M and state[n1] == U: state[n1] = M
            elif state[n2] == A and state[n1] == U: state[n1] = A
    else:
        r = np.random.rand()
        if state[n1] == A and r < 1/3: state[n1] = U
        elif state[n1] == U:
            if r < 1/3: state[n1] = A
            elif r < 2/3: state[n1] = M
        elif state[n1] == M and r < 1/3: state[n1] = U
    return state

def step_noncoop_B(state, F):
    """ Non-cooperative version of Case B """
    α = F / (F + 1)
    n1 = np.random.randint(N)
    if np.random.rand() < α:
        n2 = np.random.randint(N)
        if state[n2] != U:
            if state[n2] == M and state[n1] == U: state[n1] = M
            elif state[n2] == A and state[n1] == U: state[n1] = A
    else:
        return step_coop_B(state, 0)
    return state

# --------- Case C: Only demodifying enzymes recruited ---------

def step_coop_C(state, F):
    """ Cooperative recruitment for demodification only """
    α = F / (F + 1)
    n1 = np.random.randint(N)
    if np.random.rand() < α:
        n2, n3 = np.random.randint(N, size=2)
        if state[n2] == state[n3]:
            if state[n2] == M and state[n1] == A: state[n1] = U
            elif state[n2] == A and state[n1] == M: state[n1] = U
    else:
        r = np.random.rand()
        if state[n1] == A and r < 1/3: state[n1] = U
        elif state[n1] == U:
            if r < 1/3: state[n1] = A
            elif r < 2/3: state[n1] = M
        elif state[n1] == M and r < 1/3: state[n1] = U
    return state

def step_noncoop_C(state, F):
    """ Non-cooperative version of Case C """
    α = F / (F + 1)
    n1 = np.random.randint(N)
    if np.random.rand() < α:
        n2 = np.random.randint(N)
        if state[n2] == M and state[n1] == A: state[n1] = U
        elif state[n2] == A and state[n1] == M: state[n1] = U
    else:
        return step_coop_C(state, 0)
    return state

# ------------------------------------------------
# Simulation: Run dynamics and record histograms
# ------------------------------------------------

def simulate_distribution(step_fn, F, runs):
    """ Simulate and collect distribution of (M - A) across multiple runs """
    hist_total = Counter()
    for _ in range(runs):
        state = np.random.choice([A, U, M], size=N)
        for _ in range(Teq):
            state = step_fn(state.copy(), F)
        for t in range(Tmax):
            state = step_fn(state, F)
            if t % N == 0:
                m = np.sum(state == M)
                a = np.sum(state == A)
                delta = m - a
                hist_total[delta] += 1
    return hist_total

def get_normalized(hist):
    """ Normalize histogram to probability distribution """
    x = sorted(hist.keys())
    y = np.array([hist[k] for k in x], dtype=float)
    y /= y.sum()
    return x, y

def plot_distribution(case_name, hist_coop, hist_noncoop):
    """ Plot cooperative vs non-cooperative distributions """
    x1, y1 = get_normalized(hist_coop)
    x2, y2 = get_normalized(hist_noncoop)

    plt.figure(figsize=(6, 4))
    plt.plot(x1, y1, label='Cooperative', linestyle='-', color='blue')
    plt.plot(x2, y2, label='Non-Cooperative', linestyle='--', color='blue')
    plt.yscale('log')
    plt.xlabel('M - A')
    plt.ylabel('P(M - A)')
    plt.title(f'{case_name} (F = {F})')
    plt.xlim(-60, 60)
    plt.ylim(1e-5, 1e0)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# ------------------------------------------------
# Run All Three Simulation Cases
# ------------------------------------------------

print("Running Case A...")
plot_distribution("Case A", simulate_distribution(step_coop_A, F, runs), simulate_distribution(step_noncoop_A, F, runs))

print("Running Case B...")
plot_distribution("Case B", simulate_distribution(step_coop_B, F, runs), simulate_distribution(step_noncoop_B, F, runs))

print("Running Case C...")
plot_distribution("Case C", simulate_distribution(step_coop_C, F, runs), simulate_distribution(step_noncoop_C, F, runs))
