import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
N = 60                              # Number of nucleosomes
F_values = np.logspace(-1, 2, 30)   # Feedback strength values (log scale from 0.1 to 100)
Tmax = 20000 * N                    # Total simulation time
Teq = 10 * N                        # Equilibration time
num_runs = 10                       # Number of independent runs for averaging

# Nucleosome states: A (acetylated), U (unmodified), M (methylated)
A, U, M = 0, 1, 2

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

# --- Simulation function with time-averaged G calculation ---
def run_simulation_avg(step_fn, F, Tmax, num_runs=10):
    G_vals_all_runs = []
    for _ in range(num_runs):
        state = np.random.choice([A, U, M], size=N) # Initial random state
        for _ in range(Teq): # Equilibration phase
            state = step_fn(state, F)
        G_vals = []
        for i in range(Tmax): # Main simulation loop
            state = step_fn(state, F)
            if i % N == 0:
                m = np.sum(state == M)
                a = np.sum(state == A)
                if (m + a) > 0:
                    G_vals.append(abs(m - a) / (m + a)) # Order parameter G
        if G_vals:
            G_vals_all_runs.append(np.mean(G_vals))
    return np.mean(G_vals_all_runs)

# --- Run simulations for all cases and F values ---
Gc_A, Gn_A = [], [] # Case A
Gc_B, Gn_B = [], [] # Case B
Gc_C, Gn_C = [], [] # Case C

for F in F_values:
    print(f"Running F={F:.2f}")
    Gc_A.append(run_simulation_avg(step_coop_A, F, Tmax, num_runs))
    Gn_A.append(run_simulation_avg(step_noncoop_A, F, Tmax, num_runs))
    Gc_B.append(run_simulation_avg(step_coop_B, F, Tmax, num_runs))
    Gn_B.append(run_simulation_avg(step_noncoop_B, F, Tmax, num_runs))
    Gc_C.append(run_simulation_avg(step_coop_C, F, Tmax, num_runs))
    Gn_C.append(run_simulation_avg(step_noncoop_C, F, Tmax, num_runs))

# --- Plotting the G vs F curves for all cases ---
plt.figure(figsize=(12, 4))
x_min, x_max = 0.1, 100

for idx, (Gc, Gn, case) in enumerate(zip([Gc_A, Gc_B, Gc_C], [Gn_A, Gn_B, Gn_C], ['A', 'B', 'C']), start=1):
    ax = plt.subplot(1, 3, idx)
    ax.plot(F_values, Gc, '-', label='Cooperative') # Solid line
    ax.plot(F_values, Gn, '--', label='Non-cooperative') # Dashed line
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel('F') # Feedback strength
    ax.set_ylabel('G') # Order parameter
    ax.set_title(f'Case {case}')
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.show()
