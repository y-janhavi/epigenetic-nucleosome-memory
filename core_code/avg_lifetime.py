import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Simulation Parameters
# -------------------------------
N = 60                  # Number of nucleosomes
Tmax = 800000           # Simulation steps per run
n_runs = 10             # Number of runs for averaging
Teq = 10*N              # Equilibrium time steps
F_values = np.concatenate([
    np.logspace(-1, np.log10(1.0), 10),
    np.logspace(np.log10(1.2), np.log10(4.0), 8)]) # Feedback-to-noise ratios on log scale

# Nucleosome States
A, U, M = 0, 1, 2 # Acetylated, Unmodified, Methylated

# -------------------------------
# Step Function
# -------------------------------
def step(state, F):
    """Performs one simulation update step using feedback or noise."""
    alpha = F / (1 + F)
    n1 = np.random.randint(N) # Target nucleosome

    if np.random.rand() < alpha:
        # Feedback-based recruitment step
        n2 = np.random.randint(N)
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
        # Noisy random conversion (non-feedback)
        r = np.random.rand()
        if state[n1] == A and r < 1/3:
            state[n1] = U
        elif state[n1] == U:
            if r < 1/3: state[n1] = A
            elif r < 2/3: state[n1] = M
            else: state[n1] = U
        elif state[n1] == M and r < 1/3:
            state[n1] = U
    return state

# -------------------------------
# Lifetime Calculation Function
# -------------------------------
def calculate_lifetime(F):
    """
    Measures average duration for which the system stays in a single
    stable chromatin state (high-M or high-A) before switching.
    """
    state = np.random.choice([A, M, U], size=N)
    current_state = None
    current_duration = 0
    durations = []

    for _ in range(Tmax):
        state = step(state, F)
        M1 = np.sum(state == M)
        A1 = np.sum(state == A)

        # Apply 1.5× rule to detect dominant state
        if M1 > 1.5 * A1:
            label = 'M' # High-M state
        elif A1 > 1.5 * M1:
            label = 'A' # High-A state
        else:
            label = None # Mixed or unclassified

        # Track how long the system remains in the same state
        if label == current_state:
            current_duration += 1
        else:
            if current_state is not None:
                durations.append(current_duration)
            current_state = label
            current_duration = 1

    # Add final duration to the list
    if current_state is not None:
        durations.append(current_duration)

    # Return average lifetime (in simulation steps)
    return np.mean(durations) if durations else 0

# -------------------------------
# Gap Score Calculation
# -------------------------------
def compute_gap_measure(F):

    """Computes the average gap score G = ⟨|M - A| / (M + A)⟩"""
    
    state = np.random.choice([A,M,U], size=N)
    gap_sum = 0

    for t in range(Teq):
        state = step(state, F)

    for _ in range(Tmax):
        state = step(state, F)
        M1 = np.sum(state == M)
        A1 = np.sum(state == A)
        if M1 + A1 > 0:
            gap_sum += np.abs(M1 - A1) / (M1 + A1)

    return gap_sum / Tmax

# -------------------------------
# Run Simulations Across F Values
# -------------------------------
lifetimes = []
gap_measures = []

for F in F_values:
    # Repeat simulations and average results
    lifetime_runs = [calculate_lifetime(F) for _ in range(n_runs)]
    gap_runs = [compute_gap_measure(F) for _ in range(n_runs)]

    lifetimes.append(np.mean(lifetime_runs))
    gap_measures.append(np.mean(gap_runs))

    print(f"F = {F:.2f} | Lifetime = {lifetimes[-1]:.1f} | Gap G = {gap_measures[-1]:.3f}")

# -------------------------------
# Plot: Lifetime vs F (Figure 2E)
# -------------------------------
plt.figure()
plt.plot(F_values, lifetimes, marker='o', color='navy')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1, 4.0)
plt.ylim(1e0, 1e6)
plt.xlabel('F (log scale)')
plt.ylabel('Average State Lifetime (log scale)')
plt.title('Lifetime vs Feedback-to-Noise Ratio (F)')
plt.grid(True)
plt.xticks([0.1, 0.2, 0.5, 1, 2, 4], labels=['0.1', '0.2', '0.5', '1', '2', '4'])
plt.tight_layout()
plt.show()

# -------------------------------
# Plot: Gap Score vs F (Figure 2F)
# -------------------------------
plt.figure()
plt.plot(F_values, gap_measures, marker='o', color='navy')
plt.xscale('log')
plt.xlim(0.1, 4.0)
plt.ylim(0, 1.0)
plt.xlabel('F (log scale)')
plt.ylabel('Gap Score (G)')
plt.title('Gap Score vs Feedback-to-Noise Ratio (F)')
plt.grid(True)
plt.xticks([0.1, 0.2, 0.5, 1, 2, 4], labels=['0.1', '0.2', '0.5', '1', '2', '4'])
plt.tight_layout()
plt.show()

