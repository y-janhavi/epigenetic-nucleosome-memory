import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
N = 60                             # Total number of nucleosomes
F_values = [0.4, 1.0, 1.4, 2.0]    # Feedback-to-noise ratios (F)
Tmax = 5000 * N                    # Total number of simulation steps
Teq = 10                           # Initial equilibrium steps before recording begins

# Nucleosome states
A, U, M = 0, 1, 2 # A = acetylated, U = unmodified, M = methylated

# ------------------------------------------------------------
# Step Function: Implements state change per time step
# ------------------------------------------------------------
def step(state, F):
    alpha = F / (F + 1) # Probability of recruited (feedback-driven) conversion
    n1 = np.random.randint(N) # Randomly selected nucleosome

    if np.random.rand() < alpha:
        # Recruited conversion
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
        # Noisy conversion (random flipping between states)
        r = np.random.rand()
        if state[n1] == A:
            if r < 1/3:
                state[n1] = U
        elif state[n1] == U:
            if r < 1/3:
                state[n1] = A
            elif r < 2/3:
                state[n1] = M
            else:
                state[n1] = U
        elif state[n1] == M:
            if r < 1/3:
                state[n1] = U
    return state

# ------------------------------------------------------------
# Simulation and Plotting
# ------------------------------------------------------------
fig, axs = plt.subplots(len(F_values), 2, figsize=(12, 8))
G_value = []

for idx, F in enumerate(F_values):
    # Initialize nucleosome states randomly
    state = np.random.choice([A, M, U], size=N)

    # Initial equilibration steps (optional but realistic)
    for t in range(Teq):
        state = step(state, F)

    # Initialize record-keeping variables
    M_counts = []
    A_counts = []
    U_counts = []
    G_record = [] # G = |M - A| / (M + A)
    time = []
    t = 0

    # Main simulation loop
    for i in range(Tmax):
        state = step(state, F)
        t += 1 / N # Normalize time by number of nucleosomes

        # Record data every N steps
        if i % N == 0:
            M1 = np.sum(state == M)
            A_count = np.sum(state == A)
            U_count = np.sum(state == U)

            if M1 + A_count > 0:
                G_record.append(np.abs(M1 - A_count) / (M1 + A_count))

            M_counts.append(M1)
            A_counts.append(A_count)
            U_counts.append(U_count)
            time.append(t)

    # Compute averages and asymmetry G
    avg_M = np.mean(M_counts)
    avg_A = np.mean(A_counts)
    avg_U = np.mean(U_counts)
    G = np.mean(G_record)
    G_value.append(G)

    # Display summary for current F
    print(f"F = {F}, G = {G:.2f}, M = {avg_M:.2f}, A = {avg_A:.2f}, U = {avg_U:.2f}")

    # -----------------------------
    # Plot 1: Methylated count vs Time
    # -----------------------------
    axs[idx, 0].plot(time, M_counts)
    axs[idx, 0].set_title(f"F = {F}, G = {G:.2f}")
    axs[idx, 0].set_xlabel("Time (avg. conversions per nucleosome)")
    axs[idx, 0].set_ylabel("Number of M nucleosomes")

    # -----------------------------
    # Plot 2: Histogram of P(M)
    # -----------------------------
    axs[idx, 1].hist(M_counts, bins=np.arange(0, N + 1), density=True, align='left', rwidth=0.8)
    axs[idx, 1].set_title(f"F = {F}")
    axs[idx, 1].set_xlabel("Number of M nucleosomes")
    axs[idx, 1].set_ylabel("Probability")

# Format all plots nicely
plt.tight_layout()
plt.show()

# Display alpha values for clarity
for F in F_values:
    alpha = F / (F + 1)
    print(f"F = {F}, alpha = {alpha:.2f}")