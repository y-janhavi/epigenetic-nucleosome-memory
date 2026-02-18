# epigenetic-nucleosome-memory
Stochastic simulation of the Dodd et al. (2007) model for epigenetic cell memory, exploring nucleosome modification feedback and bistability.
# Epigenetic Cell Memory: Theoretical Analysis of Nucleosome Modification

This repository contains a numerical reproduction of the stochastic model from the paper:
**Dodd et al. (2007). "Theoretical Analysis of Epigenetic Cell Memory by Nucleosome Modification." Cell.**

## 🔬 Project Overview
This project simulates how histone modifications (M, U, A states) and positive feedback loops create a stable "epigenetic memory." It explores the transition from noisy random modifications to stable bistable states.

## 📂 Code Structure & Figures
Each script is designed to reproduce a specific thermodynamic or biological property discussed in the paper:

1. **Bistability Analysis**: Script to reproduce the primary bistability plots showing M and A state stability.
2. **Average Lifetime**: Calculates how long a specific epigenetic state persists before a stochastic flip.
3. **Cooperativity vs. Non-Cooperative**: Compares how different recruitment mechanisms affect memory.
4. **Probability Distributions**: Density plots for nucleosome occupancy in various recruitment modes.
5. **Spatial Constraints**: Implements the model where recruitment depends on the physical distance between nucleosomes.
6. **Spatial Probability Density**: Statistical analysis of the effects of 1D/3D spatial constraints.

## 🛠️ Requirements
- Python 3.x
- NumPy
- Matplotlib
- SciPy
