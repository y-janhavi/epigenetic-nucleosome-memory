# epigenetic-nucleosome-memory
Stochastic simulation of the Dodd et al. (2007) model for epigenetic cell memory, exploring nucleosome modification feedback and bistability.
# Epigenetic Cell Memory: Theoretical Analysis of Nucleosome Modification

This repository contains a numerical reproduction of the stochastic model from the paper:
**Dodd et al. (2007). "Theoretical Analysis of Epigenetic Cell Memory by Nucleosome Modification." Cell.**

## Project Overview
Projects in this discipline revolve around the reproduction and computation of a stochastic model - more specifically, the one outlined by Dodd et al. (Cell, 2007) - which sought to understand the underlying principles of memory in terms of epigenetics within the context of nucleosome modification. The model is based on an array of nu-cleosomes with each of them being in one of the following three states: methylated (M),acetylated (A) and unmodified (U). To gain insight on the system’s cooperativity and spatial constraints, we administered the model noise-and-feedback controlled transitions and studied the system under varying conditions with respect to bistability, lifetime of lasting memory states, and formation of memory states. Every memory state of the model in this case was systematically tested and motioned to gauge how the gap scores and lifetimes of the said states dictate the level of feedback versus noise present in the state. With the feedback having dual cooperativity, the computations proved that final stability is indeed possible, but only in single-feedback systems. Testing other constraints such as neighbor-only spatial interaction constraints and power-law recruitment provided the values necessary to allow the said memory states to remain stable, showcasing the necessity of long-range interactions within the formation of stable epigenetic memories. Key findings from the study conducted by Dodd et al. were successfully reproduced with the use of other simulation results, further broadening the findings on interaction

## Code Structure & Figures
Each script is designed to reproduce a specific thermodynamic or biological property discussed in the paper:

1. **Bistability Analysis**: Script to reproduce the primary bistability plots showing M and A state stability.
2. **Average Lifetime**: Calculates how long a specific epigenetic state persists before a stochastic flip.
3. **Cooperativity vs. Non-Cooperative**: Compares how different recruitment mechanisms affect memory.
4. **Probability Distributions**: Density plots for nucleosome occupancy in various recruitment modes.
5. **Spatial Constraints**: Implements the model where recruitment depends on the physical distance between nucleosomes.
6. **Spatial Probability Density**: Statistical analysis of the effects of 1D/3D spatial constraints.

## Requirements
- Python 3.x
- NumPy
- Matplotlib
- SciPy
