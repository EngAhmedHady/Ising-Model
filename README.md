# 2D Ising Model Simulation

This repository contains the implementation of a Monte Carlo simulation for the two-dimensional Ising model, utilizing the Metropolis-Hastings algorithm. The project investigates ferromagnetism properties, focusing on magnetization and energy at varying temperatures.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [References](#references)

## Introduction

The Ising model is a simplified theoretical framework for understanding ferromagnetism, invented by Wilhelm Lenz and named after his student, Ernst Ising. This project simulates the 2D Ising model on a square lattice using periodic boundary conditions. The Hamiltonian for the system is given by:

$$H(\sigma) = -J \sum_{\langle i, j \rangle} \sigma_i \sigma_j, \quad J > 0$$

where:
- $J$: interaction strength (positive for ferromagnetic materials).
- $\sigma_i = \pm 1$: spin states at lattice sites.

The simulation analyzes the system at various temperatures, focusing on magnetization, energy, susceptibility, and specific heat.

## Features

- Simulates the 2D Ising model using Monte Carlo methods.
- Implements the Metropolis algorithm to update spin configurations.
- Computes key physical observables:
  - Magnetization
  - Energy
  - Susceptibility
  - Specific heat
- Supports different initial configurations:
  - Fully positive spins
  - Fully negative spins
  - Random spins (50% probability for $+1$ or $-1$).

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/EngAhmedHady/Ising-Model.git
   ```
2. Ensure you have a FORTRAN compiler installed (e.g., `gfortran`).
3. Compile the source code:
   ```bash
   gfortran -o ising_model Final.f95
   ```

## Usage

Run the compiled program:
```bash
./ising_model
```

Edit the input parameters (e.g., temperature range, grid size, initial configuration) directly in the `Final.f95` source file before compilation.

### Parameters
- **Grid Size:** Default is $100 \times 100$.
- **Temperature Range:** From $T = 1$ to $T = 4$ with a step size of $0.1$.
- **Iterations:** Number of Monte Carlo steps per spin.

## Results
![TemperatureResult](https://github.com/user-attachments/assets/43f6fb5f-5bb6-4b37-a15c-65b39b99d572)

### Key Observations
- **Low Temperatures:** Spins align, leading to stable magnetization and low energy.
- **Critical Temperature:** A phase transition occurs, causing abrupt changes in physical properties.
- **High Temperatures:** Random spin configurations dominate, resulting in higher energy and negligible magnetization.

### Sample Outputs
- Energy and magnetization plotted against time and temperature.
- Visualization of spin configurations under different conditions.

### Known Issues
- Instability in some random configurations may require extended simulation times to stabilize.

## References

1. Nina Kuklisova, *Solving the 2D Ising Model*, PHYS 35200, March 2013.
2. Jacques Kotze, *Introduction to Monte Carlo Methods for an Ising Model of a Ferromagnet*, March 2008.
3. Alexey Khorev, *A Monte Carlo Implementation of the Ising Model in Python*, August 2017.

For more details, refer to the [Stochastic_Report.pdf](Stochastic_Report.pdf) included in this repository.

---

Developed by Ahmed H. Hanfy as part of a project at the University of L'Aquila, supervised by Dr. Matteo Colangeli.
