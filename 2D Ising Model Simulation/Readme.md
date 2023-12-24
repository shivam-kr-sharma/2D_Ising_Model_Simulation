# Ising Model Simulation Code Files

## final_code.py

This is the main file that implements the full Ising model simulation and analysis.

- Defines functions to initialize the spin lattice, calculate energy, perform Monte Carlo steps, calculate magnetization.
- Contains the main simulation function that runs the Monte Carlo simulation across a range of temperatures.
- Stores the resulting energy and magnetization data.
- Includes a function to plot the results (energy and magnetization vs temperature).
- Can simulate different lattice sizes and interaction modes based on user input.
- Main output is plots of the phase transition behavior in energy and magnetization.
- Everything is well documented in this file. 
- The output of this file is the phase transition plots which is given one by one once you run the python file. Like, here we are just plotting energy and magnetization values against temperature, so after running this file, we first obtain the plot of energy vs temperature and then magnetization vs temperature.
- This file is designed such that the functions here can handle only one lattice at a time and then do the simulation. Here you cannot compare the observable results of different lattices by plotting it.

## metropolis.py

This implements the Metropolis Monte Carlo algorithm for simulating spin flips.

- Defines the `metropolis()` function which performs a single Monte Carlo step.
- Calculates the energy change for flipping a random spin.
- Uses Metropolis acceptance criteria to accept or reject the spin flip.
- Helper functions to pick a random site and calculate system energy.  
- This Metropolis implementation is used within `final_code.py` to actually run the Monte Carlo simulation.
- It contains the key physics of the Ising model in the energy difference calculation and Metropolis acceptance probability.
- This file can be easily customized and the user can plot the results of the simulation for the more than one lattice sizes. For example, the user can plot energy vs temperature and magnetization vs temperature on different plots but the different lattice sizes can be plotted on the same plot. So, the user can compare the plot of different lattice sizes very easily.

