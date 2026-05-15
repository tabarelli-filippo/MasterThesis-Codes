# Numerical Tools for Aircraft Engine Vibration Analysis

The design of modern aircraft engines demands significant resources in terms of both time and funding. At present, vibration is among the critical areas primarily responsible for design inefficiency, as predictions in this field still rely extensively on semi-empirical correlations.

This repository provides numerical tools to investigate engine vibrations on a global scale (shaft line assembly) and a local scale (blade rows). While Finite Elements (FE) models are the current standard for vibration assessment, their computational cost is often prohibitive for rapid design evaluations. Consequently, the introduction of low-to-mid-fidelity Reduced Order Models enables computationally efficient investigations. Leveraging the FE approach, the developed tools are consistent across both rotordynamic and aeroelastic phenomena.

For a comprehensive overview of the theoretical background and the complete work, please refer to the full thesis available at: https://thesis.unipd.it/handle/20.500.12608/107544.

## Global Scale Analysis: Rotordynamics

For the global scale, a one-dimensional shaft-line FE model was developed to evaluate the lateral vibrations of multi-disk rotor systems. By expressing the element matrices in closed form, the solution is obtained without numerical integration, significantly reducing processing time. The tool addresses both the linear and non-linear system behaviours, incorporating effects from bearings, squeeze-film dampers, seals and aerodynamic forces. 

The rotordynamic software was successfully validated against four literature cases and demonstrated consistency with modern commercial and in-house software, including DYNROT, XLRotor, and ANSYS.

**Key Capabilities:**
* **Analysis Types:** Modal analysis (natural frequencies, Campbell diagrams, critical speeds), Frequency Response (steady-state response to unbalance and shaft bow), and Time-Domain Simulation (transient and steady-state nonlinear simulations using `ode15s`).
* **Modeling Elements:** Euler-Bernoulli and Timoshenko beam formulations (including tapered elements), rigid disks, and multiple bearing/support models (linear constraints, oil film journal bearings, Thomas-Alford forces, etc.).
* **Post-Processing:** Automated mesh generation, Craig-Bampton modal reduction, and extensive visualization functions (Campbell diagrams, Root Locus, Bode plots, Poincaré maps, Waterfall diagrams).

## Local Scale Analysis: Aeroelasticity

For the local scale, the framework implements the identification of mistuning parameters, aerodynamic forcing functions, and Aerodynamic Influence Coefficients. The estimation leverages a Least-Squares Complex Frequency-Domain algorithm applied to simulated Blade Tip Timing measurements. 

The aeroelastic identification procedure was applied to a simulated rotor, based on the Purdue Rotor 2 case study, reproducing the results established by Hall et al. (2024).

## Documentation and Requirements

For the complete operation and usage of the code, a comprehensive User Manual is provided within the repository, available in both Italian and English. The manual details the theoretical background, the data structures, matrix assembly procedures, and function signatures required to run the solver.

**Requirements:** MATLAB (developed and tested on recent releases). No strict external dependencies are required beyond standard MATLAB toolboxes.
