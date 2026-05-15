# Numerical Tools for Aircraft Engine Vibration Analysis

The design of modern aircraft engines demands significant resources in terms of both time and funding. At present, vibration is among the critical areas primarily responsible for design inefficiency, as predictions in this field still rely extensively on semi-empirical correlations.

This repository provides numerical tools to investigate engine vibrations on a global scale (shaft line assembly) and a local scale (blade rows). While Finite Elements (FE) models are the current standard for vibration assessment, their computational cost is often prohibitive for rapid design evaluations. Consequently, the introduction of low-to-mid-fidelity Reduced Order Models enables computationally efficient investigations. Leveraging the FE approach, the developed tools are consistent across both rotordynamic and aeroelastic phenomena.

For a comprehensive overview of the theoretical background and the complete work, please refer to the full thesis available at: https://thesis.unipd.it/handle/20.500.12608/107544.

## Global Scale Analysis: Rotordynamics

For the global scale, a one-dimensional shaft-line FE model was developed to evaluate the lateral vibrations of multi-disk rotor systems. The theoretical formulation and implementation of the code are based on the work of **Michael I. Friswell** (*Dynamics of Rotating Machines*, Cambridge University Press, 2010).

By expressing the element matrices in closed form, the solution is obtained without numerical integration, significantly reducing processing time. The tool addresses both the linear and non-linear system behaviours, incorporating effects from bearings, squeeze-film dampers, seals and aerodynamic forces. 

The rotordynamic software was successfully validated against four literature cases and demonstrated consistency with modern commercial and in-house software, including DYNROT, XLRotor, and ANSYS.

**Key Capabilities:**
* **Analysis Types:** Modal analysis (natural frequencies, Campbell diagrams, critical speeds), Frequency Response (steady-state response to unbalance and shaft bow), and Time-Domain Simulation (transient and steady-state nonlinear simulations using `ode15s`).
* **Modeling Elements:** Euler-Bernoulli and Timoshenko beam formulations (including tapered elements), rigid disks, and multiple bearing/support models (linear constraints, oil film journal bearings, Thomas-Alford forces, etc.).
* **Post-Processing:** Automated mesh generation, Craig-Bampton modal reduction, and extensive visualization functions (Campbell diagrams, Root Locus, Bode plots, Poincaré maps, Waterfall diagrams).

## Local Scale Analysis: Aeroelasticity

For the local scale, the framework implements the identification of mistuning parameters, aerodynamic forcing functions, and Aerodynamic Influence Coefficients. The estimation leverages a Least-Squares Complex Frequency-Domain algorithm applied to simulated Blade Tip Timing measurements. The aeroelastic identification codes are implemented entirely in **Python**.

The aeroelastic identification procedure was applied to a simulated rotor, based on the Purdue Rotor 2 case study, reproducing the results established by Hall et al. (2024, [DOI: 10.1115/1.4064816](https://doi.org/10.1115/1.4064816)).

## Documentation and Requirements

For the complete operation and usage of the codes, a comprehensive User Manual is provided within the repository, available in both Italian and English. The manual details the theoretical background, the data structures, matrix assembly procedures, and function signatures required to run the analyses.

**System Requirements:**

* **Rotordynamics (MATLAB):**
  * Developed and tested using **MATLAB R2025a**. No strict external dependencies are required beyond standard MATLAB toolboxes.
* **Aeroelasticity (Python):**
  * Written in **Python**.
  * Requires the **Ipopt** (Interior Point Optimizer) solver ([GitHub repository](https://github.com/coin-or/IPOPT)) for the numerical optimization routines. The implementation of this solver is based on the work of Wächter and Biegler ([DOI: 10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y)).
  * Standard scientific Python libraries are required (e.g., `NumPy`, `SciPy`, `Matplotlib`, and a Python interface for Ipopt such as `cyipopt`).
