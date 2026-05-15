# Numerical Tools for Aircraft Engine Vibration Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2025a-blue.svg)](https://mathworks.com)
[![Python](https://img.shields.io/badge/Python-3.x-green.svg)](https://python.org)
[![Thesis](https://img.shields.io/badge/Thesis-UniPD-red.svg)](https://thesis.unipd.it/handle/20.500.12608/107544)

The design of modern aircraft engines demands significant resources in terms of both time and funding. At present, vibration is among the critical areas primarily responsible for design inefficiency, as predictions in this field still rely extensively on semi-empirical correlations.

This repository provides numerical tools to investigate engine vibrations on a **global scale** (shaft line assembly) and a **local scale** (blade rows). While Finite Elements (FE) models are the current standard for vibration assessment, their computational cost is often prohibitive for rapid design evaluations. Consequently, the introduction of low-to-mid-fidelity Reduced Order Models enables computationally efficient investigations. Leveraging the FE approach, the developed tools are consistent across both rotordynamic and aeroelastic phenomena.

For a comprehensive overview of the theoretical background and the complete work, please refer to the full thesis:
📄 [https://thesis.unipd.it/handle/20.500.12608/107544](https://thesis.unipd.it/handle/20.500.12608/107544)

---

## Repository Structure

```
MasterThesis-Codes/
├── Rotordynamics/          # MATLAB – global scale (shaft-line FE model)
├── Aeroelasticity/         # Python – local scale (blade row identification)
├── CITATION.cff            # Citation metadata
├── LICENSE                 # MIT License
└── README.md
```

---

## Quick Start

### Rotordynamics (MATLAB)

```matlab
% Navigate to the Rotordynamics folder and run the main script
cd Rotordynamics
% Open and run one of the example scripts to set up your rotor model
% See the User Manual for full data structure documentation
```

### Aeroelasticity (Python)

```bash
# Install dependencies
pip install numpy scipy matplotlib cyipopt

# Navigate to the Aeroelasticity folder and run the main script
cd Aeroelasticity
python main.py
```

> See the **User Manual** (available in Italian and English inside the repository) for detailed usage instructions, data structures, and function signatures.

---

## Global Scale Analysis: Rotordynamics

For the global scale, a one-dimensional shaft-line FE model was developed to evaluate the lateral vibrations of multi-disk rotor systems. The theoretical formulation and implementation are based on the work of **Michael I. Friswell** (*Dynamics of Rotating Machines*, Cambridge University Press, 2010).

By expressing the element matrices in closed form, the solution is obtained without numerical integration, significantly reducing processing time. The tool addresses both linear and non-linear system behaviours, incorporating effects from bearings, squeeze-film dampers, seals, and aerodynamic forces.

The rotordynamic software was successfully validated against four literature cases and demonstrated consistency with modern commercial and in-house software, including DYNROT, XLRotor, and ANSYS.

**Key Capabilities:**

| Category | Features |
|---|---|
| Analysis Types | Modal analysis (natural frequencies, Campbell diagrams, critical speeds), Frequency Response (unbalance and shaft bow), Time-Domain Simulation (transient and steady-state nonlinear via `ode15s`) |
| Modeling Elements | Euler-Bernoulli and Timoshenko beam formulations (including tapered elements), rigid disks, linear constraints, oil film journal bearings, Thomas-Alford forces |
| Post-Processing | Automated mesh generation, Craig-Bampton modal reduction, Campbell diagrams, Root Locus, Bode plots, Poincaré maps, Waterfall diagrams |

---

## Local Scale Analysis: Aeroelasticity

For the local scale, the framework implements the identification of mistuning parameters, aerodynamic forcing functions, and Aerodynamic Influence Coefficients. The estimation leverages a **Least-Squares Complex Frequency-Domain (LSCF)** algorithm applied to simulated Blade Tip Timing measurements. The aeroelastic identification codes are implemented entirely in **Python**.

The procedure was applied to a simulated rotor based on the **Purdue Rotor 2** case study, reproducing the results established by Hall et al. (2024):
📎 [DOI: 10.1115/1.4064816](https://doi.org/10.1115/1.4064816)

---

## Requirements

### Rotordynamics (MATLAB)

- **MATLAB R2025a** or later
- No external toolboxes required beyond the standard MATLAB installation

### Aeroelasticity (Python)

- **Python 3.x**
- **Ipopt** (Interior Point Optimizer) solver — [GitHub](https://github.com/coin-or/IPOPT), based on Wächter & Biegler ([DOI: 10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y))
- Standard scientific libraries: `NumPy`, `SciPy`, `Matplotlib`, `cyipopt`

Install all Python dependencies at once:

```bash
pip install -r Aeroelasticity/requirements.txt
```

---

## Documentation

A comprehensive **User Manual** is provided within the repository in both **Italian** and **English**. It covers:
- Theoretical background for both modules
- Data structures and matrix assembly procedures
- Full function signatures and parameter descriptions
- Step-by-step worked examples

---

## Citation

If you use this software in your research, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) or the following format:

> Tabarelli, F. (2025). *Numerical Tools for Aircraft Engine Vibration Analysis* [Software]. University of Padova. https://github.com/tabarelli-filippo/MasterThesis-Codes

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
