# Numerical Tools for Aircraft Engine Vibration Analysis

[cite_start]The design of modern aircraft engines demands significant resources in terms of both time and funding[cite: 1]. [cite_start]At present, vibration is among the critical areas primarily responsible for design inefficiency, as predictions in this field still rely extensively on semi-empirical correlations[cite: 2].

[cite_start]This repository provides numerical tools to investigate engine vibrations on a global scale (shaft line assembly) and a local scale (blade rows)[cite: 3]. [cite_start]While Finite Elements (FE) models are the current standard for vibration assessment, their computational cost is often prohibitive for rapid design evaluations[cite: 4]. [cite_start]Consequently, the introduction of low-to-mid-fidelity Reduced Order Models enables computationally efficient investigations[cite: 5]. [cite_start]Leveraging the FE approach, the developed tools are consistent across both rotordynamic and aeroelastic phenomena[cite: 6].

## Global Scale Analysis: Rotordynamics

[cite_start]For the global scale, a one-dimensional shaft-line FE model was developed to evaluate the lateral vibrations of multi-disk rotor systems[cite: 7]. [cite_start]By expressing the element matrices in closed form, the solution is obtained without numerical integration, significantly reducing processing time[cite: 8]. [cite_start]The tool addresses both the linear and non-linear system behaviours, incorporating effects from bearings, squeeze-film dampers, seals and aerodynamic forces[cite: 9]. 

[cite_start]The rotordynamic software was successfully validated against four literature cases and demonstrated consistency with modern commercial and in-house software, including DYNROT, XLRotor, and ANSYS[cite: 10].

**Key Capabilities:**
* **Analysis Types:** Modal analysis (natural frequencies, Campbell diagrams, critical speeds), Frequency Response (steady-state response to unbalance and shaft bow), and Time-Domain Simulation (transient and steady-state nonlinear simulations using `ode15s`).
* **Modeling Elements:** Euler-Bernoulli and Timoshenko beam formulations (including tapered elements), rigid disks, and multiple bearing/support models (linear constraints, oil film journal bearings, Thomas-Alford forces, etc.).
* **Post-Processing:** Automated mesh generation, Craig-Bampton modal reduction, and extensive visualization functions (Campbell diagrams, Root Locus, Bode plots, Poincaré maps, Waterfall diagrams).

## Local Scale Analysis: Aeroelasticity

[cite_start]For the local scale, the framework implements the identification of mistuning parameters, aerodynamic forcing functions, and Aerodynamic Influence Coefficients[cite: 9]. [cite_start]The estimation leverages a Least-Squares Complex Frequency-Domain algorithm applied to simulated Blade Tip Timing measurements[cite: 10]. 

The aeroelastic identification procedure was applied to a simulated rotor, based on the Purdue Rotor 2 case study, reproducing the results established by Hall et al. (2024) [cite_start][cite: 11].

---

## Usage Example (Rotordynamics)

The following MATLAB script demonstrates the setup for a basic linear analysis, including mesh generation, component definition, and Campbell diagram computation.

```matlab
% 1. Mesh generation
Rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, ...
            isGraded, shaftElType, shaftElProperties);

% 2. Disks setup (with unbalance)
Rotor.disk(1) = struct('type', 2, 'node', 5, ...
                       'mass', 2.3, 'Id', 3.5e-6, 'Ip', 4e-6);
Rotor.disk(1).epsilon = 1e-6;

% 3. Bearings setup
Rotor.bearing(1) = struct('type', 1, 'node', 1);     % Pinned support
Rotor.bearing(2) = struct('type', 7, 'node', 15, ... % Linear oil film
    'F', 500, 'D', 0.1, 'L', 0.03, 'c', 2e-4, 'eta', 0.03);

% 4. Forcing setup
Rotor.forcing(1).type = 1;    % Mass unbalance

% 5. Calculate natural frequencies & plot Campbell diagram
RotorSpeed = (0:100:5000) * pi/30;   % [rad/s]
[evals, evecs, kappa] = charRoots(Rotor, RotorSpeed);
[evals_s, evecs_s, kappa_s] = sortModesMAC(evals, evecs, kappa);

NX = [1 2]; isDamped = true;
plotCampbell(RotorSpeed, evals_s(1:6,:), NX, isDamped, kappa_s(:,1:6,:));

% 6. Find critical speeds
[omega_c, modes_c] = critSpeeds(Rotor, 1, isDamped, 2, "num_crit", 4);
fprintf('1st critical speed: %.1f rpm\n', omega_c(1)*30/pi);
