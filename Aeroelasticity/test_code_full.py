"""
test_code_full.py
=================
Rigid-disk aeromechanical identification on synthetic rigid-disk data.

Workflow
--------
1. Generate synthetic blade-vibration data with the rigid-disk model.
2. Set up AeromechanicalIdentification and an IPOPT optimiser with
   physics-informed parameter bounds and diagonal scaling.
3. Recover the 5N-1 parameter vector via maximum-likelihood estimation.
4. Compute Cramér-Rao bounds and produce diagnostic plots.

Dependencies
------------
    auxSolver.aero_ident      — AeromechanicalIdentification
    auxSolver.data_generation — generate_rotor_data
    auxSolver.plot_aero       — AeroPlotter
    cyipopt                   — Python interface to IPOPT
"""

import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.aero_ident import AeromechanicalIdentification
from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter


# ---------------------------------------------------------------------------
# Cached evaluator with parameter scaling
# ---------------------------------------------------------------------------

class CachedEvaluator:
    """
    Caching wrapper with diagonal parameter scaling for IPOPT.

    The optimiser operates on a scaled parameter vector x = theta / scales.
    Objective and gradient are recomputed only when x changes, avoiding
    redundant solves when the optimiser queries them separately.

    Parameters
    ----------
    ident_problem : AeromechanicalIdentification
    scale_factors : array-like, shape (5N-1,)
        Diagonal scaling; typically the expected magnitude of each parameter.
    """

    def __init__(self, ident_problem, scale_factors):
        self.problem = ident_problem
        self.scales  = np.asarray(scale_factors)

        self.last_theta_scaled = None
        self.last_J            = None
        self.last_grad_scaled  = None

    def _update(self, theta_scaled):
        if (self.last_theta_scaled is None
                or not np.array_equal(theta_scaled, self.last_theta_scaled)):
            theta_physical   = theta_scaled * self.scales
            J, grad_physical = self.problem.compute_objective_and_gradient(theta_physical)
            self.last_J           = J
            self.last_grad_scaled = grad_physical * self.scales
            self.last_theta_scaled = np.copy(theta_scaled)

    def objective(self, theta_scaled):
        """Return the cached objective value."""
        self._update(theta_scaled)
        return self.last_J

    def gradient(self, theta_scaled):
        """Return the cached scaled gradient."""
        self._update(theta_scaled)
        return self.last_grad_scaled

    def objective_and_gradient(self, theta_scaled):
        """Return (J, grad) together (convenience method)."""
        self._update(theta_scaled)
        return self.last_J, self.last_grad_scaled.copy()


# ---------------------------------------------------------------------------
# Main identification routine
# ---------------------------------------------------------------------------

def run_identification_with_simulated_data():
    """
    End-to-end identification pipeline for the rigid-disk model.

    Steps
    -----
    1.  Generate T=250 synthetic frequency-response snapshots.
    2.  Define physics-informed bounds on all 5N-1 parameters.
    3.  Solve the MLE optimisation problem with IPOPT (L-BFGS-B Hessian).
    4.  Compute Cramér-Rao bounds and generate three diagnostic plots.
    """

    # ------------------------------------------------------------------
    # 1.  Simulation parameters and data generation
    # ------------------------------------------------------------------
    N           = 33
    omega_0     = 17000.0         # nominal blade frequency [rad/s]
    T           = 250
    omega_t     = np.linspace(16250.0, 17750.0, T)
    noise_level = 1e-4

    print("Generating synthetic data (rigid-disk model)...")
    result     = generate_rotor_data('rigid', omega_t, omega_0, noise_level, N=N)
    y_measured = result['y_measured']
    theta_true = result['theta_true']

    ident_problem = AeromechanicalIdentification(y_measured, omega_t, omega_0, N)
    num_params    = 5 * N - 1

    # ------------------------------------------------------------------
    # 2.  Block offsets in theta
    # ------------------------------------------------------------------
    idx_d0 = 0
    idx_a0 = N
    idx_a1 = N + (N - 1)
    idx_fr = N + (N - 1) + N

    # ------------------------------------------------------------------
    # 3.  Physics-informed parameter bounds
    # ------------------------------------------------------------------
    lb = np.empty(num_params)
    ub = np.empty(num_params)

    lb[idx_d0 : idx_a0] = -0.10;  ub[idx_d0 : idx_a0] =  0.10   # mistuning
    lb[idx_a0 : idx_a1] = -0.005; ub[idx_a0 : idx_a1] =  0.005  # a0 stiffness
    lb[idx_a1 : idx_fr] = -0.010; ub[idx_a1 : idx_fr] =  0.010  # a1 damping
    lb[idx_a1]          =  1e-6;  ub[idx_a1]          =  0.020  # self-damping ≥ 0
    lb[idx_fr :]        = -1.0;   ub[idx_fr :]        =  1.0    # forcing

    # ------------------------------------------------------------------
    # 4.  Initial guess
    # ------------------------------------------------------------------
    np.random.seed(99)
    theta_initial = 1e-5 * np.ones(num_params)

    # ------------------------------------------------------------------
    # 5.  Diagonal scaling for numerical conditioning
    #     Scale factors approximate the expected parameter magnitudes,
    #     so the optimiser works with O(1) variables.
    # ------------------------------------------------------------------
    scale_factors = np.ones(num_params)
    scale_factors[idx_d0 : idx_a0] = 1e-2   # mistuning   O(1e-2)
    scale_factors[idx_a0 : idx_fr] = 1e-3   # a0, a1      O(1e-3)
    scale_factors[idx_fr :]        = 1e-1   # forcing     O(1e-1)

    theta_initial_scaled = theta_initial / scale_factors
    lb_scaled            = lb / scale_factors
    ub_scaled            = ub / scale_factors
    bounds               = list(zip(lb_scaled, ub_scaled))

    evaluator = CachedEvaluator(ident_problem, scale_factors)

    print(f"Initial objective : {evaluator.objective(theta_initial_scaled):.6e}")
    print(f"Objective at true : {ident_problem.objective(theta_true):.6e}")

    # ------------------------------------------------------------------
    # 6.  IPOPT optimisation (164 parameters, scaled)
    # ------------------------------------------------------------------
    print(f"\nStarting IPOPT optimisation ({num_params} parameters, scaled)...")

    ipopt_options = {
        'max_iter'                  : 5000,
        'tol'                       : 1e-9,
        'acceptable_tol'            : 1e-8,
        'print_level'               : 5,
        'hessian_approximation'     : 'limited-memory',
        'limited_memory_max_history': 20,
        'limited_memory_update_type': 'sr1',
        'bound_push'                : 1e-3,
    }

    res = minimize_ipopt(
        fun    = evaluator.objective,
        x0     = theta_initial_scaled,
        jac    = evaluator.gradient,
        bounds = bounds,
        options= ipopt_options,
    )

    # ------------------------------------------------------------------
    # 7.  Results summary
    # ------------------------------------------------------------------
    print(f"\nStatus     : {res.message}")
    print(f"Iterations : {res.nit}  |  Evaluations: {res.nfev}")
    print(f"Final J    : {res.fun:.6e}")

    theta_opt = res.x * scale_factors

    # ------------------------------------------------------------------
    # 8.  Cramér-Rao bounds and diagnostic plots
    # ------------------------------------------------------------------
    print("\nComputing Cramér-Rao bounds...")
    sigma_theta = ident_problem.compute_cramer_rao_bounds(theta_opt, noise_level ** 2)

    print("Generating plots...")
    plotter = AeroPlotter(
        sim_problem   = ident_problem,
        omega_t       = omega_t,
        y_measured    = y_measured,
        N             = N,
        theta_true    = theta_true,
        theta_opt     = theta_opt,
        theta_initial = theta_initial,
        sigma_theta   = sigma_theta,
    )

    plotter.plot_response_envelope()
    plotter.plot_parameters_2x2()
    plotter.plot_influence_coefficients_tw()


if __name__ == '__main__':
    run_identification_with_simulated_data()
