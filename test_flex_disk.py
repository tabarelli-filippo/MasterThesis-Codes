"""
test_flex_disk.py
=================
Flexible-disk aeromechanical identification on synthetic flexible-disk data.

Workflow
--------
1.  Generate synthetic blade-vibration data with the flexible-disk model
    (dense K_bd, Schur condensation).
2.  Identify with FlexibleAeromechanicalIdentification using the *correct*
    disk matrices — this should recover the smooth TW coefficient curve
    with no artefact.
3.  Diagnose residual artefact by comparing the identified TW coefficients
    against the smooth physical reference.
4.  Compute Cramér-Rao bounds and generate three diagnostic plots.

Dependencies
------------
    auxSolver.flex_aero_ident — FlexibleAeromechanicalIdentification
    auxSolver.data_generation — generate_rotor_data
    auxSolver.plot_aero       — AeroPlotter
    cyipopt                   — Python interface to IPOPT
"""

import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter
from auxSolver.flex_aero_ident import FlexibleAeromechanicalIdentification


# ---------------------------------------------------------------------------
# Cached evaluator with parameter scaling
# ---------------------------------------------------------------------------

class CachedEvaluator:
    """
    Caching wrapper with diagonal parameter scaling for IPOPT.

    Parameters
    ----------
    ident_problem : FlexibleAeromechanicalIdentification
    scale_factors : array-like, shape (5N-1,)
    """

    def __init__(self, ident_problem, scale_factors):
        self.problem = ident_problem
        self.scales  = np.asarray(scale_factors, dtype=np.float64)
        self._key  = None
        self._J    = None
        self._grad = None

    def _update(self, x_scaled):
        key = x_scaled.tobytes()
        if key != self._key:
            theta = x_scaled * self.scales
            J, g  = self.problem.compute_objective_and_gradient(theta)
            self._J    = J
            self._grad = g * self.scales
            self._key  = key

    def objective(self, x):
        """Return the cached objective value."""
        self._update(x)
        return self._J

    def gradient(self, x):
        """Return the cached scaled gradient."""
        self._update(x)
        return self._grad


# ---------------------------------------------------------------------------
# IPOPT wrapper
# ---------------------------------------------------------------------------

def _run_ipopt(evaluator, x0_scaled, label=''):
    """
    Run an unconstrained IPOPT minimisation.

    Parameters
    ----------
    evaluator : CachedEvaluator
    x0_scaled : (5N-1,) float  — scaled initial guess
    label     : str            — label printed to stdout

    Returns
    -------
    res : OptimizeResult from cyipopt
    """
    print(f"\n  Starting {label} with IPOPT...")
    ipopt_options = {
        'max_iter'                  : 2000,
        'tol'                       : 1e-12,
        'acceptable_tol'            : 1e-10,
        'print_level'               : 5,
        'hessian_approximation'     : 'limited-memory',
        'limited_memory_max_history': 30,
        'limited_memory_update_type': 'bfgs',
    }
    res = minimize_ipopt(
        fun    = evaluator.objective,
        x0     = x0_scaled,
        jac    = evaluator.gradient,
        bounds = None,
        options= ipopt_options,
    )
    return res


# ---------------------------------------------------------------------------
# Main identification routine
# ---------------------------------------------------------------------------

def run_identification():
    """
    End-to-end flexible-disk identification pipeline.

    Uses the exact disk matrices (M_dd, K_dd, K_bd) from data generation,
    so the identifier model is correctly specified.  Any residual artefact
    in the TW coefficients is due to noise, not model mismatch.
    """

    # ------------------------------------------------------------------
    # 1.  Simulation parameters and data generation
    # ------------------------------------------------------------------
    N, omega_0, noise_level = 33, 17000.0, 1e-4
    omega_t = np.linspace(16250.0, 17750.0, 200)

    print("\n[1/4] Generating data: FLEXIBLE-DISK model ...")
    result   = generate_rotor_data('flexible', omega_t, omega_0, noise_level)
    y_measured            = result['y_measured']
    theta_true            = result['theta_true']
    M_dd, K_dd, K_bd = result['M_dd'], result['K_dd'], result['K_bd']

    # ------------------------------------------------------------------
    # 2.  Parameter scaling and initial guess
    # ------------------------------------------------------------------
    num_params = 5 * N - 1
    scales = np.ones(num_params)
    scales[0 : N]        = 1e-2   # d0    O(1e-2)
    scales[N : 3*N - 1]  = 1e-3   # a0,a1 O(1e-3)
    scales[3*N - 1 :]    = 1e-1   # fr,fi O(1e-1)

    np.random.seed(99)
    theta_initial        = 1e-5 * np.ones(num_params)
    theta_initial_scaled = theta_initial / scales

    # ------------------------------------------------------------------
    # 3.  Flexible-disk identification (no bounds, correct disk model)
    # ------------------------------------------------------------------
    print("\n[2/4] Running FLEXIBLE-DISK identification ...")
    ident_flex = FlexibleAeromechanicalIdentification(
        y_measured, omega_t, omega_0, N, M_dd, K_dd, K_bd)

    eval_flex     = CachedEvaluator(ident_flex, scales)
    res_flex      = _run_ipopt(eval_flex, theta_initial_scaled.copy(),
                               label='Flexible ID')
    theta_opt     = res_flex.x * scales

    # ------------------------------------------------------------------
    # 4.  Artefact check
    #     Physical aerodynamic coefficients (Table 1) form a smooth TW
    #     curve.  Deviations larger than the noise floor indicate that
    #     disk dynamics have leaked into the aerodynamic estimate.
    # ------------------------------------------------------------------
    a0_circ_true = np.zeros(N)
    a1_circ_true = np.zeros(N)
    a1_circ_true[0]   = 0.0050
    a0_circ_true[1]   = 0.00001;  a1_circ_true[1]   = 0.0010
    a0_circ_true[2]   = 0.00000;  a1_circ_true[2]   = 0.0015
    a0_circ_true[N-2] = 0.00000;  a1_circ_true[N-2] = 0.0020
    a0_circ_true[N-1] = 0.00000;  a1_circ_true[N-1] = 0.0025
    A_tw_true_smooth = np.fft.fft(a0_circ_true + 1j * a1_circ_true)

    _, a0_opt, a1_opt, _, _ = ident_flex.unpack_theta(theta_opt)
    A_tw_opt = np.fft.fft(a0_opt + 1j * a1_opt)
    print(f"\n  Artefact check — std dev relative to smooth reference:")
    print(f"    std(Real ΔA_tw) = "
          f"{np.std(np.real(A_tw_opt) - np.real(A_tw_true_smooth)):.4e}")
    print(f"    std(Imag ΔA_tw) = "
          f"{np.std(np.imag(A_tw_opt) - np.imag(A_tw_true_smooth)):.4e}")
    print(f"  (values >> noise={noise_level:.0e} indicate a significant artefact)\n")

    # ------------------------------------------------------------------
    # 5.  Cramér-Rao bounds
    # ------------------------------------------------------------------
    print("\n[3/4] Computing Cramér-Rao bounds ...")
    sigma_theta = ident_flex.compute_cramer_rao(theta_opt, noise_level ** 2)

    # ------------------------------------------------------------------
    # 6.  Diagnostic plots
    # ------------------------------------------------------------------
    print("\n[4/4] Generating plots (flexible-disk model)...")
    plotter = AeroPlotter(
        sim_problem   = ident_flex,
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
    plotter.plot_artifact_comparison()

    return theta_opt


if __name__ == '__main__':
    run_identification()
