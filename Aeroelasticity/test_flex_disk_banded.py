"""
test_flex_disk_banded.py
========================
Flexible-disk identification with modal truncation and banded aerodynamic
influence coefficients.

Background
----------
The true data are generated with a fully coupled flexible-disk model (dense
K_bd).  Identification is performed with a *truncated* coupling matrix: nodal
diameters above ``max_nd`` are zeroed out in K_bd, deliberately introducing a
model mismatch.  Band constraints on a0 and a1 (bandwidth B = 3) prevent the
optimiser from overfitting through high-ND aerodynamic terms.

This test reproduces the "unmodelled flexibility" scenario and quantifies how
much artefact leaks into the identified TW coefficient curve.

Workflow
--------
1.  Generate data with the full flexible-disk model.
2.  Apply modal truncation (max_nd = 10) to K_bd.
3.  Set up banded constraints on a0/a1 and isolate the ND=28 forcing.
4.  Run IPOPT with the truncated flexible-disk model.
5.  Compute Cramér-Rao bounds and generate diagnostic plots.

Dependencies
------------
    auxSolver.flex_aero_ident — FlexibleAeromechanicalIdentification
    auxSolver.data_generation — generate_rotor_data
    auxSolver.plot_aero       — AeroPlotter
    cyipopt                   — Python interface to IPOPT
"""

import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.flex_aero_ident import FlexibleAeromechanicalIdentification
from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter


# ---------------------------------------------------------------------------
# Modal truncation
# ---------------------------------------------------------------------------

def apply_modal_truncation(M_dd, K_dd, K_bd, N, max_nd=10):
    """
    Zero out blade-disk coupling for nodal diameters above ``max_nd``.

    Only K_bd is modified; M_dd and K_dd are left intact so the disk
    impedance matrix remains numerically invertible.  Rows/columns of K_bd
    with ND > max_nd are set to zero, simulating unmodelled disk flexibility
    at high nodal diameters.

    Parameters
    ----------
    M_dd   : (N, N) complex  — disk modal mass (returned unchanged)
    K_dd   : (N, N) complex  — disk modal stiffness (returned unchanged)
    K_bd   : (N, N) complex  — blade-disk coupling to truncate
    N      : int             — number of blades
    max_nd : int             — maximum retained nodal diameter

    Returns
    -------
    M_dd_trunc, K_dd_trunc, K_bd_trunc : copies with truncated K_bd
    """
    M_dd_trunc = np.copy(M_dd)
    K_dd_trunc = np.copy(K_dd)
    K_bd_trunc = np.copy(K_bd)

    for i in range(N):
        nd = i if i <= N // 2 else N - i
        if nd > max_nd:
            K_bd_trunc[i, :] = 0.0

    return M_dd_trunc, K_dd_trunc, K_bd_trunc


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

def _run_ipopt(evaluator, x0_scaled, bounds, label=''):
    """
    Run a bounded IPOPT minimisation.

    Parameters
    ----------
    evaluator : CachedEvaluator
    x0_scaled : (5N-1,) float          — scaled initial guess
    bounds    : list of (lb, ub) pairs  — scaled parameter bounds
    label     : str

    Returns
    -------
    res : OptimizeResult from cyipopt
    """
    print(f"\n  Starting {label} with IPOPT...")
    ipopt_options = {
        'max_iter'                  : 5000,
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
        bounds = bounds,
        options= ipopt_options,
    )
    flag = "✓" if res.success else "~"
    print(f"  {flag} J={res.fun:.4e}  nit={res.nit}  success={res.success}")
    return res


# ---------------------------------------------------------------------------
# Main identification routine
# ---------------------------------------------------------------------------

def run_identification():
    """
    End-to-end identification with modal truncation and banded constraints.
    """

    # ------------------------------------------------------------------
    # 1.  Simulation parameters and data generation
    # ------------------------------------------------------------------
    N           = 33
    omega_0     = 17000.0
    noise_level = 1e-4
    omega_t     = np.linspace(16250.0, 17750.0, 200)

    print("\n[1/5] Generating data: FULL FLEXIBLE-DISK model ...")
    result   = generate_rotor_data('flexible', omega_t, omega_0, noise_level)
    y_measured            = result['y_measured']
    theta_true            = result['theta_true']
    M_dd_true, K_dd_true, K_bd_true = result['M_dd'], result['K_dd'], result['K_bd']

    # ------------------------------------------------------------------
    # 2.  Modal truncation of K_bd (unmodelled flexibility)
    # ------------------------------------------------------------------
    print("\n[2/5] Applying modal truncation (max_nd=10) ...")
    M_dd_trunc, K_dd_trunc, K_bd_trunc = apply_modal_truncation(
        M_dd_true, K_dd_true, K_bd_true, N, max_nd=10
    )

    # ------------------------------------------------------------------
    # 3.  Parameter scaling and block offsets
    # ------------------------------------------------------------------
    num_params = 5 * N - 1
    scales = np.ones(num_params)

    idx_d0 = 0               # mistuning     length N
    idx_a0 = N               # aero stiff.   length N-1
    idx_a1 = 2 * N - 1       # aero damping  length N
    idx_fr = 3 * N - 1       # real forcing  length N
    idx_fi = 4 * N - 1       # imag forcing  length N

    scales[idx_d0 : idx_a0] = 1e-2
    scales[idx_a0 : idx_fr] = 1e-3
    scales[idx_fr :]        = 1e-1

    # ------------------------------------------------------------------
    # 4.  Parameter bounds
    # ------------------------------------------------------------------
    print("\n[3/5] Setting bounds: banded influence matrix (B=3) + ND=28 forcing ...")
    bounds = [(-np.inf, np.inf)] * num_params

    # 4a.  Banded constraints on a0 and a1 (bandwidth B = 3).
    #      Indices outside the band are fixed to zero, preventing the
    #      optimiser from recovering disk dynamics through high-ND aero terms.
    B = 3
    for i in range(1, N):
        if i > B and i < N - B:
            bounds[idx_a0 + i - 1] = (0.0, 0.0)

    for i in range(N):
        if i > B and i < N - B:
            bounds[idx_a1 + i] = (0.0, 0.0)

    # 4b.  Single-ND forcing: only ND=28 is free, all other NDs are fixed to zero.
    target_nd = 28
    for k in range(N):
        if k != target_nd:
            bounds[idx_fr + k] = (0.0, 0.0)
            bounds[idx_fi + k] = (0.0, 0.0)

    # 4c.  Physical stability: self-damping must be positive.
    bounds[idx_a1] = (1e-6, np.inf)

    # ------------------------------------------------------------------
    # 5.  Initial guess (consistent with bounds)
    # ------------------------------------------------------------------
    np.random.seed(99)
    theta_initial = 1e-5 * np.ones(num_params)
    for idx, (lo, hi) in enumerate(bounds):
        if lo == 0.0 and hi == 0.0:
            theta_initial[idx] = 0.0

    theta_initial_scaled = theta_initial / scales

    # ------------------------------------------------------------------
    # 6.  Identification with truncated flexible-disk model + banded bounds
    # ------------------------------------------------------------------
    print("\n[4/5] Running FLEXIBLE identification (truncated model + banded bounds)...")
    ident_flex = FlexibleAeromechanicalIdentification(
        y_measured, omega_t, omega_0, N,
        M_dd_trunc, K_dd_trunc, K_bd_trunc)

    eval_flex = CachedEvaluator(ident_flex, scales)
    res_flex  = _run_ipopt(eval_flex, theta_initial_scaled.copy(), bounds,
                           label='Truncated + Banded ID')
    theta_opt = res_flex.x * scales

    # ------------------------------------------------------------------
    # 7.  Artefact diagnostics
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
    # 8.  Cramér-Rao bounds and diagnostic plots
    # ------------------------------------------------------------------
    print("\n[5/5] Computing Cramér-Rao bounds and generating plots...")
    sigma_theta = ident_flex.compute_cramer_rao(theta_opt, noise_level ** 2)

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
