"""
data_generation.py
==================
Unified module for generating simulated blade-vibration data.

Two physics models are supported via the ``mode`` parameter:

    mode='rigid'    — rigid-disk model
    mode='flexible' — flexible-disk model with a dense K_bd in TW space

Public interface
----------------
    result = generate_rotor_data(mode, omega_t, omega_0, noise_level, N=33, seed=42)

Return value (dict):
    'y_measured' : (T, N) complex  — noisy measured displacements
    'theta_true' : (5N-1,) float   — true parameters in the rigid-ID layout
    'M_dd'       : (N, N) or None  — disk modal mass   (flexible only)
    'K_dd'       : (N, N) or None  — disk modal stiffness (flexible only)
    'K_bd'       : (N, N) or None  — blade-disk coupling  (flexible only)

The keys 'M_dd', 'K_dd', 'K_bd' are always present; for mode='rigid' they are
None, so callers can unpack uniformly without branching.

Backward-compatible aliases
---------------------------
    initialize_rotor(...)                      ← data_generation_rigid
    initialize_flexible_rotor_complete(...)    ← data_generation_flexible
    initialize_flexible_rotor(...)             ← data_generation_flexible (legacy)

Rigid-disk physics
------------------
    [A_blade - gamma^2 * I] x = f

    A_blade built from (d0, a0_circ, a1_circ) following the
    AeromechanicalIdentification normalisation.  K_bd = 0 (no Schur condensation).

Flexible-disk physics
---------------------
    [A_blade - gamma^2 * I - D_disk(omega) / omega_0^2] x = E @ f_q

    gamma    = omega / omega_0
    D_disk   = K_bd @ inv(K_dd - omega^2 * M_dd) @ K_db

    K_bd is *dense* in TW space (exponential decay in ΔND), making D_disk
    non-circulant.  This is precisely the setting that reveals the mathematical
    artefact of Hall & Hall (2024) when a rigid-disk model is used for
    identification.
"""

import numpy as np
import scipy.linalg


# ---------------------------------------------------------------------------
# Shared physical parameters (Table 1 of the paper)
# ---------------------------------------------------------------------------

_D0_ARRAY = np.array([
    -0.00405, -0.01262,  0.00791,  0.00500,  0.00923,  0.00476, -0.00110, -0.00332,
    -0.00650,  0.00231, -0.00478,  0.01178,  0.00845, -0.00540, -0.00029, -0.00510,
     0.00836,  0.01254,  0.00465,  0.00314, -0.00444, -0.00786,  0.00601,  0.00515,
     0.01004,  0.00890, -0.00042,  0.00385,  0.00845,  0.01111,  0.00225,  0.00592,
     0.01173,
], dtype=np.float64)


def _build_aero_circulant(N):
    """
    Build the aerodynamic circulant first-columns (a0_circ, a1_circ) for N blades.

    Coefficients are taken from Table 1 for N=33.  For N != 33 the non-zero
    entries are mapped consistently with circular wrapping.

    Returns
    -------
    a0 : (N,) float  — first column of the stiffness circulant
    a1 : (N,) float  — first column of the damping circulant
    """
    a0 = np.zeros(N)
    a1 = np.zeros(N)
    a1[0]   = 0.0050
    a0[1]   = 0.00001;  a1[1]   = 0.0010
    a0[2]   = 0.00000;  a1[2]   = 0.0015
    a0[N-2] = 0.00000;  a1[N-2] = 0.0020
    a0[N-1] = 0.00000;  a1[N-1] = 0.0025
    return a0, a1


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def generate_rotor_data(mode, omega_t, omega_0, noise_level=1e-4, N=33, seed=42):
    """
    Generate simulated blade-vibration data for identification experiments.

    Parameters
    ----------
    mode        : {'rigid', 'flexible'}
    omega_t     : array-like, shape (T,)  — excitation frequencies [rad/s]
    omega_0     : float                   — nominal blade frequency [rad/s]
    noise_level : float                   — complex noise standard deviation
    N           : int                     — number of blades (rigid mode only;
                                            flexible mode is fixed at N=33)
    seed        : int                     — random seed for reproducibility

    Returns
    -------
    dict with keys:
        'y_measured' : (T, N) complex
        'theta_true' : (5N-1,) float
        'M_dd'       : (N, N) complex or None
        'K_dd'       : (N, N) complex or None
        'K_bd'       : (N, N) complex or None
    """
    if mode == 'rigid':
        y_measured, theta_true = _generate_rigid(N, omega_0, omega_t, noise_level, seed)
        return {
            'y_measured': y_measured,
            'theta_true': theta_true,
            'M_dd': None,
            'K_dd': None,
            'K_bd': None,
        }
    elif mode == 'flexible':
        y_measured, theta_true, M_dd, K_dd, K_bd = \
            _generate_flexible(omega_t, omega_0, noise_level, seed)
        return {
            'y_measured': y_measured,
            'theta_true': theta_true,
            'M_dd': M_dd,
            'K_dd': K_dd,
            'K_bd': K_bd,
        }
    else:
        raise ValueError(f"mode must be 'rigid' or 'flexible', got: {mode!r}")


# ---------------------------------------------------------------------------
# Rigid-disk implementation
# ---------------------------------------------------------------------------

def _generate_rigid(N, omega_0, omega_t, noise_level, seed):
    """
    Generate data from the rigid-disk model.

        [A_blade - gamma^2 * I] x = f

    Parameters match _generate_flexible; see generate_rotor_data for details.
    """
    # Local import to avoid circular dependencies
    from auxSolver.aero_ident import AeromechanicalIdentification

    np.random.seed(seed)
    T = len(omega_t)
    num_params = 5 * N - 1

    # Block offsets in theta
    idx_d0 = 0
    idx_a0 = N
    idx_a1 = N + (N - 1)
    idx_fr = N + (N - 1) + N
    idx_fi = N + (N - 1) + N + N

    theta_true = np.zeros(num_params, dtype=np.float64)

    # Stiffness mistuning
    d0 = _D0_ARRAY if N == len(_D0_ARRAY) else np.resize(_D0_ARRAY, N)
    theta_true[idx_d0: idx_d0 + N] = d0

    # Aerodynamic influence coefficients
    a0_circ, a1_circ = _build_aero_circulant(N)
    theta_true[idx_a0: idx_a0 + N - 1] = a0_circ[1:]   # a0_circ[0] = 0 by definition
    theta_true[idx_a1: idx_a1 + N]     = a1_circ

    # Single-ND forcing at ND = 28
    target_nd = 28
    if target_nd < N:
        theta_true[idx_fr + target_nd] = -0.1
        theta_true[idx_fi + target_nd] =  0.0

    # Exact response via AeromechanicalIdentification
    y_dummy     = np.zeros((T, N), dtype=np.complex128)
    sim_problem = AeromechanicalIdentification(y_dummy, omega_t, omega_0, N)

    d_0, a_0, a_1, f_r, f_i = sim_problem.unpack_theta(theta_true)
    A, f = sim_problem.build_system_matrices(d_0, a_0, a_1, f_r, f_i)

    mu_exact = np.zeros((T, N), dtype=np.complex128)
    for t in range(T):
        gamma_sq    = (omega_t[t] / omega_0) ** 2
        H_t         = A - gamma_sq * sim_problem.I_N
        mu_exact[t] = np.linalg.solve(H_t, f)

    rng = np.random.RandomState(seed)
    y_measured = (mu_exact
                  + rng.normal(0, noise_level, (T, N))
                  + 1j * rng.normal(0, noise_level, (T, N)))

    return y_measured, theta_true


# ---------------------------------------------------------------------------
# Flexible-disk implementation
# ---------------------------------------------------------------------------

def _generate_flexible(omega_t, omega_0, noise_level, seed):
    """
    Generate data from the flexible-disk model with a dense K_bd.

        [A_blade - gamma^2 * I - D_disk(omega) / omega_0^2] x = E @ f_q

    The blade-disk coupling K_bd decays exponentially with nodal-diameter
    separation, making D_disk non-circulant and introducing the model-mismatch
    artefact described in Hall & Hall (2024) when rigid-disk identification
    is applied to this data.

    N is fixed at 33 to match the physical parameter tables.
    """
    N = 33
    T = len(omega_t)

    a0_circ, a1_circ = _build_aero_circulant(N)
    d0_true          = _D0_ARRAY.copy()

    # Modal forcing at ND = 28
    f_q = np.zeros(N, dtype=np.complex128)
    f_q[28] = -0.1 + 0.0j

    # ------------------------------------------------------------------
    # Disk modal structure
    # ------------------------------------------------------------------
    M_dd = np.eye(N, dtype=np.complex128) * 50.0

    nd_array   = np.fft.fftfreq(N, d=1.0 / N)
    disk_freqs = omega_0 * (0.35 + 0.13 * np.abs(nd_array))
    # Small structural damping (1e-4) included via complex stiffness
    K_dd = np.diag(50.0 * disk_freqs ** 2 * (1.0 + 1j * 1e-4)).astype(np.complex128)

    # ------------------------------------------------------------------
    # Blade-disk coupling (dense in TW space — exponential decay in ΔND)
    # ------------------------------------------------------------------
    coupling_amplitude = 0.02 * omega_0 ** 2
    coupling_length    = 3.0          # decay length in nodal diameters

    K_bd = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            delta_nd   = min(abs(i - j), N - abs(i - j))
            K_bd[i, j] = coupling_amplitude * np.exp(-delta_nd / coupling_length)

    # Zero out coupling for ND 12..21 (structural gap in the disk design)
    for idx in range(12, 22):
        K_bd[idx, :] = 0.0
        K_bd[:, idx] = 0.0

    K_db = K_bd.conj().T

    # ------------------------------------------------------------------
    # System matrices
    # ------------------------------------------------------------------
    j_idx, k_idx = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    E = (1.0 / np.sqrt(N)) * np.exp(1j * 2 * np.pi * j_idx * k_idx / N)

    A_blade = (np.eye(N, dtype=np.complex128)
               + np.diag(d0_true)
               + scipy.linalg.circulant(a0_circ)
               + 1j * scipy.linalg.circulant(a1_circ))

    f_phys = E @ f_q

    # ------------------------------------------------------------------
    # Harmonic solution via Schur condensation
    # ------------------------------------------------------------------
    y_exact = np.zeros((T, N), dtype=np.complex128)
    for t in range(T):
        w           = omega_t[t]
        gamma_sq    = (w / omega_0) ** 2
        Z_dd        = K_dd - w ** 2 * M_dd
        D_disk_norm = K_bd @ np.linalg.solve(Z_dd, K_db) / omega_0 ** 2
        Z_t         = A_blade - gamma_sq * np.eye(N, dtype=np.complex128) - D_disk_norm
        y_exact[t]  = np.linalg.solve(Z_t, f_phys)

    rng = np.random.RandomState(seed)
    y_measured = (y_exact
                  + rng.normal(0, noise_level, (T, N))
                  + 1j * rng.normal(0, noise_level, (T, N)))

    # ------------------------------------------------------------------
    # True parameter vector in the rigid-ID layout
    # ------------------------------------------------------------------
    theta_true = np.concatenate([
        d0_true,
        a0_circ[1:],    # N-1 elements
        a1_circ,        # N   elements
        np.real(f_q),
        np.imag(f_q),
    ])

    return y_measured, theta_true, M_dd, K_dd, K_bd


# ---------------------------------------------------------------------------
# Backward-compatible aliases
# ---------------------------------------------------------------------------

def initialize_rotor(N, omega_0, omega_t, noise_level, seed=42):
    """Backward-compatible alias for data_generation_rigid.initialize_rotor."""
    res = generate_rotor_data('rigid', omega_t, omega_0, noise_level, N=N, seed=seed)
    return res['y_measured'], res['theta_true']


def initialize_flexible_rotor_complete(omega_t, omega_0, noise_level=1e-4, seed=42):
    """Backward-compatible alias for data_generation_flexible.initialize_flexible_rotor_complete."""
    res = generate_rotor_data('flexible', omega_t, omega_0, noise_level, seed=seed)
    return res['y_measured'], res['theta_true'], res['M_dd'], res['K_dd'], res['K_bd']


def initialize_flexible_rotor(omega_t, omega_0, noise_level, seed=42,
                               M_dd=None, K_dd=None, M_bd=None, K_bd=None):
    """Backward-compatible alias for data_generation_flexible.initialize_flexible_rotor (legacy)."""
    res = generate_rotor_data('flexible', omega_t, omega_0, noise_level, seed=seed)
    return res['y_measured'], res['theta_true']
