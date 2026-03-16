"""
flex_aero_ident.py
==================
Maximum-likelihood aeromechanical identification for a flexible-disk rotor.

The blade response model (Eq. 13 extended with flexible-disk condensation,
Sec. 4.3) is:

    [A(theta) - gamma_t^2 * I - D_disk(omega_t) / omega_0^2] x_t = E @ f_q

where:
    A(theta) = I + diag(d0) + circ(a0) + 1j*circ(a1)
    gamma_t  = omega_t / omega_0
    D_disk   = K_bd @ inv(K_dd - omega_t^2 * M_dd) @ K_db  (pre-computed)

The disk influence term D_disk is computed once at construction time for all
T frequencies and stored as a (T, N, N) array.

Parameter vector layout  (5N-1 real entries):
    theta = [d0(0..N-1) | a0_params(0..N-2) | a1_params(0..N-1) |
             fr(0..N-1) | fi(0..N-1)]

    a0_params[k] = a0_circ[k+1]   for k = 0..N-2
    a1_params[0] = a1_circ[0]     (self-damping)
    a1_params[k] = a1_circ[k]     for k = 1..N-1

Note on scipy.linalg.circulant(c):
    C[i, j] = c[(i - j) mod N]  — c is the *first column*.
    Consequently C[0, :] = [c[0], c[N-1], ..., c[1]]  (first row ≠ first column).
"""

import numpy as np
import scipy.linalg


class FlexibleAeromechanicalIdentification:
    """
    Flexible-disk aeromechanical identification via maximum likelihood.

    Parameters
    ----------
    y_measured : array-like, shape (T, N), complex
        Measured blade displacements at T excitation frequencies.
    omega_t : array-like, shape (T,), float
        Excitation frequencies [rad/s].
    omega_0 : float
        Nominal blade natural frequency [rad/s].
    N : int
        Number of blades.
    M_dd : array-like, shape (N, N), complex
        Modal mass matrix of the disk.
    K_dd : array-like, shape (N, N), complex
        Modal stiffness matrix of the disk.
    K_bd : array-like, shape (N, N), complex
        Blade-disk coupling stiffness matrix (dense in TW space).
    """

    def __init__(self, y_measured, omega_t, omega_0, N, M_dd, K_dd, K_bd):
        self.y = np.asarray(y_measured, dtype=np.complex128)
        self.omega_t = np.asarray(omega_t, dtype=np.float64)
        self.omega_0 = float(omega_0)
        self.N = N
        self.T = len(omega_t)
        self.I_N = np.eye(N, dtype=np.complex128)

        # Pre-compute the IDFT matrix E (Eq. 11)
        j, k = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        self.E = (1.0 / np.sqrt(N)) * np.exp(1j * 2 * np.pi * j * k / N)

        # Block index boundaries in theta (5N-1) — exposed for external scripts
        self.idx_d0 = 0
        self.idx_a0 = N
        self.idx_a1 = 2 * N - 1
        self.idx_fr = 3 * N - 1
        self.idx_fi = 4 * N - 1

        # Pre-compute the normalised disk influence matrices D_disk[t] / omega_0^2.
        # Shape: (T, N, N).  This avoids repeating the Schur complement at every
        # objective evaluation.
        K_db = K_bd.conj().T
        self.D_disk = np.zeros((self.T, N, N), dtype=np.complex128)
        for t in range(self.T):
            Z_dd = K_dd - (self.omega_t[t] ** 2) * M_dd
            self.D_disk[t] = K_bd @ np.linalg.solve(Z_dd, K_db) / (self.omega_0 ** 2)

    # ------------------------------------------------------------------
    # Parameter packing / unpacking
    # ------------------------------------------------------------------

    def unpack_theta(self, theta):
        """
        Unpack the flat parameter vector into physical sub-arrays.

        Returns
        -------
        d_0      : (N,)   diagonal stiffness mistuning
        a_0_circ : (N,)   first column of circ(a0); a_0_circ[0] = 0 by definition
        a_1_circ : (N,)   first column of circ(a1); index 0 = self-damping
        f_r      : (N,)   real part of the TW forcing vector
        f_i      : (N,)   imaginary part of the TW forcing vector
        """
        N = self.N
        d_0 = theta[0:N]

        # a0: N-1 free parameters = first column of circ(a0) at distances 1..N-1
        a_0_circ = np.zeros(N)
        a_0_circ[1:] = theta[N: 2 * N - 1]

        # a1: N parameters; index 0 = self-damping, indices 1..N-1 = distances
        a_1_circ = theta[2 * N - 1: 3 * N - 1].copy()

        f_r = theta[3 * N - 1: 4 * N - 1]
        f_i = theta[4 * N - 1: 5 * N - 1]
        return d_0, a_0_circ, a_1_circ, f_r, f_i

    def pack_theta(self, d_0, a_0_circ, a_1_circ, f_r, f_i):
        """Inverse of unpack_theta."""
        return np.concatenate([d_0, a_0_circ[1:], a_1_circ, f_r, f_i])

    # ------------------------------------------------------------------
    # System matrix construction
    # ------------------------------------------------------------------

    def build_system_matrices(self, d_0, a_0_circ, a_1_circ, f_r, f_i):
        """
        Build the blade system matrix A and the physical forcing vector f.

            A = I + diag(d0) + circ(a0) + 1j * circ(a1)
            f = E @ (f_r + 1j * f_i)

        Note: scipy.linalg.circulant(c) treats c as the first column.
        """
        A = (self.I_N
             + np.diag(d_0)
             + scipy.linalg.circulant(a_0_circ)
             + 1j * scipy.linalg.circulant(a_1_circ))
        f = self.E @ (f_r + 1j * f_i)
        return A, f

    def solve_eigendecomposition(self, A):
        """
        Bi-orthogonal eigendecomposition: A = X @ diag(lam) @ Y, Y @ X = I.

        Returns
        -------
        lam : (N,)   eigenvalues
        X   : (N, N) right eigenvectors (columns)
        Y   : (N, N) left eigenvectors (rows)
        """
        lam, X = np.linalg.eig(A)
        Y = np.linalg.inv(X)
        return lam, X, Y

    # ------------------------------------------------------------------
    # Objective and gradient
    # ------------------------------------------------------------------

    def compute_objective_and_gradient(self, theta):
        """
        Evaluate the normalised least-squares objective and its gradient.

            J(theta) = sum_t ||y_t - x_t||^2 / ||y||^2

        The gradient is computed via the adjoint method (Eqs. 33-43).

        Adjoint system:
            Z_t^H lambda_t = r_t,   r_t = y_t - x_t,
            Z_t = A - gamma_t^2 I - D_disk[t]

        Sensitivity accumulation:
            G_A = sum_t outer(x_t, conj(lambda_t))

        Gradient expressions:
            dJ/d(d0)_j      =  2 Re{ G_A[j, j] }
            dJ/d(a0_circ)_k =  2 Re{ sum_{p-q=k mod N} G_A[q, p] }
            dJ/d(a1_circ)_k = -2 Im{ sum_{p-q=k mod N} G_A[q, p] }
            dJ/d(fr)        = -2 Re{ E^H @ lambda_sum }
            dJ/d(fi)        = -2 Im{ E^H @ lambda_sum }

        Returns
        -------
        J        : float   normalised objective
        gradient : (5N-1,) gradient vector
        """
        d_0, a_0_circ, a_1_circ, f_r, f_i = self.unpack_theta(theta)
        A_blade, f = self.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)

        N = self.N
        J = 0.0
        lambda_sum = np.zeros(N, dtype=np.complex128)
        G_A = np.zeros((N, N), dtype=np.complex128)

        denom = np.sum(np.abs(self.y) ** 2) + 1e-30

        for t in range(self.T):
            gamma_sq = (self.omega_t[t] / self.omega_0) ** 2

            # Full system matrix including the flexible-disk correction
            Z_t = A_blade - gamma_sq * self.I_N - self.D_disk[t]

            x_t = np.linalg.solve(Z_t, f)
            r_t = self.y[t, :] - x_t
            J += np.real(np.vdot(r_t, r_t))

            # Adjoint solve: Z_t^H lambda_t = r_t
            lambda_t = np.linalg.solve(Z_t.conj().T, r_t)

            lambda_sum += lambda_t
            G_A += np.outer(x_t, lambda_t.conj())

        # Forcing gradient
        # Derivation: dJ/df_r = -2 Re{ (E^H lambda_sum)^H }
        S = self.E.conj().T @ lambda_sum
        grad_f_r = -2.0 * np.real(S)
        grad_f_i = -2.0 * np.imag(S)

        # Mistuning gradient (Eq. 40)
        grad_d_0 = 2.0 * np.real(np.diag(G_A))

        # Circulant gradients (Eqs. 42-43).
        # For circ(a) with first-column convention: C[p, q] = a[(p-q) mod N].
        # The anti-diagonal sum at offset k in G_A equals
        #   sum_{p-q=k mod N} G_A[q, p].
        # np.trace(G_A, offset=k) sums entries G_A[i, i+k],
        # i.e. those with q = i, p = i+k → p-q = k.
        # The wrap-around part uses offset = k - N.
        grad_a_0_circ = np.zeros(N)
        grad_a_1_circ = np.zeros(N)
        for k in range(N):
            if k == 0:
                diag_sum = np.trace(G_A)
            else:
                diag_sum = np.trace(G_A, offset=k) + np.trace(G_A, offset=k - N)
            grad_a_0_circ[k] =  2.0 * np.real(diag_sum)
            grad_a_1_circ[k] = -2.0 * np.imag(diag_sum)

        # Map circulant gradients to the theta layout:
        #   a0_params = a0_circ[1:]  →  gradient uses indices 1..N-1
        #   a1_params = a1_circ      →  gradient uses all N indices
        grad_a_0_params = grad_a_0_circ[1:]
        grad_a_1_params = grad_a_1_circ

        gradient = np.concatenate([
            grad_d_0, grad_a_0_params, grad_a_1_params, grad_f_r, grad_f_i
        ])

        return J / denom, gradient / denom

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def solve_frequency_response(self, theta, t):
        """
        Compute the blade response x_t for a single frequency index t.

        Parameters
        ----------
        theta : (5N-1,) float  — parameter vector
        t     : int            — frequency index into omega_t

        Returns
        -------
        x_t : (N,) complex
        """
        d_0, a_0_circ, a_1_circ, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)
        gamma_sq = (self.omega_t[t] / self.omega_0) ** 2
        Z_t = A - gamma_sq * self.I_N - self.D_disk[t]
        return np.linalg.solve(Z_t, f)

    def compute_cramer_rao(self, theta, sigma_noise):
        """
        Compute the Cramér-Rao lower bounds via finite-difference sensitivities.

        Fisher information matrix (Eqs. 21-22):
            I(theta) = (1/sigma^2) * sum_t (dmu_t/dtheta)^H (dmu_t/dtheta)

        Sensitivities dmu_t/dtheta are approximated by forward finite differences
        with step eps = 1e-7.

        Parameters
        ----------
        theta       : (5N-1,) float  — parameter estimate at which to evaluate
        sigma_noise : float          — measurement noise standard deviation

        Returns
        -------
        crb : (5N-1,) float  — standard-deviation lower bounds per parameter
        """
        eps = 1e-7
        num_params = len(theta)

        # Baseline responses at current theta
        d_0, a_0_circ, a_1_circ, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)

        mu = np.zeros((self.T, self.N), dtype=np.complex128)
        for t in range(self.T):
            gamma_sq = (self.omega_t[t] / self.omega_0) ** 2
            Z_t = A - gamma_sq * self.I_N - self.D_disk[t]
            mu[t] = np.linalg.solve(Z_t, f)

        # Finite-difference sensitivities
        Fisher = np.zeros((num_params, num_params))
        dmu = np.zeros((num_params, self.T, self.N), dtype=np.complex128)

        for i in range(num_params):
            theta_p = theta.copy(); theta_p[i] += eps
            d_0p, a_0p, a_1p, frp, fip = self.unpack_theta(theta_p)
            Ap, fp = self.build_system_matrices(d_0p, a_0p, a_1p, frp, fip)
            for t in range(self.T):
                gamma_sq = (self.omega_t[t] / self.omega_0) ** 2
                Z_t = Ap - gamma_sq * self.I_N - self.D_disk[t]
                mu_p = np.linalg.solve(Z_t, fp)
                dmu[i, t] = (mu_p - mu[t]) / eps

        for i in range(num_params):
            for j in range(i, num_params):
                val = np.sum(np.real(np.einsum('ti,ti->', dmu[i].conj(), dmu[j])))
                Fisher[i, j] = val / sigma_noise ** 2
                Fisher[j, i] = Fisher[i, j]

        try:
            Fisher_inv = np.linalg.inv(Fisher + 1e-20 * np.eye(num_params))
            crb = np.sqrt(np.maximum(np.diag(Fisher_inv), 0.0))
        except np.linalg.LinAlgError:
            crb = np.full(num_params, np.nan)

        return crb
