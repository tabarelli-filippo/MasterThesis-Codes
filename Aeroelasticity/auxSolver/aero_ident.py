"""
aero_ident.py
=============
Maximum-likelihood aeromechanical identification for a rigid-disk rotor.

The blade response model is (see Hall & Hall 2024, Eq. 13):

    [A(theta) - gamma_t^2 * I] x_t = f

where:
    A(theta) = I + diag(d0) + circ(a0) + 1j*circ(a1)
    gamma_t  = omega_t / omega_0
    f        = E @ (f_r + 1j*f_i)

Parameter vector layout  (5N-1 real entries):
    theta = [d0(0..N-1) | a0_params(1..N-1) | a1_params(0..N-1) |
             fr(0..N-1) | fi(0..N-1)]

Note: a0_params[k] = a0_circ[k+1] for k=0..N-2  (a0_circ[0] = 0 by definition).

Classes
-------
AeromechanicalIdentification
    Core identification class: builds system matrices, evaluates the
    least-squares objective, and computes its adjoint gradient.

CachedEvaluator
    Thin wrapper that caches the last (J, grad) evaluation to avoid
    redundant computations when an optimizer queries objective and
    gradient separately.
"""

import numpy as np
import scipy.linalg


class AeromechanicalIdentification:
    """
    Rigid-disk aeromechanical identification via maximum likelihood.

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
    """

    def __init__(self, y_measured, omega_t, omega_0, N):
        self.y = np.asarray(y_measured, dtype=np.complex128)
        self.omega_t = np.asarray(omega_t, dtype=np.float64)
        self.omega_0 = omega_0
        self.N = N
        self.T = len(omega_t)

        # Pre-compute the IDFT matrix E (Eq. 11 of the paper)
        j_idx, k_idx = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        self.E = (1.0 / np.sqrt(N)) * np.exp(1j * 2 * np.pi * j_idx * k_idx / N)
        self.I_N = np.eye(N, dtype=np.complex128)

    # ------------------------------------------------------------------
    # Parameter packing / unpacking
    # ------------------------------------------------------------------

    def unpack_theta(self, theta):
        """
        Unpack the flat parameter vector into physical sub-arrays.

        Returns
        -------
        d_0 : (N,)   diagonal stiffness mistuning
        a_0 : (N,)   first column of circ(a0); a_0[0] = 0 by definition
        a_1 : (N,)   first column of circ(a1); a_1[0] = self-damping term
        f_r : (N,)   real part of the forcing vector in TW coordinates
        f_i : (N,)   imaginary part of the forcing vector
        """
        N = self.N
        idx = 0
        d_0 = theta[idx: idx + N]; idx += N
        a_0 = np.zeros(N)              # a0[0] = 0 enforced
        a_0[1:] = theta[idx: idx + (N - 1)]; idx += (N - 1)
        a_1 = theta[idx: idx + N]; idx += N
        f_r = theta[idx: idx + N]; idx += N
        f_i = theta[idx: idx + N]
        return d_0, a_0, a_1, f_r, f_i

    # ------------------------------------------------------------------
    # System matrix construction
    # ------------------------------------------------------------------

    def build_system_matrices(self, d_0, a_0, a_1, f_r, f_i):
        """
        Build the blade system matrix A and the physical forcing vector f.

            A = I + diag(d0) + circ(a0) + 1j * circ(a1)
            f = E @ (f_r + 1j * f_i)

        Note: scipy.linalg.circulant(c) uses c as the *first column*,
        so C[i, j] = c[(i - j) mod N].
        """
        D_0 = np.diag(d_0)
        A_0_mat = scipy.linalg.circulant(a_0)
        A_1_mat = scipy.linalg.circulant(a_1)
        A = self.I_N + D_0 + A_0_mat + 1j * A_1_mat
        f = self.E @ (f_r + 1j * f_i)
        return A, f

    def solve_eigendecomposition(self, A):
        """
        Bi-orthogonal eigendecomposition A = X @ diag(lam) @ Y, Y @ X = I.

        Used for O(N^2) frequency-loop solves instead of O(N^3).

        Returns
        -------
        lam : (N,)   eigenvalues
        X   : (N, N) right eigenvectors (columns)
        Y   : (N, N) left eigenvectors (rows), normalised so that Y @ X = I
        """
        M = self.I_N
        lam, vl, X = scipy.linalg.eig(A, b=M, left=True, right=True)
        Y_raw = vl.conj().T
        diag_YMX = np.diag(Y_raw @ M @ X)
        Y = Y_raw / diag_YMX[:, np.newaxis]
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
            Z_t^H lambda_t = r_t,   r_t = y_t - x_t

        Accumulated sensitivity matrix:
            G_A = sum_t outer(x_t, lambda_t^*)

        Gradient contributions:
            dJ/d(d0)_j        =  2 Re{ G_A[j, j] }                    (Eq. 40)
            dJ/d(a0_circ)_k   =  2 Re{ sum_{p-q=k mod N} G_A[q, p] }  (Eq. 43)
            dJ/d(a1_circ)_k   = -2 Im{ sum_{p-q=k mod N} G_A[q, p] }  (Eq. 42)
            dJ/d(fr)          = -2 Re{ (E^H lambda_sum) }
            dJ/d(fi)          = -2 Im{ (E^H lambda_sum) }

        Returns
        -------
        J        : float   objective value
        gradient : (5N-1,) gradient vector
        """
        d_0, a_0, a_1, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0, a_1, f_r, f_i)

        lam, X, Y = self.solve_eigendecomposition(A)
        Y_f = Y @ f

        J = 0.0
        g_f = np.zeros(self.N, dtype=np.complex128)
        G_A = np.zeros((self.N, self.N), dtype=np.complex128)

        for t in range(self.T):
            gamma_sq = (self.omega_t[t] / self.omega_0) ** 2
            diag_inv = 1.0 / (lam - gamma_sq)

            x_t = X @ (diag_inv * Y_f)
            r_t = self.y[t, :] - x_t
            J += np.real(np.vdot(r_t, r_t))

            w_t_H = ((r_t.conj() @ X) * diag_inv) @ Y
            g_f -= w_t_H
            G_A += np.outer(x_t, w_t_H)

        # Forcing gradient (Eq. 39)
        g_f_param = g_f @ self.E
        grad_f_r  =  2 * np.real(g_f_param)
        grad_f_i  = -2 * np.imag(g_f_param)

        # Mistuning gradient (Eq. 40)
        grad_d_0 = 2 * np.real(np.diag(G_A))

        # Circulant gradients (Eqs. 42-43).
        # scipy.linalg.circulant uses the *first column* convention:
        #   C[r, c] = v[(r - c) mod N]
        # The diagonal with offset -j in G_A accumulates terms where r - c = j,
        # matching circ index j.  Due to this column-first convention, the raw
        # gradient vectors come out index-reversed and must be flipped.
        grad_a_0_raw = np.zeros(self.N - 1)
        grad_a_1_raw = np.zeros(self.N)

        for j in range(self.N):
            diag_sum = np.trace(G_A, offset=-j) + np.trace(G_A, offset=self.N - j)
            grad_a_1_raw[j] = -2 * np.imag(diag_sum)
            if j > 0:
                grad_a_0_raw[j - 1] = 2 * np.real(diag_sum)

        # Correct index reversal caused by the first-column circulant convention
        grad_a_0 = grad_a_0_raw[::-1].copy()

        grad_a_1 = np.empty(self.N)
        grad_a_1[0] = grad_a_1_raw[0]
        for j in range(1, self.N):
            grad_a_1[j] = grad_a_1_raw[self.N - j]

        gradient = np.concatenate([grad_d_0, grad_a_0, grad_a_1, grad_f_r, grad_f_i])
        scale = 1.0 / (np.sum(np.abs(self.y) ** 2) + 1e-10)

        return J * scale, gradient * scale

    def objective(self, theta):
        """Return the objective value only (convenience wrapper)."""
        return self.compute_objective_and_gradient(theta)[0]

    def gradient(self, theta):
        """Return the gradient only (convenience wrapper)."""
        return self.compute_objective_and_gradient(theta)[1]

    # ------------------------------------------------------------------
    # Cramér-Rao bounds
    # ------------------------------------------------------------------

    def compute_cramer_rao_bounds(self, theta, noise_variance):
        """
        Compute the Cramér-Rao lower bounds for the parameter estimates.

        The Fisher information matrix is:

            I(theta) = (1/sigma^2) * sum_t (dmu_t/dtheta)^H (dmu_t/dtheta)

        and CRB = sqrt(diag(I^{-1})) (Eqs. 21-22).

        Parameters
        ----------
        theta          : (5N-1,) float  — current parameter estimate
        noise_variance : float          — measurement noise variance sigma^2

        Returns
        -------
        crb : (5N-1,) float  — standard-deviation lower bounds per parameter
        """
        d_0, a_0, a_1, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0, a_1, f_r, f_i)
        N_params = len(theta)
        Fisher = np.zeros((N_params, N_params))

        # Jacobian blocks: dA/dtheta_i and df/dtheta_i
        dA = np.zeros((N_params, self.N, self.N), dtype=np.complex128)
        df = np.zeros((N_params, self.N), dtype=np.complex128)
        idx = 0

        for j in range(self.N):              # d0 block
            dA[idx, j, j] = 1.0; idx += 1

        for j in range(1, self.N):           # a0 block (N-1 params)
            e_j = np.zeros(self.N); e_j[j] = 1.0
            dA[idx] = scipy.linalg.circulant(e_j); idx += 1

        for j in range(self.N):              # a1 block
            e_j = np.zeros(self.N); e_j[j] = 1.0
            dA[idx] = 1j * scipy.linalg.circulant(e_j); idx += 1

        for j in range(self.N):              # fr block
            df[idx] = self.E[:, j]; idx += 1

        for j in range(self.N):              # fi block
            df[idx] = 1j * self.E[:, j]; idx += 1

        for t in range(self.T):
            H_t  = A - ((self.omega_t[t] / self.omega_0) ** 2) * self.I_N
            mu_t = np.linalg.solve(H_t, f)
            # Sensitivity: dmu_t/dtheta_i = H_t^{-1} (df_i - dA_i @ mu_t)
            J_mu_t = np.linalg.solve(H_t, (df - (dA @ mu_t)).T)
            Fisher += np.real(J_mu_t.conj().T @ J_mu_t)

        return np.sqrt(np.diag(noise_variance * np.linalg.pinv(Fisher, rcond=1e-6)))


# ---------------------------------------------------------------------------
# Cached evaluator
# ---------------------------------------------------------------------------

class CachedEvaluator:
    """
    Caching wrapper for AeromechanicalIdentification.

    Many optimisers (e.g. L-BFGS-B) query the objective and the gradient
    separately with the same theta.  This wrapper stores the last evaluation
    and returns the cached result on repeated calls, avoiding redundant
    full-system solves.

    Parameters
    ----------
    ident_problem : AeromechanicalIdentification
    """

    def __init__(self, ident_problem):
        self.problem    = ident_problem
        self.last_theta = None
        self.last_J     = None
        self.last_grad  = None

    def _update(self, theta):
        if self.last_theta is None or not np.array_equal(theta, self.last_theta):
            self.last_J, self.last_grad = self.problem.compute_objective_and_gradient(theta)
            self.last_theta = np.copy(theta)

    def objective(self, theta):
        """Return the cached objective, recomputing only if theta changed."""
        self._update(theta)
        return self.last_J

    def gradient(self, theta):
        """Return the cached gradient, recomputing only if theta changed."""
        self._update(theta)
        return self.last_grad
