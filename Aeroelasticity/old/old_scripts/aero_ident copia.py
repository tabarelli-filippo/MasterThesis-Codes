import numpy as np
import scipy.linalg

class AeromechanicalIdentification:
    def __init__(self, y_measured, omega_t, omega_0, N):
        self.y = np.asarray(y_measured, dtype=np.complex128)
        self.omega_t = np.asarray(omega_t, dtype=np.float64)
        self.omega_0 = omega_0
        self.N = N
        self.T = len(omega_t)

        # Operatore IDFT precalcolato (matrice E del paper, eq. 11)
        j_idx, k_idx = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        self.E = (1.0 / np.sqrt(N)) * np.exp(1j * 2 * np.pi * j_idx * k_idx / N)
        self.I_N = np.eye(N, dtype=np.complex128)

    def unpack_theta(self, theta):
        N = self.N
        idx = 0
        d_0 = theta[idx : idx+N]; idx += N
        a_0 = np.zeros(N)              # Vincolo a0[0] = 0
        a_0[1:] = theta[idx : idx+(N-1)]; idx += (N-1)
        a_1 = theta[idx : idx+N]; idx += N
        f_r = theta[idx : idx+N]; idx += N
        f_i = theta[idx : idx+N]
        return d_0, a_0, a_1, f_r, f_i

    def build_system_matrices(self, d_0, a_0, a_1, f_r, f_i):
        D_0 = np.diag(d_0)
        A_0_mat = scipy.linalg.circulant(a_0)
        A_1_mat = scipy.linalg.circulant(a_1)
        A = self.I_N + D_0 + A_0_mat + 1j * A_1_mat
        f = self.E @ (f_r + 1j * f_i)
        return A, f

    def solve_eigendecomposition(self, A):
        M = self.I_N
        lam, vl, X = scipy.linalg.eig(A, b=M, left=True, right=True)
        Y_raw = vl.conj().T
        diag_YMX = np.diag(Y_raw @ M @ X)
        Y = Y_raw / diag_YMX[:, np.newaxis]
        return lam, X, Y

    def compute_objective_and_gradient(self, theta):
        d_0, a_0, a_1, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0, a_1, f_r, f_i)

        lam, X, Y = self.solve_eigendecomposition(A)
        Y_f = Y @ f

        J = 0.0
        g_f  = np.zeros(self.N, dtype=np.complex128)
        G_A  = np.zeros((self.N, self.N), dtype=np.complex128)

        for t in range(self.T):
            gamma_sq = (self.omega_t[t] / self.omega_0)**2
            diag_inv = 1.0 / (lam - gamma_sq)

            x_t   = X @ (diag_inv * Y_f)
            r_t   = self.y[t, :] - x_t
            J    += np.real(np.vdot(r_t, r_t))

            w_t_H = ((r_t.conj() @ X) * diag_inv) @ Y
            g_f  -= w_t_H
            G_A  += np.outer(x_t, w_t_H)

        # ---- Gradiente forzante (eq. 39 del paper) ----
        g_f_param  = g_f @ self.E
        grad_f_r   = 2  * np.real(g_f_param)
        grad_f_i   = -2 * np.imag(g_f_param)

        # ---- Gradiente mistuning diagonale (eq. 40) ----
        grad_d_0 = 2 * np.real(np.diag(G_A))

        # ---- Gradiente circolanti: calcolo grezzo poi corretto ----
        #
        # Formula per il circolante scipy (prima colonna, C[r,c] = v[(r-c) mod N]):
        #   dJ/dv[j]_grezzo -> somma elementi con r-c = j mod N
        #   = trace(G_A, offset=-j) + trace(G_A, offset=N-j)
        #
        # A causa della convenzione scipy (prima colonna anziché prima riga),
        # l'assegnazione a grad_a_0[j-1] e grad_a_1[j] risulta INVERTITA
        # rispetto all'indicizzazione di theta.
        # Correzione: inversione dei vettori grezza.

        grad_a_0_raw = np.zeros(self.N - 1)
        grad_a_1_raw = np.zeros(self.N)

        for j in range(self.N):
            diag_sum = np.trace(G_A, offset=-j) + np.trace(G_A, offset=self.N - j)
            grad_a_1_raw[j] = -2 * np.imag(diag_sum)
            if j > 0:
                grad_a_0_raw[j - 1] = 2 * np.real(diag_sum)

        # FIX: correggi l'inversione di indice
        # a_0[j] corrisponde a theta[N + j - 1], ma grezzo[j-1] = dJ/da_0[N-j]
        # => inversione del vettore
        grad_a_0 = grad_a_0_raw[::-1].copy()

        # a_1: j=0 (diagonale) si mappa su se stesso; j=1..N-1 sono invertiti
        grad_a_1 = np.empty(self.N)
        grad_a_1[0] = grad_a_1_raw[0]
        for j in range(1, self.N):
            grad_a_1[j] = grad_a_1_raw[self.N - j]

        gradient = np.concatenate([grad_d_0, grad_a_0, grad_a_1, grad_f_r, grad_f_i])
        scale_factor = 1.0 / (np.sum(np.abs(self.y)**2) + 1e-10)

        return J * scale_factor, gradient * scale_factor

    def objective(self, theta):
        return self.compute_objective_and_gradient(theta)[0]

    def gradient(self, theta):
        return self.compute_objective_and_gradient(theta)[1]

    def compute_cramer_rao_bounds(self, theta, noise_variance):
        d_0, a_0, a_1, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0, a_1, f_r, f_i)
        N_params = len(theta)
        Fisher = np.zeros((N_params, N_params))
        dA = np.zeros((N_params, self.N, self.N), dtype=np.complex128)
        df = np.zeros((N_params, self.N), dtype=np.complex128)
        idx = 0
        for j in range(self.N):
            dA[idx, j, j] = 1.0; idx += 1
        for j in range(1, self.N):
            e_j = np.zeros(self.N); e_j[j] = 1.0
            dA[idx, :, :] = scipy.linalg.circulant(e_j); idx += 1
        for j in range(self.N):
            e_j = np.zeros(self.N); e_j[j] = 1.0
            dA[idx, :, :] = 1j * scipy.linalg.circulant(e_j); idx += 1
        for j in range(self.N):
            df[idx, :] = self.E[:, j]; idx += 1
        for j in range(self.N):
            df[idx, :] = 1j * self.E[:, j]; idx += 1
        for t in range(self.T):
            H_t   = A - ((self.omega_t[t]/self.omega_0)**2) * self.I_N
            mu_t  = np.linalg.solve(H_t, f)
            J_mu_t = np.linalg.solve(H_t, (df - (dA @ mu_t)).T)
            Fisher += np.real(J_mu_t.conj().T @ J_mu_t)
        return np.sqrt(np.diag(noise_variance * np.linalg.pinv(Fisher, rcond=1e-6)))


class CachedEvaluator:
    """Wrapper con cache per evitare doppi calcoli di J e del gradiente."""

    def __init__(self, ident_problem):
        self.problem   = ident_problem
        self.last_theta = None
        self.last_J     = None
        self.last_grad  = None

    def _update(self, theta):
        if self.last_theta is None or not np.array_equal(theta, self.last_theta):
            self.last_J, self.last_grad = self.problem.compute_objective_and_gradient(theta)
            self.last_theta = np.copy(theta)

    def objective(self, theta):
        self._update(theta)
        return self.last_J

    def gradient(self, theta):
        self._update(theta)
        return self.last_grad