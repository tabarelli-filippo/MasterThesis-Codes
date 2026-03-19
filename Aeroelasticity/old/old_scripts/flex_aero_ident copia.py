import numpy as np
import scipy.linalg


class FlexibleAeromechanicalIdentification:
    """
    Identificazione aeromeccanica massima verosimiglianza per rotore con disco flessibile.
    
    Modello (Eq. 13 estesa con disco flessibile, Sez. 4.3):
    
        [A(theta) - gamma_t^2 * I - D_disk(omega_t)/omega_0^2] x_t = E @ f_q
    
    dove:
        A(theta) = I + diag(d0) + circ(a0) + 1j*circ(a1)
        gamma_t  = omega_t / omega_0
        D_disk   = K_bd @ inv(K_dd - omega_t^2 * M_dd) @ K_db   (pre-calcolato)
    
    Layout di theta (5N-1 parametri reali):
        [d0(0..N-1) | a0_params(0..N-2) | a1_params(0..N-1) | fr(0..N-1) | fi(0..N-1)]
    
    dove:
        a0_params[k] = a0_circ[k+1]   per k=0..N-2  (distanze 1..N-1)
        a1_params[0] = a1_circ[0]     (auto-smorzamento)
        a1_params[k] = a1_circ[k]     per k=1..N-1
    
    NOTA su scipy.linalg.circulant(c):
        C[i,j] = c[(i-j) mod N]  — usa c come PRIMA COLONNA.
        Quindi C[0,:] = [c[0], c[N-1], c[N-2], ..., c[1]]  (prima riga != prima colonna).
    """

    def __init__(self, y_measured, omega_t, omega_0, N, M_dd, K_dd, K_bd):
        self.y = np.asarray(y_measured, dtype=np.complex128)
        self.omega_t = np.asarray(omega_t, dtype=np.float64)
        self.omega_0 = float(omega_0)
        self.N = N
        self.T = len(omega_t)
        self.I_N = np.eye(N, dtype=np.complex128)

        # Matrice DFT (Eq. 11)
        j, k = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        self.E = (1.0 / np.sqrt(N)) * np.exp(1j * 2 * np.pi * j * k / N)

        # Indici di blocco nel vettore theta (5N-1) — accessibili dagli script esterni
        self.idx_d0 = 0
        self.idx_a0 = N
        self.idx_a1 = 2 * N - 1
        self.idx_fr = 3 * N - 1
        self.idx_fi = 4 * N - 1

        # Pre-calcolo dell'influenza del disco flessibile normalizzata: D_disk/omega_0^2
        # Shape: (T, N, N)
        K_db = K_bd.conj().T
        self.D_disk = np.zeros((self.T, N, N), dtype=np.complex128)
        for t in range(self.T):
            Z_dd = K_dd - (self.omega_t[t]**2) * M_dd
            self.D_disk[t] = K_bd @ np.linalg.solve(Z_dd, K_db) / (self.omega_0**2)

    # ------------------------------------------------------------------
    # Unpacking / packing di theta
    # ------------------------------------------------------------------

    def unpack_theta(self, theta):
        """
        Estrae i blocchi di parametri da theta e ricostruisce i vettori
        della prima colonna delle matrici circolanti.

        Returns:
            d_0       : (N,)   mistuning diagonale
            a_0_circ  : (N,)   prima colonna di circ(a0); a_0_circ[0]=0 per definizione
            a_1_circ  : (N,)   prima colonna di circ(a1); a_1_circ[0]=auto-smorzamento
            f_r       : (N,)   parte reale del vettore forzante in TW coordinates
            f_i       : (N,)   parte immaginaria
        """
        N = self.N
        d_0 = theta[0:N]

        # a0: N-1 parametri = prima colonna circ(a0) alle distanze 1..N-1
        a_0_circ = np.zeros(N)
        a_0_circ[1:] = theta[N: 2*N - 1]          # a0_circ[k] = theta[N + k - 1], k=1..N-1

        # a1: N parametri; a1_circ[0]=auto-smorzamento, a1_circ[1..N-1]=distanze 1..N-1
        a_1_circ = theta[2*N - 1: 3*N - 1].copy()  # copia diretta: indice 0 = auto-smorzamento

        f_r = theta[3*N - 1: 4*N - 1]
        f_i = theta[4*N - 1: 5*N - 1]
        return d_0, a_0_circ, a_1_circ, f_r, f_i

    def pack_theta(self, d_0, a_0_circ, a_1_circ, f_r, f_i):
        """Operazione inversa di unpack_theta."""
        return np.concatenate([d_0, a_0_circ[1:], a_1_circ, f_r, f_i])

    # ------------------------------------------------------------------
    # Costruzione matrici di sistema
    # ------------------------------------------------------------------

    def build_system_matrices(self, d_0, a_0_circ, a_1_circ, f_r, f_i):
        """
        Costruisce A(theta) e il vettore forzante fisico f.

        A = I + diag(d0) + circ(a0) + 1j*circ(a1)
        f = E @ (f_r + 1j*f_i)

        scipy.linalg.circulant(c): C[i,j] = c[(i-j) mod N]
        """
        A = (self.I_N
             + np.diag(d_0)
             + scipy.linalg.circulant(a_0_circ)
             + 1j * scipy.linalg.circulant(a_1_circ))
        f = self.E @ (f_r + 1j * f_i)
        return A, f

    def solve_eigendecomposition(self, A):
        """
        Decompone A = X @ diag(lam) @ Y  con Y @ X = I (bi-ortogonalità).
        Utile per il calcolo efficiente O(N^2) invece di O(N^3) per T risoluzioni.

        Returns:
            lam : (N,)  autovalori
            X   : (N,N) matrice degli autovettori destri (colonne)
            Y   : (N,N) matrice degli autovettori sinistri (righe), Y @ X = I
        """
        lam, X = np.linalg.eig(A)
        Y = np.linalg.inv(X)
        return lam, X, Y

    # ------------------------------------------------------------------
    # Funzione obiettivo e gradiente (metodo adjoint)
    # ------------------------------------------------------------------

    def compute_objective_and_gradient(self, theta):
        """
        Calcola J(theta) = sum_t ||y_t - x_t||^2 / ||y||^2   (normalizzato)
        e il suo gradiente rispetto a theta usando il metodo adjoint.

        Derivazione del gradiente (da Eq. 33-43 dell'articolo):

            delta_J = 2 Re{ sum_t  x_t^H delta_A lambda_t
                                 + delta_f^H lambda_t }

        dove lambda_t soddisfa il sistema aggiunto:
            Z_t^H lambda_t = r_t      con r_t = y_t - x_t

        Definendo:
            G_A = sum_t  outer(x_t, lambda_t^*)   [N x N]

        i gradienti sono:
            dJ/d(d0)_j    = 2 Re{ G_A[j,j] }                         (Eq. 40)
            dJ/d(a0_circ)_k = 2 Re{ sum_{p-q=k mod N} G_A[p,q] }    (Eq. 43)
            dJ/d(a1_circ)_k = -2 Im{ sum_{p-q=k mod N} G_A[p,q] }   (Eq. 42)
            dJ/d(fr)      = 2 Re{ E^H @ sum_t lambda_t }^*
            dJ/d(fi)      = -2 Im{ E^H @ sum_t lambda_t }^*
        """
        d_0, a_0_circ, a_1_circ, f_r, f_i = self.unpack_theta(theta)
        A_blade, f = self.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)

        N = self.N
        J = 0.0
        lambda_sum = np.zeros(N, dtype=np.complex128)
        G_A = np.zeros((N, N), dtype=np.complex128)

        # Normalizzazione per condizionamento numerico
        denom = np.sum(np.abs(self.y)**2) + 1e-30

        for t in range(self.T):
            gamma_sq = (self.omega_t[t] / self.omega_0)**2

            # Matrice di sistema completa (con disco flessibile già normalizzato)
            Z_t = A_blade - gamma_sq * self.I_N - self.D_disk[t]

            # Soluzione diretta
            x_t = np.linalg.solve(Z_t, f)
            r_t = self.y[t, :] - x_t
            J += np.real(np.vdot(r_t, r_t))

            # Soluzione aggiunta: Z_t^H @ lambda_t = r_t
            lambda_t = np.linalg.solve(Z_t.conj().T, r_t)

            lambda_sum += lambda_t
            # outer(x_t, lambda_t^*) = x_t @ lambda_t^H
            G_A += np.outer(x_t, lambda_t.conj())

        # --- Gradiente della forzante ---
        # delta_J/delta_f = -2 Re{ sum_t lambda_t^H delta_f }
        # f = E @ (f_r + 1j*f_i)  =>  delta_f = E @ delta_f_q
        # => delta_J/delta_f_r = 2 Re{ E^H @ lambda_sum }   (segno: vedi Eq. 39)
        # => delta_J/delta_f_i = -2 Im{ E^H @ lambda_sum }
        # Con il segno corretto dell'adjoint (r_t = y - x, non x - y):
        #   grad_fr = -2 Re(E^H lambda_sum)  se definiamo S = lambda_sum^H E
        # Verifichiamo con la derivazione diretta:
        #   J = sum ||r_t||^2,  r_t = y_t - Z_t^{-1} E f_q
        #   dJ/df_r = -2 Re{ sum_t (Z_t^{-H} r_t)^H E }
        #           = -2 Re{ lambda_sum^H E }
        #           = -2 Re{ (E^H lambda_sum)^H }
        S = self.E.conj().T @ lambda_sum          # S = E^H @ lambda_sum, shape (N,)
        grad_f_r = -2.0 * np.real(S)
        grad_f_i = -2.0 * np.imag(S)

        # --- Gradiente mistuning (Eq. 40) ---
        # delta_A = diag(delta_d0)  =>  Tr(G_A^T delta_A) = sum_j G_A[j,j] delta_d0[j]
        grad_d_0 = 2.0 * np.real(np.diag(G_A))

        # --- Gradiente circolanti (Eq. 42-43) ---
        # Per circ(a)  con a = prima colonna: C[p,q] = a[(p-q) mod N]
        # delta_A = circ(delta_a)  =>  Tr(G_A^T delta_A) = sum_{p,q} G_A[q,p] delta_a[(p-q) mod N]
        #         = sum_k delta_a[k] * sum_{p-q=k mod N} G_A[q,p]
        # Quindi:
        #   dJ/d(a0_circ)[k] = 2 Re{ sum_{p-q=k mod N} G_A[q,p] }
        #   dJ/d(a1_circ)[k] = -2 Im{ sum_{p-q=k mod N} G_A[q,p] }
        #
        # sum_{p-q=k mod N} G_A[q,p] = sum diagonali di G_A con offset k
        # np.trace(M, offset=k) calcola sum_i M[i, i+k]  = sum_{q, p=q+k} G_A[q,p]
        # => p - q = k  ✓ per offset=k nel range 0..N-1
        # Per k=0: traccia principale (offset=0)
        # Per k=N-j (j>0): corrisponde a offset -(j) mod N → usiamo offset negativo
        #   np.trace(G_A, offset=-(N-k)) = np.trace(G_A, offset=k-N) per k=1..N-1

        grad_a_0_circ = np.zeros(N)
        grad_a_1_circ = np.zeros(N)
        for k in range(N):
            if k == 0:
                diag_sum = np.trace(G_A)
            else:
                # Diagonale con p - q = k:
                # - parte "sopra": np.trace(G_A, offset=k)  [p=q+k, q=0..N-1-k]
                # - parte "sotto" (wrap): np.trace(G_A, offset=k-N) [p=q+k-N, q=N-k..N-1]
                diag_sum = np.trace(G_A, offset=k) + np.trace(G_A, offset=k - N)
            grad_a_0_circ[k] =  2.0 * np.real(diag_sum)
            grad_a_1_circ[k] = -2.0 * np.imag(diag_sum)

        # Ricostruzione del vettore gradiente nel layout di theta
        # a0_params = a0_circ[1:]  => grad_a0_params = grad_a_0_circ[1:]
        # a1_params = a1_circ      => grad_a1_params = grad_a_1_circ
        grad_a_0_params = grad_a_0_circ[1:]    # N-1 elementi
        grad_a_1_params = grad_a_1_circ        # N elementi

        gradient = np.concatenate([grad_d_0, grad_a_0_params, grad_a_1_params, grad_f_r, grad_f_i])

        return J / denom, gradient / denom

    # ------------------------------------------------------------------
    # Utilità
    # ------------------------------------------------------------------

    def solve_frequency_response(self, theta, t):
        """Calcola la risposta x_t per un singolo indice di frequenza."""
        d_0, a_0_circ, a_1_circ, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)
        gamma_sq = (self.omega_t[t] / self.omega_0)**2
        Z_t = A - gamma_sq * self.I_N - self.D_disk[t]
        return np.linalg.solve(Z_t, f)

    def compute_cramer_rao(self, theta, sigma_noise):
        """
        Calcola la matrice di Fisher e i CRB (Eq. 21-22).
        I(theta) = (1/sigma^2) * sum_t (dmu_t/dtheta)^H (dmu_t/dtheta)
        CRB = sqrt(diag(I^{-1}))
        """
        eps = 1e-7
        num_params = len(theta)
        J0, _ = self.compute_objective_and_gradient(theta)
        denom = np.sum(np.abs(self.y)**2) + 1e-30

        # Calcola mu_t(theta) per tutti t
        d_0, a_0_circ, a_1_circ, f_r, f_i = self.unpack_theta(theta)
        A, f = self.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)

        mu = np.zeros((self.T, self.N), dtype=np.complex128)
        for t in range(self.T):
            gamma_sq = (self.omega_t[t] / self.omega_0)**2
            Z_t = A - gamma_sq * self.I_N - self.D_disk[t]
            mu[t] = np.linalg.solve(Z_t, f)

        # Differenze finite per dmu/dtheta
        Fisher = np.zeros((num_params, num_params))
        dmu = np.zeros((num_params, self.T, self.N), dtype=np.complex128)

        for i in range(num_params):
            theta_p = theta.copy(); theta_p[i] += eps
            d_0p, a_0p, a_1p, frp, fip = self.unpack_theta(theta_p)
            Ap, fp = self.build_system_matrices(d_0p, a_0p, a_1p, frp, fip)
            for t in range(self.T):
                gamma_sq = (self.omega_t[t] / self.omega_0)**2
                Z_t = Ap - gamma_sq * self.I_N - self.D_disk[t]
                mu_p = np.linalg.solve(Z_t, fp)
                dmu[i, t] = (mu_p - mu[t]) / eps

        for i in range(num_params):
            for j in range(i, num_params):
                val = np.sum(np.real(np.einsum('ti,ti->', dmu[i].conj(), dmu[j])))
                Fisher[i, j] = val / sigma_noise**2
                Fisher[j, i] = Fisher[i, j]

        try:
            Fisher_inv = np.linalg.inv(Fisher + 1e-20 * np.eye(num_params))
            crb = np.sqrt(np.maximum(np.diag(Fisher_inv), 0.0))
        except np.linalg.LinAlgError:
            crb = np.full(num_params, np.nan)

        return crb