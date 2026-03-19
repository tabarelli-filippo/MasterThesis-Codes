"""
data_generation.py
==================
Modulo unificato per la generazione di dati simulati di vibrazione delle pale.
Supporta due modalità tramite il parametro `mode`:

    mode='rigid'    — disco rigido (data_generation_rigid.py)
    mode='flexible' — disco flessibile con K_bd densa (data_generation_flexible.py)

Interfaccia unificata
---------------------
    result = generate_rotor_data(mode, omega_t, omega_0, noise_level, N=33, seed=42)

Valori di ritorno (dict):
    'y_measured'  : array (T, N) complesso — spostamenti misurati (con rumore)
    'theta_true'  : array (5N-1,) reale    — parametri veri nel layout dell'ID rigido
    'M_dd'        : array (N, N) o None    — masse modali del disco (solo flexible)
    'K_dd'        : array (N, N) o None    — rigidezze modali del disco (solo flexible)
    'K_bd'        : array (N, N) o None    — coupling pala-disco (solo flexible)

Le chiavi 'M_dd', 'K_dd', 'K_bd' sono sempre presenti nel dict; per mode='rigid'
valgono None, così il chiamante può fare unpacking generico senza if/else.

Compatibilità con i test esistenti
------------------------------------
    # test_flex_disk_banded.py  / text_flex_with_rigidity.py
    result = generate_rotor_data('flexible', omega_t, omega_0, noise_level)
    y_measured   = result['y_measured']
    theta_true   = result['theta_true']
    M_dd, K_dd, K_bd = result['M_dd'], result['K_dd'], result['K_bd']

    # test rigido
    result = generate_rotor_data('rigid', omega_t, omega_0, noise_level, N=33)
    y_measured = result['y_measured']
    theta_true = result['theta_true']

Fisica del modello flessibile
------------------------------
    [A_blade - γ²·I - D_disk(ω)/ω₀²] x = E·f_q

con γ = ω/ω₀,  D_disk = K_bd · inv(K_dd − ω²·M_dd) · K_db.
K_bd è DENSA in spazio TW (decadimento esponenziale in ΔND), il che rende
D_disk NON-circolante → l'artefatto matematico di Hall & Hall (2024) è visibile
quando si usa il modello rigido per identificare questi dati.

Fisica del modello rigido
--------------------------
    [A_blade − γ²·I] x = f

con A_blade costruito da (d0, a0_circ, a1_circ) secondo la normalizzazione di
AeromechanicalIdentification.  K_bd = 0 (nessuna condensazione di Schur).
"""

import numpy as np
import scipy.linalg


# ---------------------------------------------------------------------------
# Parametri fisici condivisi (Tabella 1 del paper)
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
    Restituisce (a0_circ, a1_circ) — prime colonne dei circolanti aerodinamici
    per N=33 (Tabella 1).  Per N≠33 i coefficienti non-zero sono mappati
    in modo consistente con il wrapping circolare.
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
# Funzione pubblica unificata
# ---------------------------------------------------------------------------

def generate_rotor_data(mode, omega_t, omega_0, noise_level=1e-4, N=33, seed=42):
    """
    Genera dati simulati di vibrazione delle pale.

    Parametri
    ----------
    mode        : str   — 'rigid' oppure 'flexible'
    omega_t     : array (T,)  — frequenze di eccitazione [rad/s]
    omega_0     : float        — frequenza nominale delle pale [rad/s]
    noise_level : float        — deviazione standard del rumore
    N           : int          — numero di pale (usato solo per mode='rigid';
                                 per 'flexible' è fissato a 33)
    seed        : int          — seme per il generatore di numeri casuali

    Restituisce
    -----------
    dict con le chiavi:
        'y_measured' : array (T, N) complesso
        'theta_true' : array (5N-1,) reale
        'M_dd'       : array (N, N) o None
        'K_dd'       : array (N, N) o None
        'K_bd'       : array (N, N) o None
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
        raise ValueError(f"mode deve essere 'rigid' o 'flexible', ricevuto: {mode!r}")


# ---------------------------------------------------------------------------
# Implementazione interna — disco RIGIDO
# ---------------------------------------------------------------------------

def _generate_rigid(N, omega_0, omega_t, noise_level, seed):
    """
    Modello a disco rigido.
    [A_blade − γ²·I] x = f

    Replica esatta di data_generation_rigid.initialize_rotor(), rifattorizzata
    per usare i parametri condivisi e accettare un seed esplicito.
    """
    from auxSolver.aero_ident import AeromechanicalIdentification  # import locale per evitare dipendenze circolari

    np.random.seed(seed)
    T = len(omega_t)
    num_params = 5 * N - 1

    # Indici nel vettore theta
    idx_d0 = 0
    idx_a0 = N
    idx_a1 = N + (N - 1)
    idx_fr = N + (N - 1) + N
    idx_fi = N + (N - 1) + N + N

    theta_true = np.zeros(num_params, dtype=np.float64)

    # Mistuning
    d0 = _D0_ARRAY if N == len(_D0_ARRAY) else np.resize(_D0_ARRAY, N)
    theta_true[idx_d0 : idx_d0 + N] = d0

    # Coefficienti aerodinamici
    a0_circ, a1_circ = _build_aero_circulant(N)
    theta_true[idx_a0 : idx_a0 + N - 1] = a0_circ[1:]   # a0_circ[0] = 0 per definizione
    theta_true[idx_a1 : idx_a1 + N]     = a1_circ

    # Forzante (ND=28)
    target_nd = 28
    if target_nd < N:
        theta_true[idx_fr + target_nd] = -0.1
        theta_true[idx_fi + target_nd] =  0.0

    # Risposta esatta tramite AeromechanicalIdentification
    y_dummy     = np.zeros((T, N), dtype=np.complex128)
    sim_problem = AeromechanicalIdentification(y_dummy, omega_t, omega_0, N)

    d_0, a_0, a_1, f_r, f_i = sim_problem.unpack_theta(theta_true)
    A, f = sim_problem.build_system_matrices(d_0, a_0, a_1, f_r, f_i)

    mu_exact = np.zeros((T, N), dtype=np.complex128)
    for t in range(T):
        gamma_sq     = (omega_t[t] / omega_0) ** 2
        H_t          = A - gamma_sq * sim_problem.I_N
        mu_exact[t]  = np.linalg.solve(H_t, f)

    rng = np.random.RandomState(seed)
    y_measured = (mu_exact
                  + rng.normal(0, noise_level, (T, N))
                  + 1j * rng.normal(0, noise_level, (T, N)))

    return y_measured, theta_true


# ---------------------------------------------------------------------------
# Implementazione interna — disco FLESSIBILE
# ---------------------------------------------------------------------------

def _generate_flexible(omega_t, omega_0, noise_level, seed):
    """
    Modello a disco flessibile con K_bd densa in spazio TW.
    [A_blade − γ²·I − D_disk(ω)/ω₀²] x = E·f_q

    Replica esatta di data_generation_flexible.initialize_flexible_rotor_complete(),
    rifattorizzata per usare i parametri condivisi e accettare un seed esplicito.
    """
    N = 33   # fissato dal modello fisico (Hall & Hall 2024)
    T = len(omega_t)

    a0_circ, a1_circ = _build_aero_circulant(N)
    d0_true          = _D0_ARRAY.copy()

    # Forzante modale
    f_q = np.zeros(N, dtype=np.complex128)
    f_q[28] = -0.1 + 0.0j

    # ------------------------------------------------------------------
    # Struttura modale del disco
    # ------------------------------------------------------------------
    M_dd = np.eye(N, dtype=np.complex128) * 50.0

    nd_array   = np.fft.fftfreq(N, d=1.0 / N)
    disk_freqs = omega_0 * (0.35 + 0.13 * np.abs(nd_array))
    K_dd       = np.diag(50.0 * disk_freqs**2 * (1.0 + 1j * 1e-4)).astype(np.complex128)

    # ------------------------------------------------------------------
    # Accoppiamento pala-disco (K_bd DENSA — decadimento esponenziale in ΔND)
    # ------------------------------------------------------------------
    coupling_amplitude = 0.02 * omega_0 ** 2
    coupling_length    = 3.0

    K_bd = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            delta_nd   = min(abs(i - j), N - abs(i - j))
            K_bd[i, j] = coupling_amplitude * np.exp(-delta_nd / coupling_length)

    # Troncamento ND 12..21
    for idx in range(12, 22):
        K_bd[idx, :] = 0.0
        K_bd[:, idx] = 0.0

    K_db = K_bd.conj().T

    # ------------------------------------------------------------------
    # Matrici del sistema
    # ------------------------------------------------------------------
    j_idx, k_idx = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    E = (1.0 / np.sqrt(N)) * np.exp(1j * 2 * np.pi * j_idx * k_idx / N)

    A_blade = (np.eye(N, dtype=np.complex128)
               + np.diag(d0_true)
               + scipy.linalg.circulant(a0_circ)
               + 1j * scipy.linalg.circulant(a1_circ))

    f_phys = E @ f_q

    # ------------------------------------------------------------------
    # Soluzione armonica con condensazione di Schur
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
    # theta_true nel layout del modello identificativo RIGIDO
    # ------------------------------------------------------------------
    theta_true = np.concatenate([
        d0_true,
        a0_circ[1:],      # N-1 elementi
        a1_circ,           # N   elementi
        np.real(f_q),
        np.imag(f_q),
    ])

    return y_measured, theta_true, M_dd, K_dd, K_bd


# ---------------------------------------------------------------------------
# Alias di compatibilità backward (per non rompere import esistenti)
# ---------------------------------------------------------------------------

def initialize_rotor(N, omega_0, omega_t, noise_level, seed=42):
    """Alias backward-compatible con data_generation_rigid.initialize_rotor."""
    res = generate_rotor_data('rigid', omega_t, omega_0, noise_level, N=N, seed=seed)
    return res['y_measured'], res['theta_true']


def initialize_flexible_rotor_complete(omega_t, omega_0, noise_level=1e-4, seed=42):
    """Alias backward-compatible con data_generation_flexible.initialize_flexible_rotor_complete."""
    res = generate_rotor_data('flexible', omega_t, omega_0, noise_level, seed=seed)
    return res['y_measured'], res['theta_true'], res['M_dd'], res['K_dd'], res['K_bd']


def initialize_flexible_rotor(omega_t, omega_0, noise_level, seed=42,
                               M_dd=None, K_dd=None, M_bd=None, K_bd=None):
    """Alias backward-compatible con data_generation_flexible.initialize_flexible_rotor (legacy)."""
    res = generate_rotor_data('flexible', omega_t, omega_0, noise_level, seed=seed)
    return res['y_measured'], res['theta_true']
