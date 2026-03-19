import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.flex_aero_ident import FlexibleAeromechanicalIdentification
from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter


def apply_modal_truncation(M_dd, K_dd, K_bd, N, max_nd=10):
    """
    Annulla l'accoppiamento tra pala e disco per diametri nodali superiori a max_nd.
    Questo simula la flessibilità del disco non modellata senza rendere singolare
    la matrice di impedenza del disco isolato.
    """
    M_dd_trunc = np.copy(M_dd)
    K_dd_trunc = np.copy(K_dd)
    K_bd_trunc = np.copy(K_bd)

    for i in range(N):
        # Calcolo del diametro nodale effettivo
        nd = i if i <= N // 2 else N - i
        
        if nd > max_nd:
            # Azzera SOLO l'accoppiamento. 
            # M_dd e K_dd restano intatti per garantire l'invertibilità numerica,
            # ma l'influenza dinamica del disco sulle pale diventa nulla.
            K_bd_trunc[i] = 0.0
            
    return M_dd_trunc, K_dd_trunc, K_bd_trunc


class CachedEvaluator:
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

    def objective(self, x): self._update(x); return self._J
    def gradient(self,  x): self._update(x); return self._grad


def optimize(evaluator, x0_scaled, bounds, label=''):
    print(f"\n  Avvio {label} con IPOPT...")
    ipopt_options = {
        'max_iter': 5000,
        'tol': 1e-12,
        'acceptable_tol': 1e-10,
        'print_level': 5,
        'hessian_approximation': 'limited-memory',
        'limited_memory_max_history': 30,
        'limited_memory_update_type': 'bfgs',
    }
    res = minimize_ipopt(
        fun=evaluator.objective, 
        x0=x0_scaled, 
        jac=evaluator.gradient,
        bounds=bounds, 
        options=ipopt_options
    )
    flag = "✓" if res.success else "~"
    print(f"  {flag} J={res.fun:.4e}  nit={res.nit}  success={res.success}")
    return res


def run_identification():
    # 1. Parametri e Generazione Dati (Modello Fisico Completo)
    N = 33
    omega_0 = 17000.0
    noise_level = 1e-4
    omega_t = np.linspace(16250.0, 17750.0, 200)

    print("\n[1/5] Generazione dati: modello a DISCO FLESSIBILE COMPLETO ...")
    result = generate_rotor_data('flexible', omega_t, omega_0, noise_level)
    y_measured   = result['y_measured']
    theta_true   = result['theta_true']
    M_dd_true, K_dd_true, K_bd_true = result['M_dd'], result['K_dd'], result['K_bd']

    # 2. Troncamento Modale (Generazione flessibilità non modellata)
    print("\n[2/5] Applicazione troncamento modale (max_nd=10) ...")
    M_dd_trunc, K_dd_trunc, K_bd_trunc = apply_modal_truncation(
        M_dd_true, K_dd_true, K_bd_true, N, max_nd=10
    )

    # 3. Scaling e Indici Parametri
    num_params = 5 * N - 1
    scales = np.ones(num_params)
    
    idx_d0 = 0               # Mistuning (lunghezza N)
    idx_a0 = N               # Rigidezza aerodinamica (lunghezza N-1)
    idx_a1 = 2 * N - 1       # Smorzamento aerodinamico (lunghezza N)
    idx_fr = 3 * N - 1       # Forzante reale (lunghezza N)
    idx_fi = 4 * N - 1       # Forzante immaginaria (lunghezza N)

    scales[idx_d0 : idx_a0] = 1e-2
    scales[idx_a0 : idx_fr] = 1e-3
    scales[idx_fr :]        = 1e-1

    # 4. Creazione dei Vincoli (Bounds) per Matrice Banded e Forzante Isolata
    print("\n[3/5] Creazione vincoli: matrice banded (+/- 3 pale) e forzante (ND=28) ...")
    bounds = [(-np.inf, np.inf)] * num_params

    # 4a. Vincoli Banded (B = 3) su a_0 e a_1
    B = 3
    for i in range(1, N):
        if i > B and i < N - B:
            bounds[idx_a0 + i - 1] = (0.0, 0.0)

    for i in range(N):
        if i > B and i < N - B:
            bounds[idx_a1 + i] = (0.0, 0.0)

    # 4b. Vincolo sulla Forzante (solo ND=28 è libero)
    target_nd = 28
    for k in range(N):
        if k != target_nd:
            bounds[idx_fr + k] = (0.0, 0.0)
            bounds[idx_fi + k] = (0.0, 0.0)
    
    # 4c. Vincoli Fisici di Stabilità
    bounds[idx_a1] = (1e-6, np.inf)

    np.random.seed(99)
    theta_initial_phys = 1e-5 * np.ones(num_params)
    
    for idx, (lower, upper) in enumerate(bounds):
        if lower == 0.0 and upper == 0.0:
            theta_initial_phys[idx] = 0.0
    
    theta_initial_scaled = theta_initial_phys / scales

    # 5. Identificazione FLESSIBILE CON MODELLO TRONCATO E VINCOLI MULTIPLI
    print("\n[4/5] Identificazione FLESSIBILE (Modello Troncato + Vincoli)...")
    ident_flex_trunc = FlexibleAeromechanicalIdentification(
        y_measured, omega_t, omega_0, N, M_dd_trunc, K_dd_trunc, K_bd_trunc)
    
    eval_flex = CachedEvaluator(ident_flex_trunc, scales)
    res_flex  = optimize(eval_flex, theta_initial_scaled.copy(), bounds, 'ID Troncato + Vincoli Rigorosi')
    theta_opt_flex = res_flex.x * scales

    # Coefficienti aerodinamici veri (Tabella 1) — curva smooth di riferimento
    a0_circ_true = np.zeros(N)
    a1_circ_true = np.zeros(N)
    a1_circ_true[0]   = 0.0050
    a0_circ_true[1]   = 0.00001;  a1_circ_true[1]   = 0.0010
    a0_circ_true[2]   = 0.00000;  a1_circ_true[2]   = 0.0015
    a0_circ_true[N-2] = 0.00000;  a1_circ_true[N-2] = 0.0020
    a0_circ_true[N-1] = 0.00000;  a1_circ_true[N-1] = 0.0025
    A_tw_true_smooth = np.fft.fft(a0_circ_true + 1j * a1_circ_true)

    _, a0_opt, a1_opt, _, _ = ident_flex_trunc.unpack_theta(theta_opt_flex)
    A_tw_opt = np.fft.fft(a0_opt + 1j * a1_opt)
    print(f"\n  Verifica artefatto — dev. std. rispetto alla curva smooth vera:")
    print(f"    std(Real ΔA_tw) = "
          f"{np.std(np.real(A_tw_opt) - np.real(A_tw_true_smooth)):.4e}")
    print(f"    std(Imag ΔA_tw) = "
          f"{np.std(np.imag(A_tw_opt) - np.imag(A_tw_true_smooth)):.4e}")
    print(f"  (ordini >> noise={noise_level:.0e} indicano artefatto significativo)\n")

    # 6. Analisi statistica e Grafici
    print("\n[5/5] Calcolo Cramér-Rao Bounds e Generazione grafici...")
    sigma_theta_flex = ident_flex_trunc.compute_cramer_rao(theta_opt_flex, noise_level**2)
    
    plotter_flex = AeroPlotter(
        sim_problem=ident_flex_trunc,
        omega_t=omega_t,
        y_measured=y_measured,
        N=N,
        theta_true=theta_true,
        theta_opt=theta_opt_flex,
        theta_initial=theta_initial_phys,
        sigma_theta=sigma_theta_flex
    )

    plotter_flex.plot_response_envelope()      
    plotter_flex.plot_parameters_2x2()         
    plotter_flex.plot_artifact_comparison()    

    return theta_opt_flex

if __name__ == '__main__':
    run_identification()