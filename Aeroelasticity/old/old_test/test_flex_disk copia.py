import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter
from auxSolver.flex_aero_ident import FlexibleAeromechanicalIdentification


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


def optimize(evaluator, x0_scaled, label=''):
    print(f"\n  Avvio {label} con IPOPT...")
    ipopt_options = {
        'max_iter': 2000,
        'tol': 1e-12,
        'acceptable_tol': 1e-10,
        'print_level': 5,
        'hessian_approximation': 'limited-memory',
        'limited_memory_max_history': 30,
        'limited_memory_update_type': 'bfgs',
    }
    res = minimize_ipopt(
        fun=evaluator.objective, x0=x0_scaled, jac=evaluator.gradient,
        bounds=None, options=ipopt_options
    )
    return res


def run_identification():
    # 1. Parametri e Generazione Dati
    N, omega_0, noise_level = 33, 17000.0, 1e-4
    omega_t = np.linspace(16250.0, 17750.0, 200)

    print("\n[1/4] Generazione dati: modello a DISCO FLESSIBILE ...")
    result = generate_rotor_data('flexible', omega_t, omega_0, noise_level)
    y_measured   = result['y_measured']
    theta_true   = result['theta_true']
    M_dd, K_dd, K_bd = result['M_dd'], result['K_dd'], result['K_bd']

    # 2. Scaling e Inizializzazione
    num_params = 5 * N - 1
    scales = np.ones(num_params)
    scales[0:N] = 1e-2           # d0
    scales[N:3*N-1] = 1e-3       # a0, a1
    scales[3*N-1:] = 1e-1        # fr, fi

    np.random.seed(99)
    theta_initial_phys = 1e-5 * np.ones(num_params)
    theta_initial_scaled = theta_initial_phys / scales

    # 3. Identificazione FLESSIBILE
    print("\n[2/4] Identificazione FLESSIBILE...")
    ident_flex = FlexibleAeromechanicalIdentification(
        y_measured, omega_t, omega_0, N, M_dd, K_dd, K_bd)
    
    eval_flex = CachedEvaluator(ident_flex, scales)
    res_flex  = optimize(eval_flex, theta_initial_scaled.copy(), 'ID Flessibile')
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

    _, a0_opt, a1_opt, _, _ = ident_flex.unpack_theta(theta_opt_flex)
    A_tw_opt = np.fft.fft(a0_opt + 1j * a1_opt)
    print(f"\n  Verifica artefatto — dev. std. rispetto alla curva smooth vera:")
    print(f"    std(Real ΔA_tw) = "
          f"{np.std(np.real(A_tw_opt) - np.real(A_tw_true_smooth)):.4e}")
    print(f"    std(Imag ΔA_tw) = "
          f"{np.std(np.imag(A_tw_opt) - np.imag(A_tw_true_smooth)):.4e}")
    print(f"  (ordini >> noise={noise_level:.0e} indicano artefatto significativo)\n")

    # 4. Analisi statistica (CRB) per il caso flessibile
    print("\n[3/4] Calcolo dei Cramér-Rao Bounds...")
    # Calcolo della varianza dei parametri basata sul rumore di misura
    sigma_theta = ident_flex.compute_cramer_rao(theta_opt_flex, noise_level**2)

    # 5. Generazione Plot con FlexibleAeroPlotter
    print("\n[4/4] Generazione grafici (Modello Flessibile)...")
    
    plotter_flex = AeroPlotter(
        sim_problem=ident_flex,
        omega_t=omega_t,
        y_measured=y_measured,
        N=N,
        theta_true=theta_true,
        theta_opt=theta_opt_flex,
        theta_initial=theta_initial_phys,
        sigma_theta= sigma_theta
    )
    
    # Esecuzione dei grafici specifici
    plotter_flex.plot_response_envelope()      # Verifica il fit della FRF
    plotter_flex.plot_parameters_2x2()         # Parametri identificati con barre 2sigma
    plotter_flex.plot_artifact_comparison()    # Verifica assenza di artefatti (smoothness)

    return theta_opt_flex

if __name__ == '__main__':
    run_identification()