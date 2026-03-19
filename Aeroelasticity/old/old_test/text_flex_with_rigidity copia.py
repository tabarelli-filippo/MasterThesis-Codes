"""
test_flex_with_rigidity.py
============
Dimostra l'artefatto matematico che si produce quando si usa un modello di
identificazione a DISCO RIGIDO per dati generati con un rotore a DISCO FLESSIBILE.

Flusso:
  1. Dati simulati  →  initialize_flexible_rotor_complete()
       - modello generativo coerente con il modello identificativo (stessa normalizzazione)
       - restituisce anche i coefficienti aerodinamici *veri* (smooth nel dominio TW)
         separati da theta_true, per il confronto diretto nel grafico degli artefatti
  2. Identificazione →  AeromechanicalIdentification  (modello rigido)
  3. Grafici         →  AeromechanicalPlotter
       - plot_response_envelope           : fit della risposta
       - plot_parameters_2x2              : mistuning, a0, a1, forzante
       - plot_artifact_comparison_phys    : artefatto in coordinate TW
                                           (curva smooth vera vs. identificata seghettata)
"""

import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.aero_ident import AeromechanicalIdentification
from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter

class CachedEvaluator:
    """Wrapper con scaling per l'ottimizzazione scipy."""

    def __init__(self, ident_problem, scale_factors):
        self.problem = ident_problem
        self.scales  = np.asarray(scale_factors, dtype=np.float64)
        self._cache_key  = None
        self._cache_J    = None
        self._cache_grad = None

    def _update(self, theta_scaled):
        key = theta_scaled.tobytes()
        if key != self._cache_key:
            theta_phys = theta_scaled * self.scales
            J, g_phys  = self.problem.compute_objective_and_gradient(theta_phys)
            self._cache_J    = J
            self._cache_grad = g_phys * self.scales
            self._cache_key  = key

    def objective(self, theta_scaled):
        self._update(theta_scaled)
        return self._cache_J

    def gradient(self, theta_scaled):
        self._update(theta_scaled)
        return self._cache_grad


def isolate_aerodynamic_artifact():
    # ------------------------------------------------------------------
    # 1. Parametri operativi
    # ------------------------------------------------------------------
    N           = 33
    omega_0     = 17000.0
    T           = 200
    omega_t     = np.linspace(16250.0, 17750.0, T)
    noise_level = 1e-4

    print("=" * 65)
    print("  Artefatto aerodinamico da disco flessibile — Hall & Hall 2024")
    print("=" * 65)

    # ------------------------------------------------------------------
    # 2. Generazione dati con modello a DISCO FLESSIBILE
    #    initialize_flexible_rotor_complete usa la stessa normalizzazione
    #    del modello identificativo (gamma^2 = omega^2/omega_0^2).
    # ------------------------------------------------------------------
    print("\n[1/4] Generazione dati: modello a DISCO FLESSIBILE ...")
    result = generate_rotor_data('flexible', omega_t, omega_0, noise_level)
    y_measured   = result['y_measured']
    theta_true   = result['theta_true']

    # ------------------------------------------------------------------
    # Coefficienti aerodinamici FISICI (Tabella 1 del paper).
    # Nel dominio TW producono una curva smooth (bassa energia alle alte ND),
    # che è la firma di un'influenza a corto raggio.
    # L'identificatore rigido non può distinguere questa influenza dalla
    # dinamica del disco → produce una curva seghettata (l'artefatto).
    # ------------------------------------------------------------------
    a0_circ_true = np.zeros(N)
    a1_circ_true = np.zeros(N)
    a1_circ_true[0]   = 0.0050   # auto-smorzamento
    a0_circ_true[1]   = 0.00001;  a1_circ_true[1]   = 0.0010
    a0_circ_true[2]   = 0.00000;  a1_circ_true[2]   = 0.0015
    a0_circ_true[N-2] = 0.00000;  a1_circ_true[N-2] = 0.0020
    a0_circ_true[N-1] = 0.00000;  a1_circ_true[N-1] = 0.0025

    # Trasformata TW vera: curva smooth di riferimento per il grafico artefatto
    A_tw_true_smooth = np.fft.fft(a0_circ_true + 1j * a1_circ_true)

    # ------------------------------------------------------------------
    # 3. Modello identificativo a DISCO RIGIDO
    # ------------------------------------------------------------------
    print("[2/4] Inizializzazione modello identificativo: DISCO RIGIDO ...")
    ident_problem = AeromechanicalIdentification(y_measured, omega_t, omega_0, N)

    num_params = 5 * N - 1
    idx_d0 = 0
    idx_a0 = N
    idx_fr = N + (N - 1) + N

    # Scaling: migliora il condizionamento numerico.
    scale_factors = np.ones(num_params)
    scale_factors[idx_d0 : idx_a0] = 1e-2   # mistuning  O(1e-2)
    scale_factors[idx_a0 : idx_fr] = 1e-3   # aero       O(1e-3)
    scale_factors[idx_fr :]        = 1e-1   # forzante   O(1e-1)

    np.random.seed(99)
    theta_initial        = 1e-5 * np.ones(num_params)
    theta_initial_scaled = theta_initial / scale_factors

    evaluator = CachedEvaluator(ident_problem, scale_factors)

    # ------------------------------------------------------------------
    # 4. Ottimizzazione L-BFGS-B (nessun vincolo)
    #    Il solutore sovradatta liberamente, assorbendo la dinamica del
    #    disco flessibile nei coefficienti aerodinamici identificati.
    # ------------------------------------------------------------------
    print(f"[3/4] Ottimizzazione L-BFGS-B ({num_params} parametri) ...")
    ipopt_options = {
        'max_iter': 20000, 
        'tol': 1e-9,
        'acceptable_tol': 1e-6,
        'print_level': 5,
        'hessian_approximation': 'limited-memory',
        'limited_memory_max_history': 1,
        'limited_memory_update_type': 'sr1',
    }

    res = minimize_ipopt(
        fun=evaluator.objective,
        x0=theta_initial_scaled,
        jac=evaluator.gradient,
        bounds=None,
        options=ipopt_options
    )

    theta_opt_full = res.x * scale_factors

    # ------------------------------------------------------------------
    # 5. Diagnostica overfitting
    # ------------------------------------------------------------------
    J_opt  = ident_problem.objective(theta_opt_full)
    J_true = ident_problem.objective(theta_true)

    print(f"\n{'─'*55}")
    print(f"  J finale (modello rigido ottimizzato) : {J_opt:.4e}")
    print(f"  J con parametri veri su modello rigido: {J_true:.4e}")
    print(f"  Rapporto J_true / J_opt               : {J_true / J_opt:.2f}")
    print(f"  (rapporto >> 1 conferma l'overfitting)")
    print(f"{'─'*55}")

    _, a0_opt, a1_opt, _, _ = ident_problem.unpack_theta(theta_opt_full)
    A_tw_opt = np.fft.fft(a0_opt + 1j * a1_opt)
    print(f"\n  Verifica artefatto — dev. std. rispetto alla curva smooth vera:")
    print(f"    std(Real ΔA_tw) = "
          f"{np.std(np.real(A_tw_opt) - np.real(A_tw_true_smooth)):.4e}")
    print(f"    std(Imag ΔA_tw) = "
          f"{np.std(np.imag(A_tw_opt) - np.imag(A_tw_true_smooth)):.4e}")
    print(f"  (ordini >> noise={noise_level:.0e} indicano artefatto significativo)\n")

    # ------------------------------------------------------------------
    # 6. Cramér-Rao bounds
    # ------------------------------------------------------------------
    sigma_theta = ident_problem.compute_cramer_rao_bounds(theta_opt_full, noise_level**2)

    # ------------------------------------------------------------------
    # 7. Grafici
    # ------------------------------------------------------------------
    print("[4/4] Generazione grafici ...")
    plotter = AeroPlotter(
        sim_problem=ident_problem,
        omega_t=omega_t,
        y_measured=y_measured,
        N=N,
        theta_true=theta_true,
        theta_opt=theta_opt_full,
        theta_initial=theta_initial,
        sigma_theta=sigma_theta,
        A_tw_true_smooth=A_tw_true_smooth,
    )

    plotter.plot_response_envelope()
    plotter.plot_parameters_2x2()
    plotter.plot_artifact_comparison_phys()

    print("\nDone.")


if __name__ == '__main__':
    isolate_aerodynamic_artifact()