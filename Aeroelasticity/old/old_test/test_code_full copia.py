import numpy as np
from cyipopt import minimize_ipopt

from auxSolver.aero_ident import AeromechanicalIdentification
from auxSolver.data_generation import generate_rotor_data
from auxSolver.plot_aero import AeroPlotter


class CachedEvaluator:
    """Wrapper con cache e scalatura per l'ottimizzazione a dimensionalità completa."""

    def __init__(self, ident_problem, scale_factors):
        self.problem = ident_problem
        self.scales  = np.asarray(scale_factors)

        self.last_theta_scaled = None
        self.last_J            = None
        self.last_grad_scaled  = None

    def _update(self, theta_scaled):
        if self.last_theta_scaled is None or not np.array_equal(theta_scaled, self.last_theta_scaled):
            theta_physical = theta_scaled * self.scales
            J, grad_physical = self.problem.compute_objective_and_gradient(theta_physical)
            self.last_J            = J
            self.last_grad_scaled  = grad_physical * self.scales
            self.last_theta_scaled = np.copy(theta_scaled)

    def objective(self, theta_scaled):
        self._update(theta_scaled)
        return self.last_J

    def gradient(self, theta_scaled):
        self._update(theta_scaled)
        return self.last_grad_scaled

    def objective_and_gradient(self, theta_scaled):
        self._update(theta_scaled)
        return self.last_J, self.last_grad_scaled.copy()


def run_identification_with_simulated_data():
    # 1. Parametri operativi
    N           = 33
    omega_0     = 17000.0
    T           = 250
    omega_t     = np.linspace(16250.0, 17750.0, T)
    noise_level = 1e-4

    print("Generazione dei dati sintetici (disco rigido)...")
    result     = generate_rotor_data('rigid', omega_t, omega_0, noise_level, N=N)
    y_measured = result['y_measured']
    theta_true = result['theta_true']

    ident_problem = AeromechanicalIdentification(y_measured, omega_t, omega_0, N)

    num_params = 5 * N - 1

    # 2. Offset nel vettore theta
    idx_d0 = 0
    idx_a0 = N
    idx_a1 = N + (N - 1)
    idx_fr = N + (N - 1) + N

    # 3. Bounds fisici
    lb = np.empty(num_params)
    ub = np.empty(num_params)

    lb[idx_d0 : idx_a0] = -0.10;  ub[idx_d0 : idx_a0] =  0.10   # mistuning
    lb[idx_a0 : idx_a1] = -0.005; ub[idx_a0 : idx_a1] =  0.005  # a0
    lb[idx_a1 : idx_fr] = -0.010; ub[idx_a1 : idx_fr] =  0.010  # a1
    lb[idx_a1]          =  1e-6;  ub[idx_a1]          =  0.020  # auto-smorzamento >= 0
    lb[idx_fr :]        = -1.0;   ub[idx_fr :]        =  1.0    # forzante

    # 4. Punto iniziale
    np.random.seed(99)
    theta_initial = 1e-5 * np.ones(num_params)

    # 5. Fattori di scala per condizionamento
    scale_factors = np.ones(num_params)
    scale_factors[idx_d0 : idx_a0] = 1e-2   # mistuning  O(1e-2)
    scale_factors[idx_a0 : idx_fr] = 1e-3   # a0, a1     O(1e-3)
    scale_factors[idx_fr :]        = 1e-1   # forzante   O(1e-1)

    theta_initial_scaled = theta_initial / scale_factors
    lb_scaled            = lb / scale_factors
    ub_scaled            = ub / scale_factors
    bounds_scipy         = list(zip(lb_scaled, ub_scaled))

    evaluator = CachedEvaluator(ident_problem, scale_factors)

    J_init = evaluator.objective(theta_initial_scaled)
    print(f"J iniziale : {J_init:.6e}")
    print(f"J ottimale : {ident_problem.objective(theta_true):.6e}")

    # 6. Ottimizzazione IPOPT a dimensionalità completa (164 parametri)
    print(f"\nAvvio ottimizzazione IPOPT ({num_params} parametri, scalati)...")

    ipopt_options = {
        'max_iter': 5000,
        'tol': 1e-9,
        'acceptable_tol': 1e-8,
        'print_level': 5,
        'hessian_approximation': 'limited-memory',
        'limited_memory_max_history': 20,
        'limited_memory_update_type': 'sr1',
        'bound_push': 1e-3,
    }

    res = minimize_ipopt(
        fun=evaluator.objective,
        x0=theta_initial_scaled,
        jac=evaluator.gradient,
        bounds=bounds_scipy,
        options=ipopt_options,
    )

    # 7. Risultato
    print(f"\nStato          : {res.message}")
    print(f"Iterazioni     : {res.nit}  |  Valutazioni: {res.nfev}")
    print(f"J finale       : {res.fun:.6e}")

    theta_opt_full = res.x * scale_factors

    # 8. Cramér-Rao bounds e grafici
    print("\nCalcolo Cramér-Rao bounds...")
    sigma_theta = ident_problem.compute_cramer_rao_bounds(theta_opt_full, noise_level**2)

    print("Generazione dei grafici...")
    plotter = AeroPlotter(
        sim_problem   = ident_problem,
        omega_t       = omega_t,
        y_measured    = y_measured,
        N             = N,
        theta_true    = theta_true,
        theta_opt     = theta_opt_full,
        theta_initial = theta_initial,
        sigma_theta   = sigma_theta,
    )

    plotter.plot_response_envelope()
    plotter.plot_parameters_2x2()
    plotter.plot_influence_coefficients_tw()


if __name__ == '__main__':
    run_identification_with_simulated_data()