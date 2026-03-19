"""
plot_aero.py
====================
Plotter unificato per l'identificazione aeromeccanica.
Sostituisce sia plot_aero.py (disco rigido) che plot_aero_flex.py (disco flessibile).

Utilizzo
--------
    plotter = AeroPlotter(
        sim_problem   = ident_problem,    # AeromechanicalIdentification o
                                          # FlexibleAeromechanicalIdentification
        omega_t       = omega_t,
        y_measured    = y_measured,
        N             = N,
        theta_true    = theta_true,
        theta_opt     = theta_opt,
        theta_initial = theta_initial,
        # Opzionali:
        sigma_theta      = sigma_theta,       # CRB; se None, nessuna barra d'errore
        A_tw_true_smooth = A_tw_true_smooth,  # curva smooth vera (caso rigid-on-flex)
    )

    plotter.plot_response_envelope()
    plotter.plot_parameters_2x2()
    plotter.plot_artifact_comparison()   # autodetect: smooth vs theta_true
    plotter.plot_all()                   # tutti e quattro

Differenze rispetto ai plotter originali
-----------------------------------------
- plot_response_envelope: usa sempre np.linalg.solve con D_disk[t] se disponibile
  (disc flessibile), altrimenti D_disk = 0 (disco rigido). Nessuna scorciatoia
  eigendecomposition, che ignorava il termine disco.
- plot_parameters_2x2: barre d'errore ±2σ opzionali; floor visivo del 5%
  dell'escursione del pannello per leggibilità.
- plot_artifact_comparison: se A_tw_true_smooth è fornita usa quella (caso
  rigid-on-flex, mostra l'artefatto rispetto ai coefficienti fisici veri);
  altrimenti deriva la curva vera da theta_true via unpack_theta (caso flex-on-flex).
- Tutti i metodi accettano save=True/False e restituiscono la figura matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt


class AeroPlotter:
    """
    Plotter unificato: funziona con modelli rigidi e flessibili.

    Parametri
    ----------
    sim_problem      : oggetto con i metodi unpack_theta, build_system_matrices,
                       omega_0, I_N.  Opzionalmente D_disk (array (T,N,N)) per
                       il disco flessibile.
    omega_t          : array (T,)
    y_measured       : array (T, N) complesso
    N                : int — numero di pale
    theta_true       : array (5N-1,)
    theta_opt        : array (5N-1,)
    theta_initial    : array (5N-1,)
    sigma_theta      : array (5N-1,) o None — std dev CRB
    A_tw_true_smooth : array (N,) complesso o None — curva smooth di riferimento
                       per il grafico artefatto nel caso rigid-on-flex.
                       Se None, la curva viene derivata da theta_true.
    """

    def __init__(self, sim_problem, omega_t, y_measured, N,
                 theta_true, theta_opt, theta_initial,
                 sigma_theta=None, A_tw_true_smooth=None):
        self.sim_problem      = sim_problem
        self.omega_t          = np.asarray(omega_t)
        self.y_measured       = np.asarray(y_measured)
        self.N                = N
        self.theta_true       = np.asarray(theta_true)
        self.theta_opt        = np.asarray(theta_opt)
        self.theta_initial    = np.asarray(theta_initial)
        self.sigma_theta      = np.asarray(sigma_theta) if sigma_theta is not None else None
        self.A_tw_true_smooth = (np.asarray(A_tw_true_smooth)
                                 if A_tw_true_smooth is not None else None)

    # ------------------------------------------------------------------
    # Helper interni
    # ------------------------------------------------------------------

    def _compute_max_response_identified(self):
        """
        Calcola inviluppo max(|x_t|) del modello identificato.
        Supporta sia disco rigido (D_disk assente → zero) che flessibile
        (D_disk[t] precomputed sull'oggetto sim_problem).
        """
        d_0, a_0_circ, a_1_circ, f_r, f_i = \
            self.sim_problem.unpack_theta(self.theta_opt)
        A, f   = self.sim_problem.build_system_matrices(d_0, a_0_circ, a_1_circ, f_r, f_i)
        N      = A.shape[0]
        T      = len(self.omega_t)
        has_disk = hasattr(self.sim_problem, 'D_disk')

        max_resp = np.zeros(T)
        for t in range(T):
            gamma_sq = (self.omega_t[t] / self.sim_problem.omega_0) ** 2
            D = self.sim_problem.D_disk[t] if has_disk else np.zeros((N, N), dtype=np.complex128)
            Z_t       = A - gamma_sq * self.sim_problem.I_N - D
            x_t       = np.linalg.solve(Z_t, f)
            max_resp[t] = float(np.max(np.abs(x_t)))
        return max_resp

    def _get_A_tw_true(self):
        """
        Restituisce la curva TW vera di riferimento.
        Priorità: A_tw_true_smooth passata al costruttore → derivata da theta_true.
        """
        if self.A_tw_true_smooth is not None:
            return self.A_tw_true_smooth
        _, a0_true, a1_true, _, _ = self.sim_problem.unpack_theta(self.theta_true)
        return np.fft.fft(a0_true + 1j * a1_true)

    # ------------------------------------------------------------------
    # Grafici pubblici
    # ------------------------------------------------------------------

    def plot_response_envelope(self, save=True, filename='ampiezza_vibrazione.png'):
        """
        Inviluppo della risposta in frequenza: dati misurati vs modello identificato.
        """
        max_resp_data = np.max(np.abs(self.y_measured), axis=1)
        max_resp_id   = self._compute_max_response_identified()

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(self.omega_t, np.abs(self.y_measured),
                color='lightgray', linewidth=0.5, alpha=0.4)
        ax.plot(self.omega_t, max_resp_data, 'k-',  linewidth=1.5,
                label='Max response data')
        ax.plot(self.omega_t, max_resp_id,   'r--', linewidth=1.5,
                label='Max response ID model')
        ax.set_title('Blade Vibration Response Envelope')
        ax.set_xlabel('Frequency [rad/sec]')
        ax.set_ylabel('Amplitude')
        ax.legend()
        ax.grid(True, linestyle=':', alpha=0.6)
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Salvato: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_parameters_2x2(self, tol=1e-10, save=True,
                             filename='parametri_identificati_2x2.png'):
        """
        Confronto parametri: veri (barre gialle), iniziale (blu), MLE (rosso).
        Le barre d'errore ±2σ sono mostrate solo se sigma_theta è disponibile.
        I parametri rimasti fissi (|init - opt| < tol) sono cerchiati in verde.
        """
        unpack = self.sim_problem.unpack_theta
        d0_true, a0_true, a1_true, fr_true, _ = unpack(self.theta_true)
        d0_opt,  a0_opt,  a1_opt,  fr_opt,  _ = unpack(self.theta_opt)
        d0_init, a0_init, a1_init, fr_init, _ = unpack(self.theta_initial)

        fixed_d0 = np.isclose(d0_init,     d0_opt,     atol=tol)
        fixed_a0 = np.isclose(a0_init[1:], a0_opt[1:], atol=tol)
        fixed_a1 = np.isclose(a1_init,     a1_opt,     atol=tol)
        fixed_fr = np.isclose(fr_init,     fr_opt,     atol=tol)

        has_crb = self.sigma_theta is not None
        if has_crb:
            sig_d0, sig_a0, sig_a1, sig_fr, _ = unpack(self.sigma_theta)
            print("--- Diagnostica Cramér-Rao Bounds ---")
            print(f"  CRB 2σ  d0: {2*sig_d0.max():.2e}  "
                  f"a0: {2*sig_a0[1:].max():.2e}  "
                  f"a1: {2*sig_a1.max():.2e}  "
                  f"fr: {2*sig_fr.max():.2e}")

        nd_indices = np.fft.fftfreq(self.N, 1.0 / self.N)
        sort_idx   = np.argsort(nd_indices)

        fig, axs = plt.subplots(2, 2, figsize=(15, 10))
        blade_idx = np.arange(self.N)
        a0_idx    = np.arange(1, self.N)

        def _panel(ax, x_bar, y_bar, x_pt, y_init, y_opt,
                   fixed_mask, yerr, title, xlabel):
            ax.bar(x_bar, y_bar, color='yellow', edgecolor='black',
                   alpha=0.6, label='Exact Model')
            ax.plot(x_pt, y_init, 'bo', markersize=4, label='Initial Guess')
            if has_crb and yerr is not None:
                ax.errorbar(x_pt, y_opt, yerr=yerr,
                            fmt='r.', capsize=3, markersize=8,
                            label='MLE Estimate (±2σ)')
            else:
                ax.plot(x_pt, y_opt, 'r.', markersize=8, label='MLE Estimate')
            if np.any(fixed_mask):
                ax.scatter(x_pt[fixed_mask], y_opt[fixed_mask], s=100,
                           facecolors='none', edgecolors='g',
                           linewidths=1.5, label='Fixed values')
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.legend(fontsize=8)
            ax.grid(True, linestyle=':', alpha=0.5)

        _panel(axs[0, 0],
               blade_idx, d0_true, blade_idx, d0_init, d0_opt,
               fixed_d0, (2.0 * sig_d0) if has_crb else None,
               'Stiffness Mistuning ($d_0$)', 'Blade Number')

        _panel(axs[0, 1],
               a0_idx, a0_true[1:], a0_idx, a0_init[1:], a0_opt[1:],
               fixed_a0, (2.0 * sig_a0[1:]) if has_crb else None,
               'Aerodynamic Stiffness Influence ($a_0$)', 'Influence Index $i$')

        _panel(axs[1, 0],
               blade_idx, a1_true, blade_idx, a1_init, a1_opt,
               fixed_a1, (2.0 * sig_a1) if has_crb else None,
               'Aerodynamic Damping Influence ($a_1$)', 'Influence Index $i$')

        nd_s      = nd_indices[sort_idx]
        fixed_fr_s = fixed_fr[sort_idx]
        _panel(axs[1, 1],
               nd_s, fr_true[sort_idx], nd_s, fr_init[sort_idx], fr_opt[sort_idx],
               fixed_fr_s, (2.0 * sig_fr[sort_idx]) if has_crb else None,
               'Forcing Function Magnitude (Real Part)', 'Nodal Diameter')

        fig.suptitle('Identified Aeromechanical Parameters', fontsize=13, y=1.01)
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Salvato: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_influence_coefficients_tw(self, save=True,
                                        filename='coefficienti_influenza_nodali.png'):
        """
        Coefficienti di influenza aerodinamica identificati nello spazio TW.
        """
        _, a0_opt, a1_opt, _, _ = self.sim_problem.unpack_theta(self.theta_opt)
        A_tw = np.fft.fft(a0_opt + 1j * a1_opt)
        nd_idx = np.arange(self.N)

        fig, axs = plt.subplots(1, 2, figsize=(14, 5))
        for ax, data, ylabel in [
            (axs[0], np.real(A_tw),
             'Transformed Influence Coefficients, Real($a$)'),
            (axs[1], np.imag(A_tw),
             'Transformed Influence Coefficients, Imag($a$)'),
        ]:
            ax.plot(nd_idx, data, 'ro-', markerfacecolor='none',
                    markersize=6, linewidth=1.5)
            ax.axhline(0, color='black', linewidth=0.8)
            ax.set_xlabel('Nodal Diameter, N')
            ax.set_ylabel(ylabel)
            ax.set_xticks(np.arange(0, self.N, 2))
            ax.grid(True, linestyle=':', alpha=0.6)

        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Salvato: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_artifact_comparison(self, save=True,
                                  filename='artefatto_aerodinamico.png'):
        """
        Confronto tra coefficienti TW veri e identificati.

        - Se A_tw_true_smooth è disponibile (caso rigid-on-flex): la curva
          vera è quella dei coefficienti fisici smooth → mostra l'artefatto
          introdotto dal mismatch del modello.
        - Altrimenti (caso flex-on-flex o rigido): la curva vera è derivata
          da theta_true via unpack_theta.

        Il titolo del grafico si adatta automaticamente al caso rilevato.
        """
        A_tw_true = self._get_A_tw_true()

        _, a0_opt, a1_opt, _, _ = self.sim_problem.unpack_theta(self.theta_opt)
        A_tw_opt = np.fft.fft(a0_opt + 1j * a1_opt)

        if self.A_tw_true_smooth is not None:
            label_true = 'True Aerodynamic Coefficients (Smooth)'
            label_opt  = 'Identified — Rigid Disk Model (Artifact)'
            title      = 'Aerodynamic Artifact: Rigid ID Model on Flexible Disk Data'
        else:
            label_true = 'True Coefficients (Smooth)'
            label_opt  = 'Identified (Flexible Model)'
            title      = 'Comparison: artifact coefficients vs original solution'

        nd_idx = np.arange(self.N)
        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        for ax, true_data, opt_data, ylabel in [
            (axs[0], np.real(A_tw_true), np.real(A_tw_opt),
             'Transformed Influence Coefficients, Real($a$)'),
            (axs[1], np.imag(A_tw_true), np.imag(A_tw_opt),
             'Transformed Influence Coefficients, Imag($a$)'),
        ]:
            ax.plot(nd_idx, true_data, 'b-', linewidth=2, label=label_true)
            ax.plot(nd_idx, opt_data,  'ro-', markerfacecolor='none',
                    markersize=6, linewidth=1.5, label=label_opt)
            ax.axhline(0, color='black', linewidth=0.8)
            ax.set_xlabel('Nodal Diameter, N')
            ax.set_ylabel(ylabel)
            ax.set_xticks(np.arange(0, self.N, 2))
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.legend()

        fig.suptitle(title, fontsize=14)
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Salvato: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_all(self):
        """Genera tutti e quattro i grafici in sequenza."""
        self.plot_response_envelope()
        self.plot_parameters_2x2()
        self.plot_influence_coefficients_tw()
        self.plot_artifact_comparison()


# ---------------------------------------------------------------------------
# Alias di compatibilità backward
# ---------------------------------------------------------------------------

# Mantiene i vecchi nomi importabili senza modificare i test esistenti.
AeromechanicalPlotter  = AeroPlotter   # alias per plot_aero.AeromechanicalPlotter
FlexibleAeroPlotter    = AeroPlotter   # alias per plot_aero_flex.FlexibleAeroPlotter
