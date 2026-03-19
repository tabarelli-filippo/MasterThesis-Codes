"""
plot_aero.py
============
Unified plotter for aeromechanical identification results.

Replaces the separate rigid-disk (plot_aero.py) and flexible-disk
(plot_aero_flex.py) plotters with a single class that auto-detects
which model is in use.

Usage
-----
    plotter = AeroPlotter(
        sim_problem      = ident_problem,   # AeromechanicalIdentification or
                                            # FlexibleAeromechanicalIdentification
        omega_t          = omega_t,
        y_measured       = y_measured,
        N                = N,
        theta_true       = theta_true,
        theta_opt        = theta_opt,
        theta_initial    = theta_initial,
        # Optional:
        sigma_theta      = sigma_theta,       # CRB std-devs; None → no error bars
        A_tw_true_smooth = A_tw_true_smooth,  # smooth reference curve (rigid-on-flex)
    )

    plotter.plot_response_envelope()
    plotter.plot_parameters_2x2()
    plotter.plot_artifact_comparison()   # auto-selects smooth vs theta_true reference
    plotter.plot_all()                   # all four plots

Differences from the original plotters
---------------------------------------
- plot_response_envelope: uses np.linalg.solve with the precomputed D_disk[t]
  when available (flexible disk), or D_disk = 0 (rigid disk).  The old
  eigendecomposition shortcut is not used here because it ignores the disk term.
- plot_parameters_2x2: optional ±2σ error bars from sigma_theta; a visual floor
  of 5 % of the panel range is applied for readability.
- plot_artifact_comparison: uses A_tw_true_smooth when provided (rigid-on-flex
  case, shows the identification artefact); otherwise derives the reference from
  theta_true via unpack_theta (flex-on-flex or rigid case).
- All methods accept save=True/False and return the matplotlib Figure.
"""

import numpy as np
import matplotlib.pyplot as plt


class AeroPlotter:
    """
    Unified plotter for rigid- and flexible-disk identification results.

    Parameters
    ----------
    sim_problem : object
        Identification problem instance with methods ``unpack_theta``,
        ``build_system_matrices``, and attributes ``omega_0``, ``I_N``.
        If it also has a ``D_disk`` attribute (shape (T, N, N)), the
        flexible-disk correction is included in the response plots.
    omega_t : array-like, shape (T,)
    y_measured : array-like, shape (T, N), complex
    N : int
        Number of blades.
    theta_true : array-like, shape (5N-1,)
    theta_opt : array-like, shape (5N-1,)
    theta_initial : array-like, shape (5N-1,)
    sigma_theta : array-like, shape (5N-1,) or None
        Parameter standard deviations from the Cramér-Rao bound.
        If None, error bars are omitted.
    A_tw_true_smooth : array-like, shape (N,), complex or None
        Smooth TW reference curve for the artefact comparison plot.
        Relevant when a rigid-disk model is identified on flexible-disk data.
        If None, the reference is derived from theta_true.
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
    # Private helpers
    # ------------------------------------------------------------------

    def _compute_max_response_identified(self):
        """
        Compute the max-amplitude envelope of the identified model.

        For the flexible-disk model (sim_problem.D_disk present), the
        precomputed disk term is included.  For the rigid-disk model it
        is set to zero.

        Returns
        -------
        max_resp : (T,) float
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
            Z_t        = A - gamma_sq * self.sim_problem.I_N - D
            x_t        = np.linalg.solve(Z_t, f)
            max_resp[t] = float(np.max(np.abs(x_t)))
        return max_resp

    def _get_A_tw_true(self):
        """
        Return the reference TW aerodynamic coefficient curve.

        Priority: A_tw_true_smooth (if provided) → derived from theta_true.
        """
        if self.A_tw_true_smooth is not None:
            return self.A_tw_true_smooth
        _, a0_true, a1_true, _, _ = self.sim_problem.unpack_theta(self.theta_true)
        return np.fft.fft(a0_true + 1j * a1_true)

    # ------------------------------------------------------------------
    # Public plots
    # ------------------------------------------------------------------

    def plot_response_envelope(self, save=True, filename='blade_response_envelope.png'):
        """
        Plot the frequency-response envelope: measured data vs identified model.

        Parameters
        ----------
        save     : bool   — write the figure to disk
        filename : str    — output file path (PNG)

        Returns
        -------
        fig : matplotlib.figure.Figure
        """
        max_resp_data = np.max(np.abs(self.y_measured), axis=1)
        max_resp_id   = self._compute_max_response_identified()

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(self.omega_t, np.abs(self.y_measured),
                color='lightgray', linewidth=0.5, alpha=0.4)
        ax.plot(self.omega_t, max_resp_data, 'k-',  linewidth=1.5,
                label='Max response — data')
        ax.plot(self.omega_t, max_resp_id,   'r--', linewidth=1.5,
                label='Max response — identified model')
        ax.set_title('Blade Vibration Response Envelope')
        ax.set_xlabel('Frequency [rad/s]')
        ax.set_ylabel('Amplitude')
        ax.legend()
        ax.grid(True, linestyle=':', alpha=0.6)
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_parameters_2x2(self, tol=1e-10, save=True,
                             filename='identified_parameters_2x2.png'):
        """
        2×2 parameter comparison: true (yellow bars), initial (blue), MLE (red).

        Optional ±2σ error bars are shown when sigma_theta is available.
        Parameters unchanged from initial guess (|init - opt| < tol) are
        highlighted with green circles.

        Parameters
        ----------
        tol      : float  — tolerance for detecting fixed parameters
        save     : bool
        filename : str

        Returns
        -------
        fig : matplotlib.figure.Figure
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
            print("--- Cramér-Rao Bound diagnostics ---")
            print(f"  CRB 2σ  d0: {2*sig_d0.max():.2e}  "
                  f"a0: {2*sig_a0[1:].max():.2e}  "
                  f"a1: {2*sig_a1.max():.2e}  "
                  f"fr: {2*sig_fr.max():.2e}")

        nd_indices = np.fft.fftfreq(self.N, 1.0 / self.N)
        sort_idx   = np.argsort(nd_indices)

        fig, axs   = plt.subplots(2, 2, figsize=(15, 10))
        blade_idx  = np.arange(self.N)
        a0_idx     = np.arange(1, self.N)

        def _panel(ax, x_bar, y_bar, x_pt, y_init, y_opt,
                   fixed_mask, yerr, title, xlabel):
            ax.bar(x_bar, y_bar, color='yellow', edgecolor='black',
                   alpha=0.6, label='Exact model')
            ax.plot(x_pt, y_init, 'bo', markersize=4, label='Initial guess')
            if has_crb and yerr is not None:
                ax.errorbar(x_pt, y_opt, yerr=yerr,
                            fmt='r.', capsize=3, markersize=8,
                            label='MLE estimate (±2σ)')
            else:
                ax.plot(x_pt, y_opt, 'r.', markersize=8, label='MLE estimate')
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
               'Stiffness Mistuning ($d_0$)', 'Blade number')

        _panel(axs[0, 1],
               a0_idx, a0_true[1:], a0_idx, a0_init[1:], a0_opt[1:],
               fixed_a0, (2.0 * sig_a0[1:]) if has_crb else None,
               'Aerodynamic Stiffness Influence ($a_0$)', 'Influence index $i$')

        _panel(axs[1, 0],
               blade_idx, a1_true, blade_idx, a1_init, a1_opt,
               fixed_a1, (2.0 * sig_a1) if has_crb else None,
               'Aerodynamic Damping Influence ($a_1$)', 'Influence index $i$')

        nd_s       = nd_indices[sort_idx]
        fixed_fr_s = fixed_fr[sort_idx]
        _panel(axs[1, 1],
               nd_s, fr_true[sort_idx], nd_s, fr_init[sort_idx], fr_opt[sort_idx],
               fixed_fr_s, (2.0 * sig_fr[sort_idx]) if has_crb else None,
               'Forcing function magnitude (real part)', 'Nodal diameter')

        fig.suptitle('Identified Aeromechanical Parameters', fontsize=13, y=1.01)
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_influence_coefficients_tw(self, save=True,
                                        filename='influence_coefficients_tw.png'):
        """
        Plot the identified aerodynamic influence coefficients in TW space.

        Parameters
        ----------
        save     : bool
        filename : str

        Returns
        -------
        fig : matplotlib.figure.Figure
        """
        _, a0_opt, a1_opt, _, _ = self.sim_problem.unpack_theta(self.theta_opt)
        A_tw  = np.fft.fft(a0_opt + 1j * a1_opt)
        nd_idx = np.arange(self.N)

        fig, axs = plt.subplots(1, 2, figsize=(14, 5))
        for ax, data, ylabel in [
            (axs[0], np.real(A_tw),
             'Transformed influence coefficients — Real($a$)'),
            (axs[1], np.imag(A_tw),
             'Transformed influence coefficients — Imag($a$)'),
        ]:
            ax.plot(nd_idx, data, 'ro-', markerfacecolor='none',
                    markersize=6, linewidth=1.5)
            ax.axhline(0, color='black', linewidth=0.8)
            ax.set_xlabel('Nodal diameter')
            ax.set_ylabel(ylabel)
            ax.set_xticks(np.arange(0, self.N, 2))
            ax.grid(True, linestyle=':', alpha=0.6)

        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_artifact_comparison(self, save=True,
                                  filename='aerodynamic_artifact.png'):
        """
        Compare true and identified TW aerodynamic coefficients.

        - If A_tw_true_smooth is provided (rigid-on-flex scenario): the reference
          is the smooth physical coefficients, highlighting the artefact introduced
          by the model mismatch.
        - Otherwise (flex-on-flex or rigid scenario): the reference is derived from
          theta_true via unpack_theta.

        Parameters
        ----------
        save     : bool
        filename : str

        Returns
        -------
        fig : matplotlib.figure.Figure
        """
        A_tw_true = self._get_A_tw_true()

        _, a0_opt, a1_opt, _, _ = self.sim_problem.unpack_theta(self.theta_opt)
        A_tw_opt = np.fft.fft(a0_opt + 1j * a1_opt)

        if self.A_tw_true_smooth is not None:
            label_true = 'True aerodynamic coefficients (smooth)'
            label_opt  = 'Identified — rigid-disk model (artefact)'
            title      = 'Aerodynamic artefact: rigid model identified on flexible-disk data'
        else:
            label_true = 'True coefficients (smooth)'
            label_opt  = 'Identified (flexible model)'
            title      = 'TW aerodynamic coefficients: true vs identified'

        nd_idx = np.arange(self.N)
        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        for ax, true_data, opt_data, ylabel in [
            (axs[0], np.real(A_tw_true), np.real(A_tw_opt),
             'Transformed influence coefficients — Real($a$)'),
            (axs[1], np.imag(A_tw_true), np.imag(A_tw_opt),
             'Transformed influence coefficients — Imag($a$)'),
        ]:
            ax.plot(nd_idx, true_data, 'b-', linewidth=2, label=label_true)
            ax.plot(nd_idx, opt_data,  'ro-', markerfacecolor='none',
                    markersize=6, linewidth=1.5, label=label_opt)
            ax.axhline(0, color='black', linewidth=0.8)
            ax.set_xlabel('Nodal diameter')
            ax.set_ylabel(ylabel)
            ax.set_xticks(np.arange(0, self.N, 2))
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.legend()

        fig.suptitle(title, fontsize=14)
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")
        plt.show()
        return fig

    # ------------------------------------------------------------------

    def plot_all(self):
        """Generate all four plots in sequence."""
        self.plot_response_envelope()
        self.plot_parameters_2x2()
        self.plot_influence_coefficients_tw()
        self.plot_artifact_comparison()


# ---------------------------------------------------------------------------
# Backward-compatible aliases
# ---------------------------------------------------------------------------

AeromechanicalPlotter = AeroPlotter   # replaces plot_aero.AeromechanicalPlotter
FlexibleAeroPlotter   = AeroPlotter   # replaces plot_aero_flex.FlexibleAeroPlotter
